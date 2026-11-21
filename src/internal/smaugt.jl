const SMAUGT_N = 256

struct SmaugTParameters
    mode::Union{Int,Symbol}
    k::Int
    log_q::Int
    log_p::Int
    log_p_prime::Int
    cbd_seed_bytes::Int
    message_bytes::Int
    secret_weight::Int
    public_key_bytes::Int
    pke_secret_key_bytes::Int
    secret_key_bytes::Int
    ciphertext_bytes::Int
end

function SmaugTParameters(mode, k, log_q, log_p, log_p_prime, cbd_seed_bytes,
                          message_bytes, secret_weight)
    public_key_bytes = 32 + k * log_q * SMAUGT_N ÷ 8
    pke_secret_key_bytes = k * SMAUGT_N ÷ 4
    return SmaugTParameters(
        mode, k, log_q, log_p, log_p_prime, cbd_seed_bytes, message_bytes,
        secret_weight, public_key_bytes, pke_secret_key_bytes,
        pke_secret_key_bytes + 32 + public_key_bytes,
        k * log_p * SMAUGT_N ÷ 8 + log_p_prime * SMAUGT_N ÷ 8,
    )
end

const SMAUGT_PARAMETERS = Dict(
    128 => SmaugTParameters(1, 2, 10, 8, 5, 96, 32, 70),
    192 => SmaugTParameters(3, 3, 11, 9, 4, 64, 32, 88),
    256 => SmaugTParameters(5, 4, 11, 9, 7, 128, 32, 87),
    :timer => SmaugTParameters(:timer, 2, 10, 8, 3, 96, 16, 70),
)

@inline function smaugt_set!(polynomial::Vector{Int16}, index::Int, value::Integer)
    polynomial[index + 1] = reinterpret(Int16, UInt16(value & 0xffff))
end

@inline smaugt_get(polynomial::AbstractVector, index::Int) = Int(polynomial[index + 1])
smaugt_zero_polynomial() = zeros(Int16, SMAUGT_N)

function smaugt_pack_bits(values, bits, shift)
    output = zeros(UInt8, length(values) * bits ÷ 8)
    mask = (1 << bits) - 1
    accumulator = 0
    available = 0
    position = 0
    for index in 0:(length(values) - 1)
        accumulator |= ((smaugt_get(values, index) >> shift) & mask) << available
        available += bits
        while available >= 8
            output[position + 1] = UInt8(accumulator & 0xff)
            position += 1
            accumulator >>= 8
            available -= 8
        end
    end
    return output
end

function smaugt_unpack_bits(bytes, count, bits, shift)
    output = zeros(Int16, count)
    mask = (1 << bits) - 1
    accumulator = 0
    available = 0
    position = 0
    for index in 0:(count - 1)
        while available < bits
            accumulator |= Int(bytes[position + 1]) << available
            position += 1
            available += 8
        end
        smaugt_set!(output, index, (accumulator & mask) << shift)
        accumulator >>= bits
        available -= bits
    end
    return output
end

smaugt_pack_polynomial(polynomial, bits, shift = 0) =
    smaugt_pack_bits(polynomial, bits, shift)
smaugt_unpack_polynomial(bytes, bits, shift = 0) =
    smaugt_unpack_bits(bytes, SMAUGT_N, bits, shift)

function smaugt_pack_vector(vector, bits, shift = 0)
    polynomial_bytes = SMAUGT_N * bits ÷ 8
    output = zeros(UInt8, length(vector) * polynomial_bytes)
    for index in 0:(length(vector) - 1)
        packed = smaugt_pack_polynomial(vector[index + 1], bits, shift)
        copyto!(output, index * polynomial_bytes + 1, packed, 1, polynomial_bytes)
    end
    return output
end

function smaugt_unpack_vector(bytes, count, bits, shift = 0)
    polynomial_bytes = SMAUGT_N * bits ÷ 8
    return [smaugt_unpack_polynomial(
        @view(bytes[(index * polynomial_bytes + 1):((index + 1) * polynomial_bytes)]),
        bits,
        shift,
    ) for index in 0:(count - 1)]
end

function smaugt_pack_secret(vector)
    output = zeros(UInt8, length(vector) * SMAUGT_N ÷ 4)
    for polynomial in 0:(length(vector) - 1)
        coefficients = vector[polynomial + 1]
        offset = polynomial * SMAUGT_N ÷ 4
        for index in 0:(SMAUGT_N ÷ 4 - 1)
            coefficient_index = 4index
            packed = ((1 - smaugt_get(coefficients, coefficient_index)) & 3) |
                (((1 - smaugt_get(coefficients, coefficient_index + 1)) & 3) << 2) |
                (((1 - smaugt_get(coefficients, coefficient_index + 2)) & 3) << 4) |
                (((1 - smaugt_get(coefficients, coefficient_index + 3)) & 3) << 6)
            output[offset + index + 1] = UInt8(packed)
        end
    end
    return output
end

function smaugt_unpack_secret(bytes, count)
    return [begin
        output = smaugt_zero_polynomial()
        offset = polynomial * SMAUGT_N ÷ 4
        for index in 0:(SMAUGT_N ÷ 4 - 1)
            packed = Int(bytes[offset + index + 1])
            for lane in 0:3
                smaugt_set!(output, 4index + lane, 1 - ((packed >> (2lane)) & 3))
            end
        end
        output
    end for polynomial in 0:(count - 1)]
end

function smaugt_multiply_accumulate!(accumulated, left, right, shift)
    for left_index in 0:(SMAUGT_N - 1)
        left_value = smaugt_get(left, left_index) >> shift
        for right_index in 0:(SMAUGT_N - 1)
            right_value = smaugt_get(right, right_index)
            right_value == 0 && continue
            degree = left_index + right_index
            if degree >= SMAUGT_N
                accumulated[degree - SMAUGT_N + 1] -= left_value * right_value
            else
                accumulated[degree + 1] += left_value * right_value
            end
        end
    end
    return accumulated
end

function smaugt_vector_product_add!(result, left, right, input_shift, output_shift = input_shift)
    product = zeros(Int64, SMAUGT_N)
    for index in eachindex(left)
        smaugt_multiply_accumulate!(product, left[index], right[index], input_shift)
    end
    for index in 0:(SMAUGT_N - 1)
        smaugt_set!(result, index, smaugt_get(result, index) + (product[index + 1] << output_shift))
    end
    return result
end

function smaugt_transpose_matrix_vector_product(parameters, matrix, vector)
    shift = 16 - parameters.log_q
    return [begin
        product = zeros(Int64, SMAUGT_N)
        for row in 1:parameters.k
            smaugt_multiply_accumulate!(product, matrix[row][column], vector[row], 0)
        end
        output = smaugt_zero_polynomial()
        for index in 0:(SMAUGT_N - 1)
            smaugt_set!(output, index, product[index + 1] << shift)
        end
        output
    end for column in 1:parameters.k]
end

function smaugt_matrix_vector_product_subtract(parameters, initial, matrix, vector)
    shift = 16 - parameters.log_q
    return [begin
        product = zeros(Int64, SMAUGT_N)
        for column in 1:parameters.k
            smaugt_multiply_accumulate!(product, matrix[row][column], vector[column], 0)
        end
        output = smaugt_zero_polynomial()
        for index in 0:(SMAUGT_N - 1)
            smaugt_set!(output, index,
                smaugt_get(initial[row], index) - (product[index + 1] << shift))
        end
        output
    end for row in 1:parameters.k]
end

function smaugt_sample_fixed_weight(parameters, seed)
    xof = Shake(256)
    update!(xof, seed)
    finalize!(xof)
    random_bytes_value = squeeze!(xof, Vector{UInt8}(undef, 616))
    random = Vector{UInt16}(undef, 308)
    for index in 0:307
        random[index + 1] = UInt16(random_bytes_value[2index + 1]) |
                            (UInt16(random_bytes_value[2index + 2]) << 8)
    end

    indices = zeros(Int16, SMAUGT_N)
    extra = SMAUGT_N
    for index in 0:(SMAUGT_N - 1)
        modulus = SMAUGT_N - index
        threshold = 65_536 % modulus
        product = Int(random[index + 1]) * modulus
        low = product & 0xffff
        while low < threshold
            extra < length(random) || return nothing
            product = Int(random[extra + 1]) * modulus
            extra += 1
            low = product & 0xffff
        end
        indices[index + 1] = Int16(product >> 16)
    end

    signs = squeeze!(xof, Vector{UInt8}(undef, SMAUGT_N ÷ 4))
    output = smaugt_zero_polynomial()
    cutoff = SMAUGT_N - parameters.secret_weight
    for index in 0:(SMAUGT_N - 1)
        selection = (Int(indices[index + 1]) - cutoff) >> 15
        cutoff += selection
        present = 1 + selection
        sign_index = (((index >> 4) >> 3) << 4) + (index & 15)
        sign = (Int(signs[sign_index + 1]) >> ((index >> 4) & 7)) & 1
        smaugt_set!(output, index, -present & ((sign << 1) - 1))
    end
    return output
end

function smaugt_expand_secret(parameters, seed)
    return [begin
        attempt = 0
        polynomial = nothing
        while polynomial === nothing
            polynomial = smaugt_sample_fixed_weight(
                parameters,
                concatenate_bytes(seed, UInt8[index * parameters.k, attempt]),
            )
            attempt += 1
        end
        polynomial
    end for index in 0:(parameters.k - 1)]
end

function smaugt_sample_cbd(parameters, bytes)
    output = smaugt_zero_polynomial()
    if parameters.mode in (1, :timer)
        for block in 0:(SMAUGT_N ÷ 8 - 1)
            packed = Int(bytes[3block + 1]) | (Int(bytes[3block + 2]) << 8) |
                     (Int(bytes[3block + 3]) << 16)
            magnitude = packed & 0x249249
            magnitude &= (packed >> 1) & 0x249249
            signs = (packed >> 2) & 0x249249
            for lane in 0:7
                present = (magnitude >> (3lane)) & 1
                smaugt_set!(output, 8block + lane,
                    present * (1 - (((signs >> (3lane)) & 1) << 1)))
            end
        end
    elseif parameters.mode == 3
        for index in 0:(SMAUGT_N - 1)
            byte = Int(bytes[(index >> 2) + 1])
            smaugt_set!(output, index,
                ((byte >> (2(index & 3))) & 1) -
                ((byte >> (2(index & 3) + 1)) & 1))
        end
    else
        for block in 0:(SMAUGT_N ÷ 8 - 1)
            packed = UInt32(bytes[4block + 1]) |
                     (UInt32(bytes[4block + 2]) << 8) |
                     (UInt32(bytes[4block + 3]) << 16) |
                     (UInt32(bytes[4block + 4]) << 24)
            magnitude = packed & 0x11111111
            magnitude |= (packed >> 1) & 0x11111111
            magnitude &= (packed >> 2) & 0x11111111
            signs = (packed >> 3) & 0x11111111
            for lane in 0:7
                present = Int((magnitude >> (4lane)) & 1)
                sign = Int((signs >> (4lane)) & 1)
                smaugt_set!(output, 8block + lane, present * (1 - (sign << 1)))
            end
        end
    end
    return output
end

function smaugt_expand_ephemeral(parameters, seed)
    return [smaugt_sample_cbd(
        parameters,
        shake(256, parameters.cbd_seed_bytes,
            concatenate_bytes(seed, UInt8[index])),
    ) for index in 0:(parameters.k - 1)]
end

function smaugt_load_transposed_words(bytes)
    output = zeros(UInt32, 80)
    position = 0
    for group in 0:3, word in 0:9
        target = 20group + 2word
        output[target + 1] = UInt32(bytes[position + word + 1]) |
            (UInt32(bytes[position + 10 + word + 1]) << 8) |
            (UInt32(bytes[position + 20 + word + 1]) << 16) |
            (UInt32(bytes[position + 30 + word + 1]) << 24)
        output[target + 2] = UInt32(bytes[position + 40 + word + 1]) |
            (UInt32(bytes[position + 50 + word + 1]) << 8) |
            (UInt32(bytes[position + 60 + word + 1]) << 16) |
            (UInt32(bytes[position + 70 + word + 1]) << 24)
        word == 9 && (position += 80)
    end
    return output
end

function smaugt_sample_gaussian(parameters, seed)
    shift = 16 - parameters.log_q
    return [begin
        bytes = shake(256, 320,
            concatenate_bytes(seed, UInt8[parameters.k * vector_index]))
        words = smaugt_load_transposed_words(bytes)
        output = smaugt_zero_polynomial()
        for block in 0:3, half in 0:1
            offset = 20block
            x0, x1, x2, x3, x4 = (words[offset + 2i + half + 1] for i in 0:4)
            x5, x6, x7, x8 = (words[offset + 2i + half + 1] for i in 5:8)
            signs = words[offset + 18 + half + 1]
            first = (x0 & x1 & x2 & x3 & x4 & x5 & x7 & ~x8) |
                (x0 & x3 & x4 & x5 & x6 & x8) |
                (x1 & x3 & x4 & x5 & x6 & x8) |
                (x2 & x3 & x4 & x5 & x6 & x8) |
                (~x2 & ~x3 & ~x6 & x8) |
                (~x1 & ~x3 & ~x6 & x8) |
                (x6 & x7 & ~x8) | (~x5 & ~x6 & x8) |
                (~x4 & ~x6 & x8) | (~x7 & x8)
            second = (x1 & x2 & x4 & x5 & x7 & x8) |
                     (x3 & x4 & x5 & x7 & x8) | (x6 & x7 & x8)
            for bit in 0:31
                value = Int((first >> bit) & 1) | (Int((second >> bit) & 1) << 1)
                sign = Int((signs >> bit) & 1)
                value = xor(-sign, value) + sign
                smaugt_set!(output, 64block + 32half + bit, value << shift)
            end
        end
        output
    end for vector_index in 0:(parameters.k - 1)]
end

function smaugt_expand_matrix(parameters, seed)
    polynomial_bytes = parameters.log_q * SMAUGT_N ÷ 8
    return [[smaugt_unpack_polynomial(
        shake(128, polynomial_bytes,
            concatenate_bytes(seed, UInt8[row - 1, column - 1])),
        parameters.log_q,
    ) for column in 1:parameters.k] for row in 1:parameters.k]
end

function smaugt_unpack_public_key(parameters, bytes)
    seed = copy(bytes[1:32])
    return (
        seed = seed,
        matrix = smaugt_expand_matrix(parameters, seed),
        vector = smaugt_unpack_vector(@view(bytes[33:end]), parameters.k,
            parameters.log_q, 16 - parameters.log_q),
    )
end

function smaugt_encode_message(parameters, message)
    output = smaugt_zero_polynomial()
    if parameters.mode == :timer
        for index in 0:127
            value = ((Int(message[(index >> 3) + 1]) >> (index & 7)) & 1) * 32_767
            smaugt_set!(output, index, value)
            smaugt_set!(output, index + 128, value)
        end
    else
        for index in 0:(SMAUGT_N - 1)
            value = ((Int(message[(index >> 3) + 1]) >> (index & 7)) & 1) << 15
            smaugt_set!(output, index, value)
        end
    end
    return output
end

function smaugt_flip_absolute(value)
    difference = wrapped_int16(value - 32_767)
    mask = difference >> 15
    return xor(difference + mask, mask) & 0xffff
end

function smaugt_decode_message(parameters, polynomial)
    output = zeros(UInt8, parameters.message_bytes)
    if parameters.mode == :timer
        for index in 0:127
            sum_value = (smaugt_flip_absolute(smaugt_get(polynomial, index)) +
                         smaugt_flip_absolute(smaugt_get(polynomial, index + 128))) & 0xffff
            bit = ((sum_value - 32_767) & 0xffff) >> 15
            output[(index >> 3) + 1] |= UInt8(bit << (index & 7))
        end
    else
        for index in 0:(SMAUGT_N - 1)
            rounded = wrapped_int16(smaugt_get(polynomial, index) + 0x4000)
            output[(index >> 3) + 1] |= UInt8(((rounded >> 15) & 1) << (index & 7))
        end
    end
    return output
end

function smaugt_generate_pke_keypair(parameters, seed)
    expanded = sha3(512, seed)
    secret = smaugt_expand_secret(parameters, @view(expanded[1:32]))
    matrix_seed = copy(expanded[33:end])
    matrix = smaugt_expand_matrix(parameters, matrix_seed)
    errors = smaugt_sample_gaussian(parameters, @view(expanded[1:32]))
    public_vector = smaugt_matrix_vector_product_subtract(parameters, errors, matrix, secret)
    return (
        public_key = concatenate_bytes(matrix_seed,
            smaugt_pack_vector(public_vector, parameters.log_q, 16 - parameters.log_q)),
        secret_key = smaugt_pack_secret(secret),
    )
end

function smaugt_encrypt_pke(parameters, public_key_bytes, message, seed)
    public_key = smaugt_unpack_public_key(parameters, public_key_bytes)
    ephemeral = smaugt_expand_ephemeral(parameters, seed)
    first = smaugt_transpose_matrix_vector_product(parameters, public_key.matrix, ephemeral)
    shift_p = 16 - parameters.log_p
    round_add = 1 << (shift_p - 1)
    round_mask = (0x1_0000 - (1 << shift_p)) & 0xffff
    for polynomial in first, index in 0:(SMAUGT_N - 1)
        value = ((smaugt_get(polynomial, index) + round_add) & round_mask) >> shift_p
        smaugt_set!(polynomial, index, value)
    end

    second = smaugt_encode_message(parameters, message)
    smaugt_vector_product_add!(second, public_key.vector, ephemeral, 16 - parameters.log_q)
    shift_p_prime = 16 - parameters.log_p_prime
    second_add = 1 << (shift_p_prime - 1)
    second_mask = (0x1_0000 - (1 << shift_p_prime)) & 0xffff
    for index in 0:(SMAUGT_N - 1)
        value = ((smaugt_get(second, index) + second_add) & second_mask) >> shift_p_prime
        smaugt_set!(second, index, value)
    end
    return concatenate_bytes(
        smaugt_pack_vector(first, parameters.log_p),
        smaugt_pack_polynomial(second, parameters.log_p_prime),
    )
end

function smaugt_decrypt_pke(parameters, secret_key, ciphertext)
    first_bytes = parameters.k * parameters.log_p * SMAUGT_N ÷ 8
    first = smaugt_unpack_vector(@view(ciphertext[1:first_bytes]),
        parameters.k, parameters.log_p)
    second = smaugt_unpack_polynomial(@view(ciphertext[(first_bytes + 1):end]),
        parameters.log_p_prime)
    shift_p = 16 - parameters.log_p
    shift_p_prime = 16 - parameters.log_p_prime
    for index in 0:(SMAUGT_N - 1)
        smaugt_set!(second, index, smaugt_get(second, index) << shift_p_prime)
    end
    smaugt_vector_product_add!(second, first,
        smaugt_unpack_secret(secret_key, parameters.k), 0, shift_p)
    return smaugt_decode_message(parameters, second)
end

smaugt_derive(length_value, first, second) = shake(256, length_value, first, second)

function smaugt_generate_keypair(parameters; entropy_source=random_bytes)
    fallback_seed = entropy_source(32)
    key_seed = entropy_source(32)
    pke = smaugt_generate_pke_keypair(parameters, key_seed)
    return (
        public_key = pke.public_key,
        secret_key = concatenate_bytes(pke.secret_key, fallback_seed, pke.public_key),
    )
end

function smaugt_encapsulate(parameters, public_key; entropy_source=random_bytes)
    message = entropy_source(parameters.message_bytes)
    expanded = smaugt_derive(64, message, sha3(256, public_key))
    return (
        ciphertext = smaugt_encrypt_pke(parameters, public_key, message, @view(expanded[1:32])),
        shared_secret = copy(expanded[33:64]),
    )
end

function smaugt_decapsulate(parameters, ciphertext, secret_key)
    pke_end = parameters.pke_secret_key_bytes
    pke_secret = @view(secret_key[1:pke_end])
    fallback_seed = @view(secret_key[(pke_end + 1):(pke_end + 32)])
    public_key = @view(secret_key[(pke_end + 33):end])
    message = smaugt_decrypt_pke(parameters, pke_secret, ciphertext)
    expanded = smaugt_derive(64, message, sha3(256, public_key))
    expected = smaugt_encrypt_pke(parameters, public_key, message, @view(expanded[1:32]))
    if ciphertext == expected
        return copy(expanded[33:64])
    end
    fallback = smaugt_derive(64, fallback_seed, sha3(256, ciphertext))
    return copy(fallback[33:64])
end
