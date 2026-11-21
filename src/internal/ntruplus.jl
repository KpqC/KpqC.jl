const NTRUPLUS_Q = 3457
const NTRUPLUS_Q_INVERSE = 12929
const NTRUPLUS_MONTGOMERY_R = -147
const NTRUPLUS_MONTGOMERY_R_INVERSE = -682
const NTRUPLUS_MONTGOMERY_R_SQUARED = 867
const NTRUPLUS_OMEGA = -886
const NTRUPLUS_Z_MINUS_Z5_INVERSE = -1665
const NTRUPLUS_BARRETT_MULTIPLIER = ((1 << 26) + NTRUPLUS_Q ÷ 2) ÷ NTRUPLUS_Q

@inline coefficient(value::AbstractVector, index::Int) = Int(value[index + 1])
@inline function set_coefficient!(value::AbstractVector{Int16}, index::Int, item::Integer)
    value[index + 1] = reinterpret(Int16, UInt16(item & 0xffff))
    return item
end

@inline function wrapped_int16(value::Integer)
    return Int(reinterpret(Int16, UInt16(value & 0xffff)))
end

function concatenate_bytes(values::AbstractVector{UInt8}...)
    output = UInt8[]
    sizehint!(output, sum(length, values))
    foreach(value -> append!(output, value), values)
    return output
end

function ntruplus_montgomery_reduce(value::Integer)
    input = Int(value)
    factor = wrapped_int16(input * NTRUPLUS_Q_INVERSE)
    return (input - factor * NTRUPLUS_Q) >> 16
end

function ntruplus_barrett_reduce(value::Integer)
    input = wrapped_int16(value)
    quotient = (NTRUPLUS_BARRETT_MULTIPLIER * input + (1 << 25)) >> 26
    return input - quotient * NTRUPLUS_Q
end

@inline ntruplus_multiply_field(left::Integer, right::Integer) =
    ntruplus_montgomery_reduce(left * right)

function ntruplus_invert_field(value::Integer)
    t0 = ntruplus_multiply_field(value, value)
    t0 = ntruplus_multiply_field(t0, t0)
    t0 = ntruplus_multiply_field(t0, t0)
    t0 = ntruplus_multiply_field(t0, t0)
    t1 = ntruplus_multiply_field(t0, value)
    t0 = ntruplus_multiply_field(t0, t0)
    t0 = ntruplus_multiply_field(t0, t0)
    t0 = ntruplus_multiply_field(t0, t0)
    t1 = ntruplus_multiply_field(t0, t1)
    t0 = ntruplus_multiply_field(t1, t0)
    t1 = ntruplus_multiply_field(t0, t1)
    t1 = ntruplus_multiply_field(t1, t0)
    t0 = ntruplus_multiply_field(t1, t1)
    t0 = ntruplus_multiply_field(t0, t0)
    t0 = ntruplus_multiply_field(t0, t1)
    return ntruplus_multiply_field(NTRUPLUS_MONTGOMERY_R_INVERSE, t0)
end

function ntruplus_ntt!(parameters::NtruPlusParameters, result::Vector{Int16})
    (; binary_start, degree, mode, n, ternary_end, zetas) = parameters
    twiddle_index = 1
    zeta_one = coefficient(zetas, twiddle_index)
    twiddle_index += 1

    for index in 0:(n ÷ 2 - 1)
        product = ntruplus_multiply_field(zeta_one, coefficient(result, index + n ÷ 2))
        base = coefficient(result, index)
        set_coefficient!(result, index + n ÷ 2,
            base + coefficient(result, index + n ÷ 2) - product)
        set_coefficient!(result, index, base + product)
    end

    step = n ÷ 6
    while step >= ternary_end
        for start in 0:(3step):(n - 1)
            zeta_one = coefficient(zetas, twiddle_index)
            twiddle_index += 1
            zeta_two = coefficient(zetas, twiddle_index)
            twiddle_index += 1
            for index in start:(start + step - 1)
                base = coefficient(result, index)
                first = ntruplus_multiply_field(zeta_one, coefficient(result, index + step))
                second = ntruplus_multiply_field(zeta_two, coefficient(result, index + 2step))
                third = ntruplus_multiply_field(NTRUPLUS_OMEGA, first - second)
                set_coefficient!(result, index + 2step, base - first - third)
                set_coefficient!(result, index + step, base - second + third)
                set_coefficient!(result, index, base + first + second)
            end
        end
        step ÷= 3
    end

    step = binary_start
    while step >= degree
        for start in 0:(2step):(n - 1)
            zeta_one = coefficient(zetas, twiddle_index)
            twiddle_index += 1
            for index in start:(start + step - 1)
                base = coefficient(result, index)
                product = ntruplus_multiply_field(zeta_one, coefficient(result, index + step))
                upper = base - product
                lower = base + product
                if mode != 768
                    upper = ntruplus_barrett_reduce(upper)
                    lower = ntruplus_barrett_reduce(lower)
                end
                set_coefficient!(result, index + step, upper)
                set_coefficient!(result, index, lower)
            end
        end
        step >>= 1
    end
    if mode == 768
        for index in 0:(n - 1)
            set_coefficient!(result, index, ntruplus_barrett_reduce(coefficient(result, index)))
        end
    end
    return result
end

ntruplus_ntt(parameters::NtruPlusParameters, value::Vector{Int16}) =
    ntruplus_ntt!(parameters, copy(value))

function ntruplus_inverse_ntt!(parameters::NtruPlusParameters, result::Vector{Int16})
    (; binary_start, degree, mode, n, ternary_end, zetas) = parameters
    twiddle_index = length(zetas) - 1

    step = degree
    while step <= binary_start
        for start in 0:(2step):(n - 1)
            zeta = coefficient(zetas, twiddle_index)
            twiddle_index -= 1
            for index in start:(start + step - 1)
                upper = coefficient(result, index + step)
                set_coefficient!(result, index + step,
                    ntruplus_multiply_field(zeta, upper - coefficient(result, index)))
                set_coefficient!(result, index,
                    ntruplus_barrett_reduce(coefficient(result, index) + upper))
            end
        end
        step <<= 1
    end

    step = ternary_end
    while step <= n ÷ 6
        for start in 0:(3step):(n - 1)
            zeta_two = coefficient(zetas, twiddle_index)
            twiddle_index -= 1
            zeta_one = coefficient(zetas, twiddle_index)
            twiddle_index -= 1
            for index in start:(start + step - 1)
                base = coefficient(result, index)
                middle = coefficient(result, index + step)
                upper = coefficient(result, index + 2step)
                first = ntruplus_multiply_field(NTRUPLUS_OMEGA, middle - base)
                set_coefficient!(result, index + step,
                    ntruplus_multiply_field(zeta_one, upper - base + first))
                set_coefficient!(result, index + 2step,
                    ntruplus_multiply_field(zeta_two, upper - middle - first))
                sum = base + middle + upper
                set_coefficient!(result, index,
                    mode == 768 ? sum : ntruplus_barrett_reduce(sum))
            end
        end
        step *= 3
    end

    for index in 0:(n ÷ 2 - 1)
        lower = coefficient(result, index)
        upper = coefficient(result, index + n ÷ 2)
        sum = lower + upper
        difference = ntruplus_multiply_field(
            NTRUPLUS_Z_MINUS_Z5_INVERSE,
            lower - upper,
        )
        set_coefficient!(result, index,
            ntruplus_multiply_field(parameters.inverse_n, sum - difference))
        set_coefficient!(result, index + n ÷ 2,
            ntruplus_multiply_field(parameters.inverse_two_n, difference))
    end
    return result
end

ntruplus_inverse_ntt(parameters::NtruPlusParameters, value::Vector{Int16}) =
    ntruplus_inverse_ntt!(parameters, copy(value))

function ntruplus_prepare_base3!(numerators, input, offset, zeta)
    first = ntruplus_montgomery_reduce(coefficient(input, offset + 1) * coefficient(input, offset + 2))
    second = ntruplus_montgomery_reduce(coefficient(input, offset + 2)^2)
    third = ntruplus_montgomery_reduce(
        coefficient(input, offset + 1)^2 - coefficient(input, offset) * coefficient(input, offset + 2),
    )
    first = ntruplus_montgomery_reduce(coefficient(input, offset)^2 - first * zeta)
    second = ntruplus_montgomery_reduce(second * zeta - coefficient(input, offset) * coefficient(input, offset + 1))
    denominator = ntruplus_montgomery_reduce(third * coefficient(input, offset + 1) + second * coefficient(input, offset + 2))
    denominator = ntruplus_montgomery_reduce(denominator * zeta + first * coefficient(input, offset))
    set_coefficient!(numerators, offset, first)
    set_coefficient!(numerators, offset + 1, second)
    set_coefficient!(numerators, offset + 2, third)
    return denominator
end

function ntruplus_prepare_base4!(numerators, input, offset, zeta)
    first = ntruplus_montgomery_reduce(
        coefficient(input, offset + 2)^2 -
        2coefficient(input, offset + 1) * coefficient(input, offset + 3),
    )
    second = ntruplus_montgomery_reduce(coefficient(input, offset + 3)^2)
    first = ntruplus_montgomery_reduce(coefficient(input, offset)^2 + first * zeta)
    second = ntruplus_montgomery_reduce(
        coefficient(input, offset + 1)^2 + second * zeta -
        2coefficient(input, offset) * coefficient(input, offset + 2),
    )
    third = ntruplus_montgomery_reduce(second * zeta)
    denominator = ntruplus_montgomery_reduce(first^2 - second * third)
    outputs = (
        ntruplus_montgomery_reduce(coefficient(input, offset) * first + coefficient(input, offset + 2) * third),
        -ntruplus_montgomery_reduce(coefficient(input, offset + 3) * third + coefficient(input, offset + 1) * first),
        ntruplus_montgomery_reduce(coefficient(input, offset + 2) * first + coefficient(input, offset) * second),
        -ntruplus_montgomery_reduce(coefficient(input, offset + 1) * second + coefficient(input, offset + 3) * first),
    )
    for index in 0:3
        set_coefficient!(numerators, offset + index, outputs[index + 1])
    end
    return denominator
end

function ntruplus_invert_batch(values::Vector{Int16})
    prefixes = zeros(Int16, length(values))
    product = NTRUPLUS_MONTGOMERY_R
    for index in 0:(length(values) - 1)
        value = coefficient(values, index)
        value == 0 && return nothing
        set_coefficient!(prefixes, index, product)
        product = ntruplus_multiply_field(product, value)
    end
    inverse = ntruplus_invert_field(product)
    output = zeros(Int16, length(values))
    for index in (length(values) - 1):-1:0
        set_coefficient!(output, index,
            ntruplus_multiply_field(inverse, coefficient(prefixes, index)))
        inverse = ntruplus_multiply_field(inverse, coefficient(values, index))
    end
    return output
end

function ntruplus_invert_polynomial(parameters::NtruPlusParameters, input::Vector{Int16})
    numerators = zeros(Int16, parameters.n)
    degree = parameters.degree
    segments = parameters.n ÷ degree
    denominators = zeros(Int16, segments)
    zeta_offset = length(parameters.zetas) ÷ 2
    for segment in 0:(segments - 1)
        offset = segment * degree
        base_zeta = coefficient(parameters.zetas, zeta_offset + (segment >> 1))
        zeta = iseven(segment) ? base_zeta : -base_zeta
        denominator = degree == 3 ?
            ntruplus_prepare_base3!(numerators, input, offset, zeta) :
            ntruplus_prepare_base4!(numerators, input, offset, zeta)
        set_coefficient!(denominators, segment, denominator)
    end
    inverses = ntruplus_invert_batch(denominators)
    inverses === nothing && return nothing
    output = zeros(Int16, parameters.n)
    for segment in 0:(segments - 1)
        inverse = coefficient(inverses, segment)
        offset = segment * degree
        for index in 0:(degree - 1)
            set_coefficient!(output, offset + index, ntruplus_montgomery_reduce(
                coefficient(numerators, offset + index) * inverse,
            ))
        end
    end
    return output
end

function ntruplus_multiply_base3!(output, left, right, offset, zeta, add)
    first = ntruplus_montgomery_reduce(
        coefficient(left, offset + 2) * coefficient(right, offset + 1) +
        coefficient(left, offset + 1) * coefficient(right, offset + 2),
    )
    second = ntruplus_montgomery_reduce(coefficient(left, offset + 2) * coefficient(right, offset + 2))
    first = ntruplus_montgomery_reduce(first * zeta + coefficient(left, offset) * coefficient(right, offset))
    second = ntruplus_montgomery_reduce(
        second * zeta + coefficient(left, offset) * coefficient(right, offset + 1) +
        coefficient(left, offset + 1) * coefficient(right, offset),
    )
    third = ntruplus_montgomery_reduce(
        coefficient(left, offset + 2) * coefficient(right, offset) +
        coefficient(left, offset + 1) * coefficient(right, offset + 1) +
        coefficient(left, offset) * coefficient(right, offset + 2),
    )
    for (index, value) in enumerate((first, second, third))
        result = value * NTRUPLUS_MONTGOMERY_R_SQUARED
        if add !== nothing
            result += coefficient(add, offset + index - 1) * NTRUPLUS_MONTGOMERY_R
        end
        set_coefficient!(output, offset + index - 1, ntruplus_montgomery_reduce(result))
    end
end

function ntruplus_multiply_base4!(output, left, right, offset, zeta, add)
    t0 = ntruplus_montgomery_reduce(
        coefficient(left, offset + 1) * coefficient(right, offset + 3) +
        coefficient(left, offset + 2) * coefficient(right, offset + 2) +
        coefficient(left, offset + 3) * coefficient(right, offset + 1),
    )
    t1 = ntruplus_montgomery_reduce(
        coefficient(left, offset + 2) * coefficient(right, offset + 3) +
        coefficient(left, offset + 3) * coefficient(right, offset + 2),
    )
    t2 = ntruplus_montgomery_reduce(coefficient(left, offset + 3) * coefficient(right, offset + 3))
    values = (
        ntruplus_montgomery_reduce(t0 * zeta + coefficient(left, offset) * coefficient(right, offset)),
        ntruplus_montgomery_reduce(t1 * zeta + coefficient(left, offset) * coefficient(right, offset + 1) + coefficient(left, offset + 1) * coefficient(right, offset)),
        ntruplus_montgomery_reduce(t2 * zeta + coefficient(left, offset) * coefficient(right, offset + 2) + coefficient(left, offset + 1) * coefficient(right, offset + 1) + coefficient(left, offset + 2) * coefficient(right, offset)),
        ntruplus_montgomery_reduce(coefficient(left, offset) * coefficient(right, offset + 3) + coefficient(left, offset + 1) * coefficient(right, offset + 2) + coefficient(left, offset + 2) * coefficient(right, offset + 1) + coefficient(left, offset + 3) * coefficient(right, offset)),
    )
    for (index, value) in enumerate(values)
        result = value * NTRUPLUS_MONTGOMERY_R_SQUARED
        if add !== nothing
            result += coefficient(add, offset + index - 1) * NTRUPLUS_MONTGOMERY_R
        end
        set_coefficient!(output, offset + index - 1, ntruplus_montgomery_reduce(result))
    end
end

function ntruplus_multiply_polynomials(
    parameters::NtruPlusParameters,
    left::Vector{Int16},
    right::Vector{Int16},
    add::Union{Nothing,Vector{Int16}} = nothing,
)
    output = zeros(Int16, parameters.n)
    degree = parameters.degree
    block = 2degree
    zeta_offset = length(parameters.zetas) ÷ 2
    index = 0
    for offset in 0:block:(parameters.n - 1)
        zeta = coefficient(parameters.zetas, zeta_offset + index)
        multiply = degree == 3 ? ntruplus_multiply_base3! : ntruplus_multiply_base4!
        multiply(output, left, right, offset, zeta, add)
        multiply(output, left, right, offset + degree, -zeta, add)
        index += 1
    end
    return output
end

const NTRUPLUS_HASH_F_DOMAIN = UInt8[0]
const NTRUPLUS_HASH_G_DOMAIN = UInt8[1]
const NTRUPLUS_HASH_H_DOMAIN = UInt8[2]

ntruplus_hash_f(message) = shake(256, 32, NTRUPLUS_HASH_F_DOMAIN, message)
ntruplus_hash_g(parameters, message) =
    shake(256, parameters.n ÷ 4, NTRUPLUS_HASH_G_DOMAIN, message)
ntruplus_hash_h(parameters, messages...) =
    shake(256, 32 + parameters.n ÷ 4, NTRUPLUS_HASH_H_DOMAIN, messages...)

function ntruplus_sample_centered_binomial(parameters, bytes)
    output = zeros(Int16, parameters.n)
    for index in 0:(parameters.n ÷ 8 - 1)
        first = Int(bytes[index + 1])
        second = Int(bytes[index + parameters.n ÷ 8 + 1])
        for bit in 0:7
            set_coefficient!(output, 8index + bit, (first & 1) - (second & 1))
            first >>= 1
            second >>= 1
        end
    end
    return output
end

function ntruplus_serialize(parameters, polynomial)
    output = zeros(UInt8, parameters.polynomial_bytes)
    for index in 0:(parameters.n ÷ 2 - 1)
        first = coefficient(polynomial, 2index)
        second = coefficient(polynomial, 2index + 1)
        first < 0 && (first += NTRUPLUS_Q)
        second < 0 && (second += NTRUPLUS_Q)
        output[3index + 1] = UInt8(first & 0xff)
        output[3index + 2] = UInt8(((first >> 8) | (second << 4)) & 0xff)
        output[3index + 3] = UInt8((second >> 4) & 0xff)
    end
    return output
end

function ntruplus_equal_encodings(left, right)
    difference = 0
    for index in 0:(length(left) - 1)
        left_coefficient = coefficient(left, index)
        right_coefficient = coefficient(right, index)
        left_coefficient < 0 && (left_coefficient += NTRUPLUS_Q)
        right_coefficient < 0 && (right_coefficient += NTRUPLUS_Q)
        difference |= xor(left_coefficient, right_coefficient)
    end
    return difference == 0
end

function ntruplus_deserialize(parameters, bytes)
    polynomial = zeros(Int16, parameters.n)
    canonical = true
    for index in 0:(parameters.n ÷ 2 - 1)
        first = (Int(bytes[3index + 1]) | (Int(bytes[3index + 2]) << 8)) & 0xfff
        second = ((Int(bytes[3index + 2]) >> 4) | (Int(bytes[3index + 3]) << 4)) & 0xfff
        set_coefficient!(polynomial, 2index, first)
        set_coefficient!(polynomial, 2index + 1, second)
        canonical &= first < NTRUPLUS_Q && second < NTRUPLUS_Q
    end
    return (; polynomial, canonical)
end

function ntruplus_encode_sotp(parameters, message, pad)
    bytes = copy(pad)
    for index in 0:(parameters.n ÷ 8 - 1)
        bytes[index + 1] = xor(bytes[index + 1], message[index + 1])
    end
    return ntruplus_sample_centered_binomial(parameters, bytes)
end

function ntruplus_decode_sotp(parameters, polynomial, pad)
    message = zeros(UInt8, parameters.n ÷ 8)
    invalid = 0
    for index in 0:(parameters.n ÷ 8 - 1)
        first = Int(pad[index + 1])
        second = Int(pad[index + parameters.n ÷ 8 + 1])
        byte = 0
        for bit in 0:7
            sum_value = ((second & 1) + coefficient(polynomial, 8index + bit)) & 0xffff
            invalid |= sum_value
            byte = xor(byte, (xor(sum_value, first) & 1) << bit)
            first >>= 1
            second >>= 1
        end
        message[index + 1] = UInt8(byte)
    end
    valid = (invalid >> 1) == 0
    valid || fill!(message, 0)
    return (; message, valid)
end

function ntruplus_centered_remainder_three(value)
    reduced = Int(value)
    reduced < 0 && (reduced += NTRUPLUS_Q)
    reduced -= (NTRUPLUS_Q + 1) ÷ 2
    reduced < 0 && (reduced += NTRUPLUS_Q)
    reduced -= (NTRUPLUS_Q - 1) ÷ 2
    quotient = ((10923 * reduced + (1 << 14)) >> 15) * 3
    return reduced - quotient
end

function ntruplus_centered_mod_three(polynomial)
    output = zeros(Int16, length(polynomial))
    for index in 0:(length(polynomial) - 1)
        set_coefficient!(output, index,
            ntruplus_centered_remainder_three(coefficient(polynomial, index)))
    end
    return output
end

function ntruplus_generate_secret(parameters, coins, constant)
    polynomial = ntruplus_sample_centered_binomial(
        parameters,
        shake(256, parameters.n ÷ 4, coins),
    )
    for index in 0:(length(polynomial) - 1)
        set_coefficient!(polynomial, index, 3coefficient(polynomial, index))
    end
    set_coefficient!(polynomial, 0, coefficient(polynomial, 0) + constant)
    ntruplus_ntt!(parameters, polynomial)
    inverse = ntruplus_invert_polynomial(parameters, polynomial)
    return inverse === nothing ? nothing : (; polynomial, inverse)
end

function ntruplus_generate_keypair(parameters; entropy_source=random_bytes)
    first = nothing
    second = nothing
    while first === nothing
        first = ntruplus_generate_secret(parameters, entropy_source(32), 1)
    end
    while second === nothing
        second = ntruplus_generate_secret(parameters, entropy_source(32), 0)
    end
    public_polynomial = ntruplus_multiply_polynomials(
        parameters,
        second.polynomial,
        first.inverse,
    )
    inverse_public = ntruplus_multiply_polynomials(
        parameters,
        first.polynomial,
        second.inverse,
    )
    public_key = ntruplus_serialize(parameters, public_polynomial)
    secret_key = concatenate_bytes(
        ntruplus_serialize(parameters, first.polynomial),
        ntruplus_serialize(parameters, inverse_public),
        ntruplus_hash_f(public_key),
    )
    return (; public_key, secret_key)
end

function ntruplus_encapsulate(parameters, public_key; entropy_source=random_bytes)
    decoded_public_key = ntruplus_deserialize(parameters, public_key)
    decoded_public_key.canonical || throw(KpqCError(:encapsulate, 1))
    coins = entropy_source(parameters.n ÷ 8)
    message = concatenate_bytes(coins, ntruplus_hash_f(public_key))
    expanded = ntruplus_hash_h(parameters, message)
    random_polynomial = ntruplus_sample_centered_binomial(
        parameters,
        @view(expanded[33:end]),
    )
    ntruplus_ntt!(parameters, random_polynomial)
    pad = ntruplus_hash_g(parameters, ntruplus_serialize(parameters, random_polynomial))
    encoded_message = ntruplus_encode_sotp(parameters, message, pad)
    ntruplus_ntt!(parameters, encoded_message)
    ciphertext_polynomial = ntruplus_multiply_polynomials(
        parameters,
        decoded_public_key.polynomial,
        random_polynomial,
        encoded_message,
    )
    return (
        ciphertext = ntruplus_serialize(parameters, ciphertext_polynomial),
        shared_secret = copy(expanded[1:32]),
    )
end

function ntruplus_decapsulate(parameters, ciphertext, secret_key)
    polynomial_bytes = parameters.polynomial_bytes
    decoded_ciphertext = ntruplus_deserialize(parameters, ciphertext)
    decoded_first = ntruplus_deserialize(parameters, @view(secret_key[1:polynomial_bytes]))
    decoded_inverse_public = ntruplus_deserialize(
        parameters,
        @view(secret_key[(polynomial_bytes + 1):(2polynomial_bytes)]),
    )
    (decoded_ciphertext.canonical && decoded_first.canonical &&
     decoded_inverse_public.canonical) || throw(KpqCError(:decapsulate, 1))

    product = ntruplus_multiply_polynomials(
        parameters,
        decoded_ciphertext.polynomial,
        decoded_first.polynomial,
    )
    ntruplus_inverse_ntt!(parameters, product)
    message_polynomial = ntruplus_centered_mod_three(product)
    transformed_message = ntruplus_ntt(parameters, message_polynomial)
    for index in 0:(length(transformed_message) - 1)
        set_coefficient!(transformed_message, index,
            coefficient(decoded_ciphertext.polynomial, index) -
            coefficient(transformed_message, index))
    end
    recovered_random = ntruplus_multiply_polynomials(
        parameters,
        transformed_message,
        decoded_inverse_public.polynomial,
    )
    pad = ntruplus_hash_g(parameters, ntruplus_serialize(parameters, recovered_random))
    decoded = ntruplus_decode_sotp(parameters, message_polynomial, pad)
    public_hash = @view(secret_key[(2polynomial_bytes + 1):(2polynomial_bytes + 32)])
    expanded = ntruplus_hash_h(parameters, decoded.message, public_hash)
    expected_random = ntruplus_sample_centered_binomial(
        parameters,
        @view(expanded[33:end]),
    )
    ntruplus_ntt!(parameters, expected_random)
    valid = ntruplus_equal_encodings(recovered_random, expected_random) && decoded.valid
    valid || throw(KpqCError(:decapsulate, 1))
    return copy(expanded[1:32])
end
