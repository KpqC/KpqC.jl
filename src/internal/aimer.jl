struct AimerParameters
    mode::Int
    bits::Int
    repetitions::Int
    parties::Int
    log_parties::Int
    message_prefix::UInt8
    signature_bytes::Int
    exponents::Vector{Int}
end

aimer_field_bytes(p) = p.bits ÷ 8
aimer_input_sboxes(p) = length(p.exponents) - 1
aimer_commitment_bytes(p) = 2aimer_field_bytes(p)
aimer_public_key_bytes(p) = 2aimer_field_bytes(p)
aimer_secret_key_bytes(p) = 3aimer_field_bytes(p)
aimer_tape_bytes(p) = (2aimer_input_sboxes(p) + 3) * aimer_field_bytes(p)
aimer_proof_bytes(p) = p.log_parties * aimer_field_bytes(p) +
    aimer_commitment_bytes(p) + (2aimer_input_sboxes(p) + 3) * aimer_field_bytes(p)

const AIMER_PARAMETERS = Dict(
    :aimer128f => AimerParameters(0, 128, 33, 16, 4, 0x00, 6944, [3, 7, 11]),
    :aimer128s => AimerParameters(1, 128, 17, 256, 8, 0x10, 4704, [3, 7, 11]),
    :aimer192f => AimerParameters(2, 192, 49, 16, 4, 0x20, 15408, [11, 23, 47]),
    :aimer192s => AimerParameters(3, 192, 25, 256, 8, 0x30, 10320, [11, 23, 47]),
    :aimer256f => AimerParameters(4, 256, 65, 16, 4, 0x40, 31360, [3, 7, 11, 19]),
    :aimer256s => AimerParameters(5, 256, 33, 256, 8, 0x50, 20224, [3, 7, 11, 19]),
)

struct AimerLinear
    matrices::Vector{Vector{BigInt}}
    vectors::Vector{BigInt}
end

struct AimerCheck
    x::Vector{BigInt}
    z::Vector{BigInt}
    a::Vector{BigInt}
    b::Vector{BigInt}
    c::BigInt
end

mutable struct AimerTape
    plaintext::BigInt
    y::Vector{BigInt}
    a::Vector{BigInt}
    c::BigInt
end

function AimerTape(parameters)
    return AimerTape(
        BigInt(0),
        fill(BigInt(0), aimer_input_sboxes(parameters)),
        fill(BigInt(0), aimer_input_sboxes(parameters) + 1),
        BigInt(0),
    )
end

struct AimerCommitTape
    commitment::Vector{UInt8}
    tape::AimerTape
end

mutable struct AimerHash
    shake::Shake
end

AimerHash(parameters::AimerParameters) =
    AimerHash(Shake(parameters.bits == 128 ? 128 : 256))

function aimer_hash_update!(hash::AimerHash, input, offset::Int = 0,
                            length_value::Int = length(input) - offset)
    length_value == 0 && return hash
    update!(hash.shake, @view(input[(offset + 1):(offset + length_value)]))
    return hash
end

function aimer_hash_squeeze!(hash::AimerHash, output, offset::Int = 0,
                             length_value::Int = length(output) - offset)
    temporary = Vector{UInt8}(undef, length_value)
    squeeze!(hash.shake, temporary)
    copyto!(output, offset + 1, temporary, 1, length_value)
    return output
end

function aimer_prefixed_hash(parameters, prefix)
    hash = AimerHash(parameters)
    aimer_hash_update!(hash, UInt8[prefix])
    return hash
end

function aimer_field_from_bytes(parameters, input, offset = 0)
    output = BigInt(0)
    for index in 0:(aimer_field_bytes(parameters) - 1)
        output |= BigInt(input[offset + index + 1]) << (8index)
    end
    return output
end

function aimer_field_to_bytes!(parameters, value, output, offset = 0)
    for index in 0:(aimer_field_bytes(parameters) - 1)
        output[offset + index + 1] = UInt8((value >> (8index)) & 0xff)
    end
    return output
end

function aimer_field_bytes(parameters, value)
    output = zeros(UInt8, aimer_field_bytes(parameters))
    return aimer_field_to_bytes!(parameters, value, output)
end

const AIMER_REDUCTION_TABLES = Dict{Int,Vector{BigInt}}()
const AIMER_SQUARE_BYTES = let table = Vector{UInt16}(undef, 256)
    for value in 0:255
        expanded = UInt16(0)
        for bit in 0:7
            expanded |= UInt16((value >> bit) & 1) << (2bit)
        end
        table[value + 1] = expanded
    end
    table
end

function aimer_reduction_table(parameters)
    return get!(AIMER_REDUCTION_TABLES, parameters.bits) do
        terms = parameters.bits == 256 ? (0, 2, 5, 10) : (0, 1, 2, 7)
        remainder = reduce(|, (BigInt(1) << term for term in terms))
        table = Vector{BigInt}(undef, 256)
        table[1] = BigInt(0)
        for value in 1:255
            least = trailing_zeros(value)
            table[value + 1] = xor(table[(value & (value - 1)) + 1], remainder << least)
        end
        table
    end
end

function aimer_field_reduce(parameters, value)
    result = BigInt(value)
    table = aimer_reduction_table(parameters)
    for position in (2parameters.bits - 8):-8:parameters.bits
        byte = Int((result >> position) & 0xff)
        iszero(byte) && continue
        result = xor(
            result,
            BigInt(byte) << position,
            table[byte + 1] << (position - parameters.bits),
        )
    end
    return result & ((BigInt(1) << parameters.bits) - 1)
end

function aimer_field_multiply(parameters, left, right)
    multiples = Vector{BigInt}(undef, 16)
    multiples[1] = BigInt(0)
    left_value = BigInt(left)
    for value in 1:15
        least = trailing_zeros(value)
        multiples[value + 1] = xor(
            multiples[(value & (value - 1)) + 1],
            left_value << least,
        )
    end
    product = BigInt(0)
    right_value = BigInt(right)
    for position in (parameters.bits - 4):-4:0
        product <<= 4
        product = xor(product, multiples[Int((right_value >> position) & 0x0f) + 1])
    end
    return aimer_field_reduce(parameters, product)
end

function aimer_field_square(parameters, value)
    product = BigInt(0)
    for byte_index in 0:(parameters.bits ÷ 8 - 1)
        byte = Int((value >> (8byte_index)) & 0xff)
        product |= BigInt(AIMER_SQUARE_BYTES[byte + 1]) << (16byte_index)
    end
    return aimer_field_reduce(parameters, product)
end

function aimer_field_inverse(parameters, value)
    power = BigInt(value)
    output = BigInt(1)
    for _ in 1:(parameters.bits - 1)
        power = aimer_field_square(parameters, power)
        output = aimer_field_multiply(parameters, output, power)
    end
    return output
end

function aimer_matrix_vector(matrix, vector, parameters)
    output = BigInt(0)
    for bit in 0:(parameters.bits - 1)
        !iszero((vector >> bit) & 1) && (output = xor(output, matrix[bit + 1]))
    end
    return output
end

function aimer_squeeze_invertible_matrix(parameters, hash)
    lower = fill(BigInt(0), parameters.bits)
    upper = fill(BigInt(0), parameters.bits)
    bytes = zeros(UInt8, aimer_field_bytes(parameters))
    for row in 0:(parameters.bits - 1)
        aimer_hash_squeeze!(hash, bytes)
        random = aimer_field_from_bytes(parameters, bytes)
        diagonal = BigInt(1) << row
        lower_mask = xor(
            (BigInt(1) << parameters.bits) - 1,
            (BigInt(1) << row) - 1,
        )
        upper_mask = (BigInt(1) << row) - 1
        lower[row + 1] = (random & lower_mask) | diagonal
        upper[row + 1] = (random & upper_mask) | diagonal
    end
    return [aimer_matrix_vector(lower, upper[row + 1], parameters)
            for row in 0:(parameters.bits - 1)]
end

function aimer_generate_linear(parameters, iv)
    hash = AimerHash(parameters)
    aimer_hash_update!(hash, iv)
    matrices = [aimer_squeeze_invertible_matrix(parameters, hash)
                for _ in 1:(2aimer_input_sboxes(parameters))]
    vectors = BigInt[]
    bytes = zeros(UInt8, aimer_field_bytes(parameters))
    for _ in 0:aimer_input_sboxes(parameters)
        aimer_hash_squeeze!(hash, bytes)
        push!(vectors, aimer_field_from_bytes(parameters, bytes))
    end
    return AimerLinear(matrices, vectors)
end

function aimer_encrypt(parameters, plaintext_bytes, iv)
    linear = aimer_generate_linear(parameters, iv)
    plaintext = aimer_field_from_bytes(parameters, plaintext_bytes)
    states = fill(BigInt(0), aimer_input_sboxes(parameters))
    for index in eachindex(states)
        state = xor(
            aimer_matrix_vector(linear.matrices[index], plaintext, parameters),
            linear.vectors[index],
        )
        iszero(state) && return nothing
        inverse = aimer_field_inverse(parameters, state)
        for _ in 1:parameters.exponents[index]
            state = aimer_field_square(parameters, state)
        end
        state = xor(state, 1)
        state = aimer_field_multiply(parameters, state, inverse)
        states[index] = aimer_matrix_vector(
            linear.matrices[index + aimer_input_sboxes(parameters)],
            state,
            parameters,
        )
    end
    state = xor(reduce(xor, states), linear.vectors[aimer_input_sboxes(parameters) + 1])
    iszero(state) && return nothing
    inverse = aimer_field_inverse(parameters, state)
    for _ in 1:parameters.exponents[aimer_input_sboxes(parameters) + 1]
        state = aimer_field_square(parameters, state)
    end
    state = xor(state, 1)
    state = xor(aimer_field_multiply(parameters, state, inverse), plaintext)
    return aimer_field_bytes(parameters, state)
end

function aimer_sbox_outputs(parameters, linear, plaintext, ciphertext)
    output = fill(BigInt(0), aimer_input_sboxes(parameters) + 1)
    for index in 1:aimer_input_sboxes(parameters)
        value = xor(
            aimer_matrix_vector(linear.matrices[index], plaintext, parameters),
            linear.vectors[index],
        )
        inverse = aimer_field_inverse(parameters, value)
        for _ in 1:parameters.exponents[index]
            value = aimer_field_square(parameters, value)
        end
        value = xor(value, 1)
        output[index] = aimer_field_multiply(parameters, value, inverse)
    end
    output[end] = xor(plaintext, ciphertext)
    return output
end

function aimer_tape_from_bytes(parameters, bytes)
    tape = AimerTape(parameters)
    size = aimer_field_bytes(parameters)
    offset = 0
    tape.plaintext = aimer_field_from_bytes(parameters, bytes, offset)
    offset += size
    for index in eachindex(tape.y)
        tape.y[index] = aimer_field_from_bytes(parameters, bytes, offset)
        offset += size
    end
    for index in eachindex(tape.a)
        tape.a[index] = aimer_field_from_bytes(parameters, bytes, offset)
        offset += size
    end
    tape.c = aimer_field_from_bytes(parameters, bytes, offset)
    return tape
end

function aimer_mpc(parameters, linear, tape, ciphertext, party)
    last = party == parameters.parties - 1
    count = aimer_input_sboxes(parameters) + 1
    x = fill(BigInt(0), count)
    z = fill(BigInt(0), count)
    b = fill(BigInt(0), count)
    copyto!(b, 1, tape.y, 1, length(tape.y))
    b[end] = xor(tape.plaintext, last ? ciphertext : BigInt(0))
    for index in 1:aimer_input_sboxes(parameters)
        x[index] = aimer_matrix_vector(linear.matrices[index], tape.plaintext, parameters)
        last && (x[index] = xor(x[index], linear.vectors[index]))
    end
    for index in eachindex(tape.y)
        x[end] = xor(x[end], aimer_matrix_vector(
            linear.matrices[index + aimer_input_sboxes(parameters)],
            tape.y[index],
            parameters,
        ))
    end
    last && (x[end] = xor(x[end], linear.vectors[end]))
    for index in 1:count
        z[index] = aimer_field_square(parameters, x[index])
        for _ in 2:parameters.exponents[index]
            z[index] = aimer_field_square(parameters, z[index])
        end
        last && (z[index] = xor(z[index], 1))
    end
    return AimerCheck(x, z, copy(tape.a), b, tape.c)
end

function aimer_commit_and_expand_tape(parameters, salt, repetition, party, seed)
    hash = aimer_prefixed_hash(parameters, 5)
    aimer_hash_update!(hash, salt)
    aimer_hash_update!(hash, UInt8[repetition, party])
    aimer_hash_update!(hash, seed)
    commitment = zeros(UInt8, aimer_commitment_bytes(parameters))
    aimer_hash_squeeze!(hash, commitment)
    bytes = zeros(UInt8, aimer_tape_bytes(parameters))
    aimer_hash_squeeze!(hash, bytes)
    return AimerCommitTape(commitment, aimer_tape_from_bytes(parameters, bytes))
end

function aimer_expand_tree(parameters, salt, repetition, root)
    nodes = [zeros(UInt8, aimer_field_bytes(parameters))
             for _ in 1:(2parameters.parties - 1)]
    copyto!(nodes[1], root)
    for node in 1:(parameters.parties - 1)
        hash = aimer_prefixed_hash(parameters, 4)
        aimer_hash_update!(hash, salt)
        aimer_hash_update!(hash, UInt8[repetition, node])
        aimer_hash_update!(hash, nodes[node])
        children = zeros(UInt8, 2aimer_field_bytes(parameters))
        aimer_hash_squeeze!(hash, children)
        size = aimer_field_bytes(parameters)
        copyto!(nodes[2node], 1, children, 1, size)
        copyto!(nodes[2node + 1], 1, children, size + 1, size)
    end
    return nodes
end

function aimer_reveal_path(parameters, nodes, covered)
    size = aimer_field_bytes(parameters)
    output = zeros(UInt8, parameters.log_parties * size)
    index = covered + parameters.parties
    for depth in 0:(parameters.log_parties - 1)
        source = nodes[xor(index, 1)]
        copyto!(output, depth * size + 1, source, 1, size)
        index >>= 1
    end
    return output
end

function aimer_reconstruct_tree(parameters, salt, path, repetition, covered)
    size = aimer_field_bytes(parameters)
    nodes = [zeros(UInt8, size) for _ in 1:(2parameters.parties - 2)]
    for depth in 1:(parameters.log_parties - 1)
        sibling = xor((covered + parameters.parties) >> (parameters.log_parties - depth), 1)
        copyto!(nodes[sibling - 1], 1, path,
            (parameters.log_parties - depth) * size + 1, size)
        for index in (1 << depth):((2 << depth) - 1)
            hash = aimer_prefixed_hash(parameters, 4)
            aimer_hash_update!(hash, salt)
            aimer_hash_update!(hash, UInt8[repetition, index])
            aimer_hash_update!(hash, nodes[index - 1])
            children = zeros(UInt8, 2size)
            aimer_hash_squeeze!(hash, children)
            copyto!(nodes[2index - 1], 1, children, 1, size)
            copyto!(nodes[2index], 1, children, size + 1, size)
        end
    end
    sibling = xor(covered + parameters.parties, 1)
    copyto!(nodes[sibling - 1], 1, path, 1, size)
    return nodes
end

aimer_proof_start(p, repetition) = aimer_field_bytes(p) +
    2aimer_commitment_bytes(p) + repetition * aimer_proof_bytes(p)
aimer_missing_commitment_offset(p) = p.log_parties * aimer_field_bytes(p)
aimer_delta_offset(p) = aimer_missing_commitment_offset(p) + aimer_commitment_bytes(p)
aimer_missing_alpha_offset(p) = aimer_delta_offset(p) +
    (aimer_input_sboxes(p) + 2) * aimer_field_bytes(p)

function aimer_update_fields!(hash, parameters, fields)
    for value in fields
        aimer_hash_update!(hash, aimer_field_bytes(parameters, value))
    end
end

function aimer_message_digest(parameters, public_key, prefix, message)
    hash = aimer_prefixed_hash(parameters, parameters.message_prefix)
    aimer_hash_update!(hash, public_key)
    aimer_hash_update!(hash, prefix)
    aimer_hash_update!(hash, message)
    digest = zeros(UInt8, aimer_commitment_bytes(parameters))
    aimer_hash_squeeze!(hash, digest)
    return digest
end

function aimer_challenged_indices(parameters, h2)
    hash = AimerHash(parameters)
    aimer_hash_update!(hash, h2)
    bytes = zeros(UInt8, parameters.repetitions)
    aimer_hash_squeeze!(hash, bytes)
    return [Int(byte) & (parameters.parties - 1) for byte in bytes]
end

function aimer_fields_from_bytes(parameters, input, count)
    size = aimer_field_bytes(parameters)
    return [aimer_field_from_bytes(parameters, input, index * size)
            for index in 0:(count - 1)]
end

function aimer_equal(left, right)
    length(left) == length(right) || return false
    difference = 0
    for index in eachindex(left)
        difference |= xor(left[index], right[index])
    end
    return difference == 0
end

function aimer_generate_keypair(parameters; entropy_source=random_bytes)
    size = aimer_field_bytes(parameters)
    while true
        plaintext = entropy_source(size)
        iv = entropy_source(size)
        ciphertext = aimer_encrypt(parameters, plaintext, iv)
        ciphertext === nothing && continue
        public_key = concatenate_bytes(iv, ciphertext)
        secret_key = concatenate_bytes(plaintext, public_key)
        return (; public_key, secret_key)
    end
end

function aimer_sign(parameters, message, secret_key, context;
                    entropy_source=random_bytes)
    length(context) <= 255 || throw(ArgumentError("context must not exceed 255 bytes"))
    length(secret_key) == aimer_secret_key_bytes(parameters) ||
        throw(DimensionMismatch("invalid AIMer secret key length"))

    size = aimer_field_bytes(parameters)
    signature = zeros(UInt8, parameters.signature_bytes)
    randomness = entropy_source(size)
    prefix = concatenate_bytes(UInt8[length(context)], context)
    public_key = copy(secret_key[(size + 1):end])
    digest = aimer_message_digest(parameters, public_key, prefix, message)
    plaintext = aimer_field_from_bytes(parameters, secret_key)
    ciphertext = aimer_field_from_bytes(parameters, secret_key, 2size)
    linear = aimer_generate_linear(parameters, @view(secret_key[(size + 1):(2size)]))
    outputs = aimer_sbox_outputs(parameters, linear, plaintext, ciphertext)

    seed_hash = aimer_prefixed_hash(parameters, 3)
    aimer_hash_update!(seed_hash, secret_key, 0, size)
    aimer_hash_update!(seed_hash, digest)
    aimer_hash_update!(seed_hash, randomness)
    aimer_hash_squeeze!(seed_hash, signature, 0, size)
    salt = copy(signature[1:size])
    roots = zeros(UInt8, parameters.repetitions * size)
    aimer_hash_squeeze!(seed_hash, roots)

    h1_hash = aimer_prefixed_hash(parameters, 1)
    aimer_hash_update!(h1_hash, digest)
    aimer_hash_update!(h1_hash, salt)
    all_nodes = Vector{Vector{Vector{UInt8}}}()
    all_commitments = Vector{Vector{Vector{UInt8}}}()
    checks = Vector{Vector{AimerCheck}}()

    for repetition in 0:(parameters.repetitions - 1)
        root = @view(roots[(repetition * size + 1):((repetition + 1) * size)])
        nodes = aimer_expand_tree(parameters, salt, repetition, root)
        commitments = Vector{Vector{UInt8}}()
        repetition_checks = AimerCheck[]
        delta = AimerTape(parameters)
        delta_bytes = zeros(UInt8, (aimer_input_sboxes(parameters) + 2) * size)
        for party in 0:(parameters.parties - 1)
            expanded = aimer_commit_and_expand_tape(
                parameters,
                salt,
                repetition,
                party,
                nodes[party + parameters.parties],
            )
            push!(commitments, expanded.commitment)
            tape = expanded.tape
            delta.plaintext = xor(delta.plaintext, tape.plaintext)
            for index in eachindex(delta.y)
                delta.y[index] = xor(delta.y[index], tape.y[index])
            end
            for index in eachindex(delta.a)
                delta.a[index] = xor(delta.a[index], tape.a[index])
            end
            delta.c = xor(delta.c, tape.c)
            if party == parameters.parties - 1
                delta.plaintext = xor(delta.plaintext, plaintext)
                aimer_field_to_bytes!(parameters, delta.plaintext, delta_bytes)
                tape.plaintext = xor(tape.plaintext, delta.plaintext)
                for index in 1:aimer_input_sboxes(parameters)
                    delta.y[index] = xor(delta.y[index], outputs[index])
                    aimer_field_to_bytes!(parameters, delta.y[index], delta_bytes, index * size)
                    tape.y[index] = xor(tape.y[index], delta.y[index])
                end
                for index in eachindex(delta.a)
                    delta.c = xor(
                        delta.c,
                        aimer_field_multiply(parameters, delta.a[index], outputs[index]),
                    )
                end
                aimer_field_to_bytes!(parameters, delta.c, delta_bytes,
                    (aimer_input_sboxes(parameters) + 1) * size)
                tape.c = xor(tape.c, delta.c)
            end
            push!(repetition_checks,
                aimer_mpc(parameters, linear, tape, ciphertext, party))
        end
        for commitment in commitments
            aimer_hash_update!(h1_hash, commitment)
        end
        aimer_hash_update!(h1_hash, delta_bytes)
        destination = aimer_proof_start(parameters, repetition) +
                      aimer_delta_offset(parameters)
        copyto!(signature, destination + 1, delta_bytes, 1, length(delta_bytes))
        push!(all_nodes, nodes)
        push!(all_commitments, commitments)
        push!(checks, repetition_checks)
    end

    h1_offset = size
    aimer_hash_squeeze!(h1_hash, signature, h1_offset, aimer_commitment_bytes(parameters))
    h1 = copy(signature[(h1_offset + 1):(h1_offset + aimer_commitment_bytes(parameters))])
    epsilon_hash = AimerHash(parameters)
    aimer_hash_update!(epsilon_hash, h1)
    h2_hash = aimer_prefixed_hash(parameters, 2)
    aimer_hash_update!(h2_hash, h1)
    aimer_hash_update!(h2_hash, salt)
    alpha_shares = Vector{Vector{Vector{BigInt}}}()

    for repetition in 0:(parameters.repetitions - 1)
        epsilon_bytes = zeros(UInt8, (aimer_input_sboxes(parameters) + 1) * size)
        aimer_hash_squeeze!(epsilon_hash, epsilon_bytes)
        epsilons = aimer_fields_from_bytes(
            parameters,
            epsilon_bytes,
            aimer_input_sboxes(parameters) + 1,
        )
        alpha = fill(BigInt(0), aimer_input_sboxes(parameters) + 1)
        v_shares = fill(BigInt(0), parameters.parties)
        repetition_shares = [fill(BigInt(0), aimer_input_sboxes(parameters) + 1)
                             for _ in 1:parameters.parties]
        push!(alpha_shares, repetition_shares)
        for party in 1:parameters.parties
            check = checks[repetition + 1][party]
            v_shares[party] = check.c
            for index in 1:(aimer_input_sboxes(parameters) + 1)
                value = xor(
                    check.a[index],
                    aimer_field_multiply(parameters, check.x[index], epsilons[index]),
                )
                repetition_shares[party][index] = value
                alpha[index] = xor(alpha[index], value)
                v_shares[party] = xor(
                    v_shares[party],
                    aimer_field_multiply(parameters, check.z[index], epsilons[index]),
                )
            end
        end
        for party in 1:parameters.parties, index in eachindex(alpha)
            v_shares[party] = xor(
                v_shares[party],
                aimer_field_multiply(
                    parameters,
                    checks[repetition + 1][party].b[index],
                    alpha[index],
                ),
            )
        end
        for shares in repetition_shares
            aimer_update_fields!(h2_hash, parameters, shares)
        end
        aimer_update_fields!(h2_hash, parameters, v_shares)
    end

    h2_offset = size + aimer_commitment_bytes(parameters)
    aimer_hash_squeeze!(h2_hash, signature, h2_offset, aimer_commitment_bytes(parameters))
    h2 = copy(signature[(h2_offset + 1):(h2_offset + aimer_commitment_bytes(parameters))])
    challenged = aimer_challenged_indices(parameters, h2)
    for repetition in eachindex(challenged)
        covered = challenged[repetition]
        path = aimer_reveal_path(parameters, all_nodes[repetition], covered)
        proof = aimer_proof_start(parameters, repetition - 1)
        copyto!(signature, proof + 1, path, 1, length(path))
        commitment_offset = proof + aimer_missing_commitment_offset(parameters)
        missing_commitment = all_commitments[repetition][covered + 1]
        copyto!(signature, commitment_offset + 1,
            missing_commitment, 1, length(missing_commitment))
        missing = proof + aimer_missing_alpha_offset(parameters)
        for index in 0:aimer_input_sboxes(parameters)
            aimer_field_to_bytes!(
                parameters,
                alpha_shares[repetition][covered + 1][index + 1],
                signature,
                missing + index * size,
            )
        end
    end
    return signature
end

function aimer_verify(parameters, message, signature, public_key, context)
    if length(context) > 255 || length(signature) != parameters.signature_bytes ||
       length(public_key) != aimer_public_key_bytes(parameters)
        return false
    end
    size = aimer_field_bytes(parameters)
    commitment_bytes = aimer_commitment_bytes(parameters)
    salt = copy(signature[1:size])
    h1 = copy(signature[(size + 1):(size + commitment_bytes)])
    h2 = copy(signature[(size + commitment_bytes + 1):(size + 2commitment_bytes)])
    linear = aimer_generate_linear(parameters, @view(public_key[1:size]))
    ciphertext = aimer_field_from_bytes(parameters, public_key, size)
    challenged = aimer_challenged_indices(parameters, h2)
    epsilon_hash = AimerHash(parameters)
    aimer_hash_update!(epsilon_hash, h1)
    prefix = concatenate_bytes(UInt8[length(context)], context)
    digest = aimer_message_digest(parameters, public_key, prefix, message)
    h1_hash = aimer_prefixed_hash(parameters, 1)
    aimer_hash_update!(h1_hash, digest)
    aimer_hash_update!(h1_hash, salt)
    h2_hash = aimer_prefixed_hash(parameters, 2)
    aimer_hash_update!(h2_hash, h1)
    aimer_hash_update!(h2_hash, salt)

    for repetition in eachindex(challenged)
        covered = challenged[repetition]
        proof = aimer_proof_start(parameters, repetition - 1)
        path_length = parameters.log_parties * size
        path = copy(signature[(proof + 1):(proof + path_length)])
        nodes = aimer_reconstruct_tree(parameters, salt, path, repetition - 1, covered)
        epsilon_bytes = zeros(UInt8, (aimer_input_sboxes(parameters) + 1) * size)
        aimer_hash_squeeze!(epsilon_hash, epsilon_bytes)
        epsilons = aimer_fields_from_bytes(
            parameters,
            epsilon_bytes,
            aimer_input_sboxes(parameters) + 1,
        )
        alpha = fill(BigInt(0), aimer_input_sboxes(parameters) + 1)
        alpha_shares = [fill(BigInt(0), aimer_input_sboxes(parameters) + 1)
                        for _ in 1:parameters.parties]
        b_shares = [fill(BigInt(0), aimer_input_sboxes(parameters) + 1)
                    for _ in 1:parameters.parties]
        v_shares = fill(BigInt(0), parameters.parties)

        for party in 0:(parameters.parties - 1)
            if party == covered
                commitment_offset = proof + aimer_missing_commitment_offset(parameters)
                aimer_hash_update!(h1_hash, signature, commitment_offset, commitment_bytes)
                missing = proof + aimer_missing_alpha_offset(parameters)
                for index in 0:aimer_input_sboxes(parameters)
                    value = aimer_field_from_bytes(
                        parameters,
                        signature,
                        missing + index * size,
                    )
                    alpha_shares[party + 1][index + 1] = value
                    alpha[index + 1] = xor(alpha[index + 1], value)
                end
                continue
            end
            expanded = aimer_commit_and_expand_tape(
                parameters,
                salt,
                repetition - 1,
                party,
                nodes[parameters.parties + party - 1],
            )
            aimer_hash_update!(h1_hash, expanded.commitment)
            tape = expanded.tape
            if party == parameters.parties - 1
                delta = proof + aimer_delta_offset(parameters)
                tape.plaintext = xor(
                    tape.plaintext,
                    aimer_field_from_bytes(parameters, signature, delta),
                )
                delta += size
                for index in eachindex(tape.y)
                    tape.y[index] = xor(
                        tape.y[index],
                        aimer_field_from_bytes(parameters, signature, delta),
                    )
                    delta += size
                end
                tape.c = xor(
                    tape.c,
                    aimer_field_from_bytes(parameters, signature, delta),
                )
            end
            check = aimer_mpc(parameters, linear, tape, ciphertext, party)
            v_shares[party + 1] = check.c
            for index in 1:(aimer_input_sboxes(parameters) + 1)
                alpha_value = xor(
                    check.a[index],
                    aimer_field_multiply(parameters, check.x[index], epsilons[index]),
                )
                alpha_shares[party + 1][index] = alpha_value
                alpha[index] = xor(alpha[index], alpha_value)
                b_shares[party + 1][index] = check.b[index]
                v_shares[party + 1] = xor(
                    v_shares[party + 1],
                    aimer_field_multiply(parameters, check.z[index], epsilons[index]),
                )
            end
        end

        for party in 0:(parameters.parties - 1)
            party == covered && continue
            for index in eachindex(alpha)
                v_shares[party + 1] = xor(
                    v_shares[party + 1],
                    aimer_field_multiply(
                        parameters,
                        b_shares[party + 1][index],
                        alpha[index],
                    ),
                )
            end
            v_shares[covered + 1] = xor(v_shares[covered + 1], v_shares[party + 1])
        end
        for shares in alpha_shares
            aimer_update_fields!(h2_hash, parameters, shares)
        end
        aimer_update_fields!(h2_hash, parameters, v_shares)
        delta_offset = proof + aimer_delta_offset(parameters)
        aimer_hash_update!(h1_hash, signature, delta_offset,
            (aimer_input_sboxes(parameters) + 2) * size)
    end

    computed_h1 = zeros(UInt8, commitment_bytes)
    computed_h2 = zeros(UInt8, commitment_bytes)
    aimer_hash_squeeze!(h1_hash, computed_h1)
    aimer_hash_squeeze!(h2_hash, computed_h2)
    return aimer_equal(computed_h1, h1) && aimer_equal(computed_h2, h2)
end
