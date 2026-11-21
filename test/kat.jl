using KpqC
using Test

const VECTOR_ROOT = get(ENV, "KPQC_TEST_VECTORS", "")
isempty(VECTOR_ROOT) && error("KPQC_TEST_VECTORS is not set")
isdir(VECTOR_ROOT) || error("KPQC_TEST_VECTORS does not name a directory")

function read_vectors(path::String)
    vectors = Dict{String,String}[]
    current = Dict{String,String}()
    for line in eachline(path)
        parts = split(line, " = "; limit = 2)
        length(parts) == 2 || continue
        if parts[1] == "count" && !isempty(current)
            push!(vectors, current)
            current = Dict{String,String}()
        end
        current[parts[1]] = parts[2]
    end
    isempty(current) || push!(vectors, current)
    return vectors
end

field(vector::Dict{String,String}, name::String) = hex2bytes(vector[name])

mutable struct EntropySource
    data::Vector{UInt8}
    offset::Int
end

EntropySource(data::Vector{UInt8}) = EntropySource(data, 1)

function (source::EntropySource)(length_value::Int)
    length_value >= 0 || error("entropy request must be nonnegative")
    final_offset = source.offset + length_value - 1
    final_offset <= length(source.data) || error("deterministic entropy is exhausted")
    output = copy(@view source.data[source.offset:final_offset])
    source.offset = final_offset + 1
    return output
end

entropy_remaining(source::EntropySource) = length(source.data) - source.offset + 1

function with_exact_entropy(operation, entropy::Vector{UInt8})
    source = EntropySource(entropy)
    result = operation(source)
    remaining = entropy_remaining(source)
    remaining == 0 || error("deterministic entropy has $remaining unused bytes")
    return result
end

mutable struct DRBG
    key::Vector{UInt8}
    counter::Vector{UInt8}
end

copy_drbg(state::DRBG) = DRBG(copy(state.key), copy(state.counter))

function increment_counter(counter::Vector{UInt8})
    result = copy(counter)
    for index in reverse(eachindex(result))
        result[index] += 0x01
        result[index] == 0x00 || break
    end
    return result
end

const AES_SBOX = UInt8[
    0x63,0x7c,0x77,0x7b,0xf2,0x6b,0x6f,0xc5,0x30,0x01,0x67,0x2b,0xfe,0xd7,0xab,0x76,
    0xca,0x82,0xc9,0x7d,0xfa,0x59,0x47,0xf0,0xad,0xd4,0xa2,0xaf,0x9c,0xa4,0x72,0xc0,
    0xb7,0xfd,0x93,0x26,0x36,0x3f,0xf7,0xcc,0x34,0xa5,0xe5,0xf1,0x71,0xd8,0x31,0x15,
    0x04,0xc7,0x23,0xc3,0x18,0x96,0x05,0x9a,0x07,0x12,0x80,0xe2,0xeb,0x27,0xb2,0x75,
    0x09,0x83,0x2c,0x1a,0x1b,0x6e,0x5a,0xa0,0x52,0x3b,0xd6,0xb3,0x29,0xe3,0x2f,0x84,
    0x53,0xd1,0x00,0xed,0x20,0xfc,0xb1,0x5b,0x6a,0xcb,0xbe,0x39,0x4a,0x4c,0x58,0xcf,
    0xd0,0xef,0xaa,0xfb,0x43,0x4d,0x33,0x85,0x45,0xf9,0x02,0x7f,0x50,0x3c,0x9f,0xa8,
    0x51,0xa3,0x40,0x8f,0x92,0x9d,0x38,0xf5,0xbc,0xb6,0xda,0x21,0x10,0xff,0xf3,0xd2,
    0xcd,0x0c,0x13,0xec,0x5f,0x97,0x44,0x17,0xc4,0xa7,0x7e,0x3d,0x64,0x5d,0x19,0x73,
    0x60,0x81,0x4f,0xdc,0x22,0x2a,0x90,0x88,0x46,0xee,0xb8,0x14,0xde,0x5e,0x0b,0xdb,
    0xe0,0x32,0x3a,0x0a,0x49,0x06,0x24,0x5c,0xc2,0xd3,0xac,0x62,0x91,0x95,0xe4,0x79,
    0xe7,0xc8,0x37,0x6d,0x8d,0xd5,0x4e,0xa9,0x6c,0x56,0xf4,0xea,0x65,0x7a,0xae,0x08,
    0xba,0x78,0x25,0x2e,0x1c,0xa6,0xb4,0xc6,0xe8,0xdd,0x74,0x1f,0x4b,0xbd,0x8b,0x8a,
    0x70,0x3e,0xb5,0x66,0x48,0x03,0xf6,0x0e,0x61,0x35,0x57,0xb9,0x86,0xc1,0x1d,0x9e,
    0xe1,0xf8,0x98,0x11,0x69,0xd9,0x8e,0x94,0x9b,0x1e,0x87,0xe9,0xce,0x55,0x28,0xdf,
    0x8c,0xa1,0x89,0x0d,0xbf,0xe6,0x42,0x68,0x41,0x99,0x2d,0x0f,0xb0,0x54,0xbb,0x16,
]

aes_xtime(value::UInt8) = xor(value << 1, (value >> 7) * UInt8(0x1b))

function aes_key_schedule(key::Vector{UInt8})
    length(key) == 32 || error("AES-256 key must contain 32 bytes")
    expanded = Vector{UInt8}(undef, 240)
    expanded[1:32] .= key
    generated = 32
    rcon = UInt8(1)
    while generated < 240
        temporary = copy(expanded[generated-3:generated])
        if generated % 32 == 0
            temporary = UInt8[AES_SBOX[Int(temporary[2])+1], AES_SBOX[Int(temporary[3])+1],
                AES_SBOX[Int(temporary[4])+1], AES_SBOX[Int(temporary[1])+1]]
            temporary[1] = xor(temporary[1], rcon)
            rcon = aes_xtime(rcon)
        elseif generated % 32 == 16
            temporary .= AES_SBOX[Int.(temporary) .+ 1]
        end
        for index in 1:4
            generated += 1
            expanded[generated] = xor(expanded[generated - 32], temporary[index])
        end
    end
    return expanded
end

function aes_mix_columns!(state)
    for column in 0:3
        offset = 4column + 1
        a0, a1, a2, a3 = state[offset:offset+3]
        total = xor(a0, a1, a2, a3)
        state[offset] = xor(a0, total, aes_xtime(xor(a0, a1)))
        state[offset+1] = xor(a1, total, aes_xtime(xor(a1, a2)))
        state[offset+2] = xor(a2, total, aes_xtime(xor(a2, a3)))
        state[offset+3] = xor(a3, total, aes_xtime(xor(a3, a0)))
    end
end

function aes256(key::Vector{UInt8}, input::Vector{UInt8})
    length(input) == 16 || error("AES block must contain 16 bytes")
    round_keys = aes_key_schedule(key)
    state = xor.(input, round_keys[1:16])
    for round in 1:14
        state .= AES_SBOX[Int.(state) .+ 1]
        shifted = similar(state)
        for row in 0:3, column in 0:3
            shifted[row + 4column + 1] = state[row + 4mod(column + row, 4) + 1]
        end
        state = shifted
        round < 14 && aes_mix_columns!(state)
        state .= xor.(state, round_keys[16round+1:16round+16])
    end
    return state
end

function generate_blocks!(state::DRBG, count::Int)
    output = Vector{UInt8}[]
    for _ in 1:count
        state.counter = increment_counter(state.counter)
        push!(output, aes256(state.key, state.counter))
    end
    return output
end

function update!(state::DRBG, provided::Vector{UInt8})
    length(provided) <= 48 || error("DRBG update data must not exceed 48 bytes")
    material = vcat(generate_blocks!(state, 3)...)
    padded = vcat(provided, zeros(UInt8, 48 - length(provided)))
    material .⊻= padded
    state.key = material[1:32]
    state.counter = material[33:48]
    return state
end

function DRBG(seed::Vector{UInt8})
    length(seed) == 48 || error("CTR-DRBG seed must contain 48 bytes")
    return update!(DRBG(zeros(UInt8, 32), zeros(UInt8, 16)), seed)
end

function generate!(state::DRBG, length_value::Int)
    blocks = cld(length_value, 16)
    output = blocks == 0 ? UInt8[] : vcat(generate_blocks!(state, blocks)...)
    update!(state, UInt8[])
    return output[1:length_value]
end

function with_repeated_entropy(
    operation,
    state::DRBG,
    call_bytes::Int;
    max_calls::Int = 128,
)
    probe = copy_drbg(state)
    entropy = vcat([generate!(probe, call_bytes) for _ in 1:max_calls]...)
    source = EntropySource(entropy)
    result = operation(source)
    consumed = length(entropy) - entropy_remaining(source)
    consumed > 0 || error("deterministic entropy was not consumed")
    consumed % call_bytes == 0 || error("entropy consumption was not aligned")
    for _ in 1:(consumed ÷ call_bytes)
        generate!(state, call_bytes)
    end
    return result
end

const AIMER_PARAMETER_KEYS = Dict(
    AIMer128f => :aimer128f,
    AIMer128s => :aimer128s,
    AIMer192f => :aimer192f,
    AIMer192s => :aimer192s,
    AIMer256f => :aimer256f,
    AIMer256s => :aimer256s,
)
const HAETAE_PARAMETER_MODES = Dict(HAETAE2 => 2, HAETAE3 => 3, HAETAE5 => 5)
const NTRUPLUS_PARAMETER_MODES = Dict(
    NTRUPlus768 => 768,
    NTRUPlus864 => 864,
    NTRUPlus1152 => 1152,
)
const SMAUGT_PARAMETER_MODES = Dict(
    SMAUGT128 => 128,
    SMAUGT192 => 192,
    SMAUGT256 => 256,
    TiMER => :timer,
)

function kat_generate_keypair(
    algorithm::SignatureAlgorithm,
    entropy_source::EntropySource,
)
    if haskey(AIMER_PARAMETER_KEYS, algorithm)
        generated = KpqC.aimer_generate_keypair(
            KpqC.AIMER_PARAMETERS[AIMER_PARAMETER_KEYS[algorithm]];
            entropy_source,
        )
        return KeyPair(generated.public_key, generated.secret_key)
    end

    public_key, secret_key = KpqC.haetae_generate_keypair(
        KpqC.HAETAE_PARAMETERS[HAETAE_PARAMETER_MODES[algorithm]];
        seed = entropy_source(32),
    )
    return KeyPair(public_key, secret_key)
end

function kat_sign(
    algorithm::SignatureAlgorithm,
    message::Vector{UInt8},
    secret_key::Vector{UInt8},
    context::Vector{UInt8},
    entropy_source::EntropySource,
)
    if haskey(AIMER_PARAMETER_KEYS, algorithm)
        return KpqC.aimer_sign(
            KpqC.AIMER_PARAMETERS[AIMER_PARAMETER_KEYS[algorithm]],
            message,
            secret_key,
            context;
            entropy_source,
        )
    end

    return KpqC.haetae_sign(
        KpqC.HAETAE_PARAMETERS[HAETAE_PARAMETER_MODES[algorithm]],
        message,
        secret_key,
        context;
        randomness = entropy_source(32),
    )
end

function kat_generate_keypair(
    algorithm::KEMAlgorithm,
    entropy_source::EntropySource,
)
    if haskey(NTRUPLUS_PARAMETER_MODES, algorithm)
        generated = KpqC.ntruplus_generate_keypair(
            KpqC.NTRUPLUS_PARAMETERS[NTRUPLUS_PARAMETER_MODES[algorithm]];
            entropy_source,
        )
    else
        generated = KpqC.smaugt_generate_keypair(
            KpqC.SMAUGT_PARAMETERS[SMAUGT_PARAMETER_MODES[algorithm]];
            entropy_source,
        )
    end
    return KeyPair(generated.public_key, generated.secret_key)
end

function kat_encapsulate(
    algorithm::KEMAlgorithm,
    public_key::Vector{UInt8},
    entropy_source::EntropySource,
)
    result = if haskey(NTRUPLUS_PARAMETER_MODES, algorithm)
        KpqC.ntruplus_encapsulate(
            KpqC.NTRUPLUS_PARAMETERS[NTRUPLUS_PARAMETER_MODES[algorithm]],
            public_key;
            entropy_source,
        )
    else
        KpqC.smaugt_encapsulate(
            KpqC.SMAUGT_PARAMETERS[SMAUGT_PARAMETER_MODES[algorithm]],
            public_key;
            entropy_source,
        )
    end
    return EncapsulatedSecret(result.ciphertext, result.shared_secret)
end

const SIGNATURE_CONFIGURATIONS = [
    (AIMer128f, "aimer/128f", false, 16),
    (AIMer128s, "aimer/128s", false, 16),
    (AIMer192f, "aimer/192f", false, 24),
    (AIMer192s, "aimer/192s", false, 24),
    (AIMer256f, "aimer/256f", false, 32),
    (AIMer256s, "aimer/256s", false, 32),
    (HAETAE2, "haetae/mode2", true, 32),
    (HAETAE3, "haetae/mode3", true, 32),
    (HAETAE5, "haetae/mode5", true, 32),
]

const KEM_CONFIGURATIONS = [
    (NTRUPlus768, "ntruplus/768", 0, 96),
    (NTRUPlus864, "ntruplus/864", 0, 108),
    (NTRUPlus1152, "ntruplus/1152", 0, 144),
    (SMAUGT128, "smaugt/mode1", 2, 32),
    (SMAUGT192, "smaugt/mode3", 2, 32),
    (SMAUGT256, "smaugt/mode5", 2, 32),
    (TiMER, "smaugt/modet", 2, 16),
]

const KAT_FILTER = get(ENV, "KPQC_KAT_FILTER", "all")
const KAT_RECORD_LIMIT = parse(Int, get(ENV, "KPQC_KAT_RECORD_LIMIT", "100"))
const KAT_RECORD_OFFSET = parse(Int, get(ENV, "KPQC_KAT_RECORD_OFFSET", "0"))
1 <= KAT_RECORD_LIMIT <= 100 || error("KPQC_KAT_RECORD_LIMIT must be between 1 and 100")
0 <= KAT_RECORD_OFFSET < 100 || error("KPQC_KAT_RECORD_OFFSET must be between 0 and 99")
KAT_RECORD_OFFSET + KAT_RECORD_LIMIT <= 100 ||
    error("KPQC_KAT_RECORD_OFFSET + KPQC_KAT_RECORD_LIMIT must not exceed 100")
const SELECTED_SIGNATURE_CONFIGURATIONS = if KAT_FILTER == "all"
    SIGNATURE_CONFIGURATIONS
elseif KAT_FILTER == "haetae"
    filter(configuration -> startswith(configuration[2], "haetae/"), SIGNATURE_CONFIGURATIONS)
elseif KAT_FILTER == "kem"
    eltype(SIGNATURE_CONFIGURATIONS)[]
else
    filter(configuration -> configuration[2] == KAT_FILTER, SIGNATURE_CONFIGURATIONS)
end
const SELECTED_KEM_CONFIGURATIONS = if KAT_FILTER in ("all", "kem")
    KEM_CONFIGURATIONS
elseif KAT_FILTER == "haetae" || startswith(KAT_FILTER, "aimer/")
    eltype(KEM_CONFIGURATIONS)[]
else
    filter(configuration -> configuration[2] == KAT_FILTER, KEM_CONFIGURATIONS)
end
isempty(SELECTED_SIGNATURE_CONFIGURATIONS) && isempty(SELECTED_KEM_CONFIGURATIONS) &&
    error("unknown KPQC_KAT_FILTER: $KAT_FILTER")

@testset "known-answer tests" begin
    total = 0

    for (algorithm, relative_path, is_haetae, seed_bytes) in SELECTED_SIGNATURE_CONFIGURATIONS
        available_vectors = read_vectors(joinpath(VECTOR_ROOT, relative_path, "kat.rsp"))
        @test length(available_vectors) == 100
        vectors = available_vectors[
            (KAT_RECORD_OFFSET + 1):(KAT_RECORD_OFFSET + KAT_RECORD_LIMIT)
        ]
        for vector in vectors
            state = DRBG(field(vector, "seed"))
            keys = if is_haetae
                with_exact_entropy(
                    entropy_source -> kat_generate_keypair(algorithm, entropy_source),
                    generate!(state, seed_bytes),
                )
            else
                with_repeated_entropy(
                    entropy_source -> kat_generate_keypair(algorithm, entropy_source),
                    state,
                    seed_bytes,
                )
            end
            @test keys.public_key == field(vector, "pk")
            @test keys.secret_key == field(vector, "sk")

            message = field(vector, "msg")
            expected_signature = if haskey(vector, "sig")
                field(vector, "sig")
            else
                field(vector, "sm")[(length(message) + 1):end]
            end
            signing_entropy = generate!(state, seed_bytes)
            context = if is_haetae
                context_length = Int(only(generate!(state, 1)))
                generate!(state, context_length)
            else
                UInt8[]
            end
            signature = with_exact_entropy(
                entropy_source -> kat_sign(
                    algorithm,
                    message,
                    keys.secret_key,
                    context,
                    entropy_source,
                ),
                signing_entropy,
            )
            @test signature == expected_signature
            @test verify(algorithm, message, expected_signature, keys.public_key; context)
        end
        total += length(vectors)
    end

    for (algorithm, relative_path, key_calls, encapsulation_bytes) in SELECTED_KEM_CONFIGURATIONS
        available_vectors = read_vectors(joinpath(VECTOR_ROOT, relative_path, "kat.rsp"))
        @test length(available_vectors) == 100
        vectors = available_vectors[
            (KAT_RECORD_OFFSET + 1):(KAT_RECORD_OFFSET + KAT_RECORD_LIMIT)
        ]
        for vector in vectors
            state = DRBG(field(vector, "seed"))
            keys = if key_calls == 0
                with_repeated_entropy(
                    entropy_source -> kat_generate_keypair(algorithm, entropy_source),
                    state,
                    32,
                )
            else
                key_entropy = vcat([generate!(state, 32) for _ in 1:key_calls]...)
                with_exact_entropy(
                    entropy_source -> kat_generate_keypair(algorithm, entropy_source),
                    key_entropy,
                )
            end
            @test keys.public_key == field(vector, "pk")
            @test keys.secret_key == field(vector, "sk")

            outbound = with_exact_entropy(
                entropy_source -> kat_encapsulate(
                    algorithm,
                    keys.public_key,
                    entropy_source,
                ),
                generate!(state, encapsulation_bytes),
            )
            @test outbound.ciphertext == field(vector, "ct")
            @test outbound.shared_secret == field(vector, "ss")
            @test decapsulate(algorithm, outbound.ciphertext, keys.secret_key) ==
                  outbound.shared_secret
        end
        total += length(vectors)
    end

    @test total == KAT_RECORD_LIMIT * (
        length(SELECTED_SIGNATURE_CONFIGURATIONS) + length(SELECTED_KEM_CONFIGURATIONS)
    )
end
