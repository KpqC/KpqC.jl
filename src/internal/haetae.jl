const HAETAE_N = 256
const HAETAE_Q = 64_513
const HAETAE_DOUBLE_Q = 129_026
const HAETAE_SEED_BYTES = 32
const HAETAE_CRH_BYTES = 64
const HAETAE_FIXED_SCALE = 8_192

include("haetae_cdt.jl")

struct HaetaeParameters
    mode::Int
    k::Int
    l::Int
    m::Int
    tau::Int
    b0_squared::BigInt
    b1_squared::BigInt
    b2_squared::BigInt
    singular_value_limit::Int
    base_encoded_high_z1::Int
    base_encoded_hint::Int
    alpha_hint::Int
    log_alpha_hint::Int
    poly_q_bytes::Int
    poly_2eta_bytes::Int
    public_key_bytes::Int
    secret_key_bytes::Int
    signature_bytes::Int
end

function haetae_parameters(mode::Int, k::Int, l::Int, tau::Int,
                           b0_squared::Int, b1_squared::Int, b2_squared::Int,
                           singular_value_limit::Int, base_high::Int, base_hint::Int,
                           alpha::Int, poly_q::Int, poly_eta::Int,
                           signature::Int)
    m = l - 1
    public_key = 32 + k * poly_q
    return HaetaeParameters(
        mode, k, l, m, tau,
        BigInt(b0_squared), BigInt(b1_squared), BigInt(b2_squared),
        singular_value_limit, base_high, base_hint, alpha, alpha == 256 ? 8 : 9,
        poly_q, poly_eta, public_key,
        public_key + m * 64 + k * poly_eta + 32, signature,
    )
end

const HAETAE_PARAMETERS = Dict(
    2 => haetae_parameters(2, 2, 4, 58, 96_944_109, 96_805_527, 163_265_017,
        611_098, 132, 7, 512, 480, 96, 1474),
    3 => haetae_parameters(3, 3, 6, 80, 335_438_492, 335_171_879, 479_901_314,
        852_505, 376, 127, 512, 480, 96, 2349),
    5 => haetae_parameters(5, 4, 7, 128, 499_239_142, 498_849_991, 597_386_433,
        778_065, 501, 358, 256, 512, 64, 2948),
)

haetae_poly() = zeros(Int64, HAETAE_N)
haetae_vector(length::Int) = [haetae_poly() for _ in 1:length]
haetae_clone(value) = [copy(item) for item in value]

function haetae_montgomery_reduce(value::Integer)
    low = UInt32(mod(Int128(value) * Int128(940_508_161), Int128(1) << 32))
    multiplier = Int64(reinterpret(Int32, low))
    return Int64((Int128(value) - Int128(multiplier) * HAETAE_Q) >> 32)
end

haetae_freeze(value::Integer) = Int64(mod(value, HAETAE_Q))
haetae_freeze_2q(value::Integer) = Int64(mod(value, HAETAE_DOUBLE_Q))

function haetae_centered_2q(value::Integer)
    result = Int64(rem(value, HAETAE_DOUBLE_Q))
    result < -HAETAE_Q && (result += HAETAE_DOUBLE_Q)
    result >= HAETAE_Q && (result -= HAETAE_DOUBLE_Q)
    return result
end

function haetae_decompose_z1(value::Integer)
    low = Int(value) & 255
    low >= 128 && (low -= 256)
    return (Int64((Int(value) + 128) >> 8), Int64(low))
end

function haetae_decompose_hint(value::Integer, alpha::Int, log_alpha::Int)
    high = (Int(value) + alpha ÷ 2) >> log_alpha
    maximum = (HAETAE_DOUBLE_Q - 2) ÷ alpha
    return high == maximum ? 0 : Int64(high)
end

function haetae_decompose_public_key(value::Integer)
    low = Int(value) & 1
    low -= ((Int(value) >> 1) & low) << 1
    return (Int64((Int(value) - low) >> 1), Int64(low))
end

haetae_le16(value::Integer) = UInt8[value & 0xff, (value >> 8) & 0xff]

function haetae_squared_norm(groups...)
    result = BigInt(0)
    for group in groups, coefficients in group, coefficient in coefficients
        result += BigInt(coefficient)^2
    end
    return result
end

const HAETAE_ZETAS = Int64[
  0, 26964, -16505, 22229, 30746, 20243, 19064, -31218, 9395, -30985,
  22859, -8851, 32144, 13744, 21408, 17599, -16039, -22946, 6241, -19553,
  10681, 22935, 22431, -29104, 28147, -27527, -29133, -20035, 20143,
  -11361, 30820, 25252, -22562, -6789, -10049, 9383, 16304, -12296, 16446,
  18239, -1296, -19725, -32076, 11782, -17941, 29643, -8577, 7893, -21464,
  -19646, -15130, -2391, 30608, -23970, -16608, 19616, -7941, 26533,
  -19129, 27690, 7597, -11459, 10615, -9430, 11591, 7814, 12697, 32114,
  -3761, -9604, 19813, 20353, 17456, -16267, -19555, 598, -29942, 4538,
  835, 15546, 3970, -27685, 1488, 8311, -12442, 31352, -17631, 1806,
  -5342, 9790, 29068, 16507, -29051, 22131, 6759, 15510, -14941, 28710,
  1160, -31327, 24985, 11261, -10623, -27727, 21502, 18731, -16186, -4127,
  -18832, 12050, -14501, 7929, 29563, -31064, 5913, 5322, -16405, 2844,
  29439, 5876, -9522, -18586, -9874, 23844, 30362, -21442, 9560, 17671,
  -27989, 3350, 787, -13857, 1657, -21224, -7374, -9190, 2464, 25555,
  -3529, -28772, 16588, -15739, 23475, 13666, 5764, 30980, 13633, -7401,
  -30317, 28847, 7682, -11808, -8796, 14864, -24162, -19194, 689, -1311,
  -31332, -16319, 1025, 10971, -23016, -2648, -21900, -12543, -25921,
  28254, 28521, -16160, 12380, -12882, -30332, -16630, 23439, 7742,
  17182, 17494, 5920, 13642, 7382, -18166, 21422, -30274, -28190, 13283,
  -20316, -9939, 10672, 21454, 6080, -17374, -29735, -25912, -10170, 3808,
  10639, -26985, -10865, 25636, 17261, -26851, -8253, -3304, 18282, -2202,
  -31368, -22243, 13882, 12069, -11242, -7729, -10226, 1761, -27298, -4800,
  -17737, -22805, -3528, 65, 10770, 8908, -23751, 26934, 21921, -27010,
  -21944, 8889, -1035, 23224, -9488, -5823, -994, -20206, 7655, -16251,
  -22820, -27740, 15822, 23078, 13803, -8099, 2931, 9217, -21126, -14203,
  25492, -12831, 7947, 17463, -12979, 29003, 31612, 26554, 8241, -20175,
]

function haetae_ntt!(values)
    cursor = 1
    length = 128
    while length > 0
        for start in 0:2length:HAETAE_N-1
            cursor += 1
            zeta = HAETAE_ZETAS[cursor]
            for index in start:start+length-1
                product = haetae_montgomery_reduce(zeta * values[index + length + 1])
                value = values[index + 1]
                values[index + length + 1] = value - product
                values[index + 1] = value + product
            end
        end
        length >>= 1
    end
    return values
end

function haetae_inverse_ntt!(values)
    cursor = 257
    length = 1
    while length < HAETAE_N
        for start in 0:2length:HAETAE_N-1
            cursor -= 1
            zeta = -HAETAE_ZETAS[cursor]
            for index in start:start+length-1
                value = values[index + 1]
                values[index + 1] = value + values[index + length + 1]
                values[index + length + 1] = haetae_montgomery_reduce(
                    zeta * (value - values[index + length + 1]))
            end
        end
        length <<= 1
    end
    for index in eachindex(values)
        values[index] = haetae_montgomery_reduce(-29_720 * values[index])
    end
    return values
end

haetae_transform(value) = haetae_ntt!(copy(value))
haetae_inverse_transform(value) = haetae_inverse_ntt!(copy(value))
haetae_transform_vector(value) = map(haetae_transform, value)
haetae_inverse_transform_vector(value) = map(haetae_inverse_transform, value)

const HAETAE_ROOT_QUADRANT = Int64[
    65_536, 65_531, 65_516, 65_492, 65_457, 65_413, 65_358, 65_294,
    65_220, 65_137, 65_043, 64_940, 64_827, 64_704, 64_571, 64_429,
    64_277, 64_115, 63_944, 63_763, 63_572, 63_372, 63_162, 62_943,
    62_714, 62_476, 62_228, 61_971, 61_705, 61_429, 61_145, 60_851,
    60_547, 60_235, 59_914, 59_583, 59_244, 58_896, 58_538, 58_172,
    57_798, 57_414, 57_022, 56_621, 56_212, 55_794, 55_368, 54_934,
    54_491, 54_040, 53_581, 53_114, 52_639, 52_156, 51_665, 51_166,
    50_660, 50_146, 49_624, 49_095, 48_559, 48_015, 47_464, 46_906,
    46_341, 45_769, 45_190, 44_604, 44_011, 43_412, 42_806, 42_194,
    41_576, 40_951, 40_320, 39_683, 39_040, 38_391, 37_736, 37_076,
    36_410, 35_738, 35_062, 34_380, 33_692, 33_000, 32_303, 31_600,
    30_893, 30_182, 29_466, 28_745, 28_020, 27_291, 26_558, 25_821,
    25_080, 24_335, 23_586, 22_834, 22_078, 21_320, 20_557, 19_792,
    19_024, 18_253, 17_479, 16_703, 15_924, 15_143, 14_359, 13_573,
    12_785, 11_996, 11_204, 10_411, 9_616, 8_820, 8_022, 7_224,
    6_424, 5_623, 4_821, 4_019, 3_216, 2_412, 1_608, 804, 0,
]
const HAETAE_ROOT_REAL = Int64[
    index <= 128 ? HAETAE_ROOT_QUADRANT[index + 1] :
    -HAETAE_ROOT_QUADRANT[257 - index]
    for index in 0:(HAETAE_N - 1)
]
const HAETAE_ROOT_IMAG = Int64[
    -HAETAE_ROOT_QUADRANT[abs(128 - index) + 1]
    for index in 0:(HAETAE_N - 1)
]

function haetae_reverse_bits8(value::Int)
    result = 0
    for bit in 0:7
        result = (result << 1) | ((value >> bit) & 1)
    end
    return result
end

haetae_shift16(value::Integer) = fld(value + 32_768, 65_536)

function haetae_fft(input)
    real = zeros(Int64, HAETAE_N)
    imaginary = zeros(Int64, HAETAE_N)
    for index in 0:HAETAE_N-1
        reversed = haetae_reverse_bits8(index)
        real[reversed + 1] = input[index + 1] * HAETAE_ROOT_REAL[index + 1]
        imaginary[reversed + 1] = input[index + 1] * HAETAE_ROOT_IMAG[index + 1]
    end
    for round in 1:8
        width = 1 << round
        half = width >> 1
        for start in 0:width:HAETAE_N-1, offset in 0:half-1
            even = start + offset + 1
            odd = even + half
            root_index = (offset << (9 - round)) + 1
            rr = HAETAE_ROOT_REAL[root_index]
            ri = HAETAE_ROOT_IMAG[root_index]
            odd_r = real[odd]
            odd_i = imaginary[odd]
            right_r = haetae_shift16(rr * odd_r) - haetae_shift16(ri * odd_i)
            right_i = haetae_shift16(rr * odd_i) + haetae_shift16(ri * odd_r)
            left_r = real[even]
            left_i = imaginary[even]
            real[even], imaginary[even] = left_r + right_r, left_i + right_i
            real[odd], imaginary[odd] = left_r - right_r, left_i - right_i
        end
    end
    return real, imaginary
end

function haetae_secret_singular_value(parameters, first, second)
    sums = zeros(Int64, HAETAE_N)
    for group in (first, second), item in group
        real, imaginary = haetae_fft(item)
        for index in eachindex(sums)
            sums[index] += haetae_shift16(real[index]^2) + haetae_shift16(imaginary[index]^2)
        end
    end
    selected_count = HAETAE_N ÷ parameters.tau + 1
    selected = sort(sums; rev=true)[1:selected_count]
    minimum = selected[end]
    result = 0
    for value in selected
        factor = value == minimum ? HAETAE_N % parameters.tau : parameters.tau
        result += fld(value + 0x10200, 1024) * factor
    end
    return fld(result + 32, 64)
end

haetae_add(left, right) = left .+ right
haetae_subtract(left, right) = left .- right
haetae_add_vectors(left, right) = [haetae_add(left[i], right[i]) for i in eachindex(left)]
haetae_subtract_vectors(left, right) = [haetae_subtract(left[i], right[i]) for i in eachindex(left)]
haeta_double_vector(value) = [item .* 2 for item in value]

function haetae_pointwise(left, right)
    return Int64[haetae_montgomery_reduce(left[i] * right[i]) for i in 1:HAETAE_N]
end

function haetae_pointwise_accumulate(left, right)
    output = haetae_poly()
    for item in eachindex(left), index in 1:HAETAE_N
        output[index] += haetae_montgomery_reduce(left[item][index] * right[item][index])
    end
    return output
end

haetae_freeze_poly(value) = Int64[haetae_freeze(coefficient) for coefficient in value]
haetae_freeze_poly_2q(value) = Int64[haetae_freeze_2q(coefficient) for coefficient in value]
haetae_centered_poly_2q(value) = Int64[haetae_centered_2q(coefficient) for coefficient in value]
haetae_freeze_vector(value) = map(haetae_freeze_poly, value)
haetae_freeze_vector_2q(value) = map(haetae_freeze_poly_2q, value)
haetae_reduce_vector_2q(value) = map(haetae_centered_poly_2q, value)

function haetae_decompose_poly(value)
    high = haetae_poly()
    low = haetae_poly()
    for index in 1:HAETAE_N
        high[index], low[index] = haetae_decompose_z1(value[index])
    end
    return high, low
end

haetae_compose_poly(high, low) = high .* 256 .+ low
haetae_lsb(value) = Int64[coefficient & 1 for coefficient in value]

function haetae_from_crt(mod_q, mod_2=nothing)
    output = haetae_poly()
    for index in 1:HAETAE_N
        parity = isnothing(mod_2) ? 0 : mod_2[index]
        output[index] = mod_q[index] + (HAETAE_Q & -(xor(mod_q[index], parity) & 1))
    end
    return output
end

function haetae_stream(bits::Int, seed, nonce::Int)
    instance = Shake(bits)
    update!(instance, seed)
    update!(instance, haetae_le16(nonce))
    finalize!(instance)
    return instance
end

function haetae_sample_uniform(seed, nonce::Int)
    output = haetae_poly()
    generator = haetae_stream(128, seed, nonce)
    count = 0
    bytes = Vector{UInt8}(undef, 168)
    while count < HAETAE_N
        squeeze!(generator, bytes)
        for offset in 1:2:length(bytes)-1
            candidate = Int(bytes[offset]) | (Int(bytes[offset + 1]) << 8)
            if candidate < HAETAE_Q
                count += 1
                output[count] = candidate
                count == HAETAE_N && break
            end
        end
    end
    return output
end

function haetae_sample_eta(seed, nonce::Int)
    output = haetae_poly()
    generator = haetae_stream(256, seed, nonce)
    count = 0
    bytes = Vector{UInt8}(undef, 136)
    while count < HAETAE_N
        squeeze!(generator, bytes)
        for byte in bytes
            byte >= 243 && continue
            value = Int(byte)
            for _ in 1:5
                count += 1
                remainder = value % 3
                output[count] = remainder == 2 ? -1 : remainder
                value ÷= 3
                count == HAETAE_N && break
            end
            count == HAETAE_N && break
        end
    end
    return output
end

function haetae_expand_matrix(parameters, seed, include_first_column::Bool)
    matrix = Vector{Vector{Vector{Int64}}}(undef, parameters.k)
    for row in 0:parameters.k-1
        result = haetae_vector(include_first_column ? parameters.l : parameters.m)
        for column in 0:parameters.m-1
            result[column + (include_first_column ? 2 : 1)] =
                haetae_sample_uniform(seed, (row << 8) + column)
        end
        matrix[row + 1] = result
    end
    return matrix
end

function haetae_expand_public_vector(parameters, seed)
    nonce = (parameters.k << 8) + parameters.m
    return [haetae_sample_uniform(seed, nonce + index) for index in 0:parameters.k-1]
end

function haetae_expand_secret(parameters, seed, nonce::Int)
    current = nonce
    first = Vector{Vector{Int64}}()
    second = Vector{Vector{Int64}}()
    for _ in 1:parameters.m
        push!(first, haetae_sample_eta(seed, current)); current += 1
    end
    for _ in 1:parameters.k
        push!(second, haetae_sample_eta(seed, current)); current += 1
    end
    return first, second
end

haetae_multiply_matrix(matrix, input) = [haetae_pointwise_accumulate(row, input) for row in matrix]

function haetae_public_key_decompose(value)
    high = haetae_vector(length(value))
    low = haetae_vector(length(value))
    for item in eachindex(value), index in 1:HAETAE_N
        high[item][index], low[item][index] = haetae_decompose_public_key(value[item][index])
    end
    return high, low
end

function haetae_hint_high_bits(parameters, value)
    return [Int64[haetae_decompose_hint(coefficient, parameters.alpha_hint,
        parameters.log_alpha_hint) for coefficient in item] for item in value]
end

function haetae_challenge(parameters, packed, digest)
    output = haetae_poly()
    generator = Shake(256)
    update!(generator, packed)
    update!(generator, digest)
    finalize!(generator)
    if parameters.mode == 5
        bytes = squeeze!(generator, Vector{UInt8}(undef, HAETAE_SEED_BYTES))
        weight = sum(count_ones, bytes)
        mask = weight < 128 ? 0x00 : weight > 128 ? 0xff : (-(Int(bytes[1]) & 1) & 0xff)
        for index in 0:HAETAE_N-1
            output[index + 1] = (xor(Int(bytes[(index >> 3) + 1]), mask) >> (index & 7)) & 1
        end
        return output
    end
    random = squeeze!(generator, Vector{UInt8}(undef, 136))
    position = 1
    for index in HAETAE_N-parameters.tau:HAETAE_N-1
        selected = 0
        while true
            if position > length(random)
                squeeze!(generator, random)
                position = 1
            end
            selected = Int(random[position]); position += 1
            selected <= index && break
        end
        output[index + 1] = output[selected + 1]
        output[selected + 1] = 1
    end
    return output
end

function haetae_recover_crt_vector(mod_q, first_mod_2)
    return [haetae_from_crt(mod_q[index], index == 1 ? first_mod_2 : nothing)
        for index in eachindex(mod_q)]
end

haetae_multiply_by_alpha(parameters, value) = [item .* parameters.alpha_hint for item in value]
haetae_divide_vector_2(value) = [item .>> 1 for item in value]

function haetae_normalize_hint(parameters, value)
    modulus = (HAETAE_DOUBLE_Q - 2) ÷ parameters.alpha_hint
    return [Int64[mod(coefficient, modulus) for coefficient in item] for item in value]
end

haetae_multiply_poly_vector(value, factor) = [haetae_pointwise(item, factor) for item in value]
haetae_from_montgomery(value) = [Int64[haetae_montgomery_reduce(c * 4_214) for c in item] for item in value]
haetae_double_negate_montgomery(value) = [Int64[haetae_montgomery_reduce(c * -28_642) for c in item] for item in value]
haetae_conditional_add_q(value) = [Int64[c < 0 ? c + HAETAE_Q : c for c in item] for item in value]

const HAETAE_DISTRIBUTIONS = Dict(
    2 => (
        hint=[(0,382),(382,246),(628,65),(693,7),(700,1),(701,1),(702,1),(703,1),(704,1),(705,1),(706,7),(713,65),(778,246)],
        hint_cut=6, hint_offset=239,
        high=[(0,1),(1,1),(2,1),(3,5),(8,58),(66,246),(312,398),(710,247),(957,59),(1016,5),(1021,1),(1022,1),(1023,1)],
        high_offset=6),
    3 => (
        hint=[(0,266),(266,213),(479,113),(592,40),(632,9),(641,1),(642,1),(643,1),(644,1),(645,1),(646,1),(647,1),(648,1),(649,9),(658,40),(698,113),(811,213)],
        hint_cut=8, hint_offset=235,
        high=[(0,1),(1,1),(2,1),(3,1),(4,8),(12,37),(49,112),(161,216),(377,269),(646,216),(862,112),(974,38),(1012,8),(1020,1),(1021,1),(1022,1),(1023,1)],
        high_offset=8),
    5 => (
        hint=[(0,122),(122,117),(239,102),(341,81),(422,59),(481,39),(520,23),(543,13),(556,6),(562,3),(565,1),(566,1),(567,1),(568,1),(569,1),(570,1),(571,1),(572,1),(573,1),(574,1),(575,1),(576,1),(577,1),(578,1),(579,3),(582,6),(588,13),(601,24),(625,39),(664,59),(723,81),(804,102),(906,118)],
        hint_cut=16, hint_offset=471,
        high=[(0,1),(1,1),(2,1),(3,1),(4,2),(6,13),(19,48),(67,118),(185,204),(389,245),(634,204),(838,118),(956,48),(1004,14),(1018,2),(1020,1),(1021,1),(1022,1),(1023,1)],
        high_offset=9),
)

function haetae_flatten(value)
    output = Vector{Int64}(undef, length(value) * HAETAE_N)
    for index in eachindex(value)
        copyto!(output, (index - 1) * HAETAE_N + 1, value[index], 1, HAETAE_N)
    end
    return output
end

function haetae_encode_symbols(symbols, distribution)
    state = UInt32(1 << 23)
    emitted = UInt8[]
    for index in length(symbols):-1:1
        symbol = Int(symbols[index]) + 1
        1 <= symbol <= length(distribution) || return nothing
        start, frequency = distribution[symbol]
        maximum = UInt64(((1 << 23) >> 10) << 8) * UInt64(frequency)
        while UInt64(state) >= maximum
            push!(emitted, UInt8(state & 0xff))
            state >>= 8
        end
        state = UInt32((UInt64(state) ÷ UInt64(frequency)) * 1024 +
            (UInt64(state) % UInt64(frequency)) + UInt64(start))
    end
    output = Vector{UInt8}(undef, 4 + length(emitted))
    for i in 0:3
        output[i + 1] = UInt8((state >> (8i)) & 0xff)
    end
    for index in eachindex(emitted)
        output[4 + index] = emitted[end - index + 1]
    end
    return output
end

function haetae_decode_symbols(input, count::Int, distribution)
    length(input) >= 4 || return nothing
    state = UInt32(input[1]) | (UInt32(input[2]) << 8) |
        (UInt32(input[3]) << 16) | (UInt32(input[4]) << 24)
    UInt32(1 << 23) <= state < UInt32(1 << 31) || return nothing
    offset = 5
    output = Vector{Int64}(undef, count)
    for index in 1:count
        cumulative = Int(state & 1023)
        symbol = 0
        for candidate in eachindex(distribution)
            start, frequency = distribution[candidate]
            if start <= cumulative < start + frequency
                symbol = candidate
                state = UInt32(UInt64(frequency) * UInt64(state >> 10) +
                    UInt64(cumulative - start))
                break
            end
        end
        symbol == 0 && return nothing
        output[index] = symbol - 1
        while state < UInt32(1 << 23) && offset <= length(input)
            state = (state << 8) | UInt32(input[offset])
            offset += 1
        end
    end
    return state == UInt32(1 << 23) && offset == length(input) + 1 ? output : nothing
end

function haetae_reshape(flat, count::Int)
    output = haetae_vector(count)
    for index in 1:count
        copyto!(output[index], 1, flat, (index - 1) * HAETAE_N + 1, HAETAE_N)
    end
    return output
end

function haetae_encode_hint(parameters, hint)
    specification = HAETAE_DISTRIBUTIONS[parameters.mode]
    values = haetae_flatten(hint)
    for index in eachindex(values)
        value = values[index]
        specification.hint_cut < value <= specification.hint_cut + specification.hint_offset && return nothing
        values[index] = value > specification.hint_cut + specification.hint_offset ?
            value - specification.hint_offset : value
    end
    return haetae_encode_symbols(values, specification.hint)
end

function haetae_decode_hint(parameters, input)
    specification = HAETAE_DISTRIBUTIONS[parameters.mode]
    values = haetae_decode_symbols(input, parameters.k * HAETAE_N, specification.hint)
    isnothing(values) && return nothing
    for index in eachindex(values)
        values[index] > specification.hint_cut && (values[index] += specification.hint_offset)
    end
    return haetae_reshape(values, parameters.k)
end

function haetae_encode_high(parameters, high)
    specification = HAETAE_DISTRIBUTIONS[parameters.mode]
    values = haetae_flatten(high)
    values .+= specification.high_offset
    return haetae_encode_symbols(values, specification.high)
end

function haetae_decode_high(parameters, input)
    specification = HAETAE_DISTRIBUTIONS[parameters.mode]
    values = haetae_decode_symbols(input, parameters.l * HAETAE_N, specification.high)
    isnothing(values) && return nothing
    values .-= specification.high_offset
    return haetae_reshape(values, parameters.l)
end

function haetae_pack_bits(values, bits::Int)
    output = zeros(UInt8, HAETAE_N * bits ÷ 8)
    accumulator = UInt64(0)
    available = 0
    offset = 1
    for coefficient in values
        accumulator |= UInt64(coefficient) << available
        available += bits
        while available >= 8
            output[offset] = UInt8(accumulator & 0xff); offset += 1
            accumulator >>= 8
            available -= 8
        end
    end
    return output
end

function haetae_unpack_bits(bytes, bits::Int)
    output = haetae_poly()
    accumulator = UInt64(0)
    available = 0
    offset = 1
    mask = (UInt64(1) << bits) - 1
    for index in 1:HAETAE_N
        while available < bits
            accumulator |= UInt64(bytes[offset]) << available
            offset += 1
            available += 8
        end
        output[index] = Int64(accumulator & mask)
        accumulator >>= bits
        available -= bits
    end
    return output
end

haetae_pack_q(p, value) = haetae_pack_bits(value, p.mode == 5 ? 16 : 15)
haetae_unpack_q(p, value) = haetae_unpack_bits(value, p.mode == 5 ? 16 : 15)
haetae_pack_eta(value) = haetae_pack_bits(1 .- value, 2)
haetae_unpack_eta(value) = 1 .- haetae_unpack_bits(value, 2)
haetae_pack_double_eta(value) = haetae_pack_bits(2 .- value, 3)
haetae_unpack_double_eta(value) = 2 .- haetae_unpack_bits(value, 3)

function haetae_pack_public_key(parameters, seed, public_vector)
    return vcat(seed, (haetae_pack_q(parameters, item) for item in public_vector)...)
end

function haetae_unpack_public_key(parameters, public_key)
    seed = public_key[1:HAETAE_SEED_BYTES]
    offset = HAETAE_SEED_BYTES + 1
    packed = haetae_vector(parameters.k)
    for row in 1:parameters.k
        packed[row] = haetae_unpack_q(parameters,
            public_key[offset:offset + parameters.poly_q_bytes - 1])
        offset += parameters.poly_q_bytes
    end
    matrix = haetae_expand_matrix(parameters, seed, true)
    for row in 1:parameters.k, column in 2:parameters.l
        matrix[row][column] .*= 2
    end
    if parameters.mode == 5
        for row in 1:parameters.k
            matrix[row][1] = packed[row]
        end
    else
        public_random = haetae_expand_public_vector(parameters, seed)
        first_column = haetae_transform_vector(
            haeta_double_vector(haetae_subtract_vectors(public_random, haeta_double_vector(packed))))
        for row in 1:parameters.k
            matrix[row][1] = first_column[row]
        end
    end
    return matrix
end

function haetae_pack_secret_key(parameters, public_key, first, second, key)
    first_packed = [haetae_pack_eta(item) for item in first]
    second_packed = parameters.mode == 5 ? [haetae_pack_eta(item) for item in second] :
        [haetae_pack_double_eta(item) for item in second]
    return vcat(public_key, first_packed..., second_packed..., key)
end

function haetae_unpack_secret_key(parameters, secret_key)
    public_key = secret_key[1:parameters.public_key_bytes]
    offset = parameters.public_key_bytes + 1
    first = haetae_vector(parameters.m)
    for index in 1:parameters.m
        first[index] = haetae_unpack_eta(secret_key[offset:offset + 63]); offset += 64
    end
    second = haetae_vector(parameters.k)
    for index in 1:parameters.k
        bytes = parameters.poly_2eta_bytes
        packed = secret_key[offset:offset + bytes - 1]
        second[index] = parameters.mode == 5 ? haetae_unpack_eta(packed) : haetae_unpack_double_eta(packed)
        offset += bytes
    end
    return haetae_unpack_public_key(parameters, public_key), first, second,
        copy(secret_key[offset:offset + HAETAE_SEED_BYTES - 1])
end

function haetae_pack_high(parameters, value)
    if parameters.mode != 5
        output = UInt8[]
        for item in value
            packed = zeros(UInt8, HAETAE_N * 9 ÷ 8)
            packed[1:HAETAE_N] .= UInt8.(item)
            append!(output, packed)
        end
        return output
    end
    return vcat((haetae_pack_bits(item, 9) for item in value)...)
end

haetae_pack_lsb(value) = haetae_pack_bits(value, 1)

function haetae_pack_signature(parameters, challenge, low, high, hint)
    encoded_high = haetae_encode_high(parameters, high)
    encoded_hint = haetae_encode_hint(parameters, hint)
    (isnothing(encoded_high) || isnothing(encoded_hint)) && return nothing
    parameters.base_encoded_high_z1 <= length(encoded_high) <= parameters.base_encoded_high_z1 + 255 || return nothing
    parameters.base_encoded_hint <= length(encoded_hint) <= parameters.base_encoded_hint + 255 || return nothing
    used = HAETAE_SEED_BYTES + parameters.l * HAETAE_N + 2 + length(encoded_high) + length(encoded_hint)
    used <= parameters.signature_bytes || return nothing
    signature = zeros(UInt8, parameters.signature_bytes)
    for index in 0:HAETAE_N-1
        signature[(index >> 3) + 1] |= UInt8(challenge[index + 1] << (index & 7))
    end
    offset = HAETAE_SEED_BYTES + 1
    for item in low
        for index in 1:HAETAE_N
            signature[offset + index - 1] = reinterpret(UInt8, Int8(item[index]))
        end
        offset += HAETAE_N
    end
    signature[offset] = length(encoded_high) - parameters.base_encoded_high_z1; offset += 1
    signature[offset] = length(encoded_hint) - parameters.base_encoded_hint; offset += 1
    copyto!(signature, offset, encoded_high, 1, length(encoded_high)); offset += length(encoded_high)
    copyto!(signature, offset, encoded_hint, 1, length(encoded_hint))
    return signature
end

function haetae_unpack_signature(parameters, signature)
    length(signature) == parameters.signature_bytes || return nothing
    challenge = haetae_poly()
    for index in 0:HAETAE_N-1
        challenge[index + 1] = (signature[(index >> 3) + 1] >> (index & 7)) & 1
    end
    offset = HAETAE_SEED_BYTES + 1
    low = haetae_vector(parameters.l)
    for item in low
        for index in 1:HAETAE_N
            item[index] = Int(reinterpret(Int8, signature[offset + index - 1]))
        end
        offset += HAETAE_N
    end
    high_length = Int(signature[offset]) + parameters.base_encoded_high_z1; offset += 1
    hint_length = Int(signature[offset]) + parameters.base_encoded_hint; offset += 1
    offset + high_length + hint_length - 1 <= length(signature) || return nothing
    high = haetae_decode_high(parameters, signature[offset:offset + high_length - 1]); offset += high_length
    hint = haetae_decode_hint(parameters, signature[offset:offset + hint_length - 1]); offset += hint_length
    (isnothing(high) || isnothing(hint)) && return nothing
    all(iszero, @view signature[offset:end]) || return nothing
    return challenge, low, high, hint
end

mutable struct HaetaeFixed96
    low::UInt64
    high::UInt64
end

const HAETAE_U48 = UInt64((1 << 48) - 1)

function haetae_normalize_fixed(value::HaetaeFixed96)
    return HaetaeFixed96(value.low & HAETAE_U48, value.high + (value.low >> 48))
end

function haetae_mul64(left::UInt64, right::UInt64)
    product = UInt128(left) * UInt128(right)
    return UInt64(product & typemax(UInt64)), UInt64(product >> 64)
end

function haetae_mul48(left::UInt64, right::UInt64)
    low64, high64 = haetae_mul64(left, right)
    return low64 & HAETAE_U48, xor(high64 << 16, low64 >> 48)
end

function haetae_conditional_negate(value::HaetaeFixed96, sign::Int)
    mask = sign == 0 ? UInt64(0) : typemax(UInt64)
    return haetae_normalize_fixed(HaetaeFixed96(
        (xor(value.low, mask & HAETAE_U48) + UInt64(sign)),
        xor(value.high, mask),
    ))
end

haetae_add_fixed(left::HaetaeFixed96, right::HaetaeFixed96) =
    HaetaeFixed96(left.low + right.low, left.high + right.high)

function haetae_multiply_fixed(left::HaetaeFixed96, right::HaetaeFixed96)
    low, high = haetae_mul48(left.low, right.low)
    low = high + (((low >> 47) + 1) >> 1)
    cross_low, cross_high = haetae_mul48(left.low, right.high)
    low += cross_low
    high = cross_high
    cross_low, cross_high = haetae_mul48(left.high, right.low)
    low += cross_low
    high += cross_high
    low += UInt64(1) << 27
    low >>= 28
    low += (high << 20) & HAETAE_U48
    high >>= 28
    upper_low, upper_high = haetae_mul64(left.high, right.high)
    low += (upper_low << 20) & HAETAE_U48
    high += (upper_low >> 28) + (upper_high << 36)
    return haetae_normalize_fixed(HaetaeFixed96(low, high))
end

function haetae_square_fixed(value::HaetaeFixed96)
    low, high = haetae_mul48(value.low, value.low)
    low = (low >> 48) + high
    cross_low, cross_high = haetae_mul48(value.low, value.high)
    low += cross_low << 1
    high = cross_high << 1
    low >>= 28
    low += (high << 20) & HAETAE_U48
    high >>= 28
    upper_low, upper_high = haetae_mul64(value.high, value.high)
    low += (upper_low << 20) & HAETAE_U48
    high += (upper_low >> 28) + (upper_high << 36)
    return haetae_normalize_fixed(HaetaeFixed96(low, high))
end

haetae_subtract_fixed(left, right) = haetae_add_fixed(left, haetae_conditional_negate(right, 1))

function haetae_subtract_three_halves(value)
    result = haetae_conditional_negate(value, 1)
    result.high += UInt64(3) << 27
    return haetae_normalize_fixed(result)
end

function haetae_multiply_unsigned_signed(value, signed)
    sign = Int((signed.high >> 63) & 1)
    absolute = haetae_conditional_negate(signed, sign)
    return haetae_conditional_negate(haetae_multiply_fixed(absolute, value), sign)
end

const HAETAE_HYPERBALL = Dict(
    2 => (cube=HaetaeFixed96(0x770077e2e41a, 0x1162),
        three=HaetaeFixed96(0x693861ad937b, 0x9caa56),
        scalar=UInt64(2_643_021_496_320)),
    3 => (cube=HaetaeFixed96(0x1a2935cfae68, 0x978),
        three=HaetaeFixed96(0x7ad215218533, 0x7ff1c9),
        scalar=UInt64(4_916_390_789_120)),
    5 => (cube=HaetaeFixed96(0x700ff3e8890d, 0x702),
        three=HaetaeFixed96(0x5768588eed31, 0x73bd40),
        scalar=UInt64(5_997_831_421_952)),
)

function haetae_inverse_sqrt(mode::Int, half::HaetaeFixed96)
    specification = HAETAE_HYPERBALL[mode]
    result = haetae_subtract_fixed(specification.three,
        haetae_multiply_fixed(half, specification.cube))
    for _ in 1:6
        factor = haetae_subtract_three_halves(
            haetae_multiply_fixed(half, haetae_square_fixed(result)))
        result = haetae_multiply_unsigned_signed(result, factor)
    end
    return result
end

function haetae_multiply_fixed_high(value::HaetaeFixed96, scalar::UInt64)
    low, high = haetae_mul48(value.low, scalar)
    upper_low, upper_high = haetae_mul48(value.high, scalar)
    high += upper_low
    low += UInt64(1) << 27
    low >>= 28
    low += (high << 20) & HAETAE_U48
    high >>= 28
    high += upper_high << 20
    return haetae_normalize_fixed(HaetaeFixed96(low, high))
end

function haetae_multiply_round13(value::UInt64, factor::HaetaeFixed96, sign::Int)
    fixed = HaetaeFixed96((value & 0xffff_ffff) << 16, value >> 32)
    product = haetae_multiply_fixed(fixed, factor)
    rounded = Int64((product.high + (UInt64(1) << 14)) >> 15)
    return sign == 0 ? rounded : -rounded
end

function haetae_signed_mul_high48(left::Int64, right::UInt64)
    return Int64((Int128(left) * Int128(right) + ((Int128(1) << 48) - 1)) >> 48)
end

function haetae_approx_exp(value::UInt64)
    signed = reinterpret(Int64, value)
    result = Int64(55_868_746)
    result = haetae_signed_mul_high48(result, value) - 743_564_434
    result = haetae_signed_mul_high48(result, value) + 6_953_427_278
    result = haetae_signed_mul_high48(result, value) - 55_833_338_892
    result = haetae_signed_mul_high48(result, value) + 390_932_311_155
    result = haetae_signed_mul_high48(result, value) - 2_345_623_661_771
    result = haetae_signed_mul_high48(result, value) + 11_728_123_872_951
    result = haetae_signed_mul_high48(result, value) - 46_912_496_106_200
    result = haetae_signed_mul_high48(result, value) + 140_737_488_354_861
    result = haetae_signed_mul_high48(result, value) - 281_474_976_710_650
    result = haetae_signed_mul_high48(result, value) + 281_474_976_710_657
    return reinterpret(UInt64, result)
end

function haetae_sample_gaussian83(lower::UInt32, middle::UInt32, upper::UInt32)
    start = 1
    stop = length(HAETAE_CDT_UPPER) + 1
    while start < stop
        index = (start + stop) >> 1
        threshold_upper = HAETAE_CDT_UPPER[index]
        threshold_middle = HAETAE_CDT_MIDDLE[index]
        is_lower = threshold_upper < upper || (threshold_upper == upper &&
            (threshold_middle < middle || (threshold_middle == middle && HAETAE_CDT_LOWER[index] < lower)))
        is_lower ? (start = index + 1) : (stop = index)
    end
    return start - 1
end

haetae_le32(bytes, offset) = UInt32(bytes[offset]) | (UInt32(bytes[offset+1]) << 8) |
    (UInt32(bytes[offset+2]) << 16) | (UInt32(bytes[offset+3]) << 24)
haetae_le24(bytes, offset) = UInt64(bytes[offset]) | (UInt64(bytes[offset+1]) << 8) |
    (UInt64(bytes[offset+2]) << 16)
haetae_le48(bytes, offset) = UInt64(haetae_le32(bytes, offset)) |
    (UInt64(bytes[offset+4]) | (UInt64(bytes[offset+5]) << 8)) << 32

function haetae_gaussian_candidate(random)
    gaussian_low = haetae_le32(random, 1)
    gaussian_middle = haetae_le32(random, 5)
    gaussian_high = UInt32(random[9]) | (UInt32(random[10]) << 8) | (UInt32(random[11] & 7) << 16)
    rejection = haetae_le48(random, 12)
    x = haetae_sample_gaussian83(gaussian_low, gaussian_middle, gaussian_high)
    value = HaetaeFixed96(haetae_le48(random, 18),
        haetae_le24(random, 24) | (UInt64(x) << 24))
    sample = ((value.low >> 15) + 1) >> 1
    sample += value.high << 32
    square = haetae_square_fixed(value)
    exponent = square.high - (UInt64(x * x) << 20)
    exponent = (exponent << 20) | (square.low >> 28)
    exponent = (exponent + 1) >> 1
    even_rejection = xor(rejection, rejection & 1)
    accepted = even_rejection < haetae_approx_exp(exponent) && (sample != 0 || (rejection & 1) != 0)
    return accepted, sample, square
end

function haetae_sample_gaussian_n(seed, nonce::Int, sample_count::Int, square_sum::HaetaeFixed96)
    stream = haetae_stream(256, seed, nonce)
    signs = squeeze!(stream, Vector{UInt8}(undef, sample_count ÷ 8))
    samples = UInt64[]
    random = Vector{UInt8}(undef, 26)
    while length(samples) < sample_count
        squeeze!(stream, random)
        accepted, sample, square = haetae_gaussian_candidate(random)
        accepted || continue
        square_sum.low += square.low
        square_sum.high += square.high
        push!(samples, length(samples) < HAETAE_N || sample_count == HAETAE_N ? sample : 0)
    end
    normalized = haetae_normalize_fixed(square_sum)
    square_sum.low, square_sum.high = normalized.low, normalized.high
    return samples[1:HAETAE_N], signs
end

function haetae_sample_hyperball(parameters, seed, nonce::Int)
    dimension = parameters.l + parameters.k
    specification = HAETAE_HYPERBALL[parameters.mode]
    limit = parameters.b0_squared * BigInt(HAETAE_FIXED_SCALE^2)
    next_nonce = nonce
    while true
        square_sum = HaetaeFixed96(0, 0)
        samples = Vector{Vector{UInt64}}()
        signs = Vector{Vector{UInt8}}()
        for index in 1:dimension
            sampled, sampled_signs = haetae_sample_gaussian_n(
                seed, next_nonce, index <= 2 ? HAETAE_N + 1 : HAETAE_N, square_sum)
            next_nonce += 1
            push!(samples, sampled); push!(signs, sampled_signs)
        end
        square_sum.low = (square_sum.low + 1) >> 1
        square_sum.low += (square_sum.high & 1) << 47
        square_sum.high >>= 1
        square_sum = haetae_normalize_fixed(square_sum)
        inverse = haetae_inverse_sqrt(parameters.mode, square_sum)
        factor = haetae_multiply_fixed_high(inverse, specification.scalar)
        fixed = Vector{Vector{Int64}}()
        for item_index in eachindex(samples)
            item = samples[item_index]
            item_signs = signs[item_index]
            output = Vector{Int64}(undef, length(item))
            for coefficient in eachindex(item)
                sign = (item_signs[((coefficient - 1) >> 3) + 1] >> ((coefficient - 1) & 7)) & 1
                output[coefficient] = haetae_multiply_round13(item[coefficient], factor, sign)
            end
            push!(fixed, output)
        end
        first = fixed[1:parameters.l]
        second = fixed[parameters.l+1:end]
        if haetae_squared_norm(first, second) <= limit
            bit = shake(256, 1, seed, haetae_le16(next_nonce))[1]
            return first, second, Int(bit), next_nonce
        end
    end
end

haetae_round_fixed(value) = [Int64[(coefficient + HAETAE_FIXED_SCALE ÷ 2) >> 13 for coefficient in item] for item in value]

function haetae_add_fixed_integer(fixed, integer, sign::Int)
    factor = sign == 0 ? HAETAE_FIXED_SCALE : -HAETAE_FIXED_SCALE
    return [fixed[i] .+ factor .* integer[i] for i in eachindex(fixed)]
end

haetae_double_minus(value, subtrahend) = [value[i] .* 2 .- subtrahend[i] for i in eachindex(value)]

function haetae_generate_keypair(parameters; seed=random_bytes(HAETAE_SEED_BYTES))
    expanded = shake(256, 2 * HAETAE_SEED_BYTES + HAETAE_CRH_BYTES, seed)
    matrix_seed = expanded[1:HAETAE_SEED_BYTES]
    secret_seed = expanded[HAETAE_SEED_BYTES+1:HAETAE_SEED_BYTES+HAETAE_CRH_BYTES]
    key = expanded[HAETAE_SEED_BYTES+HAETAE_CRH_BYTES+1:end]
    matrix = haetae_expand_matrix(parameters, matrix_seed, false)
    public_random = parameters.mode == 5 ? nothing : haetae_expand_public_vector(parameters, matrix_seed)
    counter = 0
    while true
        first, second = haetae_expand_secret(parameters, secret_seed, counter)
        counter += parameters.m + parameters.k
        if parameters.mode == 5 && haetae_secret_singular_value(parameters, first, second) >
                parameters.singular_value_limit
            continue
        end
        transformed_first = haetae_transform_vector(first)
        public_vector = haetae_multiply_matrix(matrix, transformed_first)
        if parameters.mode == 5
            transformed_second = haetae_transform_vector(second)
            public_vector = haetae_conditional_add_q(haetae_double_negate_montgomery(
                haetae_add_vectors(haetae_from_montgomery(public_vector), transformed_second)))
        else
            public_vector = haetae_freeze_vector(haetae_add_vectors(
                haetae_add_vectors(haetae_inverse_transform_vector(public_vector), second), public_random))
            high, low = haetae_public_key_decompose(public_vector)
            public_vector = high
            second = [second[i] .- low[i] for i in eachindex(second)]
            haetae_secret_singular_value(parameters, first, second) >
                parameters.singular_value_limit && continue
        end
        public_key = haetae_pack_public_key(parameters, matrix_seed, public_vector)
        secret_key = haetae_pack_secret_key(parameters, public_key, first, second, key)
        return public_key, secret_key
    end
end

function haetae_context(context)
    length(context) <= 255 || return nothing
    return vcat(UInt8(length(context)), context)
end

function haetae_decompose_vector(value)
    high = Vector{Vector{Int64}}()
    low = Vector{Vector{Int64}}()
    for item in value
        h, l = haetae_decompose_poly(item)
        push!(high, h); push!(low, l)
    end
    return high, low
end

function haetae_add_negative_hint(parameters, value)
    modulus = (HAETAE_DOUBLE_Q - 2) ÷ parameters.alpha_hint
    return [Int64[coefficient < 0 ? coefficient + modulus : coefficient for coefficient in item]
        for item in value]
end

function haetae_sign(parameters, message, secret_key, context=UInt8[];
                     randomness=random_bytes(HAETAE_SEED_BYTES))
    length(randomness) == HAETAE_SEED_BYTES || throw(ArgumentError("randomness must be 32 bytes"))
    prefix = haetae_context(context)
    isnothing(prefix) && throw(ArgumentError("context cannot exceed 255 bytes"))
    matrix, secret_first, secret_second, key = haetae_unpack_secret_key(parameters, secret_key)
    public_key = secret_key[1:parameters.public_key_bytes]
    digest = shake(256, HAETAE_CRH_BYTES, public_key, prefix, message)
    seed = shake(256, HAETAE_CRH_BYTES, key, randomness, digest)
    transformed_first = haetae_transform_vector(secret_first)
    transformed_second = haetae_transform_vector(secret_second)
    b0_limit = parameters.b0_squared * BigInt(HAETAE_FIXED_SCALE^2)
    b1_limit = parameters.b1_squared * BigInt(HAETAE_FIXED_SCALE^2)
    counter = 0
    while true
        sampled_first, sampled_second, sampled_bit, next_nonce =
            haetae_sample_hyperball(parameters, seed, counter)
        counter = next_nonce
        rounded_first = haetae_round_fixed(sampled_first)
        rounded_second = haetae_round_fixed(sampled_second)
        rounded_first_zero = rounded_first[1]
        product = haetae_inverse_transform_vector(haetae_multiply_matrix(
            matrix, haetae_transform_vector(rounded_first)))
        product = haetae_add_vectors(product, haeta_double_vector(rounded_second))
        product = haetae_freeze_vector_2q(haetae_recover_crt_vector(product, rounded_first_zero))
        high_bits = haetae_hint_high_bits(parameters, product)
        parity = haetae_lsb(rounded_first_zero)
        packed = vcat(haetae_pack_high(parameters, high_bits), haetae_pack_lsb(parity))
        challenge_value = haetae_challenge(parameters, packed, digest[1:HAETAE_SEED_BYTES])
        transformed_challenge = haetae_transform(challenge_value)
        challenge_first = Vector{Vector{Int64}}([challenge_value])
        for index in 2:parameters.l
            multiplied = haetae_multiply_poly_vector([transformed_first[index - 1]], transformed_challenge)
            push!(challenge_first, haetae_inverse_transform_vector(multiplied)[1])
        end
        challenge_second = haetae_inverse_transform_vector(
            haetae_multiply_poly_vector(transformed_second, transformed_challenge))
        sign_bit = sampled_bit & 1
        first = haetae_add_fixed_integer(sampled_first, challenge_first, sign_bit)
        second = haetae_add_fixed_integer(sampled_second, challenge_second, sign_bit)
        haetae_squared_norm(first, second) > b1_limit && continue
        doubled_first = haetae_double_minus(first, sampled_first)
        doubled_second = haetae_double_minus(second, sampled_second)
        if (sampled_bit & 2) != 0 && haetae_squared_norm(doubled_first, doubled_second) < b0_limit
            continue
        end
        rounded_z1 = haetae_round_fixed(first)
        rounded_z2 = haetae_round_fixed(second)
        hint_base = haetae_subtract_vectors(product, haeta_double_vector(rounded_z2))
        hint_base = haetae_hint_high_bits(parameters, haetae_freeze_vector_2q(hint_base))
        hint = haetae_add_negative_hint(parameters, haetae_subtract_vectors(high_bits, hint_base))
        high_z1, low_z1 = haetae_decompose_vector(rounded_z1)
        signature = haetae_pack_signature(parameters, challenge_value, low_z1, high_z1, hint)
        isnothing(signature) || return signature
    end
end

function haetae_verify(parameters, public_key, message, signature, context=UInt8[])
    length(public_key) == parameters.public_key_bytes || return false
    prefix = haetae_context(context)
    isnothing(prefix) && return false
    decoded = haetae_unpack_signature(parameters, signature)
    isnothing(decoded) && return false
    challenge_value, low_bits, high_bits, hint = decoded
    matrix = haetae_unpack_public_key(parameters, public_key)
    first = [haetae_compose_poly(high_bits[index], low_bits[index]) for index in eachindex(high_bits)]
    norm = haetae_squared_norm(first)
    parity = haetae_lsb(haetae_subtract(first[1], challenge_value))
    product = haetae_inverse_transform_vector(
        haetae_multiply_matrix(matrix, haetae_transform_vector(first)))
    product = haetae_freeze_vector_2q(haetae_recover_crt_vector(product, parity))
    recovered_high = haetae_hint_high_bits(parameters, product)
    recovered_high = haetae_normalize_hint(parameters, haetae_add_vectors(recovered_high, hint))
    second = haetae_subtract_vectors(haetae_multiply_by_alpha(parameters, recovered_high), product)
    second[1] = haetae_add(second[1], parity)
    second = haetae_divide_vector_2(haetae_reduce_vector_2q(second))
    norm += haetae_squared_norm(second)
    norm <= parameters.b2_squared || return false
    packed = vcat(haetae_pack_high(parameters, recovered_high), haetae_pack_lsb(parity))
    digest = shake(256, HAETAE_SEED_BYTES, public_key, prefix, message)
    expected = haetae_challenge(parameters, packed, digest)
    return expected == challenge_value
end
