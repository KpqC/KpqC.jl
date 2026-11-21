const KECCAK_ROUND_CONSTANTS = UInt64[
    0x0000000000000001, 0x0000000000008082, 0x800000000000808a,
    0x8000000080008000, 0x000000000000808b, 0x0000000080000001,
    0x8000000080008081, 0x8000000000008009, 0x000000000000008a,
    0x0000000000000088, 0x0000000080008009, 0x000000008000000a,
    0x000000008000808b, 0x800000000000008b, 0x8000000000008089,
    0x8000000000008003, 0x8000000000008002, 0x8000000000000080,
    0x000000000000800a, 0x800000008000000a, 0x8000000080008081,
    0x8000000000008080, 0x0000000080000001, 0x8000000080008008,
]

const KECCAK_ROTATIONS = UInt8[
     0,  1, 62, 28, 27,
    36, 44,  6, 55, 20,
     3, 10, 43, 25, 39,
    41, 45, 15, 21,  8,
    18,  2, 61, 56, 14,
]

function keccak_f1600!(state::Vector{UInt64})
    columns = Vector{UInt64}(undef, 5)
    deltas = Vector{UInt64}(undef, 5)
    transformed = Vector{UInt64}(undef, 25)

    for round_constant in KECCAK_ROUND_CONSTANTS
        for x in 0:4
            columns[x + 1] = xor((state[x + 5y + 1] for y in 0:4)...)
        end
        for x in 0:4
            deltas[x + 1] = columns[mod(x - 1, 5) + 1] ⊻
                            bitrotate(columns[mod(x + 1, 5) + 1], 1)
        end
        for y in 0:4, x in 0:4
            state[x + 5y + 1] ⊻= deltas[x + 1]
        end

        for y in 0:4, x in 0:4
            destination_x = y
            destination_y = mod(2x + 3y, 5)
            transformed[destination_x + 5destination_y + 1] =
                bitrotate(state[x + 5y + 1], Int(KECCAK_ROTATIONS[x + 5y + 1]))
        end

        for y in 0:4, x in 0:4
            state[x + 5y + 1] = transformed[x + 5y + 1] ⊻
                ((~transformed[mod(x + 1, 5) + 5y + 1]) &
                 transformed[mod(x + 2, 5) + 5y + 1])
        end
        state[1] ⊻= round_constant
    end
    return state
end

mutable struct Shake
    state::Vector{UInt64}
    rate::Int
    suffix::UInt8
    position::Int
    finalized::Bool
end

function Shake(security_bits::Int; suffix::UInt8 = 0x1f)
    rate = security_bits == 128 ? 168 : security_bits == 512 ? 72 : 136
    security_bits in (128, 192, 256, 512) ||
        throw(ArgumentError("unsupported SHAKE security level: $security_bits"))
    return Shake(zeros(UInt64, 25), rate, suffix, 0, false)
end

function xor_byte!(instance::Shake, position::Int, byte::UInt8)
    lane = position >> 3
    shift = (position & 7) << 3
    instance.state[lane + 1] ⊻= UInt64(byte) << shift
end

function update!(instance::Shake, input::AbstractVector{UInt8})
    instance.finalized && throw(ArgumentError("cannot absorb after finalization"))
    for byte in input
        xor_byte!(instance, instance.position, byte)
        instance.position += 1
        if instance.position == instance.rate
            keccak_f1600!(instance.state)
            instance.position = 0
        end
    end
    return instance
end

function update_byte!(instance::Shake, byte::Integer)
    update!(instance, UInt8[byte])
end

function finalize!(instance::Shake)
    instance.finalized && return instance
    xor_byte!(instance, instance.position, instance.suffix)
    xor_byte!(instance, instance.rate - 1, 0x80)
    keccak_f1600!(instance.state)
    instance.position = 0
    instance.finalized = true
    return instance
end

function squeeze!(instance::Shake, output::AbstractVector{UInt8})
    finalize!(instance)
    for index in eachindex(output)
        if instance.position == instance.rate
            keccak_f1600!(instance.state)
            instance.position = 0
        end
        lane = instance.position >> 3
        shift = (instance.position & 7) << 3
        output[index] = UInt8((instance.state[lane + 1] >> shift) & 0xff)
        instance.position += 1
    end
    return output
end

function shake(security_bits::Int, length::Int, inputs::AbstractVector{UInt8}...)
    instance = Shake(security_bits)
    foreach(input -> update!(instance, input), inputs)
    return squeeze!(instance, Vector{UInt8}(undef, length))
end

function sha3(security_bits::Int, inputs::AbstractVector{UInt8}...)
    security_bits in (256, 512) ||
        throw(ArgumentError("unsupported SHA-3 output size: $security_bits"))
    instance = Shake(security_bits; suffix = 0x06)
    foreach(input -> update!(instance, input), inputs)
    return squeeze!(instance, Vector{UInt8}(undef, security_bits ÷ 8))
end
