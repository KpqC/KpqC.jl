function random_bytes(length_value::Int)
    length_value >= 0 || throw(ArgumentError("random byte length must be nonnegative"))
    output = Vector{UInt8}(undef, length_value)
    rand!(Random.RandomDevice(), output)
    return output
end
