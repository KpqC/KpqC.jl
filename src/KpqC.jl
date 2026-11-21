module KpqC

using Random
import Base: show, showerror, sign

export Algorithm,
    SignatureAlgorithm,
    KEMAlgorithm,
    SignatureSizes,
    KEMSizes,
    KeyPair,
    EncapsulatedSecret,
    KpqCError,
    AIMer128f,
    AIMer128s,
    AIMer192f,
    AIMer192s,
    AIMer256f,
    AIMer256s,
    HAETAE2,
    HAETAE3,
    HAETAE5,
    NTRUPlus768,
    NTRUPlus864,
    NTRUPlus1152,
    SMAUGT128,
    SMAUGT192,
    SMAUGT256,
    TiMER,
    algorithm_id,
    signature_algorithms,
    kem_algorithms,
    signature_sizes,
    kem_sizes,
    generate_keypair,
    sign,
    verify,
    encapsulate,
    decapsulate

"""A supported KpqC digital-signature parameter set."""
@enum SignatureAlgorithm::UInt8 begin
    AIMer128f = 1
    AIMer128s
    AIMer192f
    AIMer192s
    AIMer256f
    AIMer256s
    HAETAE2
    HAETAE3
    HAETAE5
end

"""A supported KpqC key-encapsulation parameter set."""
@enum KEMAlgorithm::UInt8 begin
    NTRUPlus768 = 1
    NTRUPlus864
    NTRUPlus1152
    SMAUGT128
    SMAUGT192
    SMAUGT256
    TiMER
end

"""Union of the supported signature and key-encapsulation algorithms."""
const Algorithm = Union{SignatureAlgorithm, KEMAlgorithm}

"""Byte lengths used by a signature algorithm."""
struct SignatureSizes
    public_key::Int
    secret_key::Int
    signature::Int
end

"""Byte lengths used by a key-encapsulation algorithm."""
struct KEMSizes
    public_key::Int
    secret_key::Int
    ciphertext::Int
    shared_secret::Int
end

"""
    KeyPair

A generated public and secret key. The default display redacts the secret key.
"""
struct KeyPair
    public_key::Vector{UInt8}
    secret_key::Vector{UInt8}
end

"""
    EncapsulatedSecret

A KEM ciphertext and its shared secret. The default display redacts the secret.
"""
struct EncapsulatedSecret
    ciphertext::Vector{UInt8}
    shared_secret::Vector{UInt8}
end

"""Error returned when a KpqC operation rejects its input."""
struct KpqCError <: Exception
    operation::Symbol
    status::Int
end

include("internal/keccak.jl")
include("internal/entropy.jl")
include("internal/ntruplus_parameters.jl")
include("internal/ntruplus.jl")
include("internal/smaugt.jl")
include("internal/aimer.jl")
include("internal/haetae.jl")

showerror(io::IO, error::KpqCError) =
    print(io, "KpqC ", error.operation, " failed with status ", error.status)

show(io::IO, keys::KeyPair) = print(
    io,
    "KeyPair(public_key=<",
    length(keys.public_key),
    " bytes>, secret_key=<redacted>)",
)

show(io::IO, value::EncapsulatedSecret) = print(
    io,
    "EncapsulatedSecret(ciphertext=<",
    length(value.ciphertext),
    " bytes>, shared_secret=<redacted>)",
)

const SIGNATURE_METADATA = Dict(
    AIMer128f => (id = "aimer-128f", sizes = SignatureSizes(32, 48, 6944)),
    AIMer128s => (id = "aimer-128s", sizes = SignatureSizes(32, 48, 4704)),
    AIMer192f => (id = "aimer-192f", sizes = SignatureSizes(48, 72, 15408)),
    AIMer192s => (id = "aimer-192s", sizes = SignatureSizes(48, 72, 10320)),
    AIMer256f => (id = "aimer-256f", sizes = SignatureSizes(64, 96, 31360)),
    AIMer256s => (id = "aimer-256s", sizes = SignatureSizes(64, 96, 20224)),
    HAETAE2 => (id = "haetae-mode2", sizes = SignatureSizes(992, 1408, 1474)),
    HAETAE3 => (id = "haetae-mode3", sizes = SignatureSizes(1472, 2112, 2349)),
    HAETAE5 => (id = "haetae-mode5", sizes = SignatureSizes(2080, 2752, 2948)),
)

const KEM_METADATA = Dict(
    NTRUPlus768 => (id = "NTRU+768", sizes = KEMSizes(1152, 2336, 1152, 32)),
    NTRUPlus864 => (id = "NTRU+864", sizes = KEMSizes(1296, 2624, 1296, 32)),
    NTRUPlus1152 => (id = "NTRU+1152", sizes = KEMSizes(1728, 3488, 1728, 32)),
    SMAUGT128 => (id = "SMAUG-T128", sizes = KEMSizes(672, 832, 672, 32)),
    SMAUGT192 => (id = "SMAUG-T192", sizes = KEMSizes(1088, 1312, 992, 32)),
    SMAUGT256 => (id = "SMAUG-T256", sizes = KEMSizes(1440, 1728, 1376, 32)),
    TiMER => (id = "TiMER", sizes = KEMSizes(672, 832, 608, 32)),
)

"""Return every supported signature algorithm."""
signature_algorithms() = collect(instances(SignatureAlgorithm))

"""Return every supported key-encapsulation algorithm."""
kem_algorithms() = collect(instances(KEMAlgorithm))

"""Return the stable textual identifier for `algorithm`."""
algorithm_id(algorithm::SignatureAlgorithm) = SIGNATURE_METADATA[algorithm].id
algorithm_id(algorithm::KEMAlgorithm) = KEM_METADATA[algorithm].id

"""Return the public-key, secret-key, and signature sizes for `algorithm`."""
signature_sizes(algorithm::SignatureAlgorithm) = SIGNATURE_METADATA[algorithm].sizes

"""Return the public-key, secret-key, ciphertext, and shared-secret sizes."""
kem_sizes(algorithm::KEMAlgorithm) = KEM_METADATA[algorithm].sizes

bytes(value::AbstractVector{UInt8}) = value isa Vector{UInt8} ? value : collect(value)

function require_size(value::AbstractVector{UInt8}, expected::Int, label::String)
    length(value) == expected || throw(DimensionMismatch(
        "$label must contain $expected bytes, received $(length(value))",
    ))
end

function require_context(context::AbstractVector{UInt8})
    length(context) <= 255 || throw(ArgumentError("context must not exceed 255 bytes"))
end

"""
    generate_keypair(algorithm) -> KeyPair

Generate a public and secret key pair for a signature or KEM algorithm.
"""
function generate_keypair(algorithm::SignatureAlgorithm)
    if algorithm in (
        AIMer128f,
        AIMer128s,
        AIMer192f,
        AIMer192s,
        AIMer256f,
        AIMer256s,
    )
        key = algorithm == AIMer128f ? :aimer128f :
              algorithm == AIMer128s ? :aimer128s :
              algorithm == AIMer192f ? :aimer192f :
              algorithm == AIMer192s ? :aimer192s :
              algorithm == AIMer256f ? :aimer256f : :aimer256s
        generated = aimer_generate_keypair(AIMER_PARAMETERS[key])
        return KeyPair(generated.public_key, generated.secret_key)
    end
    parameters = HAETAE_PARAMETERS[
        algorithm == HAETAE2 ? 2 : algorithm == HAETAE3 ? 3 : 5
    ]
    public_key, secret_key = haetae_generate_keypair(parameters)
    return KeyPair(public_key, secret_key)
end

function generate_keypair(algorithm::KEMAlgorithm)
    if algorithm in (NTRUPlus768, NTRUPlus864, NTRUPlus1152)
        parameters = NTRUPLUS_PARAMETERS[
            algorithm == NTRUPlus768 ? 768 : algorithm == NTRUPlus864 ? 864 : 1152
        ]
        generated = ntruplus_generate_keypair(parameters)
        return KeyPair(generated.public_key, generated.secret_key)
    end
    parameters = SMAUGT_PARAMETERS[
        algorithm == SMAUGT128 ? 128 :
        algorithm == SMAUGT192 ? 192 :
        algorithm == SMAUGT256 ? 256 : :timer
    ]
    generated = smaugt_generate_keypair(parameters)
    return KeyPair(generated.public_key, generated.secret_key)
end

"""
    sign(algorithm, message, secret_key; context=UInt8[]) -> Vector{UInt8}

Create a detached signature. `context` may contain at most 255 bytes.
"""
function sign(
    algorithm::SignatureAlgorithm,
    message_value::AbstractVector{UInt8},
    secret_key_value::AbstractVector{UInt8};
    context::AbstractVector{UInt8} = UInt8[],
)
    sizes = signature_sizes(algorithm)
    require_size(secret_key_value, sizes.secret_key, "secret key")
    require_context(context)

    if algorithm in (
        AIMer128f,
        AIMer128s,
        AIMer192f,
        AIMer192s,
        AIMer256f,
        AIMer256s,
    )
        key = algorithm == AIMer128f ? :aimer128f :
              algorithm == AIMer128s ? :aimer128s :
              algorithm == AIMer192f ? :aimer192f :
              algorithm == AIMer192s ? :aimer192s :
              algorithm == AIMer256f ? :aimer256f : :aimer256s
        return aimer_sign(
            AIMER_PARAMETERS[key],
            bytes(message_value),
            bytes(secret_key_value),
            bytes(context),
        )
    end
    parameters = HAETAE_PARAMETERS[
        algorithm == HAETAE2 ? 2 : algorithm == HAETAE3 ? 3 : 5
    ]
    return haetae_sign(
        parameters,
        bytes(message_value),
        bytes(secret_key_value),
        bytes(context),
    )
end

"""
    verify(algorithm, message, signature, public_key; context=UInt8[]) -> Bool

Verify a detached signature. Invalid signatures return `false`.
"""
function verify(
    algorithm::SignatureAlgorithm,
    message_value::AbstractVector{UInt8},
    signature_value::AbstractVector{UInt8},
    public_key_value::AbstractVector{UInt8};
    context::AbstractVector{UInt8} = UInt8[],
)
    sizes = signature_sizes(algorithm)
    require_size(public_key_value, sizes.public_key, "public key")
    require_context(context)
    length(signature_value) == sizes.signature || return false

    if algorithm in (
        AIMer128f,
        AIMer128s,
        AIMer192f,
        AIMer192s,
        AIMer256f,
        AIMer256s,
    )
        key = algorithm == AIMer128f ? :aimer128f :
              algorithm == AIMer128s ? :aimer128s :
              algorithm == AIMer192f ? :aimer192f :
              algorithm == AIMer192s ? :aimer192s :
              algorithm == AIMer256f ? :aimer256f : :aimer256s
        return aimer_verify(
            AIMER_PARAMETERS[key],
            bytes(message_value),
            bytes(signature_value),
            bytes(public_key_value),
            bytes(context),
        )
    end
    parameters = HAETAE_PARAMETERS[
        algorithm == HAETAE2 ? 2 : algorithm == HAETAE3 ? 3 : 5
    ]
    return haetae_verify(
        parameters,
        bytes(public_key_value),
        bytes(message_value),
        bytes(signature_value),
        bytes(context),
    )
end

"""
    encapsulate(algorithm, public_key) -> EncapsulatedSecret

Create a ciphertext and shared secret for the holder of `public_key`.
"""
function encapsulate(algorithm::KEMAlgorithm, public_key_value::AbstractVector{UInt8})
    sizes = kem_sizes(algorithm)
    require_size(public_key_value, sizes.public_key, "public key")

    if algorithm in (NTRUPlus768, NTRUPlus864, NTRUPlus1152)
        parameters = NTRUPLUS_PARAMETERS[
            algorithm == NTRUPlus768 ? 768 : algorithm == NTRUPlus864 ? 864 : 1152
        ]
        result = ntruplus_encapsulate(parameters, bytes(public_key_value))
        return EncapsulatedSecret(result.ciphertext, result.shared_secret)
    end
    parameters = SMAUGT_PARAMETERS[
        algorithm == SMAUGT128 ? 128 :
        algorithm == SMAUGT192 ? 192 :
        algorithm == SMAUGT256 ? 256 : :timer
    ]
    result = smaugt_encapsulate(parameters, bytes(public_key_value))
    return EncapsulatedSecret(result.ciphertext, result.shared_secret)
end

"""
    decapsulate(algorithm, ciphertext, secret_key) -> Vector{UInt8}

Recover the shared secret using the recipient's secret key.
"""
function decapsulate(
    algorithm::KEMAlgorithm,
    ciphertext_value::AbstractVector{UInt8},
    secret_key_value::AbstractVector{UInt8},
)
    sizes = kem_sizes(algorithm)
    require_size(ciphertext_value, sizes.ciphertext, "ciphertext")
    require_size(secret_key_value, sizes.secret_key, "secret key")

    if algorithm in (NTRUPlus768, NTRUPlus864, NTRUPlus1152)
        parameters = NTRUPLUS_PARAMETERS[
            algorithm == NTRUPlus768 ? 768 : algorithm == NTRUPlus864 ? 864 : 1152
        ]
        return ntruplus_decapsulate(
            parameters,
            bytes(ciphertext_value),
            bytes(secret_key_value),
        )
    end
    parameters = SMAUGT_PARAMETERS[
        algorithm == SMAUGT128 ? 128 :
        algorithm == SMAUGT192 ? 192 :
        algorithm == SMAUGT256 ? 256 : :timer
    ]
    return smaugt_decapsulate(
        parameters,
        bytes(ciphertext_value),
        bytes(secret_key_value),
    )
end

end
