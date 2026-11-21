using KpqC
using Test

const SMOKE_MESSAGE = Vector{UInt8}(codeunits("release-manifest:v3"))
const SMOKE_CONTEXT = Vector{UInt8}(codeunits("test-suite"))

@testset "public API smoke tests" begin
    @test length(signature_algorithms()) == 9
    @test length(kem_algorithms()) == 7

    @testset "signature round trip" begin
        keys = generate_keypair(HAETAE2)
        signature = sign(
            HAETAE2,
            SMOKE_MESSAGE,
            keys.secret_key;
            context = SMOKE_CONTEXT,
        )

        @test length(signature) == signature_sizes(HAETAE2).signature
        @test verify(
            HAETAE2,
            SMOKE_MESSAGE,
            signature,
            keys.public_key;
            context = SMOKE_CONTEXT,
        )
        @test !verify(HAETAE2, SMOKE_MESSAGE, signature, keys.public_key)
    end

    @testset "key encapsulation round trip" begin
        keys = generate_keypair(NTRUPlus768)
        outbound = encapsulate(NTRUPlus768, keys.public_key)

        @test length(outbound.ciphertext) == kem_sizes(NTRUPlus768).ciphertext
        @test decapsulate(
            NTRUPlus768,
            outbound.ciphertext,
            keys.secret_key,
        ) == outbound.shared_secret
    end
end
