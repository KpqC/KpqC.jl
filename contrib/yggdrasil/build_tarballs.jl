using BinaryBuilder

name = "libkpqc"
version = v"0.1.0"

# Build the exact vendored sources used by KpqC.jl.  The commit is immutable so
# that the resulting JLL artifacts remain reproducible.
sources = [
    GitSource(
        "https://github.com/KpqC/KpqC.jl.git",
        "59e97b45b8e863c0430f609ea13f43d9b6dc9e4d",
    ),
]

script = raw"""
cd "${WORKSPACE}/srcdir/KpqC.jl"

mkdir -p build "${libdir}"

for source in native/core_*.c; do
    include_flags=()

    case "$(basename "${source}")" in
        core_haetae*)
            include_flags+=("-Inative/third_party/HAETAE/include")
            ;;
        core_ntruplus768.c)
            include_flags+=("-Inative/third_party/NTRUplus/NTRU+768")
            ;;
        core_ntruplus864.c)
            include_flags+=("-Inative/third_party/NTRUplus/NTRU+864")
            ;;
        core_ntruplus1152.c)
            include_flags+=("-Inative/third_party/NTRUplus/NTRU+1152")
            ;;
        core_smaugt*|core_timer.c)
            include_flags+=("-Inative/third_party/SMAUG-T/include")
            ;;
    esac

    object="build/$(basename "${source}" .c).o"
    "${CC}" -c -O2 -std=c11 -fPIC "${include_flags[@]}" \
        -o "${object}" "${source}"
done

link_flags=(-lm)

if [[ "${target}" == *-apple-* ]]; then
    shared_flag="-dynamiclib"
else
    shared_flag="-shared"
fi

if [[ "${target}" == *-w64-mingw32 ]]; then
    link_flags+=("-lbcrypt")
fi

"${CC}" "${shared_flag}" build/*.o "${link_flags[@]}" \
    -o "${libdir}/libkpqc_julia.${dlext}"

install_license LICENSE
license_dir="${prefix}/share/licenses/libkpqc"
install -Dvm 0644 THIRD_PARTY_NOTICES.md \
    "${license_dir}/THIRD_PARTY_NOTICES.md"
install -Dvm 0644 native/third_party/AIMer/LICENSE \
    "${license_dir}/AIMer-LICENSE"
install -Dvm 0644 native/third_party/HAETAE/LICENSE \
    "${license_dir}/HAETAE-LICENSE"
install -Dvm 0644 native/third_party/NTRUplus/LICENSE \
    "${license_dir}/NTRUplus-LICENSE"
install -Dvm 0644 native/third_party/SMAUG-T/LICENSE \
    "${license_dir}/SMAUG-T-LICENSE"
"""

platforms = supported_platforms()

products = [
    LibraryProduct("libkpqc_julia", :libkpqc_julia),
]

dependencies = Dependency[]

build_tarballs(
    ARGS,
    name,
    version,
    sources,
    script,
    platforms,
    products,
    dependencies;
    julia_compat = "1.10",
)
