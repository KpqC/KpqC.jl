#ifndef KPQC_NATIVE_EXPORT_H
#define KPQC_NATIVE_EXPORT_H

#if defined(_WIN32)
#define KPQC_EXPORT __declspec(dllexport)
#elif defined(__GNUC__) || defined(__clang__)
#define KPQC_EXPORT __attribute__((visibility("default")))
#else
#define KPQC_EXPORT
#endif

#endif
