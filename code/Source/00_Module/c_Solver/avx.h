!
#ifdef AVX512
#define avx 64
#define avi 32
#elif AVX2 || AVX
#define avx 32
#define avi 16
#elif SSE
#define avx 16
#define avi 8
#else
#define avx 8
#define avi 4
#endif
