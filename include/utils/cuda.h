#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "cuda_runtime.h"
#include "cuda_fp16.h"

#ifndef CUDA_H
#define CUDA_H

// https://stackoverflow.com/a/27992604
#ifdef __INTELLISENSE__
dim3 blockIdx;
dim3 blockDim;
dim3 threadIdx;
dim3 gridDim;
#define CU_DIM(grid, block)
#define CU_DIM_MEM(grid, block, bytes)
#define CU_DIM_MEM_STREAM(grid, block, bytes, stream)
extern void __syncthreads();
extern int atomicAdd(int* address, int val);
#else
#define CU_DIM(grid, block) <<<grid, block>>>
#define CU_DIM_MEM(grid, block, bytes) <<<grid, block, bytes>>>
#define CU_DIM_MEM_STREAM(grid, block, bytes, stream) <<<grid, block, bytes, stream>>>
#endif

#define CU_ASSERT(err) do { _cuda_assert((err), __FILE__, __LINE__); } while (0);
inline void _cuda_assert(cudaError_t err, const char *file, int line, bool abort=true)
{
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA FATAL: %s at %s:%d\n", cudaGetErrorString(err), file, line);
        if (abort) {
            exit(1);
        }
    }
}

#define cuIdx(D) (threadIdx.D + blockIdx.D * blockDim.D)
#define cuDim(D) (blockDim.D * gridDim.D)

template <typename T> T *malloc_device(size_t sz) {
    T* out = NULL;
    if (sz > 0) {
        CU_ASSERT(cudaMalloc(&out, sizeof(T) * sz));
    }
    return out;
}

template <typename T> T *to_device(const T *in, size_t sz, T *out = NULL) {
    if (out == NULL) {
        CU_ASSERT(cudaMalloc(&out, sizeof(T) * sz));
    }
    if (sz > 0) {
        CU_ASSERT(cudaMemcpy(out, in, sizeof(T) * sz, cudaMemcpyDefault));
    }
    return out;
}

template <typename T> T *from_device(const T *in, size_t sz, T *out = NULL) {
    if (out == NULL) {
        CU_ASSERT(cudaHostAlloc(&out, sizeof(T) * sz, 0));
    }
    if (sz > 0) {
        CU_ASSERT(cudaMemcpy(out, in, sizeof(T) * sz, cudaMemcpyDefault));
    }
    return out;
}

// math

static int3 operator- (const int3 &a, const int3 &b) {
    return int3 { a.x - b.x, a.y - b.y, a.z - b.z };
}

static int3 operator+ (const int3 &a, const int3 &b) {
    return int3 { a.x + b.x, a.y + b.y, a.z + b.z };
}

static int3 operator- (const int3 &a, const int b) {
    return int3 { a.x - b, a.y - b, a.z - b };
}

static int3 operator+ (const int3 &a, const int b) {
    return int3 { a.x + b, a.y + b, a.z + b };
}

static int3 abs(int3 &a) {
    return int3 { ::abs(a.x), ::abs(a.y), ::abs(a.z) };
}

static int sum(int3 &a) {
    return a.x + a.y + a.z;
}

static float3 operator- (const float3 &a, const float3 &b) {
    return float3 { a.x - b.x, a.y - b.y, a.z - b.z };
}

static float length(const float3 &a) {
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

inline __host__ __device__ double fmin(double a, double b) {
    return a < b ? a : b;
}

inline __host__ __device__ double fmax(double a, double b) {
    return a > b ? a : b;
}

inline __host__ __device__ double3 operator+(double3 a, double3 b) {
    return double3 { a.x + b.x, a.y + b.y, a.z + b.z };
}

inline __host__ __device__ double3 operator+(double3 a, double b) {
    return double3 { a.x + b, a.y + b, a.z + b };
}

inline __host__ __device__ double3 operator-(double3 a, double3 b) {
    return double3 { a.x - b.x, a.y - b.y, a.z - b.z };
}

inline __host__ __device__ double3 operator-(double3 a, double b) {
    return double3 { a.x - b, a.y - b, a.z - b };
}

inline __host__ __device__ double3 operator*(double3 a, double b) {
    return double3 { a.x * b, a.y * b, a.z * b };
}

inline __host__ __device__ double3 operator/(double3 a, double b) {
    return double3 { a.x / b, a.y / b, a.z / b };
}

inline __host__ __device__ double3 lerp(double3 a, double3 b, double f) {
    return a * (1 - f) + b * f;
}

inline __host__ __device__ double lerp(double a, double b, double f) {
    return a * (1 - f) + b * f;
}

inline __host__ __device__ double dot(double3 a, double3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline __host__ __device__ double length(double3 a) {
    return sqrt(dot(a, a));
}

inline __host__ __device__ double3 cross(double3 a, double3 b) { 
    return double3 { a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x }; 
}

inline __host__ __device__ double3 normalize(double3 v) {
    return v / length(v);
}

inline __host__ __device__ double3 fmin(double3 a, double3 b) { 
    return double3 { fmin(a.x, b.x), fmin(a.y, b.y), fmin(a.z, b.z) }; 
}

inline __host__ __device__ double3 fmax(double3 a, double3 b) { 
    return double3 { fmax(a.x, b.x), fmax(a.y, b.y), fmax(a.z, b.z) }; 
}

inline __host__ __device__ double round_by(double a, double tol) {
    return round(a / tol) * tol;
}

inline __host__ __device__ double3 round_by(double3 a, double tol) {
    return double3 { round_by(a.x, tol), round_by(a.y, tol), round_by(a.z, tol) };
}

inline bool operator<(double3 a, double3 b) {
    return (a.x < b.x || (a.x == b.x && a.y < b.y) || (a.x == b.x && a.y == b.y && a.z < b.z));
}

#endif
