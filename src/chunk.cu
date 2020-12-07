#include "utils/cuda_utils.h"
#include "chunk.h"

#include "stdio.h"

#ifdef _WIN32
#define DLL_EXPORT __declspec(dllexport)
#else
#define DLL_EXPORT
#endif

#ifdef __INTELLISENSE__
#define $nx 1
#define $ny 1
#define $nz 1
#define $sg -1
#define $sd 0
#endif

template <int NX, int NY, int NZ, int SG, int SD, int NXY = NX * NY, int NXYZ = NXY * NZ>
class Chunk_$i {
    __device__ __forceinline__ int get_idx(int i, int j, int k) {
        return i + j * NX + k * NXY;
    }
    __device__ __forceinline__ void get_ijk(int g, int &i, int &j, int &k) {
        k = g / NXY;
        g -= k * NXY;
        j = g / NX;
        g -= j * NX;
        i = g;
    }
public:
    // core
    float Ex[NXYZ], Ey[NXYZ], Ez[NXYZ], Hx[NXYZ], Hy[NXYZ], Hz[NXYZ];
    float LEx[NXYZ], LEy[NXYZ], LEz[NXYZ], LHx[NXYZ], LHy[NXYZ], LHz[NXYZ];
    float REx[NXYZ], REy[NXYZ], REz[NXYZ], RHx[NXYZ], RHy[NXYZ], RHz[NXYZ];

    float *sig;

    __device__ __forceinline__ void init(float *le, float *re, float *lh, float *rh) {
        for (int i = cuIdx(x); i < NXYZ; i += cuDim(x)) {
            LEx[i] = le[i]; LEy[i] = le[i + NXYZ]; LEz[i] = le[i + NXYZ * 2];
            REx[i] = re[i]; REy[i] = re[i + NXYZ]; REz[i] = re[i + NXYZ * 2];
            LHx[i] = lh[i]; LHy[i] = lh[i + NXYZ]; LHz[i] = lh[i + NXYZ * 2];
            RHx[i] = rh[i]; RHy[i] = rh[i + NXYZ]; RHz[i] = rh[i + NXYZ * 2];
        }
    }
    __device__ __forceinline__ void step_h() {
        register int i, j, k;
        for (auto g = cuIdx(x); g < NXYZ; g += cuDim(x)) {
            get_ijk(g, i, j, k);
            if (i > 0 && j > 0 && k > 0 && i < NX - 1 && j < NY - 1 && k < NZ - 1) {
                Hx[g] = LHx[g] * Hx[g] + RHx[g] * (Ey[get_idx(i+1, j, k)] - Ey[get_idx(i+1, j, k+1)] - Ez[get_idx(i+1, j, k)] + Ez[get_idx(i+1, j+1, k)]);
                Hy[g] = LHy[g] * Hy[g] + RHy[g] * (Ez[get_idx(i, j+1, k)] - Ez[get_idx(i+1, j+1, k)] - Ex[get_idx(i, j+1, k)] + Ex[get_idx(i, j+1, k+1)]);
                Hz[g] = LHz[g] * Hz[g] + RHz[g] * (Ex[get_idx(i, j, k+1)] - Ex[get_idx(i, j+1, k+1)] - Ey[get_idx(i, j, k+1)] + Ey[get_idx(i+1, j, k+1)]);
            }
        }
    }
    __device__ __forceinline__ void step_e(float s, float *p) {
        register int i, j, k;
        for (auto g = cuIdx(x); g < NXYZ; g += cuDim(x)) {
            get_ijk(g, i, j, k);
            if (i > 0 && j > 0 && k > 0 && i < NX - 1 && j < NY - 1 && k < NZ - 1) {
                register float sx = 0, sy = 0, sz = 0;
                if (g == SG) {
                    SD == 0 ? (sx = s) : SD == 1 ? (sy = s) : (sz = s);
                }
                Ex[g] = LEx[g] * Ex[g] + REx[g] * (Hy[get_idx(i, j-1, k-1)] - Hy[get_idx(i, j-1, k)] - Hz[get_idx(i, j-1, k-1)] + Hz[get_idx(i, j, k-1)] + sx);
                Ey[g] = LEy[g] * Ey[g] + REy[g] * (Hz[get_idx(i-1, j, k-1)] - Hz[get_idx(i, j, k-1)] - Hx[get_idx(i-1, j, k-1)] + Hx[get_idx(i-1, j, k)] + sy);
                Ez[g] = LEz[g] * Ez[g] + REz[g] * (Hx[get_idx(i-1, j-1, k)] - Hx[get_idx(i-1, j, k)] - Hy[get_idx(i-1, j-1, k)] + Hy[get_idx(i, j-1, k)] + sz);
                if (g == SG) {
                    *p = SD == 0 ? Ex[g] : SD == 1 ? Ey[g] : Ez[g];
                }
            }
        }
    }
};

__device__ Chunk_$i<$nx, $ny, $nz, $sg, $sd> chunk_$i;
__global__ void kernel_init_$i(float *le, float *re, float *lh, float *rh) {
    chunk_$i.init(le, re, lh, rh);
}
__global__ void kernel_step_h_$i() {
    chunk_$i.step_h();
}
__global__ void kernel_step_e_$i(float s, float *p) {
    chunk_$i.step_e(s, p);
}

extern "C" DLL_EXPORT int init_$i(float *le, float *re, float *lh, float *rh) {
    constexpr int NXYZ = $nx * $ny * $nz, NVAR = NXYZ * 3;
    kernel_init_$i CU_DIM(2048, 256) (
        to_device(le, NVAR), to_device(re, NVAR),
        to_device(lh, NVAR), to_device(rh, NVAR));
    CU_ASSERT(cudaDeviceSynchronize());
    chunk_$i.sig = malloc_device<float>(1);
    return 0;
}

extern "C" DLL_EXPORT float step_$i(float s) {
    kernel_step_h_$i CU_DIM(2048, 256) ();
    kernel_step_e_$i CU_DIM(2048, 256) (s, chunk_$i.sig);
    CU_ASSERT(cudaGetLastError());
    CU_ASSERT(cudaMemcpy(&s, chunk_$i.sig, sizeof(float), cudaMemcpyDefault));
    return s;
}

extern "C" DLL_EXPORT int quit_$i() {
    return 0;
}
