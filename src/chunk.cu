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
    // pml
    float EP1[NXYZ], EP2[NXYZ], HP1[NXYZ], HP2[NXYZ];
    __device__ __forceinline__ void step_h() {
        register int i, j, k;
        for (auto g = cuIdx(x); g < NXYZ; g += cuDim(x)) {
            get_ijk(g, i, j, k);
            if (i > 0 && j > 0 && k > 0 && i < NX - 1 && j < NY - 1 && k < NZ - 1) {
                Hx[g] = LHx[g] * Hx[g] + RHx[g] * (Ey[get_idx(i, j, k+1)] - Ey[g] - Ez[get_idx(i, j+1, k)] + Ez[g]);
                Hy[g] = LHy[g] * Hy[g] + RHy[g] * (Ez[get_idx(i+1, j, k)] - Ez[g] - Ex[get_idx(i, j, k+1)] + Ex[g]);
                Hz[g] = LHz[g] * Hz[g] + RHz[g] * (Ex[get_idx(i, j+1, k)] - Ex[g] - Ey[get_idx(i+1, j, k)] + Ey[g]);
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
                Ex[g] = LEx[g] * Ex[g] + REx[g] * (Hy[get_idx(i, j, k-1)] - Hy[g] - Hz[get_idx(i, j-1, k)] + Hz[g] + sx);
                Ey[g] = LEy[g] * Ey[g] + REy[g] * (Hz[get_idx(i-1, j, k)] - Hz[g] - Hx[get_idx(i, j, k-1)] + Hx[g] + sy);
                Ez[g] = LEz[g] * Ez[g] + REz[g] * (Hx[get_idx(i, j-1, k)] - Hx[g] - Hy[get_idx(i-1, j, k)] + Hy[g] + sz);
                if (g == SG) {
                    *p = SD == 0 ? Ex[g] : SD == 1 ? Ey[g] : Ez[g];
                }
            }
        }
    }
};

float *sig_$i;
__device__ Chunk_$i<$nx, $ny, $nz, $sg, $sd> chunk_$i;
__global__ void kernel_init_$i(float *le, float *re, float *lh, float *rh) {
    constexpr int NXYZ = $nx * $ny * $nz;
    auto &c = chunk_$i;
    for (int i = cuIdx(x); i < NXYZ; i += cuDim(x)) {
        c.LEx[i] = le[i]; c.LEy[i] = le[i + NXYZ]; c.LEz[i] = le[i + NXYZ * 2];
        c.REx[i] = re[i]; c.REy[i] = re[i + NXYZ]; c.REz[i] = re[i + NXYZ * 2];
        c.LHx[i] = lh[i]; c.LHy[i] = lh[i + NXYZ]; c.LHz[i] = lh[i + NXYZ * 2];
        c.RHx[i] = rh[i]; c.RHy[i] = rh[i + NXYZ]; c.RHz[i] = rh[i + NXYZ * 2];
    }
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
    CU_ASSERT(cudaGetLastError());
    sig_$i = malloc_device<float>(1);
    return 0;
}

extern "C" DLL_EXPORT float step_$i(float s) {
    kernel_step_h_$i CU_DIM(2048, 256) ();
    kernel_step_e_$i CU_DIM(2048, 256) (s, sig_$i);
    CU_ASSERT(cudaGetLastError());
    CU_ASSERT(cudaMemcpy(&s, sig_$i, sizeof(float), cudaMemcpyDefault));
    return s;
}

extern "C" DLL_EXPORT int quit_$i() {
    return 0;
}
