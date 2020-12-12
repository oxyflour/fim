#include <vector>

#include "cuda.h"

#ifndef GRID_H
#define GRID_H

namespace grid {
    typedef struct Grid {
        std::vector<double> xs, ys, zs;
        float3 At(int3 idx);
        int GetFlatIndex(int3 idx, int dir);
        int3 FindIndex(float3 pos, float epsi);
        std::vector<int3> ParsePort(float3 src, float3 dst, float epsi);
    } Grid;
}

#endif
