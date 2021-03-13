#include <vector>

#include "utils/cuda.h"

#ifndef GRID_H
#define GRID_H

namespace grid {
    struct Grid {
        Grid(std::vector<double> xs, std::vector<double> ys, std::vector<double> zs);
        std::vector<double> xs, ys, zs;
        int nx, ny, nz, nxy, nxyz, nvar;

        float3 At(int3 idx);
        int GetIndex(int3 idx, int dir);
        int GetIndex(int i, int j, int k, int d = 0);
        int3 FindIndex(float3 pos, float epsi);
        std::vector<int3> ParsePort(float3 src, float3 dst, float epsi);
    };
}

#endif
