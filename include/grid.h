#include <vector>

#include "utils/cuda.h"

#ifndef GRID_H
#define GRID_H

namespace grid {
    using namespace std;

    struct Grid {
        Grid(vector<double> xs, vector<double> ys, vector<double> zs);
        vector<double> xs, ys, zs;
        int nx, ny, nz, nxy, nxyz, nvar;

        float3 At(int3 idx);
        int GetIndex(int4 idx);
        int GetIndex(int3 idx, int dir = 0);
        int GetIndex(int i, int j, int k, int d = 0);
        int4 GetIndex(int n);
        int3 FindIndex(float3 pos, float epsi);
        vector<int3> ParsePort(float3 src, float3 dst, float epsi);
    };
}

#endif
