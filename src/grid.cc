#include <math.h>

#include "cuda.h"
#include "grid.h"

using namespace std;
using namespace grid;

float3 Grid::At(int3 idx) {
    return float3 { (float) xs[idx.x], (float) ys[idx.y], (float) zs[idx.z] };
}

int Grid::GetFlatIndex(int3 idx, int dir) {
    int nx = xs.size(), ny = ys.size(), nz = zs.size();
    return idx.x + idx.y * nx + idx.z * nx * ny + dir * nx * ny * nz;
}

int3 Grid::FindIndex(float3 pos, float epsi) {
    int nx = xs.size(), ny = ys.size(), nz = zs.size();
    for (int i = 0; i < nx; i ++) {
        for (int j = 0; j < ny; j ++) {
            for (int k = 0; k < nz; k ++) {
                float3 pt = { (float) xs[i], (float) ys[j], (float) zs[k] };
                if (length(pos - pt) < epsi) {
                    return int3 { i, j, k };
                }
            }
        }
    }
    return int3 { -1, -1, -1 };
}

vector<int3> Grid::ParsePort(float3 src, float3 dst, float epsi) {
    auto i = FindIndex(src, epsi),
        j = FindIndex(dst, epsi);
    auto len = sum(abs(i - j));
    vector<int3> ret;
    if (len == 0) {
        return ret;
    }

    ret.push_back(i);
    int3 dir[6] = { 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1 };
    for (int s = 1; s < len + 1; s ++) {
        float dist = 0xffffffff;
        int3 idx;
        for (int m = 0; m < 6; m ++) {
            auto n = i + dir[m];
            auto r = length(At(n) - dst);
            if (r < dist) {
                dist = r;
                idx = n;
            }
        }
        ret.push_back(i = idx);
    }

    return ret;
}
