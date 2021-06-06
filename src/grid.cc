#include <math.h>

#include "cuda.h"
#include "grid.h"

using namespace std;
using namespace grid;

Grid::Grid(vector<double> xs, vector<double> ys, vector<double> zs) {
    this->xs = xs;
    this->ys = ys;
    this->zs = zs;
    this->nx = xs.size();
    this->ny = ys.size();
    this->nz = zs.size();
    this->nxy = this->nx * this->ny;
    this->nxyz = this->nx * this->ny * this->nz;
    this->nvar = this->nxyz * 3;
}

double3 Grid::At(int3 idx) {
    return double3 { xs[idx.x], ys[idx.y], zs[idx.z] };
}

int Grid::GetIndex(int4 idx) {
    return GetIndex(idx.x, idx.y, idx.z, idx.w);
}

int Grid::GetIndex(int3 idx, int dir) {
    return GetIndex(idx.x, idx.y, idx.z, dir);
}

int Grid::GetIndex(int i, int j, int k, int d) {
    return i + j * nx + k * nxy + d * nxyz;
}

int4 Grid::GetIndex(int i) {
    auto d = i / nxyz;
    i -= d * nxyz;
    auto k = i / nxy;
    i -= k * nxy;
    auto j = i / nx;
    i -= j * nx;
    return int4 { i, j, k, d };
}

int3 Grid::FindIndex(double3 pos, float epsi) {
    int nx = xs.size(), ny = ys.size(), nz = zs.size();
    for (int i = 0; i < nx; i ++) {
        for (int j = 0; j < ny; j ++) {
            for (int k = 0; k < nz; k ++) {
                auto pt = At(int3 { i, j, k });
                if (length(pos - pt) < epsi) {
                    return int3 { i, j, k };
                }
            }
        }
    }
    return int3 { -1, -1, -1 };
}

vector<int3> Grid::ParsePort(double3 src, double3 dst, float epsi) {
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
