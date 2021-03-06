#include "utils/cuda.h"

#include "utils.h"
#include "grid.h"
#include "cst.h"

#ifndef FIT_H
#define FIT_H

namespace fit {
    using namespace std;

    struct Port {
        double3 src, dst;
        vector<int3> pos;
        int idx;
        int dir;
        float imp = 50;
        float power = 1;
        Port(grid::Grid &grid, cst::port_type &port, float epsi = 1e-6);
    };

    struct Matrix {
        grid::Grid *grid;
        float *eps, *mue;
        Matrix(grid::Grid &grid, float *eps, float *mue);
    };

    struct Coefficient {
        Coefficient(Matrix &mat, float dt);
        ~Coefficient();
        grid::Grid *grid;
        float *le, *re, *lh, *rh;
        void Add(Port &port);
        vector<Port> ports;
    };
}

#endif
