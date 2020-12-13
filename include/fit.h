#include "cuda.h"
#include "utils.h"
#include "grid.h"
#include "cst.h"

#ifndef FIT_H
#define FIT_H

namespace fit {
    typedef struct Port {
        float3 src, dst;
        std::vector<int3> pos;
        int idx;
        int dir;
        float imp = 50;
        float power = 1;
        Port(grid::Grid &grid, cst::port_type &port, float epsi = 1e-6);
    } Port;

    typedef struct Matrix {
        grid::Grid *grid;
        float *eps, *mue;
        Matrix(grid::Grid &grid, float *eps, float *mue);
    } Matrix;

    typedef struct Coefficient {
        Coefficient(Matrix &mat, float dt);
        ~Coefficient();
        grid::Grid *grid;
        float *le, *re, *lh, *rh;
        void Add(Port &port);
        std::vector<Port> ports;
    } Coefficient;
}

#endif
