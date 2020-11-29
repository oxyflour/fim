#include "utils.h"
#include "grid.h"

#ifndef FIT_H
#define FIT_H

extern "C" {
    typedef int (__cdecl *InitPTR)(float *le, float *re, float *lh, float *rh);
    typedef float (__cdecl *StepPTR)(float s);
    typedef int (__cdecl *QuitPTR)();
}

namespace fit {
    typedef struct Port {
        float3 src, dst;
        std::vector<int3> pos;
        int idx;
        int dir;
        Port(Grid &grid, float3 src, float3 dst, float epsi = 1e-6);
    } Port;

    typedef struct Matrix {
        Grid *grid;
        float *eps, *mue;
        Matrix(Grid &grid, float *eps, float *mue) : eps{eps}, mue{mue}, grid{&grid} { };
    } Matrix;

    typedef struct Coefficient {
        Coefficient(Matrix &mat, float dt);
        ~Coefficient();
        Grid *grid;
        float *le, *re, *lh, *rh;
        void UpdateFromPort(Port &port);
    } Coefficient;

    typedef struct Solver {
        Solver(Matrix &mat, float dt, std::vector<Port> &ports);
        ~Solver();
        float Step(float s);
        utils::DLL *dll;
        InitPTR FnInit;
        StepPTR FnStep;
        QuitPTR FnQuit;
    } Solver;
}

#endif
