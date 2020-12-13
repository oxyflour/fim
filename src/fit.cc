#include <iostream>
#include <fstream>

#include "fit.h"
#include "utils.h"

using namespace std;
using namespace grid;
using namespace fit;

Port::Port(Grid &grid, cst::port_type &port, float epsi) {
    src = port.src;
    dst = port.dst;
    pos = grid.ParsePort(src, dst, epsi);
    auto c = pos.size() / 2;
    auto d = pos[c] - pos[c - 1];
    idx = grid.GetIndex(pos[c - 1], 0);
    dir = d.x ? 0 : d.y ? 1 : 2;
    power = sqrt(4 / imp);
};

Matrix::Matrix(grid::Grid &grid, float *eps, float *mue) {
    this->grid = &grid;
    this->eps = eps;
    this->mue = mue;
}

Coefficient::Coefficient(Matrix &mat, float dt) {
    grid = mat.grid;
    auto nvar = grid->xs.size() * grid->ys.size() * grid->zs.size() * 3;
    le = new float[nvar];
    re = new float[nvar];
    lh = new float[nvar];
    rh = new float[nvar];

    ifstream e1("E:\\e1.txt"), e2("E:\\e2.txt"), h1("E:\\h1.txt"), h2("E:\\h2.txt");

    float kap = 0, rho = 0;
    for (int i = 0; i < nvar; i ++) {
        e1 >> le[i];
        e2 >> re[i];
        h1 >> lh[i];
        h2 >> rh[i];
        /*
        auto ep = mat.eps[i], mu = mat.mue[i];
        le[i] = (1 - kap * ep * dt / 2) / (1 + kap * ep * dt / 2);
        re[i] = dt * ep / (1 + kap * ep * dt / 2);
        lh[i] = (1 - rho * mu * dt / 2) / (1 + rho * mu * dt / 2);
        rh[i] = dt * mu / (1 + rho * mu * dt / 2);
         */
    }
}

Coefficient::~Coefficient() {
    delete le;
    delete re;
    delete lh;
    delete rh;
}

void Coefficient::Add(Port &port) {
    auto &pos = port.pos;
    for (int i = 0, len = pos.size(), c = len / 2 - 1; i < len - 1; i ++) {
        auto d = pos[i + 1] - pos[i];
        auto g = grid->GetIndex(pos[i], d.x ? 0 : d.y ? 1 : 2);
        if (i != c) {
            //le[g] = 1;
            //re[g] = 0;
        } else {
            //auto scale = 1 + re[g] / port.imp;
            //le[g] /= scale;
            //re[g] /= scale;
        }
    }
    ports.push_back(port);
}
