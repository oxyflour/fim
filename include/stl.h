#include <string>
#include <vector>

#include "utils/cuda.h"
#include "grid.h"

#ifndef STL_H
#define STL_H

namespace stl {
    struct Bound {
        double3 min = { 1e9, 1e9, 1e9 }, max = { -1e9, -1e9, -1e9 };
    };

    struct Joint {
        int2 e;
        double u;
    };

    struct Mesh {
        std::vector<double3> vertices;
        std::vector<int3> faces;
        std::vector<double3> normals;
        std::vector<Bound> bounds;
        double3 min = { 1e9, 1e9, 1e9 }, max = { -1e9, -1e9, -1e9 };
        double3 get(Joint &item);
    };

    void save(std::string file, Mesh &mesh);
    Mesh load(std::string file, double tol = 1e-6);

    struct Mesher {
        Mesher(grid::Grid &grid, Mesh &mesh);
    private:
        void SplitX(Mesh &mesh, int i, double x);
        std::vector<double> xs, ys, zs;

        double3 *vertices, *normals;
        int3 *faces;
        Bound *bounds;
        int faceNum, vertNum;

        double tol;
    };
}

#endif
