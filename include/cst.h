#include <map>
#include <vector>

#include <windows.h>

#include "vendor/cst.h"
#include "utils/cuda.h"
#include "utils.h"
#include "grid.h"

#ifndef CST_H
#define CST_H

namespace cst {
    struct port_type {
        float3 src, dst;
    };

    struct excitation_type {
        std::vector<float> x, y;
    };

    struct solid_type {
        std::string name, material, stl;
    };

    struct units_type {
        float geometry, time, frequency;
    };

    struct Project {
        Project(std::string path, std::string version, bool useCache = true, bool keepCache = false);
        ~Project();

        // will be filled by MakeCacheAndLoadSettings
        float dt = -1;
        std::vector<port_type> ports;
        std::vector<solid_type> solids;
        excitation_type excitation;
        units_type units;
        std::vector<double> xs, ys, zs;

#ifdef USE_CST_DLL
        grid::Grid GetHexGrid();
        float *GetMatrix(int mat);
        double *Get1DResult(std::string tree, int num, int type);
#endif

private:
        std::string path;
        std::string version;
        CSTProjHandle handle;

        bool keepCache = false;
        std::string cachePath;
        std::string MakeCacheAndLoadSettings(std::string cst, bool useCache);

#ifdef USE_CST_DLL
        utils::DLL *dll;
#endif
    };
}

#endif
