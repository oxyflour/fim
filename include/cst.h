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
    using namespace std;

    struct port_type {
        float3 src, dst;
    };

    struct excitation_type {
        vector<float> x, y;
    };

    struct solid_type {
        string name, material, stl;
    };

    struct units_type {
        float geometry, time, frequency;
    };

    struct Project {
        Project(string path, string version, bool useCache = true, bool keepCache = false);
        ~Project();

        // will be filled by MakeCacheAndLoadSettings
        float dt = -1;
        vector<port_type> ports;
        vector<solid_type> solids;
        excitation_type excitation;
        units_type units;
        vector<double> xs, ys, zs;

#ifdef USE_CST_DLL
        grid::Grid GetHexGrid();
        float *GetMatrix(int mat);
        double *Get1DResult(string tree, int num, int type);
#endif

private:
        string path;
        string version;
        CSTProjHandle handle;

        bool keepCache = false;
        string cachePath;
        string MakeCacheAndLoadSettings(string cst, bool useCache);

#ifdef USE_CST_DLL
        utils::DLL *dll;
#endif
    };
}

#endif
