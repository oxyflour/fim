#include <regex>
#include <fstream>
#include <filesystem>
#include <iostream>

#include <nlohmann/json.hpp>
#include <picosha2.h>

#include "generated/cst_run_bas.h"

#include "cuda.h"
#include "cst.h"
#include "utils.h"
#include "utils/check.h"

using json = nlohmann::json;
using namespace std;
using namespace grid;
using namespace cst;

static map<string, utils::DLL*> dllCache;
static auto getCstDir(string version) {
    HKEY handle;
    auto key = wstring(L"SOFTWARE\\Wow6432Node\\CST AG\\CST DESIGN ENVIRONMENT\\") + utils::utf8ToWstring(version);
    auto ret = RegOpenKeyExW(HKEY_LOCAL_MACHINE, key.c_str(), 0, KEY_READ, &handle);
    CHECK(ret == ERROR_SUCCESS, "RegQuery Error: " + to_string(ret));

    wchar_t buf[1024] = { 0 };
    DWORD len = sizeof(buf) / sizeof(wchar_t);
    ret = RegQueryValueExW(handle, L"INSTALLPATH", 0, NULL, (LPBYTE) buf, &len);
    CHECK(ret == ERROR_SUCCESS, "RegQuery INSTALLPATH Error: " + to_string(ret));
    return utils::wstringToUtf8(buf);
}
static auto getCstExe(string cstDir) {
    return cstDir + "AMD64\\CST DESIGN ENVIRONMENT_AMD64.exe";
}
static auto getCstDll(string cstDir) {
    wchar_t cwd[1024] = { 0 };
    DWORD len = sizeof(cwd) / sizeof(wchar_t);
    GetCurrentDirectoryW(len, cwd);

    auto dir = cstDir + "AMD64";
    SetCurrentDirectoryW(utils::utf8ToWstring(dir).c_str());
    auto dll = new utils::DLL(dir + "\\CSTResultReader_AMD64.dll");
    SetCurrentDirectoryW(cwd);

    return dll;
}

map<string, float> UNITS = {
    { "ns", 1e-9f },
};

static auto hashOfFile(string path) {
    auto cstFile = ifstream(path, ios::binary);
    vector<unsigned char> hash(picosha2::k_digest_size);
    picosha2::hash256(cstFile, hash.begin(), hash.end());
    cstFile.close();
    return picosha2::hash256_hex_string(hash);
}

string Project::MakeCacheAndLoadSettings(string cstDir, bool useCache) {
    auto fileHash = hashOfFile(path).substr(0, 16) + "-" + picosha2::hash256_hex_string(path).substr(0, 16);

    auto tmpPath = filesystem::temp_directory_path() / ("cst-parser-" + fileHash),
        basPath = tmpPath / "run.bas",
        logPath = tmpPath / "run.log",
        jsonPath = tmpPath / "run.json",
        cstPath = tmpPath / "input.cst",
        retPath = tmpPath / "input" / "Result" / "Model.log";
    if (!filesystem::exists(tmpPath / "ok") || !useCache) {
        filesystem::create_directories(tmpPath);
        filesystem::copy_file(path, cstPath, filesystem::copy_options::overwrite_existing);
        utils::writeFile(basPath.u8string(), SRC_CST_RUN_BAS);

        _putenv(("CST_PATH=" + cstPath.u8string()).c_str());
        _putenv(("JSON_PATH=" + jsonPath.u8string()).c_str());
        _putenv("EXPORT_SOLIDS=TRUE");
        // _putenv("BUILD_MATRIX=TRUE");

        // https://stackoverflow.com/questions/9964865/c-system-not-working-when-there-are-spaces-in-two-different-parameters
        auto cmd = "\"\"" + getCstExe(cstDir) + "\" -i -m \"" + basPath.u8string() + "\" > \"" + logPath.u8string() + "\" 2>&1\"";
        auto code = system(cmd.c_str());
        CHECK(code == 0, "Parse CST Project failed: " + cmd + ", see " + logPath.u8string());
    }

    ifstream retInput(retPath);
    string line;
    regex re("\\s*without subcycles:\\s+(\\S+) (\\w+).*");
    smatch dtMatch;
    while (std::getline(retInput, line)) {
        if (regex_match(line, dtMatch, re)) {
            auto unit = dtMatch[2];
            if (UNITS.count(unit) > 0) {
                dt = stof(dtMatch[1]) * UNITS[unit];
            } else {
                cerr << "WARN: unknown unit '" << unit << "' in " << retPath << endl;
            }
        }
    }
    retInput.close();
    // CHECK(dt > 0, "cannot get dt from '" + logPath.u8string() + "'");

    auto meta = json::parse(utils::readFile(jsonPath.u8string()));
    for (auto item : meta["ports"]) {
        auto jSrc = item["src"], jDst = item["dst"];
        auto src = make_float3(jSrc[0].get<float>(), jSrc[1].get<float>(), jSrc[2].get<float>());
        auto dst = make_float3(jDst[0].get<float>(), jDst[1].get<float>(), jDst[2].get<float>());
        ports.push_back({ src, dst });
    }
    for (auto item : meta["solids"]) {
        auto name = item["name"].get<string>(),
            material = item["material"].get<string>(),
            stl = item["stl"].get<string>();
        solids.push_back(solid_type { name, material, stl });
    }
    auto jUnits = meta["units"];
    units.geometry = jUnits["geometry"].get<float>();
    units.time = jUnits["time"].get<float>();
    units.frequency = jUnits["frequency"].get<float>();

    for (auto x : meta["grid"]["xs"]) {
        xs.push_back(x.get<double>());
    }
    for (auto y : meta["grid"]["ys"]) {
        ys.push_back(y.get<double>());
    }
    for (auto z : meta["grid"]["zs"]) {
        zs.push_back(z.get<double>());
    }

    ifstream sourceInput(jsonPath.u8string() + ".excitation.txt");
    string header;
    std::getline(sourceInput, header);
    std::getline(sourceInput, header);
    float x, y;
    while (sourceInput >> x >> y) {
        excitation.x.push_back(x * units.time);
        excitation.y.push_back(y);
    }
    sourceInput.close();

    utils::writeFile((tmpPath / "ok").u8string(), "");
    return cstPath.u8string();
}

Project::Project(string path, string version, bool useCache, bool keepCache) {
    this->path = path;
    this->version = version;
    this->keepCache = keepCache;

    auto cstDir = getCstDir(version);
    cachePath = MakeCacheAndLoadSettings(cstDir, useCache);

#ifdef USE_CST_DLL
    if (dllCache.count(version) == 0) {
        dllCache[version] = getCstDll(cstDir);
    }
    dll = dllCache[version];
    auto OpenProject = (CST_OpenProject_PTR) dll->getProc("CST_OpenProject");
    auto ret = OpenProject ? OpenProject(cachePath.c_str(), &handle) : -1;
    CHECK(ret == 0, "Open CST project '" + path + "' failed with code " + to_string(ret) + ", cache " + cachePath);
#endif
}

Project::~Project() {
    if (!keepCache) {
        filesystem::remove_all(filesystem::path(cachePath).parent_path());
    }
#ifdef USE_CST_DLL
    auto CloseProject = (CST_CloseProject_PTR) dll->getProc("CST_CloseProject");
    CloseProject(&handle);
    dll = NULL;
#endif
}

#ifdef USE_CST_DLL
Grid Project::GetHexGrid() {
    CHECK(dll != NULL, "project " + path + " already destroyed");

    int nxyz[3] = { 0 };
    auto getMeshInfo = (CST_GetHexMeshInfo_PTR) dll->getProc("CST_GetHexMeshInfo");
    auto ret = getMeshInfo ? getMeshInfo(&handle, nxyz) : -1;
    CHECK(ret == 0, "GetHexMeshInfo Error: " + to_string(ret));

    int nx = nxyz[0], ny = nxyz[1], nz = nxyz[2], sz = nx + ny + nz;
    auto array = vector<double>(sz);
    auto getHexMesh = (CST_GetHexMesh_PTR) dll->getProc("CST_GetHexMesh");
    ret = getHexMesh ? getHexMesh(&handle, array.data()) : -1;
    CHECK(ret == 0, "GetHexMesh Error: " + to_string(ret));

    auto xs = vector<double>(nx);
    for (int i = 0; i < nx; i ++) {
        xs[i] = array[i];
    }
    auto ys = vector<double>(ny);
    for (int i = 0; i < ny; i ++) {
        ys[i] = array[i + nx];
    }
    auto zs = vector<double>(nz);
    for (int i = 0; i < nz; i ++) {
        zs[i] = array[i + nx + ny];
    }
    return Grid(xs, ys, zs);
}

float *Project::GetMatrix(int mat) {
    CHECK(dll != NULL, "project " + path + " already destroyed");

    int nxyz[3] = { 0 };
    auto getMeshInfo = (CST_GetHexMeshInfo_PTR) dll->getProc("CST_GetHexMeshInfo");
    auto ret = getMeshInfo ? getMeshInfo(&handle, nxyz) : -1;
    CHECK(ret == 0, "GetHexMeshInfo Error: " + to_string(ret));

    int nx = nxyz[0], ny = nxyz[1], nz = nxyz[2], sz = nx * ny * nz * 3;
    auto array = new float[sz];
    auto getMaterialMatrix = (CST_GetMaterialMatrixHexMesh_PTR) dll->getProc("CST_GetMaterialMatrixHexMesh");
    ret = getMaterialMatrix ? getMaterialMatrix(&handle, mat, array) : -1;
    CHECK(ret == 0, "GetMaterialMatrix Error: " + to_string(ret));

    if (mat == 101) {
        // reorder for mue
        auto clone = new float[sz];
        int nxy = nx * ny, nxyz = nxy * nz;
        auto idx = [=](int i, int j, int k, int d){ return (i % nx) + (j % ny) * nx + (k % nz) * nxy + d * nxyz; };
        for (int i = 0; i < nx; i ++) {
            for (int j = 0; j < ny; j ++) {
                for (int k = 0; k < nz; k ++) {
                    clone[idx(i, j, k, 0)] = array[idx(i + 1, j, k, 0)];
                    clone[idx(i, j, k, 1)] = array[idx(i, j + 1, k, 1)];
                    clone[idx(i, j, k, 2)] = array[idx(i, j, k + 1, 2)];
                }
            }
        }
        return clone;
    } else {
        return array;
    }
}

double *Project::Get1DResult(string tree, int num, int type) {
    CHECK(dll != NULL, "project " + path + " already destroyed");

    auto getResultSize = (CST_Get1DResultSize_PTR) dll->getProc("CST_Get1DResultSize");
    int size = 0;
    auto ret = getResultSize ? getResultSize(&handle, tree.c_str(), num, &size) : -1;
    CHECK(ret == 0, "Get1DResultSize Error: " + to_string(ret));

    auto getResultData = type == 0 ?
        (CST_Get1DRealDataAbszissa_PTR) dll->getProc("CST_Get1DRealDataAbszissa") :
        (CST_Get1DRealDataOrdinate_PTR) dll->getProc("CST_Get1DRealDataOrdinate");
    auto array = new double[size];
    ret = getResultData ? getResultData(&handle, tree.c_str(), num, array) : -1;
    CHECK(ret == 0, "Get1DResultData Error: " + to_string(ret));

    return array;
}
#endif
