#include "cst.h"
#include "utils.h"
#include "generated/cst_run_bas.h"

#include <fstream>
#include <filesystem>
using namespace std;

#include <nlohmann/json.hpp>
using json = nlohmann::json;

static map<string, utils::DLL*> dllCache;
static auto getCstPath(string &version) {
    HKEY handle;
    auto key = wstring(L"SOFTWARE\\Wow6432Node\\CST AG\\CST DESIGN ENVIRONMENT\\") + utils::utf8ToWstring(version);
    auto ret = RegOpenKeyExW(HKEY_LOCAL_MACHINE, key.c_str(), 0, KEY_READ, &handle);
    ASSERT(ret == ERROR_SUCCESS, "RegQuery Error: " + to_string(ret));

    wchar_t buf[1024] = { 0 };
    DWORD len = sizeof(buf) / sizeof(wchar_t);
    ret = RegQueryValueExW(handle, L"INSTALLPATH", 0, NULL, (LPBYTE) buf, &len);
    ASSERT(ret == ERROR_SUCCESS, "RegQuery INSTALLPATH Error: " + to_string(ret));
    return utils::wstringToUtf8(buf);
}

void cst::Project::ForkAndExportSettings(string cstPath) {
    auto tmp = filesystem::temp_directory_path() / ("cst-parse-" + utils::random(8));
    filesystem::create_directories(tmp);

    auto exe = cstPath + "AMD64\\CST DESIGN ENVIRONMENT_AMD64.exe",
        bas = (tmp / "run.bas").u8string(),
        log = (tmp / "run.log").u8string(),
        json = (tmp / "run.json").u8string(),
        // https://stackoverflow.com/questions/9964865/c-system-not-working-when-there-are-spaces-in-two-different-parameters
        cmd = "\"\"" + exe + "\" -i -m \"" + bas + "\" > " + log + " 2>&1\"";
    utils::writeFile(bas, SRC_CST_RUN_BAS);
    _putenv(("CST_PATH=" + path).c_str());
    _putenv(("JSON_PATH=" + json).c_str());
    _putenv("EXPORT_SOLIDS=TRUE");
    _putenv("BUILD_MATRIX=TRUE");
    auto code = system(cmd.c_str());
    ASSERT(code == 0, "Parse CST Project failed: " + cmd + ", log " + log);
    auto meta = json::parse(utils::readFile(json));

    for (auto item : meta["ports"]) {
        auto jSrc = item["src"], jDst = item["dst"];
        auto src = float3 { jSrc[0].get<float>(), jSrc[1].get<float>(), jSrc[2].get<float>() };
        auto dst = float3 { jDst[0].get<float>(), jDst[1].get<float>(), jDst[2].get<float>() };
        ports.push_back(port_type { src, dst });
    }

    for (auto item : meta["solids"]) {
        auto name = item["name"].get<string>(),
            material = item["material"].get<string>();
        solids.push_back(solid_type { name, material });
    }

    ifstream input(json + ".excitation.txt");
    string header;
    getline(input, header);
    getline(input, header);
    float x, y;
    while (input >> x >> y) {
        excitation.x.push_back(x);
        excitation.y.push_back(y);
    }
    input.close();

    auto jUnits = meta["units"];
    units.geometry = jUnits["geometry"].get<float>();
    units.time = jUnits["time"].get<float>();
    units.frequency = jUnits["frequency"].get<float>();

    filesystem::remove_all(tmp);
}

cst::Project::Project(string &path, string &version) {
    this->path = path;
    this->version = version;

    auto cstPath = getCstPath(version);
    if (dllCache.count(version) == 0) {
        wchar_t cwd[1024] = { 0 };
        DWORD len = sizeof(cwd) / sizeof(wchar_t);
        GetCurrentDirectoryW(len, cwd);

        auto dir = cstPath + "AMD64";
        SetCurrentDirectoryW(utils::utf8ToWstring(dir).c_str());

        dllCache[version] = new utils::DLL(dir + "\\CSTResultReader_AMD64.dll");
        SetCurrentDirectoryW(cwd);
    }

    ForkAndExportSettings(cstPath);

    ASSERT(dllCache.count(version) == 1, "Load Dll Version " + version + " Failed");
    dll = dllCache[version];
    auto OpenProject = (CST_OpenProject_PTR) dll->getProc("CST_OpenProject");
    auto ret = OpenProject ? OpenProject(path.c_str(), &handle) : -1;
    ASSERT(ret == 0, "Open CST Project '" + path + "' Failed: " + to_string(ret));
}

cst::Project::~Project() {
    ASSERT(dll != NULL, "project " + path + " already destroyed");
    auto CloseProject = (CST_CloseProject_PTR) dll->getProc("CST_CloseProject");
    CloseProject(&handle);
    dll = NULL;
}

Grid cst::Project::GetHexGrid() {
    ASSERT(dll != NULL, "project " + path + " already destroyed");

    int nxyz[3] = { 0 };
    auto getMeshInfo = (CST_GetHexMeshInfo_PTR) dll->getProc("CST_GetHexMeshInfo");
    auto ret = getMeshInfo ? getMeshInfo(&handle, nxyz) : -1;
    ASSERT(ret == 0, "GetHexMeshInfo Error: " + to_string(ret));

    int nx = nxyz[0], ny = nxyz[1], nz = nxyz[2], sz = nx + ny + nz;
    auto array = vector<double>(sz);
    auto getHexMesh = (CST_GetHexMesh_PTR) dll->getProc("CST_GetHexMesh");
    ret = getHexMesh ? getHexMesh(&handle, array.data()) : -1;
    ASSERT(ret == 0, "GetHexMesh Error: " + to_string(ret));

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
    return Grid { xs, ys, zs };
}

vector<float> cst::Project::GetMatrix(int mat) {
    ASSERT(dll != NULL, "project " + path + " already destroyed");

    int nxyz[3] = { 0 };
    auto getMeshInfo = (CST_GetHexMeshInfo_PTR) dll->getProc("CST_GetHexMeshInfo");
    auto ret = getMeshInfo ? getMeshInfo(&handle, nxyz) : -1;
    ASSERT(ret == 0, "GetHexMeshInfo Error: " + to_string(ret));

    int nx = nxyz[0], ny = nxyz[1], nz = nxyz[2], sz = nx * ny * nz * 3;
    auto array = vector<float>(sz);
    auto getMaterialMatrix = (CST_GetMaterialMatrixHexMesh_PTR) dll->getProc("CST_GetMaterialMatrixHexMesh");
    ret = getMaterialMatrix ? getMaterialMatrix(&handle, mat, array.data()) : -1;
    ASSERT(ret == 0, "GetMaterialMatrix Error: " + to_string(ret));

    if (mat == 101) {
        // reorder for mue
        auto arr = vector<float>(sz);
        int nxy = nx * ny, nxyz = nxy * nz;
        for (int i = 0; i < nx; i ++) {
            for (int j = 0; j < ny; j ++) {
                for (int k = 0; k < nz; k ++) {
                    arr[i + j * nx + k * nxy         ] = array[((i+1) % nx) + j * nx + k * nxy         ];
                    arr[i + j * nx + k * nxy + nxyz  ] = array[i + ((j+1) % ny) * nx + k * nxy + nxyz  ];
                    arr[i + j * nx + k * nxy + nxyz*2] = array[i + j * nx + ((k+1) % nz) * nxy + nxyz*2];
                }
            }
        }
        return arr;
    } else {
        return array;
    }
}

vector<double> cst::Project::Get1DResult(string tree, int num, int type) {
    ASSERT(dll != NULL, "project " + path + " already destroyed");

    auto getResultSize = (CST_Get1DResultSize_PTR) dll->getProc("CST_Get1DResultSize");
    int size = 0;
    auto ret = getResultSize ? getResultSize(&handle, tree.c_str(), num, &size) : -1;
    ASSERT(ret == 0, "Get1DResultSize Error: " + to_string(ret));

    auto getResultData = type == 0 ?
        (CST_Get1DRealDataAbszissa_PTR) dll->getProc("CST_Get1DRealDataAbszissa") :
        (CST_Get1DRealDataOrdinate_PTR) dll->getProc("CST_Get1DRealDataOrdinate");
    auto array = vector<double>(size);
    ret = getResultData ? getResultData(&handle, tree.c_str(), num, array.data()) : -1;
    ASSERT(ret == 0, "Get1DResultData Error: " + to_string(ret));

    return array;
}
