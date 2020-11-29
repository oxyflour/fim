#include <string>
#include <map>
#include <vector>

#ifndef UTILS_H
#define UTILS_H

namespace utils {
    std::wstring utf8ToWstring(const std::string& str);
    std::string wstringToUtf8(const std::wstring& str);
}

#ifdef _WIN32
#include <Windows.h>
typedef HINSTANCE lib_type;
typedef FARPROC lib_func;
inline lib_type dlopen(const char *libname) {
    return LoadLibraryW(utils::utf8ToWstring(libname).c_str());
}
inline lib_func dlsym(lib_type mod, const char *symname) {
    return GetProcAddress(mod, symname);
}
#else
// TODO
#endif

namespace utils {
    std::string random(const int len);

    std::wstring utf8ToWstring(const std::string& str);
    std::string wstringToUtf8(const std::wstring& str);

    std::string dirname(const std::string& fname);

    std::string readFile(const std::string& fname);
    void writeFile(const std::string& fname, const std::string& content);

    float interp1(const std::vector<float> &xs, const std::vector<float> &ys, const float x);

    typedef struct DLL {
        DLL(std::string path);
        lib_type module;
        std::string path;
        std::map<std::string, FARPROC> procs;
        lib_func __stdcall getProc(std::string name);
    } DLL;
}

#endif
