#include <string>
#include <map>
#include <vector>
#include <chrono>

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
inline int dlclose(lib_type mod) {
    return FreeLibrary(mod);
}
#else
typedef void* lib_type;
typedef void* lib_func;
#endif

// https://stackoverflow.com/a/36315819
const auto PBSTR = "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||";
const auto PBWIDTH = 60;

namespace utils {
    void outputProgress(double percentage);

    std::string random(const int len);

    std::wstring utf8ToWstring(const std::string& str);
    std::string wstringToUtf8(const std::wstring& str);

    std::string readFile(const std::string& fname);
    void writeFile(const std::string& fname, const std::string& content);

    float interp1(const std::vector<float> &xs, const std::vector<float> &ys, const float x);

    template <typename T> std::vector<T> range(T from, T to, T delta) {
        auto ret = std::vector<T>();
        for (auto val = from; val < to; val += delta) {
            ret.push_back(val);
        }
        return ret;
    }
    template <typename T> std::vector<T> inline range(T from, T to) {
        return range(from, to, 1);
    }
    template <typename T> std::vector<T> inline range(T to) {
        return range(0, to, 1);
    }

    typedef struct DLL {
        DLL(std::string path);
        ~DLL();
        lib_type module;
        std::string path;
        std::map<std::string, FARPROC> procs;
        lib_func __stdcall getProc(std::string name);
    } DLL;

    // well, it's impossible to declare auto return type of a function now
    template <class T, class F> auto map(std::vector<T> &arr, F fun) {
        auto ret = std::vector<typename std::result_of<F(const typename T&)>::type>(arr.size());
        for (const auto &item : arr) {
            ret.push_back(fun(item));
        }
        return ret;
    };

    std::chrono::system_clock::time_point clockNow();
    float secondsSince(std::chrono::system_clock::time_point &start);
}

#endif
