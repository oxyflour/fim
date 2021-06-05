#include <string>
#include <map>
#include <vector>
#include <chrono>

#ifndef UTILS_H
#define UTILS_H

namespace utils {
    using namespace std;
    wstring utf8ToWstring(const string& str);
    string wstringToUtf8(const wstring& str);
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

    string random(const int len);

    wstring utf8ToWstring(const string& str);
    string wstringToUtf8(const wstring& str);

    string readFile(const string& fname);
    void writeFile(const string& fname, const string& content);

    float interp1(const vector<float> &xs, const vector<float> &ys, const float x);

    template <typename T> vector<T> range(T from, T to, T delta) {
        auto ret = vector<T>();
        for (auto val = from; val < to; val += delta) {
            ret.push_back(val);
        }
        return ret;
    }
    template <typename T> vector<T> inline range(T from, T to) {
        return range(from, to, 1);
    }
    template <typename T> vector<T> inline range(T to) {
        return range(0, to, 1);
    }

    typedef struct DLL {
        DLL(string path);
        ~DLL();
        lib_type module;
        string path;
        map<string, FARPROC> procs;
        lib_func __stdcall getProc(string name);
    } DLL;

    template <class T, class F> auto fmap(vector<T> &arr, F fun) {
        auto ret = vector<typename result_of<F(const typename T&)>::type>(arr.size());
        for (const auto &item : arr) {
            ret.push_back(fun(item));
        }
        return ret;
    };

    chrono::system_clock::time_point clockNow();
    float secondsSince(chrono::system_clock::time_point &start);
}

#endif
