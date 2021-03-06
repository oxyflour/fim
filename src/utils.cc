#include <fstream>
#include <sstream>
#include <codecvt>
#include <vector>

#include "cuda.h"
#include "utils.h"
#include "utils/check.h"

using namespace std;

string utils::random(const int len) {
    string tmp_s;
    static const char alphanum[] =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";
    srand((unsigned) time(NULL));
    tmp_s.reserve(len);
    for (int i = 0; i < len; ++i) 
        tmp_s += alphanum[rand() % (sizeof(alphanum) - 1)];
    return tmp_s;
}

void utils::outputProgress(double percentage) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}

// https://stackoverflow.com/questions/2573834/c-convert-string-or-char-to-wstring-or-wchar-t
wstring utils::utf8ToWstring(const string& str) {
    wstring_convert<codecvt_utf8<wchar_t>> myconv;
    return myconv.from_bytes(str);
}

string utils::wstringToUtf8(const wstring& str) {
    wstring_convert<codecvt_utf8<wchar_t>> myconv;
    return myconv.to_bytes(str);
}

string utils::readFile(const string& fname) {
    ifstream stream(fname);
    stringstream content;
    content << stream.rdbuf();
    stream.close();
    return content.str();
}

void utils::writeFile(const string& fname, const string& content) {
    ofstream stream(fname);
    stream << content;
    stream.close();
}

int binSearch(const vector<float> &xs, const float x) {
    if (x < xs[0]) {
        return -1;
    }
    int i = 0, j = xs.size();
    while (i < j - 1) {
        int k = (i + j) / 2;
        if (xs[k] < x) {
            i = k;
        } else {
            j = k;
        }
    }
    return i;
}

float utils::interp1(const vector<float> &xs, const vector<float> &ys, const float x) {
    auto i = binSearch(xs, x);
    if (i < 0) {
        return ys[0];
    } else if (i >= xs.size() - 1) {
        return ys[ys.size() - 1];
    } else {
        auto dx = xs[i + 1] - xs[i];
        return ys[i] * (xs[i + 1] - x) / dx + ys[i + 1] * (x - xs[i]) / dx;
    }
}

utils::DLL::DLL(string path) {
    this->path = path;
    module = dlopen(path.c_str());
#ifdef _WIN32
    auto code = GetLastError();
#else
    auto code = -1;
#endif
    CHECK(module != NULL, "LoadLibrary Error: " + path + " (" + to_string(code) + ")");
}

utils::DLL::~DLL() {
    dlclose(module);
}

lib_func __stdcall utils::DLL::getProc(string name) {
    if (procs.count(name) == 1) {
        return procs[name];
    }
    return procs[name] = dlsym(module, name.c_str());
}

chrono::system_clock::time_point utils::clockNow() {
    return chrono::system_clock::now();
}

float utils::secondsSince(chrono::system_clock::time_point &start) {
    return chrono::duration<float>(chrono::system_clock::now() - start).count();
}
