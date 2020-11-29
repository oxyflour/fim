#include "utils.h"
#include "utils/assert.h"

#include <fstream>
#include <sstream>
#include <codecvt>
#include <vector>

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

// https://stackoverflow.com/questions/2573834/c-convert-string-or-char-to-wstring-or-wchar-t
wstring utils::utf8ToWstring(const string& str) {
    wstring_convert<codecvt_utf8<wchar_t>> myconv;
    return myconv.from_bytes(str);
}

string utils::wstringToUtf8(const wstring& str) {
    wstring_convert<codecvt_utf8<wchar_t>> myconv;
    return myconv.to_bytes(str);
}

string utils::dirname(const string& fname) {
     size_t pos = fname.find_last_of("\\/");
     return string::npos == pos ? "" : fname.substr(0, pos);
}

string utils::readFile(const string& fname) {
    ifstream stream(fname);
    stringstream content;
    content << stream.rdbuf();
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
    } else if (i > xs.size() - 1) {
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
    ASSERT(module != NULL, "LoadLibrary Error: " + path + " (" + to_string(code) + ")");
}

lib_func __stdcall utils::DLL::getProc(string name) {
    if (procs.count(name) == 1) {
        return procs[name];
    }
    return procs[name] = dlsym(module, name.c_str());
}
