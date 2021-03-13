// https://github.com/lava/matplotlib-cpp
#define WITHOUT_NUMPY
// https://bugs.python.org/issue38728
#ifdef _WIN32
    #ifdef _DEBUG
        #undef _DEBUG
        #include "matplotlibcpp.h"
        #define _DEBUG
    #else
        #include "matplotlibcpp.h"
    #endif
#endif

namespace plt = matplotlibcpp;
