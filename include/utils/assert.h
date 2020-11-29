#include <string>
#include <stdio.h>

#ifndef ASSERT_H
#define ASSERT_H

#define ASSERT(ok, msg) do { if (!(ok)) _panic(msg, __FILE__, __LINE__); } while (0);
inline void _panic(std::string msg, const char *file, int line, bool abort=true)
{
    fprintf(stderr, "FATAL: %s at %s:%d\n", msg.c_str(), file, line);
    if (abort) {
        exit(1);
    }
}

#endif
