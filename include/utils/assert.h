#include <string>
#include <exception>
#include <iostream>

#ifndef ASSERT_H
#define ASSERT_H

#define ASSERT(ok, msg) do { if (!(ok)) _panic(msg, __FILE__, __LINE__); } while (0);
inline void _panic(std::string msg, const char *file, int line)
{
    auto m = "FATAL: " + msg + " at " + file + ":" + std::to_string(line) + "\n";
    throw std::exception(m.c_str());
}

#endif
