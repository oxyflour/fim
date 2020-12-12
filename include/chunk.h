#include "fit.h"

#ifndef CHUNK_H
#define CHUNK_H

extern "C" int init_$i(fit::Coefficient *coe);
extern "C" float step_$i(float s);
extern "C" int quit_$i();

#endif
