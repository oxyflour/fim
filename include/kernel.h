#include "fit.h"

#ifndef KERNEL_H
#define KERNEL_H

extern "C" int init_$i(fit::Coefficient *coe);
extern "C" float step_$i(float s);
extern "C" int quit_$i();

#endif
