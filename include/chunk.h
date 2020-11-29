#ifndef CHUNK_H
#define CHUNK_H

extern "C" int init_$i(float *le, float *re, float *lh, float *rh);
extern "C" float step_$i(float s);
extern "C" int quit_$i();

#endif
