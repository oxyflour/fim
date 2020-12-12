#include "fit.h"
#include "chunk.h"

namespace solver {
    typedef struct Solver {
        Solver(fit::Matrix &mat, float dt, std::vector<fit::Port> &ports);
        ~Solver();
        float Step(float s);
        utils::DLL *dll;
        decltype(&init_$i) FnInit;
        decltype(&step_$i) FnStep;
        decltype(&quit_$i) FnQuit;

        // cpu solver
        fit::Coefficient *coe;
        float *E, *H;
    } Solver;
}
