#include "fit.h"
#include "kernel.h"

namespace solver {
    using namespace std;

    class Solver {
    public:
        Solver(fit::Matrix &mat, float dt, vector<fit::Port> &ports);
        ~Solver();
        fit::Coefficient *coe;
        virtual float Step(float s);
    };

    class CPUSolver : public Solver {
        grid::Grid *grid;
        int idx, dir;
        float *E, *H,
            *Ex, *Ey, *Ez, *Hx, *Hy, *Hz,
            *LEx, *LEy, *LEz, *LHx, *LHy, *LHz,
            *REx, *REy, *REz, *RHx, *RHy, *RHz;
    public:
        CPUSolver(fit::Matrix &mats, float dt, vector<fit::Port> &ports);
        ~CPUSolver();
        virtual float Step(float s);
    };

    class GPUSolver : public Solver {
        utils::DLL *dll;
        decltype(&init_$i) FnInit;
        decltype(&step_$i) FnStep;
        decltype(&quit_$i) FnQuit;
    public:
        GPUSolver(fit::Matrix &mats, float dt, vector<fit::Port> &ports);
        ~GPUSolver();
        virtual float Step(float s);
    };
}
