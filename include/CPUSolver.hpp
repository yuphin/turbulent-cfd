#pragma once
#include "Solver.hpp"

struct CPUSolver : public Solver {
    CPUSolver() = default;
    virtual void solve_pre_pressure(Real &dt) override;
    virtual void solve_pressure(Real &res, uint32_t  &it) override;
    virtual void solve_post_pressure() override;
};