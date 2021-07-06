#pragma once
#include "Solver.hpp"
#include <lib/pcgsolver.h>
struct CPUSolver : public Solver {
    CPUSolver() = default;
    virtual void initialize() override;
    virtual void solve_pre_pressure(Real &dt) override;
    virtual void solve_pressure(Real &res, uint32_t  &it) override;
    virtual void solve_post_pressure() override;
private:
    Real solve_pcg(uint32_t &it);
    Real solve_sor();
    void build_pcg_matrix();
    SparseMatrix<Real> A;

};