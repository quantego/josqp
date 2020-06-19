package com.quantego.josqp;

final class SolveLinSysTestSolsData {
    final double[] test_solve_KKT_rhs;
    final double[] test_solve_KKT_x;
    final int test_solve_KKT_n;
    final int test_solve_KKT_m;
    final CSCMatrix test_solve_KKT_Pu;
    final CSCMatrix test_solve_KKT_A;
    final double test_solve_KKT_rho;
    final CSCMatrix test_solve_KKT_KKT;
    final double test_solve_KKT_sigma;

    SolveLinSysTestSolsData(double[] test_solve_KKT_rhs, double[] test_solve_KKT_x,
                            int test_solve_KKT_n, int test_solve_KKT_m, CSCMatrix test_solve_KKT_Pu,
                            CSCMatrix test_solve_KKT_A, double test_solve_KKT_rho, CSCMatrix test_solve_KKT_KKT,
                            double test_solve_KKT_sigma) {
        this.test_solve_KKT_rhs = test_solve_KKT_rhs;
        this.test_solve_KKT_x = test_solve_KKT_x;
        this.test_solve_KKT_n = test_solve_KKT_n;
        this.test_solve_KKT_m = test_solve_KKT_m;
        this.test_solve_KKT_Pu = test_solve_KKT_Pu;
        this.test_solve_KKT_A = test_solve_KKT_A;
        this.test_solve_KKT_rho = test_solve_KKT_rho;
        this.test_solve_KKT_KKT = test_solve_KKT_KKT;
        this.test_solve_KKT_sigma = test_solve_KKT_sigma;
    }
}
