package com.quantego.josqp;

final class SolveLinSysTestDataGenerator {
    static SolveLinSysTestSolsData generateData() {

        final double[] test_solve_KKT_rhs = new double[] { -1.05436800185014400988,
                -1.00889150318170695009, -0.06752199193321305193, -0.08142504953145184021,
                2.83521598124578755318, -0.24617516933298103088, 0.70551086502222482011, };
        final double[] test_solve_KKT_x = new double[] { 1.59146823014194582768,
                -0.39979419212238193060, -0.51903227424323850059, 0.55635674324897865795,
                0.86927159669526288255, -0.38765340434320433305, 0.0 };
        final int test_solve_KKT_n = 3;
        final int test_solve_KKT_m = 4;


        final int test_solve_KKT_Pu_m = 3;
        final int test_solve_KKT_Pu_n = 3;
        final int test_solve_KKT_Pu_nzmax = 2;
        final double[] test_solve_KKT_Pu_x = new double[] { 0.38349652976488485256,
                0.15100215738121633424 };
        final int[] test_solve_KKT_Pu_i = new int[] { 0, 1 };
        final int[] test_solve_KKT_Pu_p = new int[] { 0, 1, 2, 2 };
        final CSCMatrix test_solve_KKT_Pu = new CSCMatrix(test_solve_KKT_Pu_m, test_solve_KKT_Pu_n,
                test_solve_KKT_Pu_nzmax, test_solve_KKT_Pu_p, test_solve_KKT_Pu_i,
                test_solve_KKT_Pu_x);


        final int test_solve_KKT_A_m = 4;
        final int test_solve_KKT_A_n = 3;
        final int test_solve_KKT_A_nzmax = 4;
        final double[] test_solve_KKT_A_x = new double[] { 0.40730783228994515976,
                0.54620731990215798390, 0.96963240582679277590, 0.17698462366793654699 };
        final int[] test_solve_KKT_A_i = new int[] { 0, 1, 2, 0, };
        final int[] test_solve_KKT_A_p = new int[] { 0, 2, 3, 4, };
        final CSCMatrix test_solve_KKT_A = new CSCMatrix(test_solve_KKT_A_m, test_solve_KKT_A_n,
                test_solve_KKT_A_nzmax, test_solve_KKT_A_p, test_solve_KKT_A_i, test_solve_KKT_A_x);

        final double test_solve_KKT_rho = 4.00000000000000000000;


        final int test_solve_KKT_KKT_m = 7;
        final int test_solve_KKT_KKT_n = 7;
        final int test_solve_KKT_KKT_nzmax = 15;
        final double[] test_solve_KKT_KKT_x = new double[] { 1.38349652976488490808,
                0.40730783228994515976, 0.54620731990215798390, 1.15100215738121636200,
                0.96963240582679277590, 1.0, 0.17698462366793654699, 0.40730783228994515976,
                0.17698462366793654699, -0.25, 0.54620731990215798390, -0.25,
                0.96963240582679277590, -0.25, -0.25 };
        final int[] test_solve_KKT_KKT_i = new int[] { 0, 3, 4, 1, 5, 2, 3, 0, 2, 3, 0, 4, 1, 5,
                6 };
        final int[] test_solve_KKT_KKT_p = new int[] { 0, 3, 5, 7, 10, 12, 14, 15 };
        final CSCMatrix test_solve_KKT_KKT = new CSCMatrix(test_solve_KKT_KKT_m,
                test_solve_KKT_KKT_n, test_solve_KKT_KKT_nzmax, test_solve_KKT_KKT_p,
                test_solve_KKT_KKT_i, test_solve_KKT_KKT_x);

        final double test_solve_KKT_sigma = 1.0;

        return new SolveLinSysTestSolsData(test_solve_KKT_rhs, test_solve_KKT_x, test_solve_KKT_n,
                test_solve_KKT_m, test_solve_KKT_Pu, test_solve_KKT_A, test_solve_KKT_rho,
                test_solve_KKT_KKT, test_solve_KKT_sigma);
    }
}
