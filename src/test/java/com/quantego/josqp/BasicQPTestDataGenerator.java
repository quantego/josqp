package com.quantego.josqp;

final class BasicQPTestDataGenerator {

    static BasicQPTestSolsData generateData() {

        final double[] x_test = new double[] { 0.29999999999999998890, 0.69999999999999995559 };
        final double obj_value_test = 1.87999999999999989342;
        final double[] q_new = new double[] { 2.50000000000000000000, 3.20000000000000017764 };
        final OSQP.Status status_test = OSQP.Status.SOLVED;
        final double[] l_new = new double[] { 0.80000000000000004441, -3.39999999999999991118,
                -OSQP.OSQP_INFTY, 0.5 };
        final double[] y_test = new double[] { -2.89999999999999991118, 0.0, 0.20000000000000001110,
                0.0 };
        final double[] u_new = new double[] { 1.60000000000000008882, 1.0, OSQP.OSQP_INFTY, 0.5 };
        return new BasicQPTestSolsData(x_test, obj_value_test, q_new, status_test, l_new, y_test,
                u_new);
    }

    static OSQP.Data generateProblem() {

        // Problem dimensions
        final int n = 2;
        final int m = 4;
        // Problem vectors
        final double[] l = new double[] { 1.0, 0.0, 0.0, -OSQP.OSQP_INFTY };
        final double[] u = new double[] { 1.0, 0.69999999999999995559, 0.69999999999999995559,
                OSQP.OSQP_INFTY };
        final double[] q = new double[] { 1.0, 1.0 };

        // Matrix A
        // ---------
        final int A_m = 4;
        final int A_n = 2;
        final int A_nzmax = 5;
        final double[] A_x = new double[] { 1.0, 1.0, 1.0, 1.0, 1.0 };
        final int[] A_i = new int[] { 0, 1, 0, 2, 3 };
        final int[] A_p = new int[] { 0, 2, 5 };
        final CSCMatrix A = new CSCMatrix(A_m, A_n, A_nzmax, A_p, A_i, A_x);

        // Matrix P
        // ---------
        final int P_m = 2;
        final int P_n = 2;
        final int P_nzmax = 3;
        final double[] P_x = new double[] { 4.0, 1.0, 2.0 };
        final int[] P_i = new int[] { 0, 0, 1 };
        final int[] P_p = new int[] { 0, 1, 3 };
        final CSCMatrix P = new CSCMatrix(P_m, P_n, P_nzmax, P_p, P_i, P_x);

        return new OSQP.Data(n, m, P, A, q, l, u);
    }

}
