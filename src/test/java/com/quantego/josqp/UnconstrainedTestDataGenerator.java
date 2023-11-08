package com.quantego.josqp;

final class UnconstrainedTestDataGenerator {
    static OSQP.Data generateProblem() {
        final int n = 5;
        final int m = 0;

        final double[] l = new double[0];
        final double[] u = new double[0];
        final double[] q = new double[] { -1.10593508000000007030, -1.65451545000000010965,
                -2.36346860000000003055, 1.13534534999999991989, -1.01701413999999989990, };


        final int A_m = 0;
        final int A_n = 5;
        final int A_nzmax = 0;
        final double[] A_x = null;
        final int[] A_i = null;
        final int[] A_p = new int[] { 0, 0, 0, 0, 0, 0, 0 };
        final CSCMatrix A = new CSCMatrix(A_m, A_n, A_nzmax, A_p, A_i, A_x);


        final int P_m = 5;
        final int P_n = 5;
        final int P_nzmax = 5;
        final double[] P_x = { 0.61702199999999995939, 0.92032449000000005057,
                0.20011437000000001363, 0.50233256999999997827, 0.34675589000000001105, };
        final int[] P_i = new int[] { 0, 1, 2, 3, 4 };
        final int[] P_p = new int[] { 0, 1, 2, 3, 4, 5 };
        final CSCMatrix P = new CSCMatrix(P_m, P_n, P_nzmax, P_p, P_i, P_x);

        return new OSQP.Data(n, m, P, A, q, l, u,0);
    }

    static UnconstrainedTestSolsData generateData() {
        final double[] x_test = new double[] { 1.79237541999999994147, 1.79775228000000009132,
                11.81058885000000024945, -2.26014677999999991087, 2.93293975000000006759, };
        final double obj_value_test = -19.20975202681327687060;
        final OSQP.Status status_test = OSQP.Status.SOLVED;
        return new UnconstrainedTestSolsData(x_test, obj_value_test, status_test);
    }

}
