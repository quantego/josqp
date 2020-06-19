package com.quantego.josqp;

final class PrimalDualInfeasibilityTestDataGenerator {

    static PrimalDualInfeasibilityTestSolsData generateData() {
        final OSQP.Status status4 = OSQP.Status.PRIMAL_INFEASIBLE;
        final OSQP.Status status1 = OSQP.Status.SOLVED;
        final OSQP.Status status3 = OSQP.Status.DUAL_INFEASIBLE;


        final int A12_m = 3;
        final int A12_n = 2;
        final int A12_nzmax = 4;
        final double[] A12_x = new double[] { 1.0, 1.0, 1.0, 1.0 };
        final int[] A12_i = new int[] { 0, 1, 0, 2 };
        final int[] A12_p = new int[] { 0, 2, 4 };
        final CSCMatrix A12 = new CSCMatrix(A12_m, A12_n, A12_nzmax, A12_p, A12_i, A12_x);

        final double[] u4 = new double[] { 0.0, 3.0, OSQP.OSQP_INFTY };


        final int A34_m = 3;
        final int A34_n = 2;
        final int A34_nzmax = 3;
        final double[] A34_x = new double[] { 1.0, 1.0, 1.0 };
        final int[] A34_i = new int[] { 0, 1, 2 };
        final int[] A34_p = new int[] { 0, 2, 3 };
        final CSCMatrix A34 = new CSCMatrix(A34_m, A34_n, A34_nzmax, A34_p, A34_i, A34_x);

        final double[] u1 = new double[] { 5.0, 3.0, 3.0 };
        final double[] l = new double[] { 0.0, 1.0, 1.0 };
        final double[] u3 = new double[] { 2.0, 3.0, OSQP.OSQP_INFTY };
        final double[] u2 = new double[] { 0.0, 3.0, 3.0 };
        final double[] q = new double[] { 1.0, -1.0 };


        final int P_m = 2;
        final int P_n = 2;
        final int P_nzmax = 1;
        final double[] P_x = new double[] { 1.0 };
        final int[] P_i = new int[] { 0 };
        final int[] P_p = { 0, 1, 1 };
        final CSCMatrix P = new CSCMatrix(P_m, P_n, P_nzmax, P_p, P_i, P_x);

        final double obj_value1 = -1.5;
        final double[] y1 = new double[] { 0.0, -2.0, 1.0 };
        final double[] x1 = new double[] { 1.0, 3.0 };
        final OSQP.Status status2 = OSQP.Status.PRIMAL_INFEASIBLE;

        return new PrimalDualInfeasibilityTestSolsData(status4, status1, status3, A12, u4, A34, u1,
                l, u3, u2, q, P, obj_value1, y1, x1, status2);
    }
}
