package com.quantego.josqp;

final class PrimalDualInfeasibilityTestSolsData {
    final OSQP.Status status4;
    final OSQP.Status status1;
    final OSQP.Status status3;
    final CSCMatrix A12;
    final double[] u4;
    final CSCMatrix A34;
    final double[] u1;
    final double[] l;
    final double[] u3;
    final double[] u2;
    final double[] q;
    final CSCMatrix P;
    final double obj_value1;
    final double[] y1;
    final double[] x1;
    final OSQP.Status status2;

    PrimalDualInfeasibilityTestSolsData(OSQP.Status status4, OSQP.Status status1,
                                        OSQP.Status status3, CSCMatrix A12, double[] u4, CSCMatrix A34, double[] u1, double[] l,
                                        double[] u3, double[] u2, double[] q, CSCMatrix P, double obj_value1, double[] y1,
                                        double[] x1, OSQP.Status status2) {
        this.status4 = status4;
        this.status1 = status1;
        this.status3 = status3;
        this.A12 = A12;
        this.u4 = u4;
        this.A34 = A34;
        this.u1 = u1;
        this.l = l;
        this.u3 = u3;
        this.u2 = u2;
        this.q = q;
        this.P = P;
        this.obj_value1 = obj_value1;
        this.y1 = y1;
        this.x1 = x1;
        this.status2 = status2;
    }
}
