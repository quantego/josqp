package com.quantego.josqp;

final class UnconstrainedTestSolsData {
    final double[] x_test;
    final double obj_value_test;
    final OSQP.Status status_test;

    UnconstrainedTestSolsData(double[] x_test, double obj_value_test, OSQP.Status status_test) {
        this.x_test = x_test;
        this.obj_value_test = obj_value_test;
        this.status_test = status_test;
    }
}
