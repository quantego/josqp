package com.quantego.josqp;

final class BasicQPTestSolsData {
    final double[] x_test;
    final double obj_value_test;
    final double[] q_new;
    final OSQP.Status status_test;
    final double[] l_new;
    final double[] y_test;
    final double[] u_new;

    BasicQPTestSolsData(double[] x_test, double obj_value_test, double[] q_new,
                        OSQP.Status status_test, double[] l_new, double[] y_test, double[] u_new) {
        this.x_test = x_test;
        this.obj_value_test = obj_value_test;
        this.q_new = q_new;
        this.status_test = status_test;
        this.l_new = l_new;
        this.y_test = y_test;
        this.u_new = u_new;
    }
}