package com.quantego.josqp;

import com.quantego.josqp.OSQP.Status;

public class DataGenerator {


    public static class basic_qp_sols_data{
        final double []x_test;
        final double []y_test;
        final double obj_value_test;
        final Status status_test;
        final double []q_new;
        final double []l_new;
        final double []u_new;

        public double[] getX_test() {
            return x_test;
        }

        public double[] getY_test() {
            return y_test;
        }

        public double getObj_value_test() {
            return obj_value_test;
        }

        public Status getStatus_test() {
            return status_test;
        }

        public double[] getQ_new() {
            return q_new;
        }

        public double[] getL_new() {
            return l_new;
        }

        public double[] getU_new() {
            return u_new;
        }

        public basic_qp_sols_data(double[] x_test, double[] y_test, double obj_value_test, Status status_test, double[] q_new,
                                  double[] l_new, double[] u_new) {
            super();
            this.x_test = x_test;
            this.y_test = y_test;
            this.obj_value_test = obj_value_test;
            this.status_test = status_test;
            this.q_new = q_new;
            this.l_new = l_new;
            this.u_new = u_new;
        }
    }


    public static OSQP.Data generate_problem_basic_qp(){
        CSCMatrix A = new CSCMatrix(4,2,5,new int[3],new int[5], new double[5]);

        A.Ax[0] = 1.00000000000000000000;
        A.Ax[1] = 1.00000000000000000000;
        A.Ax[2] = 1.00000000000000000000;
        A.Ax[3] = 1.00000000000000000000;
        A.Ax[4] = 1.00000000000000000000;

        A.Ai[0] = 0;
        A.Ai[1] = 1;
        A.Ai[2] = 0;
        A.Ai[3] = 2;
        A.Ai[4] = 3;

        A.Ap[0] = 0;
        A.Ap[1] = 2;
        A.Ap[2] = 5;

        CSCMatrix P = new CSCMatrix(2,2,3,new int[3],new int[3], new double[3]);

        P.Ax[0] = 4.00000000000000000000;
        P.Ax[1] = 1.00000000000000000000;
        P.Ax[2] = 2.00000000000000000000;

        P.Ai[0] = 0;
        P.Ai[1] = 0;
        P.Ai[2] = 1;

        P.Ap[0] = 0;
        P.Ap[1] = 1;
        P.Ap[2] = 3;

        OSQP.Data data = new OSQP.Data(2,4,P,A,new double[2],new double[4],new double[4]);

        data.l[0] = 1.00000000000000000000;
        data.l[1] = 0.00000000000000000000;
        data.l[2] = 0.00000000000000000000;
        data.l[3] = -OSQP.OSQP_INFTY;

        data.u[0] = 1.00000000000000000000;
        data.u[1] = 0.69999999999999995559;
        data.u[2] = 0.69999999999999995559;
        data.u[3] = OSQP.OSQP_INFTY;

        data.q[0] = 1.00000000000000000000;
        data.q[1] = 1.00000000000000000000;

        return data;
    }

    public static basic_qp_sols_data generate_problem_basic_qp_sols_data() {
        basic_qp_sols_data qp_sols_data = new basic_qp_sols_data(new double[2], new double[4], 1.87999999999999989342,OSQP.Status.SOLVED, new double[2], new double[4], new double[4]);

        qp_sols_data.x_test[0] = 0.29999999999999998890;
        qp_sols_data.x_test[1] = 0.69999999999999995559;

        qp_sols_data.y_test[0] = -2.89999999999999991118;
        qp_sols_data.y_test[1] = 0.00000000000000000000;
        qp_sols_data.y_test[2] = 0.20000000000000001110;
        qp_sols_data.y_test[3] = 0.00000000000000000000;

        qp_sols_data.q_new[0] = 2.50000000000000000000;
        qp_sols_data.q_new[1] = 3.20000000000000017764;

        qp_sols_data.l_new[0] = 0.80000000000000004441;
        qp_sols_data.l_new[1] = -3.39999999999999991118;
        qp_sols_data.l_new[2] = -OSQP.OSQP_INFTY;
        qp_sols_data.l_new[3] = 0.50000000000000000000;

        qp_sols_data.u_new[0] = 1.60000000000000008882;
        qp_sols_data.u_new[1] = 1.00000000000000000000;
        qp_sols_data.u_new[2] = OSQP.OSQP_INFTY;
        qp_sols_data.u_new[3] = 0.50000000000000000000;

        return qp_sols_data;

    }

}
