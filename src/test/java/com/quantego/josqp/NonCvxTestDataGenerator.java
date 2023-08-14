package com.quantego.josqp;

import com.quantego.josqp.OSQP.Status;

public class NonCvxTestDataGenerator {


    public static class non_cvx_sols_data{

        public int sigma_new;
        public non_cvx_sols_data(int sigma){
            this.sigma_new = sigma;
        };
    }

    public static OSQP.Data  generate_problem_non_cvx(){


        // Matrix A
        //---------
        CSCMatrix A = new CSCMatrix(5,2,8,new int[3],new int[8], new double[8]);
        A.Ax[0] = -1.00000000000000000000;
        A.Ax[1] = -1.00000000000000000000;
        A.Ax[2] = 2.00000000000000000000;
        A.Ax[3] = 3.00000000000000000000;
        A.Ax[4] = -1.00000000000000000000;
        A.Ax[5] = 3.00000000000000000000;
        A.Ax[6] = 5.00000000000000000000;
        A.Ax[7] = 4.00000000000000000000;
        //A.i = (c_int*) c_malloc(8 * sizeof(c_int));
        A.Ai[0] = 0;
        A.Ai[1] = 2;
        A.Ai[2] = 3;
        A.Ai[3] = 4;
        A.Ai[4] = 1;
        A.Ai[5] = 2;
        A.Ai[6] = 3;
        A.Ai[7] = 4;
        //A.p = (c_int*) c_malloc((2 + 1) * sizeof(c_int));
        A.Ap[0] = 0;
        A.Ap[1] = 4;
        A.Ap[2] = 8;


        // Matrix P
        //---------
        CSCMatrix P = new CSCMatrix(2,2,2,new int[3],new int[3], new double[3]);
        P.Ax[0] = 2.00000000000000000000;
        P.Ax[1] = 5.00000000000000000000;
        P.Ax[2] = 1.00000000000000000000;
        P.Ai[0] = 0;
        P.Ai[1] = 0;
        P.Ai[2] = 1;
        P.Ap[0] = 0;
        P.Ap[1] = 1;
        P.Ap[2] = 3;

        OSQP.Data data = new OSQP.Data(2,5,P,A,new double[2],new double[5],new double[5], 0);

        // Problem vectors
        data.l[0] = -OSQP.OSQP_INFTY;
        data.l[1] = -OSQP.OSQP_INFTY;
        data.l[2] = -OSQP.OSQP_INFTY;
        data.l[3] = -OSQP.OSQP_INFTY;
        data.l[4] = -OSQP.OSQP_INFTY;
        data.u[0] = 0.00000000000000000000;
        data.u[1] = 0.00000000000000000000;
        data.u[2] = -15.00000000000000000000;
        data.u[3] = 100.00000000000000000000;
        data.u[4] = 80.00000000000000000000;
        data.q[0] = 3.00000000000000000000;
        data.q[1] = 4.00000000000000000000;

        return data;

    }


    public static non_cvx_sols_data   generate_problem_non_cvx_sols_data(){

        non_cvx_sols_data  data = new non_cvx_sols_data(5);

        return data;

    };
}