package com.quantego.josqp;

final class UpdateMatricesTestSolsData {
    final CSCMatrix test_form_KKT_A;
    final CSCMatrix test_form_KKT_KKT_new;
    final double[] test_solve_y;
    final double test_form_KKT_sigma;
    final CSCMatrix test_form_KKT_A_new;
    final double[] test_solve_u;
    final double test_solve_obj_value;
    final double[] test_solve_q;
    final CSCMatrix test_form_KKT_KKT;
    final double[] test_solve_l;
    final CSCMatrix test_solve_A_new;
    final CSCMatrix test_solve_Pu;
    final double[] test_solve_P_new_y;
    final double[] test_solve_P_new_x;
    final double[] test_solve_x;
    final double test_solve_P_A_new_obj_value;
    final double test_form_KKT_rho;
    final CSCMatrix test_form_KKT_KKTu;
    final double test_solve_A_new_obj_value;
    final CSCMatrix test_solve_Pu_new;
    final CSCMatrix test_form_KKT_KKTu_new;
    final int m;
    final int test_form_KKT_m;
    final int test_form_KKT_n;
    final int n;
    final OSQP.Status test_solve_status;
    final double test_solve_P_new_obj_value;
    final OSQP.Status test_solve_A_new_status;
    final OSQP.Status test_solve_P_A_new_status;
    final double[] test_solve_A_new_x;
    final CSCMatrix test_form_KKT_Pu;
    final double[] test_solve_A_new_y;
    final OSQP.Status test_solve_P_new_status;
    final double[] test_solve_P_A_new_x;
    final CSCMatrix test_solve_A;
    final double[] test_solve_P_A_new_y;
    final CSCMatrix test_form_KKT_Pu_new;

    UpdateMatricesTestSolsData(
            // @formatter:off
            CSCMatrix test_form_KKT_A,
            CSCMatrix test_form_KKT_KKT_new,
            double[] test_solve_y,
            double test_form_KKT_sigma,
            CSCMatrix test_form_KKT_A_new,
            double[] test_solve_u,
            double test_solve_obj_value,
            double[] test_solve_q,
            CSCMatrix test_form_KKT_KKT,
            double[] test_solve_l,
            CSCMatrix test_solve_A_new,
            CSCMatrix test_solve_Pu,
            double[] test_solve_P_new_y,
            double[] test_solve_P_new_x,
            double[] test_solve_x,
            double test_solve_P_A_new_obj_value,
            double test_form_KKT_rho,
            CSCMatrix test_form_KKT_KKTu,
            double test_solve_A_new_obj_value,
            CSCMatrix test_solve_Pu_new,
            CSCMatrix test_form_KKT_KKTu_new,
            int m,
            int test_form_KKT_m,
            int test_form_KKT_n,
            int n,
            OSQP.Status test_solve_status,
            double test_solve_P_new_obj_value,
            OSQP.Status test_solve_A_new_status,
            OSQP.Status test_solve_P_A_new_status,
            double[] test_solve_A_new_x,
            CSCMatrix test_form_KKT_Pu,
            double[] test_solve_A_new_y,
            OSQP.Status test_solve_P_new_status,
            double[] test_solve_P_A_new_x,
            CSCMatrix test_solve_A,
            double[] test_solve_P_A_new_y,
            CSCMatrix test_form_KKT_Pu_new
    ) {
        this.test_form_KKT_A = test_form_KKT_A;
        this.test_form_KKT_KKT_new = test_form_KKT_KKT_new;
        this.test_solve_y = test_solve_y;
        this.test_form_KKT_sigma = test_form_KKT_sigma;
        this.test_form_KKT_A_new = test_form_KKT_A_new;
        this.test_solve_u = test_solve_u;
        this.test_solve_obj_value = test_solve_obj_value;
        this.test_solve_q = test_solve_q;
        this.test_form_KKT_KKT = test_form_KKT_KKT;
        this.test_solve_l = test_solve_l;
        this.test_solve_A_new = test_solve_A_new;
        this.test_solve_Pu = test_solve_Pu;
        this.test_solve_P_new_y = test_solve_P_new_y;
        this.test_solve_P_new_x = test_solve_P_new_x;
        this.test_solve_x = test_solve_x;
        this.test_solve_P_A_new_obj_value = test_solve_P_A_new_obj_value;
        this.test_form_KKT_rho = test_form_KKT_rho;
        this.test_form_KKT_KKTu = test_form_KKT_KKTu;
        this.test_solve_A_new_obj_value = test_solve_A_new_obj_value;
        this.test_solve_Pu_new = test_solve_Pu_new;
        this.test_form_KKT_KKTu_new = test_form_KKT_KKTu_new;
        this.m = m;
        this.test_form_KKT_m = test_form_KKT_m;
        this.test_form_KKT_n = test_form_KKT_n;
        this.n = n;
        this.test_solve_status = test_solve_status;
        this.test_solve_P_new_obj_value = test_solve_P_new_obj_value;
        this.test_solve_A_new_status = test_solve_A_new_status;
        this.test_solve_P_A_new_status = test_solve_P_A_new_status;
        this.test_solve_A_new_x = test_solve_A_new_x;
        this.test_form_KKT_Pu = test_form_KKT_Pu;
        this.test_solve_A_new_y = test_solve_A_new_y;
        this.test_solve_P_new_status = test_solve_P_new_status;
        this.test_solve_P_A_new_x = test_solve_P_A_new_x;
        this.test_solve_A = test_solve_A;
        this.test_solve_P_A_new_y = test_solve_P_A_new_y;
        this.test_form_KKT_Pu_new = test_form_KKT_Pu_new;
    }
}
