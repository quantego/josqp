package com.quantego.josqp;
import static org.junit.jupiter.api.Assertions.assertTrue;
import org.junit.jupiter.api.Test;


import com.quantego.josqp.LinAlgDataGenerator.lin_alg_sols_data;
import java.lang.Math;

public class LinAlgTest {
    final static double TESTS_TOL = 0.0001;
    final static int OSQP_TIME_LIMIT_REACHED = -6;
    public double[] csc_to_dns(CSCMatrix matrix)
    {
        int i = 0;
        int j = 0;
        int index = 0;
        double[]A = new double[matrix.m * matrix.n];
        for ( index = 0; index < matrix.Ap[matrix.n]; ++index)
        {
            i = matrix.Ai[index];
            while(matrix.Ap[j+1] <= index)
            {
                j = j+1;
            }
            A[j* (matrix.m) + i] = matrix.Ax[index];
        }
        return A;


    };
    public static Boolean is_eq_csc (CSCMatrix A, CSCMatrix B, double tol)
    {
        int j;
        int i;
        if(A.n != B.n)
        {
            return false;
        }
        for ( j = 0; j < A.n; ++j)
        {
            if(A.Ap[j] != B.Ap[j])
            {
                return false;
            }
            for ( i = A.Ap[j]; i < A.Ap[j+1]; ++i)
            {
                if((A.Ai[i] != B.Ai[i]) || ((Math.abs(A.Ax[i] - B.Ax[i]) > tol)))
                {
                    return false;
                }
            }
        }
        return true;
    };

    @Test
    public void test_constr_sparse_mat() {

        double []Adns; // Conversion to dense matrix
        lin_alg_sols_data data =  LinAlgDataGenerator.generate_problem_lin_alg_sols_data();

        // Convert sparse to dense
        Adns = csc_to_dns(data.test_sp_matrix_A);

        // Compute norm of the elementwise difference with
        assertTrue((LinAlg.vec_norm_inf_diff(Adns, data.test_sp_matrix_Adns,data.test_sp_matrix_A.m*data.test_sp_matrix_A.n) < TESTS_TOL),"Linear algebra tests: error in constructing sparse/dense matrix!");


    };
    @Test
    public void test_vec_operations() {
        double  norm_inf, vecprod; // normInf;
        double []ew_reciprocal;
        double []add_scaled;
        double []vec_ew_max_vec_test, vec_ew_min_vec_test;

        lin_alg_sols_data data = LinAlgDataGenerator.generate_problem_lin_alg_sols_data();


        // Add scaled
        add_scaled = LinAlg.vec_copy(data.test_vec_ops_v1, data.test_vec_ops_n);
        LinAlg.vec_add_scaled(add_scaled,
                add_scaled,
                data.test_vec_ops_v2,
                data.test_vec_ops_n,
                data.test_vec_ops_sc);
        assertTrue(
                LinAlg.vec_norm_inf_diff(add_scaled, data.test_vec_ops_add_scaled,
                        data.test_vec_ops_n) < TESTS_TOL ,"Linear algebra tests: error in vector operation, adding scaled vector"
        );

        // Norm_inf of the difference
        assertTrue(
                Math.abs(LinAlg.vec_norm_inf_diff(data.test_vec_ops_v1,
                        data.test_vec_ops_v2,
                        data.test_vec_ops_n) -
                        data.test_vec_ops_norm_inf_diff) <
                        TESTS_TOL, "Linear algebra tests: error in vector operation, norm_inf of difference");

        // norm_inf
        norm_inf = LinAlg.vec_norm_inf(data.test_vec_ops_v1, data.test_vec_ops_n);
        assertTrue(Math.abs(norm_inf - data.test_vec_ops_norm_inf) < TESTS_TOL,"Linear algebra tests: error in vector operation, norm_inf"
        );

        // Elementwise reciprocal
        ew_reciprocal = new double[data.test_vec_ops_n];
        LinAlg.vec_ew_recipr(data.test_vec_ops_v1, ew_reciprocal, data.test_vec_ops_n);
        assertTrue(
                LinAlg.vec_norm_inf_diff(ew_reciprocal, data.test_vec_ops_ew_reciprocal,data.test_vec_ops_n) < TESTS_TOL, "Linear algebra tests: error in vector operation, elementwise reciprocal"
        );


        // Vector product
        vecprod = LinAlg.vec_prod(data.test_vec_ops_v1,
                data.test_vec_ops_v2,
                data.test_vec_ops_n);
        assertTrue(Math.abs(vecprod - data.test_vec_ops_vec_prod) < TESTS_TOL,"Linear algebra tests: error in vector operation, vector product"
        );

        // Elementwise maximum between two vectors
        vec_ew_max_vec_test = new double[data.test_vec_ops_n];
        LinAlg.vec_ew_max_vec(data.test_vec_ops_v1,
                data.test_vec_ops_v2,
                vec_ew_max_vec_test,
                data.test_vec_ops_n);
        assertTrue(Math.abs(vecprod - data.test_vec_ops_vec_prod) < TESTS_TOL,
                "Linear algebra tests: error in vector operation, elementwise maximum between vectors"
        );
        // Elementwise minimum between two vectors
        vec_ew_min_vec_test = new double[data.test_vec_ops_n ];
        //c_float *)c_malloc(data.test_vec_ops_n * sizeof(c_float));
        LinAlg.vec_ew_min_vec(data.test_vec_ops_v1,
                data.test_vec_ops_v2,
                vec_ew_min_vec_test,
                data.test_vec_ops_n);
        assertTrue(
                LinAlg.vec_norm_inf_diff(vec_ew_min_vec_test, data.test_vec_ops_ew_min_vec,
                        data.test_vec_ops_n) < TESTS_TOL,"Linear algebra tests: error in vector operation, elementwise minimum between vectors"
        );
        return ;
    };
    @Test
    public void test_mat_operations() {
        CSCMatrix Ad, dA; // Matrices used for tests
        // csc *A_ewsq, *A_ewabs;     // Matrices used for tests
        int exitflag = 0;

        // c_float trace, fro_sq;
        double []inf_norm_cols_rows_test;


        lin_alg_sols_data data = LinAlgDataGenerator.generate_problem_lin_alg_sols_data();


        // Copy matrices
        Ad = CSCMatrix.copy_csc_mat(data.test_mat_ops_A);
        dA = CSCMatrix.copy_csc_mat(data.test_mat_ops_A);



        // Premultiply matrix A
        LinAlg.mat_premult_diag(dA, data.test_mat_ops_d);
        assertTrue(
                is_eq_csc(dA, data.test_mat_ops_prem_diag, TESTS_TOL), "Linear algebra tests: error in matrix operation, premultiply diagonal");


        // Postmultiply matrix A
        LinAlg.mat_postmult_diag(Ad, data.test_mat_ops_d);
        assertTrue(

                is_eq_csc(Ad, data.test_mat_ops_postm_diag, TESTS_TOL),"Linear algebra tests: error in matrix operation, postmultiply diagonal");

        // Maximum norm over columns
        inf_norm_cols_rows_test = new double[data.test_mat_ops_n];
        LinAlg.mat_inf_norm_cols(data.test_mat_ops_A, inf_norm_cols_rows_test);
        assertTrue(
                LinAlg.vec_norm_inf_diff(inf_norm_cols_rows_test, data.test_mat_ops_inf_norm_cols,
                        data.test_mat_ops_n) < TESTS_TOL,"Linear algebra tests: error in matrix operation, max norm over columns");

        // Maximum norm over rows
        LinAlg.mat_inf_norm_rows(data.test_mat_ops_A, inf_norm_cols_rows_test);
        assertTrue(
                LinAlg.vec_norm_inf_diff(inf_norm_cols_rows_test,
                        data.test_mat_ops_inf_norm_rows,
                        data.test_mat_ops_n) < TESTS_TOL,"Linear algebra tests: error in matrix operation, max norm over rows");

        return ;
    };
    @Test
    public void test_mat_vec_multiplication() {
        double[]Ax, ATy, Px, Ax_cum, ATy_cum, Px_cum;

        lin_alg_sols_data data = LinAlgDataGenerator.generate_problem_lin_alg_sols_data();


        // Allocate vectors
        Ax  = new double[data.test_mat_vec_m];
        ATy = new double[data.test_mat_vec_n];
        Px  = new double[data.test_mat_vec_n ];


        // Matrix-vector multiplication:  y = Ax
        LinAlg.mat_vec(data.test_mat_vec_A, data.test_mat_vec_x, Ax, 0,0,1);
        assertTrue(
                LinAlg.vec_norm_inf_diff(Ax, data.test_mat_vec_Ax,
                        data.test_mat_vec_m) < TESTS_TOL,"Linear algebra tests: error in matrix-vector operation, matrix-vector multiplication");

        // Cumulative matrix-vector multiplication:  y += Ax
        Ax_cum = LinAlg.vec_copy(data.test_mat_vec_y, data.test_mat_vec_m);
        LinAlg.mat_vec(data.test_mat_vec_A, data.test_mat_vec_x, Ax_cum,0,0,1);
        assertTrue(
                LinAlg.vec_norm_inf_diff(Ax_cum, data.test_mat_vec_Ax_cum,
                        data.test_mat_vec_m) < TESTS_TOL,"Linear algebra tests: error in matrix-vector operation, cumulative matrix-vector multiplication");

        // Matrix-transpose-vector multiplication:  x = A'*y
        LinAlg.mat_tpose_vec(data.test_mat_vec_A, data.test_mat_vec_y, ATy, 0,0,0, false);
        assertTrue(
                LinAlg.vec_norm_inf_diff(ATy, data.test_mat_vec_ATy,
                        data.test_mat_vec_n) < TESTS_TOL, "Linear algebra tests: error in matrix-vector operation, matrix-transpose-vector multiplication");

        // Cumulative matrix-transpose-vector multiplication:  x += A'*y
        ATy_cum = LinAlg.vec_copy(data.test_mat_vec_x, data.test_mat_vec_n);
        LinAlg.mat_tpose_vec(data.test_mat_vec_A, data.test_mat_vec_y, ATy_cum,0,0 ,1, false);
        assertTrue(
                LinAlg.vec_norm_inf_diff(ATy_cum, data.test_mat_vec_ATy_cum,
                        data.test_mat_vec_n) < TESTS_TOL,"Linear algebra tests: error in matrix-vector operation, cumulative matrix-transpose-vector multiplication");

        // Symmetric-matrix-vector multiplication (only upper part is stored)
        LinAlg.mat_vec(data.test_mat_vec_Pu, data.test_mat_vec_x, Px, 0,0,1);          // upper
        // traingular
        // part
        LinAlg.mat_tpose_vec(data.test_mat_vec_Pu, data.test_mat_vec_x, Px,0,0, 1, true); // lower
        // traingular
        // part
        // (without
        // diagonal)
        assertTrue(
                LinAlg.vec_norm_inf_diff(Px, data.test_mat_vec_Px,
                        data.test_mat_vec_n) < TESTS_TOL,"Linear algebra tests: error in matrix-vector operation, symmetric matrix-vector multiplication");


        // Cumulative symmetric-matrix-vector multiplication
        Px_cum = LinAlg.vec_copy(data.test_mat_vec_x, data.test_mat_vec_n);
        LinAlg.mat_vec(data.test_mat_vec_Pu, data.test_mat_vec_x, Px_cum,0,0, 1);          // upper
        // traingular
        // part
        LinAlg.mat_tpose_vec(data.test_mat_vec_Pu, data.test_mat_vec_x, Px_cum, 0,0,1, true); // lower
        // traingular
        // part
        // (without
        // diagonal)
        assertTrue(
                LinAlg.vec_norm_inf_diff(Px_cum, data.test_mat_vec_Px_cum,
                        data.test_mat_vec_n) < TESTS_TOL,"Linear algebra tests: error in matrix-vector operation, cumulative symmetric matrix-vector multiplication");
        return ;
    };
    @Test
    public static void test_extract_upper_triangular() {
        double []inf_norm_cols_test;
        lin_alg_sols_data data = LinAlgDataGenerator.generate_problem_lin_alg_sols_data();

        // Extract upper triangular part
        CSCMatrix Ptriu =  CSCMatrix.csc_to_triu(data.test_mat_extr_triu_P);

        assertTrue(
                is_eq_csc(data.test_mat_extr_triu_Pu, Ptriu, TESTS_TOL),"Linear algebra tests: error in forming upper triangular matrix!");

        // Compute infinity norm over columns of the original matrix by using the
        // upper triangular part only
        inf_norm_cols_test = new double[data.test_mat_extr_triu_n];
        LinAlg.mat_inf_norm_cols_sym_triu(Ptriu, inf_norm_cols_test);
        assertTrue(LinAlg.vec_norm_inf_diff(inf_norm_cols_test,
                data.test_mat_extr_triu_P_inf_norm_cols,
                data.test_mat_extr_triu_n) < TESTS_TOL, "Linear algebra tests: error in forming upper triangular matrix, infinity norm over columns");

        return ;
    }
    @Test
    public void test_quad_form_upper_triang() {
        double quad_form_t;

        lin_alg_sols_data data =  LinAlgDataGenerator.generate_problem_lin_alg_sols_data();

        // Compute quadratic form
        quad_form_t = LinAlg.quad_form(data.test_qpform_Pu, data.test_qpform_x);

        assertTrue(
                (Math.abs(quad_form_t - data.test_qpform_value) < TESTS_TOL),"Linear algebra tests: error in computing quadratic form using upper triangular matrix!");

        // cleanup
        //clean_problem_lin_alg_sols_data(data);

        return ;
    }

};


