package com.quantego.josqp;

import com.quantego.josqp.QDLDL;

public class LinAlg {
    /* VECTOR FUNCTIONS ----------------------------------------------------------*/
    public static void vec_add_scaled(float[] a, float[] b, float[] c, float sc){
        if(!(a.length == b.length && b.length == c.length)){
            throw new IllegalArgumentException("a,b,c need to be of the same length");
        }
        for (int i=0; i<a.length;i++){
            c[i] = a[i] + sc*b[i];
        }
    }

    public static float vec_scaled_norm_inf(float[] S, float[] v){
        if(S.length != v.length){
            throw new IllegalArgumentException("S and v need to be of the same length");
        }
        float max = 0;
        for (int i =0;i<v.length;i++){
            float abs_Sv_i = Math.abs(S[i] * v[i]);
            if(abs_Sv_i > max){
                max = abs_Sv_i;
            }
        }
        return max;
    }

    public static float vec_norm_inf(float [] v){
        float max = 0;
        for (int i = 0; i<v.length; i++){
            float abs_v_i = Math.abs(v[i]);
            if(abs_v_i > max){
                max = abs_v_i;
            }
        }
        return max;
    }

    public static float vec_norm_inf_diff(float[] a, float[] b){
        if(!(a.length == b.length)){
            throw new IllegalArgumentException("a,b need to be of the same length");
        }
        float nmDiff = 0;
        float tmp;
        for (int i = 0; i<a.length; i++){
            tmp = Math.abs(a[i]-b[i]);
            if(tmp > nmDiff){
                nmDiff = tmp;
            }
        }
        return nmDiff;
    }

    public static float vec_mean(float[] a){
        float mean = 0;
        for(int i = 0; i<a.length; i++) {
            mean += a[i];
        }
        return mean/a.length;
    }

    public static void int_vec_set_scalar(float[] a , int sc) {
        for (int i = 0; i < a.length; i++) {
            a[i] = sc;
        }
    }

    public static void vec_set_scalar(float[] a, float sc) {
        for (int i = 0; i < a.length; i++) {
            a[i] = sc;
        }
    }


    public static void vec_add_scalar(float[] a, float sc) {
        for (int i = 0; i < a.length; i++) {
            a[i] += sc;
        }
    }

    public static void vec_mult_scalar(float[] a, float sc) {
        for (int i = 0; i < a.length; i++) {
            a[i] *= sc;
        }
    }
    public static float[] vec_copy(float[] a){
    float[] b = new float[a.length];
    for (int i = 0; i<a.length; i++) b[i] = a[i];
        return b;
    }

    public static void prea_int_vec_copy(int[] a) {
        int[] b = new int[a.length];
        for (int i = 0; i < a.length; i++) b[i] = a[i];
    }

    public static void prea_vec_copy(float[] a) {
        float[] b = new float[a.length];
        for (int i = 0; i < a.length; i++) b[i] = a[i];
    }

    public static void vec_ew_recipr(float[] a , float[] b) {
        for (int i = 0; i < a.length; i++) {
            b[i] = (float) 1.0 / a[i];
        }
    }

    public static float vec_prod(float[] a, float[] b){
        float prod = 0;
        for(int i = 0; i<a.length;i++){
            prod+=a[i]*b[i];
        }

        return prod;
    }

    public static void vec_ew_prod(float[] a , float[] b, float[] c) {
        for (int i = 0; i < a.length; i++) {
            c[i] = b[i] * a[i];
        }
    }

    public static void vec_ew_sqrt(float[] a) {
        for (int i = 0; i < a.length; i++) {
            a[i] = (float) Math.sqrt(a[i]);
        }
    }

    public static void vec_ew_max(float[] a, float max_val) {
        for (int i = 0; i < a.length; i++) {
            a[i] = Math.max(a[i], max_val);
        }
    }

    public static void vec_ew_min(float[] a, float min_val) {
        for (int i = 0; i < a.length; i++) {
            a[i] = Math.min(a[i], min_val);
        }
    }

    public static void vec_ew_max_vec(float[] a, float[] b, float[] c) {
        for (int i = 0; i < a.length; i++) {
            c[i] = Math.max(a[i], b[i]);
        }
    }

    public static void vec_ew_min_vec(float[] a, float[] b, float[] c) {
        for (int i = 0; i < a.length; i++) {
            c[i] = Math.min(a[i], b[i]);
        }
    }
    /* VECTOR FUNCTIONS ----------------------------------------------------------*/



    /* MATRIX FUNCTIONS ----------------------------------------------------------*/
    public static void mat_mult_scalar(QDLDL.SparseMatrix A, float sc) {
        int nnzA = A.Ap[A.n];
        for (int i = 0; i < nnzA; i++) {
            A.Ax[i] *= sc;
        }
    }

    public static void mat_premult_diag(QDLDL.SparseMatrix A, float[] d) {
        for (int j = 0; j < A.n; j++) {                // Cycle over columns
            for (int i = A.Ap[j]; i < A.Ap[j + 1]; i++) { // Cycle every row in the column
                A.Ax[i] *= d[A.Ai[i]];                  // Scale by corresponding element
                // of d for row i
            }
        }
    }

    public static void mat_postmult_diag(QDLDL.SparseMatrix A, float[] d) {
        for (int j = 0; j < A.n; j++) {                // Cycle over columns j
            for (int i = A.Ap[j]; i < A.Ap[j + 1]; i++) { // Cycle every row i in column j
                A.Ax[i] *= d[j];                        // Scale by corresponding element
                // of d for column j
            }
        }
    }

    void mat_vec(QDLDL.SparseMatrix A, float[] x, float[] y, boolean plus_eq) {
        //TODO: Check m,n
        if (!plus_eq) {
            // y = 0
            for (int i = 0; i < A.n; i++) {
                y[i] = 0;
            }
        }

        // if A is empty
        if (A.Ap[A.n] == 0) {
            return;
        }

        if (plus_eq) {
            // y -=  A*x
            for (int j = 0; j < A.n; j++) {
                for (int i = A.Ap[j]; i < A.Ap[j + 1]; i++) {
                    y[A.Ai[i]] -= A.Ax[i] * x[j];
                }
            }
        } else {
            // y +=  A*x
            for (int j = 0; j < A.n; j++) {
                for (int i = A.Ap[j]; i < A.Ap[j + 1]; i++) {
                    y[A.Ai[i]] += A.Ax[i] * x[j];
                }
            }
        }
    }

    public static void mat_tpose_vec(QDLDL.SparseMatrix A, float[] x, float[] y, boolean plus_eq, boolean skip_diag) {

        if (!plus_eq) {
            // y = 0
            for (int i = 0; i < A.n; i++) {
                y[i] = 0;
            }
        }

        // if A is empty
        if (A.Ap[A.n] == 0) {
            return;
        }

        if (plus_eq) {
            // y -=  A*x
            if (skip_diag) {
                for (int j = 0; j < A.n; j++) {
                    for (int k = A.Ap[j]; k < A.Ap[j + 1]; k++) {
                        int i  = A.Ai[k];
                        y[j] -= i == j ? 0 : A.Ax[k] * x[i];
                    }
                }
            } else {
                for (int j = 0; j < A.n; j++) {
                    for (int k = A.Ap[j]; k < A.Ap[j + 1]; k++) {
                        y[j] -= A.Ax[k] * x[A.Ai[k]];
                    }
                }
            }
        } else {
            // y +=  A*x
            if (skip_diag) {
                for (int j = 0; j < A.n; j++) {
                    for (int k = A.Ap[j]; k < A.Ap[j + 1]; k++) {
                        int i     = A.Ai[k];
                        y[j] += i == j ? 0 : A.Ax[k] * x[i];
                    }
                }
            } else {
                for (int j = 0; j < A.n; j++) {
                    for (int k = A.Ap[j]; k < A.Ap[j + 1]; k++) {
                        y[j] += A.Ax[k] * x[A.Ai[k]];
                    }
                }
            }
        }
    }

    public static void mat_inf_norm_cols(QDLDL.SparseMatrix M, float[] E) {

        // Initialize zero max elements
        for (int j = 0; j < M.n; j++) {
            E[j] = 0;
        }

        // Compute maximum across columns
        for (int j = 0; j < M.n; j++) {
            for (int ptr = M.Ap[j]; ptr < M.Ap[j + 1]; ptr++) {
                E[j] = (float) Math.max(Math.abs(M.Ax[ptr]), E[j]);
            }
        }
    }

    public static void mat_inf_norm_rows(QDLDL.SparseMatrix M, float[] E) {
        // Initialize zero max elements
        for (int j = 0; j < M.n; j++) {
            E[j] = 0;
        }

        // Compute maximum across rows
        for (int j = 0; j < M.n; j++) {
            for (int ptr = M.Ap[j]; ptr < M.Ap[j + 1]; ptr++) {
                int i    = M.Ai[ptr];
                E[i] = (float) Math.max(Math.abs(M.Ax[ptr]), E[i]);
            }
        }
    }

    public static void mat_inf_norm_cols_sym_triu(QDLDL.SparseMatrix M, float[] E) {

        // Initialize zero max elements
        for (int j = 0; j < M.n; j++) {
            E[j] = 0;
        }

        // Compute maximum across columns
        // Note that element (i, j) contributes to
        // -> Column j (as expected in any matrices)
        // -> Column i (which is equal to row i for symmetric matrices)
        for (int j = 0; j < M.n; j++) {
            for (int ptr = M.Ap[j]; ptr < M.Ap[j + 1]; ptr++) {
                int i     = M.Ai[ptr];
                int abs_x = Math.abs(M.Ai[ptr]);
                E[j]  = Math.max(abs_x, E[j]);

                if (i != j) {
                    E[i] = Math.max(abs_x, E[i]);
                }
            }
        }
    }


    public static float quad_form(QDLDL.SparseMatrix P, float[] x) {
        float quad_form = 0;

        for (int j = 0; j < P.n; j++) {                      // Iterate over columns
            for (int ptr = P.Ap[j]; ptr < P.Ap[j + 1]; ptr++) { // Iterate over rows
                int i = P.Ai[ptr];                                // Row index

                if (i == j) {                                 // Diagonal element
                    quad_form += (float) .5 * P.Ax[ptr] * x[i] * x[i];
                }
                else if (i < j) {                             // Off-diagonal element
                    quad_form += P.Ax[ptr] * x[i] * x[j];
                }
                else {
                    throw new IllegalArgumentException("quad_form matrix is not upper triangular");
                }
            }
        }
        return quad_form;
    }


    /* MATRIX FUNCTIONS ----------------------------------------------------------*/

}