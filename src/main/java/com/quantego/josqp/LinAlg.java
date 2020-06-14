package com.quantego.josqp;


public class LinAlg {
    /* VECTOR FUNCTIONS ----------------------------------------------------------*/
    public static void vec_add_scaled(double[] c, double[] a, double[] b, int n, double sc){

        for (int i=0; i<n;i++){
            c[i] = a[i] + sc*b[i];
        }
    }

    public static double vec_scaled_norm_inf(double[] S, double[] v, int l){

        double max = 0;
        for (int i =0;i<l;i++){
            double abs_Sv_i = Math.abs(S[i] * v[i]);
            if(abs_Sv_i > max){
                max = abs_Sv_i;
            }
        }
        return max;
    }

    public static double vec_norm_inf(double [] v, int l){
        double max = 0;
        for (int i = 0; i<l; i++){
            double abs_v_i = Math.abs(v[i]);
            if(abs_v_i > max){
                max = abs_v_i;
            }
        }
        return max;
    }

    public static double vec_norm_inf_diff(double[] a, double[] b, int l){
        double nmDiff = 0;
        double tmp;
        for (int i = 0; i<l; i++){
            tmp = Math.abs(a[i]-b[i]);
            if(tmp > nmDiff){
                nmDiff = tmp;
            }
        }
        return nmDiff;
    }

    public static double vec_mean(double[] a, int n){
        double mean = 0;
        for(int i = 0; i<n; i++) {
            mean += a[i];
        }
        return mean/n;
    }

    public static void int_vec_set_scalar(int[] a , int sc, int n) {
        for (int i = 0; i < n; i++) {
            a[i] = sc;
        }
    }

    public static void vec_set_scalar(double[] a, double sc, int n) {
        for (int i = 0; i < n; i++) {
            a[i] = sc;
        }
    }


    public static void vec_add_scalar(double[] a, double sc, int n) {
        for (int i = 0; i < n; i++) {
            a[i] += sc;
        }
    }

    public static void vec_mult_scalar(double[] a, double sc, int n) {
        for (int i = 0; i < n; i++) {
            a[i] *= sc;
        }
    }
    /*public static double[] vec_copy(double[] a){
    double[] b = new double[a.length];
    for (int i = 0; i<a.length; i++) b[i] = a[i];
        return b;
    }*/

    public static void prea_int_vec_copy(int[] a, int[] b, int n) {
        for (int i = 0; i < n; i++) b[i] = a[i];
    }

    public static void prea_vec_copy(double[] a, double[] b, int n) {
        for (int i = 0; i < n; i++) b[i] = a[i];
    }

    public static void vec_ew_recipr(double[] a , double[] b, int n) {
        for (int i = 0; i < n; i++) {
            b[i] = 1.0 / a[i];
        }
    }

    public static double vec_prod(double[] a, double[] b, int n){
        double prod = 0;
        for(int i = 0; i<n;i++){
            prod+=a[i]*b[i];
        }

        return prod;
    }

    public static void vec_ew_prod(double[] a , double[] b, double[] c, int n) {
    	try {
	        for (int i = 0; i < n; i++) {
	            c[i] = b[i] * a[i];
	        }
    	} catch(Exception e) {
    		System.out.println();
    	}
    }

    public static void vec_ew_sqrt(double[] a, int n) {
        for (int i = 0; i < n; i++) {
            a[i] = Math.sqrt(a[i]);
        }
    }

    public static void vec_ew_max(double[] a, int n, double max_val) {
        for (int i = 0; i < n; i++) {
            a[i] = Math.max(a[i], max_val);
        }
    }

    public static void vec_ew_min(double[] a, int n, double min_val) {
        for (int i = 0; i < n; i++) {
            a[i] = Math.min(a[i], min_val);
        }
    }

    public static void vec_ew_max_vec(double[] a, double[] b, double[] c, int n) {
        for (int i = 0; i < n; i++) {
            c[i] = Math.max(a[i], b[i]);
        }
    }

    public static void vec_ew_min_vec(double[] a, double[] b, double[] c, int n) {
        for (int i = 0; i < n; i++) {
            c[i] = Math.min(a[i], b[i]);
        }
    }
    /* VECTOR FUNCTIONS ----------------------------------------------------------*/



    /* MATRIX FUNCTIONS ----------------------------------------------------------*/
    public static void mat_mult_scalar(CSCMatrix A, double sc) {
        int nnzA = A.Ap[A.n];
        for (int i = 0; i < nnzA; i++) {
            A.Ax[i] *= sc;
        }
    }

    public static void mat_premult_diag(CSCMatrix A, double[] d) {
        for (int j = 0; j < A.n; j++) {                // Cycle over columns
            for (int i = A.Ap[j]; i < A.Ap[j + 1]; i++) { // Cycle every row in the column
                A.Ax[i] *= d[A.Ai[i]];                  // Scale by corresponding element
                // of d for row i
            }
        }
    }

    public static void mat_postmult_diag(CSCMatrix A, double[] d) {
        for (int j = 0; j < A.n; j++) {                // Cycle over columns j
            for (int i = A.Ap[j]; i < A.Ap[j + 1]; i++) { // Cycle every row i in column j
                A.Ax[i] *= d[j];                        // Scale by corresponding element
                // of d for column j
            }
        }
    }

    public static void mat_vec(CSCMatrix A, double[] x, double[] y, int startx, int starty, int plus_eq) {
        //TODO: Check m,n
        if (plus_eq==0) {
            // y = 0
            for (int i = 0; i < A.m; i++) {
                y[i+starty] = 0;
            }
        }

        // if A is empty
        if (A.Ap[A.n] == 0) {
            return;
        }

        if (plus_eq==-1) {
            // y -=  A*x
            for (int j = 0; j < A.n; j++) {
                for (int i = A.Ap[j]; i < A.Ap[j + 1]; i++) {
                    y[A.Ai[i]+starty] -= A.Ax[i] * x[j+startx];
                }
            }
        } else {
            // y +=  A*x
            for (int j = 0; j < A.n; j++) {
                for (int i = A.Ap[j]; i < A.Ap[j + 1]; i++) {
                    y[A.Ai[i]+starty] += A.Ax[i] * x[j+startx];
                }
            }
        }
    }

    public static void mat_tpose_vec(CSCMatrix A, double[] x, double[] y, int startx, int starty, int plus_eq, boolean skip_diag) {

        if (plus_eq==0) {
            // y = 0
            for (int i = 0; i < A.n; i++) {
                y[i+starty] = 0;
            }
        }

        // if A is empty
        if (A.Ap[A.n] == 0) {
            return;
        }

        if (plus_eq==-1) {
            // y -=  A*x
            if (skip_diag) {
                for (int j = 0; j < A.n; j++) {
                    for (int k = A.Ap[j]; k < A.Ap[j + 1]; k++) {
                        int i  = A.Ai[k];
                        y[j+starty] -= i == j ? 0 : A.Ax[k] * x[i+startx];
                    }
                }
            } else {
                for (int j = 0; j < A.n; j++) {
                    for (int k = A.Ap[j]; k < A.Ap[j + 1]; k++) {
                        y[j+starty] -= A.Ax[k] * x[A.Ai[k]+startx];
                    }
                }
            }
        } else {
            // y +=  A*x
            if (skip_diag) {
                for (int j = 0; j < A.n; j++) {
                    for (int k = A.Ap[j]; k < A.Ap[j + 1]; k++) {
                        int i     = A.Ai[k];
                        y[j+starty] += i == j ? 0 : A.Ax[k] * x[i+startx];
                    }
                }
            } else {
                for (int j = 0; j < A.n; j++) {
                    for (int k = A.Ap[j]; k < A.Ap[j + 1]; k++) {
                        y[j+starty] += A.Ax[k] * x[A.Ai[k]+startx];
                    }
                }
            }
        }
    }

    public static void mat_inf_norm_cols(CSCMatrix M, double[] E) {

        // Initialize zero max elements
        for (int j = 0; j < M.n; j++) {
            E[j] = 0;
        }

        // Compute maximum across columns
        for (int j = 0; j < M.n; j++) {
            for (int ptr = M.Ap[j]; ptr < M.Ap[j + 1]; ptr++) {
                E[j] =  Math.max(Math.abs(M.Ax[ptr]), E[j]);
            }
        }
    }

    public static void mat_inf_norm_rows(CSCMatrix M, double[] E) {
        // Initialize zero max elements
        for (int j = 0; j < M.m; j++) {
            E[j] = 0;
        }

        // Compute maximum across rows
        for (int j = 0; j < M.n; j++) {
            for (int ptr = M.Ap[j]; ptr < M.Ap[j + 1]; ptr++) {
                int i    = M.Ai[ptr];
                E[i] = Math.max(Math.abs(M.Ax[ptr]), E[i]);
            }
        }
    }

    public static void mat_inf_norm_cols_sym_triu(CSCMatrix M, double[] E) {
    	

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
                double abs_x = Math.abs(M.Ax[ptr]);
                E[j]  = Math.max(abs_x, E[j]);

                if (i != j) {
                    E[i] = Math.max(abs_x, E[i]);
                }
            }
        }
    }


    public static double quad_form(CSCMatrix P, double[] x) {
        double quad_form = 0;

        for (int j = 0; j < P.n; j++) {                      // Iterate over columns
            for (int ptr = P.Ap[j]; ptr < P.Ap[j + 1]; ptr++) { // Iterate over rows
                int i = P.Ai[ptr];                                // Row index

                if (i == j) {                                 // Diagonal element
                    quad_form += 0.5 * P.Ax[ptr] * x[i] * x[i];
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