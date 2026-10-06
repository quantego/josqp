package com.quantego.josqp;

/** Scalar kernels keep the default artifact compatible with Java 8. */
class NumericKernels {
    static final NumericKernels INSTANCE = create();

    private static NumericKernels create() {
        if (!Boolean.getBoolean("josqp.vector")) return new NumericKernels();
        try {
            return (NumericKernels) Class.forName("com.quantego.josqp.VectorKernels")
            .getDeclaredConstructor().newInstance();
        } catch (ReflectiveOperationException | LinkageError e) {
            throw new IllegalStateException("Vector backend requires the vector build and "
            + "--add-modules jdk.incubator.vector", e);
        }
    }

    public void vec_add_scaled(double[] c, double[] a, double[] b, int n, double sc){

        for (int i=0; i<n;i++){
            c[i] = a[i] + sc*b[i];
        }
    }

    public double vec_scaled_norm_inf(double[] S, double[] v, int l){

        double max = 0;
        for (int i =0;i<l;i++){
            double abs_Sv_i = Math.abs(S[i] * v[i]);
            if(abs_Sv_i > max){
                max = abs_Sv_i;
            }
        }
        return max;
    }

    public double vec_norm_inf(double [] v, int l){
        double max = 0;
        for (int i = 0; i<l; i++){
            double abs_v_i = Math.abs(v[i]);
            if(abs_v_i > max){
                max = abs_v_i;
            }
        }
        return max;
    }

    public double vec_norm_inf_diff(double[] a, double[] b, int l){
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

    public void vec_mult_scalar(double[] a, double sc, int n) {
        for (int i = 0; i < n; i++) {
            a[i] *= sc;
        }
    }

    public void vec_ew_prod(double[] a , double[] b, double[] c, int n) {
        for (int i = 0; i < n; i++) {
            c[i] = b[i] * a[i];
        }
    }

    public void vec_ew_recipr(double[] a , double[] b, int n) {
        for (int i = 0; i < n; i++) {
            b[i] = 1.0 / a[i];
        }
    }

    public void vec_ew_sqrt(double[] a, int n) {
        for (int i = 0; i < n; i++) {
            a[i] = Math.sqrt(a[i]);
        }
    }

    public void mat_mult_scalar(CSCMatrix A, double sc) {
        int nnzA = A.Ap[A.n];
        for (int i = 0; i < nnzA; i++) {
            A.Ax[i] *= sc;
        }
    }

    public void mat_postmult_diag(CSCMatrix A, double[] d) {
        for (int j = 0; j < A.n; j++) {                // Cycle over columns j
            for (int i = A.Ap[j]; i < A.Ap[j + 1]; i++) { // Cycle every row i in column j
                A.Ax[i] *= d[j];                        // Scale by corresponding element
                // of d for column j
            }
        }
    }

    public void compute_rhs(OSQP.Workspace work) {
        int i; // Index

        for (i = 0; i < work.data.n; i++) {
            // Cycle over part related to x variables
            work.xz_tilde[i] = work.settings.sigma * work.x_prev[i] -
            work.data.q[i];
        }

        for (i = 0; i < work.data.m; i++) {
            // Cycle over dual variable in the first step (nu)
            work.xz_tilde[i + work.data.n] = work.z_prev[i] - work.rho_inv_vec[i] *
            work.y[i];
        }
    }

    public void update_x(OSQP.Workspace work) {
        int i;

        // update x
        for (i = 0; i < work.data.n; i++) {
            work.x[i] = work.settings.alpha * work.xz_tilde[i] +
            (1.0 - work.settings.alpha) * work.x_prev[i];
        }

        // update delta_x
        for (i = 0; i < work.data.n; i++) {
            work.delta_x[i] = work.x[i] - work.x_prev[i];
        }
    }

    public void update_z(OSQP.Workspace work) {
        int i;

        // update z
        for (i = 0; i < work.data.m; i++) {
            work.z[i] = work.settings.alpha * work.xz_tilde[i + work.data.n] +
            (1.0 - work.settings.alpha) * work.z_prev[i] +
            work.rho_inv_vec[i] * work.y[i];
        }

        // project z
        project(work, work.z);
    }

    public void update_y(OSQP.Workspace work) {
        int i; // Index

        for (i = 0; i < work.data.m; i++) {
            work.delta_y[i] = work.rho_vec[i] *
            (work.settings.alpha *
            work.xz_tilde[i + work.data.n] +
            (1.0 - work.settings.alpha) * work.z_prev[i] -
            work.z[i]);
            work.y[i] += work.delta_y[i];
        }
    }

    public void project(OSQP.Workspace work, double[] z) {
        int i, m;

        m = work.data.m;

        for (i = 0; i < m; i++) {
            z[i] = Math.min(Math.max(z[i],
            work.data.l[i]), // Between lower
            work.data.u[i]);       // and upper bounds
        }
    }

    public void project_normalcone(OSQP.Workspace work, double[] z, double[] y) {
        int i, m;

        // NB: Use z_prev as temporary vector

        m = work.data.m;

        for (i = 0; i < m; i++) {
            work.z_prev[i] = z[i] + y[i];
            z[i]            = Math.min(Math.max(work.z_prev[i], work.data.l[i]),
            work.data.u[i]);
            y[i] = work.z_prev[i] - z[i];
        }
    }
}
