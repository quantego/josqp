package com.quantego.josqp;

import java.util.Random;
import org.junit.Test;
import static org.junit.Assert.*;

public class VectorKernelsTest {
    private final NumericKernels scalar = new NumericKernels();
    private final NumericKernels vector = new VectorKernels();

    private void same(double[] a, double[] b) {
        assertEquals(a.length, b.length);
        for (int i = 0; i < a.length; i++)
            assertEquals("element " + i, Double.doubleToLongBits(a[i]), Double.doubleToLongBits(b[i]));
    }

    @Test public void denseKernelsHandleTailsAliasingAndSpecialValues() {
        Random rng = new Random(482);
        for (int n = 0; n < 130; n++) {
            double[] a = new double[n + 5], b = new double[n + 5];
            for (int i = 0; i < a.length; i++) { a[i] = rng.nextDouble() - 0.5; b[i] = rng.nextDouble() - 0.5; }
            if (n > 5) { a[0] = Double.NaN; a[1] = -0.0; a[2] = Double.POSITIVE_INFINITY; b[2] = 0.0; b[3] = Double.NEGATIVE_INFINITY; }
            assertEquals(scalar.vec_norm_inf(a, n), vector.vec_norm_inf(a, n), 0);
            assertEquals(scalar.vec_norm_inf_diff(a, b, n), vector.vec_norm_inf_diff(a, b, n), 0);
            assertEquals(scalar.vec_scaled_norm_inf(a, b, n), vector.vec_scaled_norm_inf(a, b, n), 0);
            for (int alias = 0; alias < 3; alias++) {
                double[] sa = a.clone(), sb = b.clone(), va = a.clone(), vb = b.clone();
                double[] sc = alias == 0 ? new double[a.length] : alias == 1 ? sa : sb;
                double[] vc = alias == 0 ? new double[a.length] : alias == 1 ? va : vb;
                scalar.vec_add_scaled(sc, sa, sb, n, 0.3); vector.vec_add_scaled(vc, va, vb, n, 0.3); same(sc, vc);
                sa = a.clone(); sb = b.clone(); va = a.clone(); vb = b.clone();
                sc = alias == 0 ? new double[a.length] : alias == 1 ? sa : sb;
                vc = alias == 0 ? new double[a.length] : alias == 1 ? va : vb;
                scalar.vec_ew_prod(sa, sb, sc, n); vector.vec_ew_prod(va, vb, vc, n); same(sc, vc);
            }
            double[] s = a.clone(), v = a.clone();
            scalar.vec_ew_recipr(s, s, n); vector.vec_ew_recipr(v, v, n); same(s, v);
            s = a.clone(); v = a.clone(); scalar.vec_ew_sqrt(s, n); vector.vec_ew_sqrt(v, n); same(s, v);
            s = a.clone(); v = a.clone(); scalar.vec_mult_scalar(s, 0.3, n); vector.vec_mult_scalar(v, 0.3, n); same(s, v);
        }
    }

    private OSQP.Workspace workspace(int n, int m, long seed) {
        OSQP.Workspace w = new OSQP.Workspace();
        w.data = new OSQP.Data(n, m, new CSCMatrix(n, n, 0, new int[n + 1], new int[0], new double[0]),
                new CSCMatrix(m, n, 0, new int[n + 1], new int[0], new double[0]),
                new double[n], new double[m], new double[m], 0);
        w.settings = new OSQP.Settings();
        w.x = new double[n]; w.x_prev = new double[n]; w.delta_x = new double[n]; w.xz_tilde = new double[n + m];
        w.z = new double[m]; w.z_prev = new double[m]; w.y = new double[m]; w.delta_y = new double[m];
        w.rho_vec = new double[m]; w.rho_inv_vec = new double[m];
        Random rng = new Random(seed);
        for (double[] arr : new double[][]{w.x, w.x_prev, w.xz_tilde, w.z, w.z_prev, w.y, w.rho_vec, w.rho_inv_vec, w.data.q})
            for (int i = 0; i < arr.length; i++) arr[i] = rng.nextDouble() - 0.5;
        for (int i = 0; i < m; i++) { w.data.l[i] = -0.2; w.data.u[i] = 0.4; }
        if (m > 3) { w.data.l[0] = Double.NEGATIVE_INFINITY; w.data.u[0] = Double.POSITIVE_INFINITY;
            w.data.l[1] = -0.0; w.data.u[1] = 0.0; w.z[2] = Double.NaN; w.data.l[3] = Double.NaN; }
        return w;
    }

    @Test public void admmAndProjectionPreserveArithmeticBits() {
        for (int n = 1; n < 70; n++) for (int m : new int[]{0, 1, n, n + 1}) {
            OSQP.Workspace s = workspace(n, m, n + m), v = workspace(n, m, n + m);
            scalar.compute_rhs(s); vector.compute_rhs(v); same(s.xz_tilde, v.xz_tilde);
            scalar.update_x(s); vector.update_x(v); same(s.x, v.x); same(s.delta_x, v.delta_x);
            scalar.update_z(s); vector.update_z(v); same(s.z, v.z);
            scalar.update_y(s); vector.update_y(v); same(s.y, v.y); same(s.delta_y, v.delta_y);
            scalar.project(s, s.z); vector.project(v, v.z); same(s.z, v.z);
            scalar.project_normalcone(s, s.z, s.y); vector.project_normalcone(v, v.z, v.y);
            same(s.z_prev, v.z_prev); same(s.z, v.z); same(s.y, v.y);
        }
    }

    @Test public void sparseColumnScalingHandlesEmptyColumnsAndOffsets() {
        int[] ap = {0, 0, 3, 4, 133, 133};
        int[] ai = new int[133]; double[] ax = new double[135];
        for (int i = 0; i < ax.length; i++) ax[i] = i - 3.1;
        CSCMatrix s = new CSCMatrix(5, 5, 135, ap.clone(), ai.clone(), ax.clone());
        CSCMatrix v = new CSCMatrix(5, 5, 135, ap.clone(), ai.clone(), ax.clone());
        double[] d = {1, -2.0, 0.0, 0.3, Double.NaN};
        scalar.mat_postmult_diag(s, d); vector.mat_postmult_diag(v, d); same(s.Ax, v.Ax);
        scalar.mat_mult_scalar(s, -0.0); vector.mat_mult_scalar(v, -0.0); same(s.Ax, v.Ax);
    }
}
