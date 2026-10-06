package com.quantego.josqp;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

/** Optional SIMD kernels. Keep arithmetic order; deliberately do not use FMA. */
final class VectorKernels extends NumericKernels {
    static final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_PREFERRED;

    private static DoubleVector maxIgnoringNaN(DoubleVector max, DoubleVector value) {
        return max.blend(value, value.compare(VectorOperators.GT, max));
    }

    @Override public void vec_add_scaled(double[] c, double[] a, double[] b, int n, double sc) {
        if (n < 64) { super.vec_add_scaled(c, a, b, n, sc); return; }
        int i = 0;
        for (int end = SPECIES.loopBound(n); i < end; i += SPECIES.length()) {
            DoubleVector.fromArray(SPECIES, a, i).add(DoubleVector.fromArray(SPECIES, b, i).mul(sc)).intoArray(c, i);
        }
        for (; i < n; i++) c[i] = a[i] + sc * b[i];
    }

    @Override public double vec_norm_inf(double[] v, int n) {
        if (n < 64) return super.vec_norm_inf(v, n);
        DoubleVector m0 = DoubleVector.zero(SPECIES), m1 = m0, m2 = m0, m3 = m0;
        int i = 0, lanes = SPECIES.length(), end = SPECIES.loopBound(n);
        for (; i <= end - 4 * lanes; i += 4 * lanes) {
            m0 = maxIgnoringNaN(m0, DoubleVector.fromArray(SPECIES, v, i).abs());
            m1 = maxIgnoringNaN(m1, DoubleVector.fromArray(SPECIES, v, i + 1 * lanes).abs());
            m2 = maxIgnoringNaN(m2, DoubleVector.fromArray(SPECIES, v, i + 2 * lanes).abs());
            m3 = maxIgnoringNaN(m3, DoubleVector.fromArray(SPECIES, v, i + 3 * lanes).abs());
        }
        DoubleVector max = m0.max(m1).max(m2).max(m3);
        for (; i < end; i += lanes) max = maxIgnoringNaN(max, DoubleVector.fromArray(SPECIES, v, i).abs());
        double result = max.reduceLanes(VectorOperators.MAX);
        for (; i < n; i++) if (Math.abs(v[i]) > result) result = Math.abs(v[i]);
        return result;
    }

    @Override public double vec_scaled_norm_inf(double[] s, double[] v, int n) {
        if (n < 64) return super.vec_scaled_norm_inf(s, v, n);
        DoubleVector m0 = DoubleVector.zero(SPECIES), m1 = m0, m2 = m0, m3 = m0;
        int i = 0, lanes = SPECIES.length(), end = SPECIES.loopBound(n);
        for (; i <= end - 4 * lanes; i += 4 * lanes) {
            m0 = maxIgnoringNaN(m0, DoubleVector.fromArray(SPECIES, s, i).mul(DoubleVector.fromArray(SPECIES, v, i)).abs());
            m1 = maxIgnoringNaN(m1, DoubleVector.fromArray(SPECIES, s, i + 1 * lanes).mul(DoubleVector.fromArray(SPECIES, v, i + 1 * lanes)).abs());
            m2 = maxIgnoringNaN(m2, DoubleVector.fromArray(SPECIES, s, i + 2 * lanes).mul(DoubleVector.fromArray(SPECIES, v, i + 2 * lanes)).abs());
            m3 = maxIgnoringNaN(m3, DoubleVector.fromArray(SPECIES, s, i + 3 * lanes).mul(DoubleVector.fromArray(SPECIES, v, i + 3 * lanes)).abs());
        }
        DoubleVector max = m0.max(m1).max(m2).max(m3);
        for (; i < end; i += lanes) max = maxIgnoringNaN(max, DoubleVector.fromArray(SPECIES, s, i).mul(DoubleVector.fromArray(SPECIES, v, i)).abs());
        double result = max.reduceLanes(VectorOperators.MAX);
        for (; i < n; i++) if (Math.abs(s[i] * v[i]) > result) result = Math.abs(s[i] * v[i]);
        return result;
    }

    @Override public double vec_norm_inf_diff(double[] a, double[] b, int n) {
        if (n < 64) return super.vec_norm_inf_diff(a, b, n);
        DoubleVector m0 = DoubleVector.zero(SPECIES), m1 = m0, m2 = m0, m3 = m0;
        int i = 0, lanes = SPECIES.length(), end = SPECIES.loopBound(n);
        for (; i <= end - 4 * lanes; i += 4 * lanes) {
            m0 = maxIgnoringNaN(m0, DoubleVector.fromArray(SPECIES, a, i).sub(DoubleVector.fromArray(SPECIES, b, i)).abs());
            m1 = maxIgnoringNaN(m1, DoubleVector.fromArray(SPECIES, a, i + 1 * lanes).sub(DoubleVector.fromArray(SPECIES, b, i + 1 * lanes)).abs());
            m2 = maxIgnoringNaN(m2, DoubleVector.fromArray(SPECIES, a, i + 2 * lanes).sub(DoubleVector.fromArray(SPECIES, b, i + 2 * lanes)).abs());
            m3 = maxIgnoringNaN(m3, DoubleVector.fromArray(SPECIES, a, i + 3 * lanes).sub(DoubleVector.fromArray(SPECIES, b, i + 3 * lanes)).abs());
        }
        DoubleVector max = m0.max(m1).max(m2).max(m3);
        for (; i < end; i += lanes) max = maxIgnoringNaN(max, DoubleVector.fromArray(SPECIES, a, i).sub(DoubleVector.fromArray(SPECIES, b, i)).abs());
        double result = max.reduceLanes(VectorOperators.MAX);
        for (; i < n; i++) if (Math.abs(a[i] - b[i]) > result) result = Math.abs(a[i] - b[i]);
        return result;
    }

    @Override public void vec_ew_prod(double[] a, double[] b, double[] c, int n) {
        if (n < 64) { super.vec_ew_prod(a, b, c, n); return; }
        int i = 0;
        for (int end = SPECIES.loopBound(n); i < end; i += SPECIES.length())
            DoubleVector.fromArray(SPECIES, b, i).mul(DoubleVector.fromArray(SPECIES, a, i)).intoArray(c, i);
        for (; i < n; i++) c[i] = b[i] * a[i];
    }

    @Override public void vec_ew_recipr(double[] a, double[] b, int n) {
        if (n < 64) { super.vec_ew_recipr(a, b, n); return; }
        int i = 0;
        DoubleVector one = DoubleVector.broadcast(SPECIES, 1.0);
        for (int end = SPECIES.loopBound(n); i < end; i += SPECIES.length())
            one.div(DoubleVector.fromArray(SPECIES, a, i)).intoArray(b, i);
        for (; i < n; i++) b[i] = 1.0 / a[i];
    }

    @Override public void vec_ew_sqrt(double[] a, int n) {
        if (n < 64) { super.vec_ew_sqrt(a, n); return; }
        int i = 0;
        for (int end = SPECIES.loopBound(n); i < end; i += SPECIES.length())
            DoubleVector.fromArray(SPECIES, a, i).sqrt().intoArray(a, i);
        for (; i < n; i++) a[i] = Math.sqrt(a[i]);
    }

    private void scale(double[] a, int start, int end, double sc) {
        int i = start;
        int bound = start + SPECIES.loopBound(end - start);
        for (; i < bound; i += SPECIES.length())
            DoubleVector.fromArray(SPECIES, a, i).mul(sc).intoArray(a, i);
        for (; i < end; i++) a[i] *= sc;
    }

    @Override public void vec_mult_scalar(double[] a, double sc, int n) { scale(a, 0, n, sc); }
    @Override public void mat_mult_scalar(CSCMatrix a, double sc) { scale(a.Ax, 0, a.Ap[a.n], sc); }
    @Override public void mat_postmult_diag(CSCMatrix a, double[] d) {
        for (int j = 0; j < a.n; j++) {
            int start = a.Ap[j], end = a.Ap[j + 1];
            if (end - start < 32) {
                for (int p = start; p < end; p++) a.Ax[p] *= d[j];
            } else scale(a.Ax, start, end, d[j]);
        }
    }

    @Override public void project(OSQP.Workspace w, double[] z) {
        if (w.data.m < 64) { super.project(w, z); return; }
        int i = 0, n = w.data.m;
        for (int end = SPECIES.loopBound(n); i < end; i += SPECIES.length()) {
            DoubleVector.fromArray(SPECIES, z, i).max(DoubleVector.fromArray(SPECIES, w.data.l, i))
                    .min(DoubleVector.fromArray(SPECIES, w.data.u, i)).intoArray(z, i);
        }
        for (; i < n; i++) z[i] = Math.min(Math.max(z[i], w.data.l[i]), w.data.u[i]);
    }

    @Override public void project_normalcone(OSQP.Workspace w, double[] z, double[] y) {
        if (w.data.m < 64) { super.project_normalcone(w, z, y); return; }
        int i = 0, n = w.data.m;
        for (int end = SPECIES.loopBound(n); i < end; i += SPECIES.length()) {
            DoubleVector sum = DoubleVector.fromArray(SPECIES, z, i).add(DoubleVector.fromArray(SPECIES, y, i));
            DoubleVector projected = sum.max(DoubleVector.fromArray(SPECIES, w.data.l, i))
                    .min(DoubleVector.fromArray(SPECIES, w.data.u, i));
            sum.intoArray(w.z_prev, i);
            projected.intoArray(z, i);
            sum.sub(projected).intoArray(y, i);
        }
        for (; i < n; i++) {
            w.z_prev[i] = z[i] + y[i];
            z[i] = Math.min(Math.max(w.z_prev[i], w.data.l[i]), w.data.u[i]);
            y[i] = w.z_prev[i] - z[i];
        }
    }

    @Override public void compute_rhs(OSQP.Workspace w) {
        if (Math.max(w.data.n, w.data.m) < 64) { super.compute_rhs(w); return; }
        int i = 0, n = w.data.n, m = w.data.m;
        for (int end = SPECIES.loopBound(n); i < end; i += SPECIES.length()) {
            DoubleVector.fromArray(SPECIES, w.x_prev, i).mul(w.settings.sigma)
                    .sub(DoubleVector.fromArray(SPECIES, w.data.q, i)).intoArray(w.xz_tilde, i);
        }
        for (; i < n; i++) w.xz_tilde[i] = w.settings.sigma * w.x_prev[i] - w.data.q[i];
        i = 0;
        for (int end = SPECIES.loopBound(m); i < end; i += SPECIES.length()) {
            DoubleVector.fromArray(SPECIES, w.z_prev, i).sub(DoubleVector.fromArray(SPECIES, w.rho_inv_vec, i)
                    .mul(DoubleVector.fromArray(SPECIES, w.y, i))).intoArray(w.xz_tilde, n + i);
        }
        for (; i < m; i++) w.xz_tilde[n + i] = w.z_prev[i] - w.rho_inv_vec[i] * w.y[i];
    }



    @Override public void update_z(OSQP.Workspace w) {
        if (w.data.m < 64) { super.update_z(w); return; }
        int i = 0, n = w.data.m, offset = w.data.n;
        double alpha = w.settings.alpha, beta = 1.0 - alpha;
        for (int end = SPECIES.loopBound(n); i < end; i += SPECIES.length()) {
            DoubleVector.fromArray(SPECIES, w.xz_tilde, offset + i).mul(alpha)
                    .add(DoubleVector.fromArray(SPECIES, w.z_prev, i).mul(beta))
                    .add(DoubleVector.fromArray(SPECIES, w.rho_inv_vec, i).mul(DoubleVector.fromArray(SPECIES, w.y, i)))
                    .max(DoubleVector.fromArray(SPECIES, w.data.l, i))
                    .min(DoubleVector.fromArray(SPECIES, w.data.u, i)).intoArray(w.z, i);
        }
        for (; i < n; i++) {
            w.z[i] = Math.min(Math.max(alpha * w.xz_tilde[offset + i] + beta * w.z_prev[i]
                    + w.rho_inv_vec[i] * w.y[i], w.data.l[i]), w.data.u[i]);
        }
    }

    @Override public void update_y(OSQP.Workspace w) {
        if (w.data.m < 64) { super.update_y(w); return; }
        int i = 0, n = w.data.m, offset = w.data.n;
        double alpha = w.settings.alpha, beta = 1.0 - alpha;
        for (int end = SPECIES.loopBound(n); i < end; i += SPECIES.length()) {
            DoubleVector dy = DoubleVector.fromArray(SPECIES, w.rho_vec, i)
                    .mul(DoubleVector.fromArray(SPECIES, w.xz_tilde, offset + i).mul(alpha)
                    .add(DoubleVector.fromArray(SPECIES, w.z_prev, i).mul(beta))
                    .sub(DoubleVector.fromArray(SPECIES, w.z, i)));
            dy.intoArray(w.delta_y, i);
            DoubleVector.fromArray(SPECIES, w.y, i).add(dy).intoArray(w.y, i);
        }
        for (; i < n; i++) {
            w.delta_y[i] = w.rho_vec[i] * (alpha * w.xz_tilde[offset + i] + beta * w.z_prev[i] - w.z[i]);
            w.y[i] += w.delta_y[i];
        }
    }
}
