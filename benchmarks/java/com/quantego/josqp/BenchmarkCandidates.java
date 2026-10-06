package com.quantego.josqp;

import jdk.incubator.vector.DoubleVector;

/** Experimental kernels and scalar controls for reproducible comparisons. */
final class BenchmarkCandidates {
    private static final jdk.incubator.vector.VectorSpecies<Double> SPECIES = VectorKernels.SPECIES;
    static double unrolledScalarNorm(double[] a, int n) {
        double m0 = 0, m1 = 0, m2 = 0, m3 = 0;
        int i = 0;
        for (; i <= n - 4; i += 4) {
            double v0 = Math.abs(a[i]), v1 = Math.abs(a[i + 1]);
            double v2 = Math.abs(a[i + 2]), v3 = Math.abs(a[i + 3]);
            if (v0 > m0) m0 = v0;
            if (v1 > m1) m1 = v1;
            if (v2 > m2) m2 = v2;
            if (v3 > m3) m3 = v3;
        }
        double result = Math.max(Math.max(m0, m1), Math.max(m2, m3));
        for (; i < n; i++) if (Math.abs(a[i]) > result) result = Math.abs(a[i]);
        return result;
    }
    static void vectorUpdateX(OSQP.Workspace w) {
        int i = 0, n = w.data.n;
        double alpha = w.settings.alpha, beta = 1.0 - alpha;
        for (int end = SPECIES.loopBound(n); i < end; i += SPECIES.length()) {
            DoubleVector previous = DoubleVector.fromArray(SPECIES, w.x_prev, i);
            DoubleVector x = DoubleVector.fromArray(SPECIES, w.xz_tilde, i).mul(alpha).add(previous.mul(beta));
            x.intoArray(w.x, i);
            x.sub(previous).intoArray(w.delta_x, i);
        }
        for (; i < n; i++) {
            w.x[i] = alpha * w.xz_tilde[i] + beta * w.x_prev[i];
            w.delta_x[i] = w.x[i] - w.x_prev[i];
        }
    }
}
