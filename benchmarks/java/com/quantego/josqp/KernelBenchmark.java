package com.quantego.josqp;

import java.util.Arrays;
import java.util.Locale;
import java.util.Random;
import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorOperators;

/** Exploratory paired kernel benchmark; solver measurements remain the deciding metric. */
public final class KernelBenchmark {
    private static volatile double sink;
    private static final NumericKernels SCALAR = new NumericKernels();
    private static final NumericKernels VECTOR = new VectorKernels();
    private static double[] a, b, c;
    private static CSCMatrix matrix;
    private static OSQP.Workspace work;
    private static int size;

    private static double sparseGather(CSCMatrix m, double[] x, double[] y) {
        int lanes = VectorKernels.SPECIES.length();
        for (int j = 0; j < m.n; j++) {
            int p = m.Ap[j], end = m.Ap[j + 1];
            DoubleVector sum = DoubleVector.zero(VectorKernels.SPECIES);
            int bound = p + VectorKernels.SPECIES.loopBound(end - p);
            for (; p < bound; p += lanes) {
                sum = sum.add(DoubleVector.fromArray(VectorKernels.SPECIES, m.Ax, p)
                        .mul(DoubleVector.fromArray(VectorKernels.SPECIES, x, 0, m.Ai, p)));
            }
            double value = sum.reduceLanes(VectorOperators.ADD);
            for (; p < end; p++) value += m.Ax[p] * x[m.Ai[p]];
            y[j] = value;
        }
        return y[0];
    }

    private static double operation(String name, boolean vector) {
        NumericKernels k = vector ? VECTOR : SCALAR;
        switch (name) {
            case "norm": return k.vec_norm_inf(a, size);
            case "norm_vs_unrolled_scalar":
                return vector ? VECTOR.vec_norm_inf(a, size) : BenchmarkCandidates.unrolledScalarNorm(a, size);
            case "norm_diff": return k.vec_norm_inf_diff(a, b, size);
            case "scaled_norm": return k.vec_scaled_norm_inf(a, b, size);
            case "add_scaled": k.vec_add_scaled(c, a, b, size, 0.3); return c[0];
            case "ew_prod": k.vec_ew_prod(a, b, c, size); return c[0];
            case "reciprocal": k.vec_ew_recipr(a, c, size); return c[0];
            case "admm_x": k.update_x(work); return work.x[0];
            case "admm_x_candidate":
                if (vector) BenchmarkCandidates.vectorUpdateX(work); else SCALAR.update_x(work);
                return work.x[0];
            case "admm_z": k.update_z(work); return work.z[0];
            case "admm_y": k.update_y(work); return work.y[0];
            case "projection": k.project(work, work.z); return work.z[0];
            case "sparse_column_scale": k.mat_postmult_diag(matrix, b); return matrix.Ax[0];
            case "sparse_gather":
                if (vector) return sparseGather(matrix, a, c);
                LinAlg.mat_tpose_vec(matrix, a, c, 0, 0, 0, false); return c[0];
            default: throw new IllegalArgumentException(name);
        }
    }

    private static double batch(String op, boolean vector, int reps) {
        double checksum = 0;
        long start = System.nanoTime();
        for (int i = 0; i < reps; i++) {
            // Mutable input prevents a loop-invariant result being hoisted.
            a[0] = 0.1 + (i & 15) * 1e-5;
            checksum += operation(op, vector);
        }
        double ns = (double) (System.nanoTime() - start) / reps;
        sink = checksum;
        return ns;
    }

    public static void main(String[] args) {
        Locale.setDefault(Locale.ROOT);
        Random rng = new Random(20261006);
        System.err.println("species=" + VectorKernels.SPECIES + "; java=" + System.getProperty("java.version")
                + "; arch=" + System.getProperty("os.arch"));
        System.out.println("kernel,n,entries_per_column,sample,scalar_ns,vector_ns,speedup");
        for (int n : new int[]{8, 128, 4096, 65536}) {
            size = n;
            a = new double[n]; b = new double[n]; c = new double[n];
            for (int i = 0; i < n; i++) { a[i] = rng.nextDouble() + 0.1; b[i] = 1.0; }
            work = new OSQP.Workspace();
            work.data = new OSQP.Data(n, n, new CSCMatrix(n, n, 0, new int[n + 1], new int[0], new double[0]),
                    new CSCMatrix(n, n, 0, new int[n + 1], new int[0], new double[0]),
                    new double[n], new double[n], new double[n], 0);
            Arrays.fill(work.data.l, -0.2); Arrays.fill(work.data.u, 0.7);
            work.settings = new OSQP.Settings();
            work.x = new double[n]; work.x_prev = a; work.xz_tilde = new double[2 * n]; work.delta_x = new double[n];
            work.z = new double[n]; work.z_prev = b; work.y = new double[n]; work.delta_y = new double[n];
            work.rho_vec = new double[n]; work.rho_inv_vec = new double[n];
            Arrays.fill(work.rho_vec, 0.1); Arrays.fill(work.rho_inv_vec, 10.0);
            for (int i = 0; i < work.xz_tilde.length; i++) work.xz_tilde[i] = rng.nextDouble();
            for (int entries : new int[]{2, 8, 32, 128}) {
                int cols = Math.max(1, n / entries);
                int[] ap = new int[cols + 1], ai = new int[cols * entries];
                double[] ax = new double[ai.length];
                for (int j = 0; j < cols; j++) {
                    ap[j] = j * entries;
                    for (int p = ap[j]; p < ap[j] + entries; p++) { ai[p] = rng.nextInt(n); ax[p] = rng.nextDouble(); }
                }
                ap[cols] = ai.length;
                matrix = new CSCMatrix(n, cols, ai.length, ap, ai, ax);
                String[] ops = entries == 2 ? new String[]{"norm", "norm_vs_unrolled_scalar", "norm_diff", "scaled_norm", "add_scaled", "ew_prod", "reciprocal", "admm_x", "admm_x_candidate", "admm_z", "admm_y", "projection", "sparse_column_scale", "sparse_gather"}
                        : new String[]{"sparse_column_scale", "sparse_gather"};
                for (String op : ops) {
                    if (op.equals("sparse_gather")) {
                        LinAlg.mat_tpose_vec(matrix, a, c, 0, 0, 0, false);
                        double[] expected = c.clone();
                        sparseGather(matrix, a, c);
                        for (int i = 0; i < cols; i++) if (Math.abs(expected[i] - c[i]) > 1e-10)
                            throw new AssertionError("Gather mismatch");
                    }
                    int reps = Math.max(64, Math.min(100000, 1000000 / n));
                    long deadline = System.nanoTime() + 300000000L;
                    do { batch(op, false, reps); batch(op, true, reps); } while (System.nanoTime() < deadline);
                    for (int sample = 0; sample < 7; sample++) {
                        double scalar, vector;
                        if ((sample & 1) == 0) { scalar = batch(op, false, reps); vector = batch(op, true, reps); }
                        else { vector = batch(op, true, reps); scalar = batch(op, false, reps); }
                        System.out.printf("%s,%d,%d,%d,%.3f,%.3f,%.4f%n", op, n, entries, sample, scalar, vector, scalar / vector);
                    }
                }
            }
        }
    }
}
