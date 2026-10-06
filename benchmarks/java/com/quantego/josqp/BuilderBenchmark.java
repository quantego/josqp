package com.quantego.josqp;

import java.lang.management.ManagementFactory;
import java.lang.management.MemoryPoolMXBean;
import java.lang.management.MemoryType;
import java.lang.management.MemoryUsage;
import com.sun.management.ThreadMXBean;

/**
 * Measures the cost of building a model algebraically: wall time, allocated bytes and peak heap.
 * All numbers are for the build phase only, the solve is not started.
 *
 * Usage: BuilderBenchmark <variables> <rows> <nonzerosPerRow> [algebraic|bulk]
 */
public final class BuilderBenchmark {

    public static void main(String[] args) {
        int n = args.length > 0 ? Integer.parseInt(args[0]) : 100_000;
        int rows = args.length > 1 ? Integer.parseInt(args[1]) : 1_000;
        int perRow = args.length > 2 ? Integer.parseInt(args[2]) : 200;
        boolean bulk = args.length > 3 && args[3].equals("bulk");

        ThreadMXBean threads = (ThreadMXBean) ManagementFactory.getThreadMXBean();
        long allocBefore = threads.getCurrentThreadAllocatedBytes();
        long t0 = System.nanoTime();

        Model.Builder builder = Model.getBuilder();
        Model model;
        if (bulk) {
            // The index arrays are prepared outside the measured region, as a user would.
            int[] all = new int[n];
            double[] linear = new double[n];
            double[] quadratic = new double[n];
            double[] lower = new double[n];
            double[] upper = new double[n];
            for (int i = 0; i < n; i++) {
                all[i] = i;
                linear[i] = 1e-3;
                quadratic[i] = 1.0;
                lower[i] = -1.0;
                upper[i] = 1.0;
            }
            builder.addVariables(n, lower, upper);
            builder.linearObjective(all, linear);
            builder.quadraticObjective(all, all, quadratic);
            builder.offset(42.0);
            builder.minimize();
            for (int r = 0; r < rows; r++) {
                int[] indices = new int[perRow];
                double[] coefficients = new double[perRow];
                for (int k = 0; k < perRow; k++) {
                    indices[k] = (r * 31 + k * 7) % n;
                    coefficients[k] = 1.0 / (1 + k);
                }
                builder.constraint(indices, coefficients, Double.NEGATIVE_INFINITY, 1.0);
            }
            model = builder.build();
        } else {
            Model.Variable[] x = builder.addVariables(n);
            Model.Objective objective = builder.setObjective();
            for (int i = 0; i < n; i++) {
                objective.add(1.0, x[i], x[i]);
                objective.add(1e-3, x[i]);
            }
            objective.add(42.0);
            objective.minimize();
            for (int i = 0; i < n; i++) {
                x[i].lb(-1.0).ub(1.0);
            }
            for (int r = 0; r < rows; r++) {
                Model.Constraint c = builder.addConstraint();
                for (int k = 0; k < perRow; k++) {
                    c.add(1.0 / (1 + k), x[(r * 31 + k * 7) % n]);
                }
                c.leq(1.0);
            }
            model = builder.build();
        }

        long tBuild = System.nanoTime() - t0;
        long allocated = threads.getCurrentThreadAllocatedBytes() - allocBefore;

        long t1 = System.nanoTime();
        OSQP.Status status = model.solve();
        long tSolve = System.nanoTime() - t1;

        long peak = 0;
        for (MemoryPoolMXBean pool : ManagementFactory.getMemoryPoolMXBeans())
            if (pool.getType() == MemoryType.HEAP) peak += pool.getPeakUsage().getUsed();
        System.out.printf("mode=%s n=%d structural_rows=%d user_rows=%d user_nnz=%d%n",
                bulk ? "bulk" : "algebraic", n, n + rows, rows, rows * perRow);
        System.out.printf("build_ms=%.1f%n", tBuild / 1e6);
        System.out.printf("allocated_bytes_build=%d%n", allocated);
        System.out.printf("solve_ms=%.1f%n", tSolve / 1e6);
        System.out.printf("peak_heap_mb=%.1f%n", peak / 1048576.0);
        System.out.printf("objective=%.6f status=%s%n", model.getObj(), status);
    }
}
