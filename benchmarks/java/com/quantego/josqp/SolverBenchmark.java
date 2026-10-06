package com.quantego.josqp;

import java.lang.management.ManagementFactory;
import java.util.Locale;
import java.util.logging.Level;

/** Fresh cold-start solves with separate parse/setup/solve timing and allocation counts. */
public final class SolverBenchmark {
    private static volatile double sink;

    public static void main(String[] args) throws java.io.IOException {
        Locale.setDefault(Locale.ROOT);
        OSQP.LOG.setLevel(Level.OFF);
        String file = args[0];
        int warmup = Integer.parseInt(args[1]), samples = Integer.parseInt(args[2]);
        int maxIter = Integer.parseInt(args[3]);
        long parseStart = System.nanoTime();
        OSQP.Data original = Parser.readQmps(file).getData();
        double parseMs = (System.nanoTime() - parseStart) / 1e6;
        System.err.printf("DATA n=%d m=%d nnz=%d parse_ms=%.3f%n", original.n, original.m,
                original.A.Ap[original.n] + original.P.Ap[original.n], parseMs);
        com.sun.management.ThreadMXBean bean = (com.sun.management.ThreadMXBean) ManagementFactory.getThreadMXBean();
        boolean allocations = bean.isThreadAllocatedMemorySupported();
        if (allocations) bean.setThreadAllocatedMemoryEnabled(true);
        long thread = Thread.currentThread().getId();
        long warmupDeadline = System.nanoTime() + (long) (Double.parseDouble(
                System.getProperty("josqp.bench.warmupSeconds", "0")) * 1e9);
        System.out.println("sample,file,n,m,nnz,parse_ms,setup_ms,solve_ms,allocated_bytes,status,iterations,objective,primal_residual,dual_residual,x_hash,y_hash");
        for (int sample = -warmup; sample < samples; sample++) {
            OSQP.Settings settings = new OSQP.Settings();
            settings.verbose = false;
            settings.warm_start = false;
            settings.adaptive_rho_interval = 100;
            settings.max_iter = maxIter;
            settings.time_limit = Double.parseDouble(System.getProperty("josqp.bench.timeLimit", "Infinity"));
            settings.polish = false;
            long start = System.nanoTime();
            OSQP solver = new OSQP(new OSQP.Data(original), settings);
            double setupMs = (System.nanoTime() - start) / 1e6;
            long before = allocations ? bean.getThreadAllocatedBytes(thread) : -1;
            start = System.nanoTime();
            OSQP.Status status = solver.solve();
            double solveMs = (System.nanoTime() - start) / 1e6;
            long allocated = allocations ? bean.getThreadAllocatedBytes(thread) - before : -1;
            OSQP.Info info = solver.work.info;
            sink = info.obj_val;
            if (sample >= 0) {
                System.out.printf("%d,%s,%d,%d,%d,%.6f,%.6f,%.6f,%d,%s,%d,%.17g,%.17g,%.17g,%d,%d%n",
                        sample, file, original.n, original.m, original.A.Ap[original.n] + original.P.Ap[original.n],
                        parseMs, setupMs, solveMs, allocated, status, info.iter, info.obj_val, info.pri_res, info.dua_res,
                        java.util.Arrays.hashCode(solver.work.solution.x), java.util.Arrays.hashCode(solver.work.solution.y));
            }
            String dump = System.getProperty("josqp.bench.solutionFile");
            if (sample == 0 && dump != null) {
                try (java.io.DataOutputStream out = new java.io.DataOutputStream(new java.io.FileOutputStream(dump))) {
                    out.writeInt(original.n); out.writeInt(original.m);
                    for (double x : solver.work.solution.x) out.writeLong(Double.doubleToLongBits(x));
                    for (double y : solver.work.solution.y) out.writeLong(Double.doubleToLongBits(y));
                }
            }
            if (sample == -1 && System.nanoTime() < warmupDeadline) sample--;
        }
    }
}
