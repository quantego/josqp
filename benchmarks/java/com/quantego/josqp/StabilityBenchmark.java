package com.quantego.josqp;

import java.util.Locale;
import java.util.logging.Level;

/** Independent checks on original, unscaled input after a timed cold solve. */
public final class StabilityBenchmark {
    static double norm(double[] a) {
        double n = 0;
        for (double v : a) { if (!Double.isFinite(v)) return Double.NaN; n = Math.max(n, Math.abs(v)); }
        return n;
    }
    public static void main(String[] args) {
        Locale.setDefault(Locale.ROOT);
        OSQP.LOG.setLevel(Level.OFF);
        OSQP.Data original = Parser.readQmps(args[0]).getData();
        OSQP.Settings s = new OSQP.Settings();
        s.verbose = false; s.warm_start = false;
        // 100 keeps earlier audit runs reproducible; -1 leaves the solver default in place.
        int rhoInterval = Integer.getInteger("audit.rhoInterval", 100);
        if (rhoInterval >= 0) s.adaptive_rho_interval = rhoInterval;
        s.max_iter = Integer.parseInt(System.getProperty("audit.maxIter", "100000"));
        s.time_limit = Double.parseDouble(System.getProperty("audit.seconds", "10"));
        s.eps_abs = s.eps_rel = Double.parseDouble(System.getProperty("audit.eps", "0.001"));
        s.eps_rel = Double.parseDouble(System.getProperty("audit.epsRel", Double.toString(s.eps_rel)));
        s.eps_prim_inf = s.eps_dual_inf = Double.parseDouble(System.getProperty("audit.infeasEps", "0.0001"));
        s.polish = Boolean.parseBoolean(System.getProperty("audit.polish", "false"));
        s.scaling = Integer.parseInt(System.getProperty("audit.scaling", "10"));
        // Reflection keeps the harness compilable against frozen older solvers.
        try { OSQP.Settings.class.getField("check_dualgap").setBoolean(s,
                Boolean.parseBoolean(System.getProperty("audit.gap", "true"))); }
        catch (NoSuchFieldException ignored) { }
        catch (IllegalAccessException ex) { throw new IllegalStateException(ex); }
        try {
            OSQP.Settings.class.getField("check_linear_residual").setBoolean(s, Boolean.getBoolean("audit.linear"));
            OSQP.Settings.class.getField("linear_refine_iter").setInt(s, Integer.getInteger("audit.refine",0));
            OSQP.Settings.class.getField("adaptive_rho_kkt").setBoolean(s, Boolean.getBoolean("audit.kktRho"));
            // 2 and 0.1 are the solver defaults; 1 / 1.0 reproduce upstream OSQP's verdicts exactly.
            OSQP.Settings.class.getField("infeasibility_streak").setInt(s, Integer.getInteger("audit.infeasStreak",2));
            OSQP.Settings.class.getField("infeasibility_margin").setDouble(s,
                    Double.parseDouble(System.getProperty("audit.infeasMargin", "0.1")));
        } catch (NoSuchFieldException ignored) { }
        catch (IllegalAccessException ex) { throw new IllegalStateException(ex); }
        OSQP.Data input = new OSQP.Data(original);
        double costScale = Double.parseDouble(System.getProperty("audit.costScale", "1"));
        for (int i=0;i<input.P.Ap[input.n];i++) input.P.Ax[i] *= costScale;
        for (int i=0;i<input.n;i++) input.q[i] *= costScale;
        if (Boolean.getBoolean("audit.warmup")) {
            // One untimed solve on a throwaway solver so the JIT is warm; the measured
            // solver below is still built and solved from a fresh state.
            OSQP warm = new OSQP(new OSQP.Data(original), s);
            warm.solve();
        }
        long t = System.nanoTime();
        OSQP solver = new OSQP(input,s);
        double setup = (System.nanoTime()-t)/1e6;
        t=System.nanoTime(); solver.solve(); double solve=(System.nanoTime()-t)/1e6;
        double[] x=solver.getPrimalSolution(), y=solver.getDualSolution();
        double[] px=new double[original.n], ax=new double[original.m], aty=new double[original.n];
        double xpx=0,qx=0,support=0,pri=0,comp=0,cone=0;
        // Explicit sparse multiplication, independent of LinAlg/NumericKernels.
        for(int c=0;c<original.n;c++) {
            qx += original.q[c]*x[c];
            for(int k=original.P.Ap[c];k<original.P.Ap[c+1];k++) {
                int r=original.P.Ai[k]; double v=original.P.Ax[k];
                px[r]+=v*x[c]; if(r!=c) px[c]+=v*x[r];
            }
            for(int k=original.A.Ap[c];k<original.A.Ap[c+1];k++) {
                int r=original.A.Ai[k]; double v=original.A.Ax[k];
                ax[r]+=v*x[c]; aty[c]+=v*y[r]/costScale;
            }
        }
        for(int i=0;i<original.n;i++) xpx+=x[i]*px[i];
        for(int i=0;i<original.m;i++) {
            double yi=y[i]/costScale;
            if(original.l[i]>-OSQP.OSQP_INFTY) pri=Math.max(pri,original.l[i]-ax[i]);
            if(original.u[i]<OSQP.OSQP_INFTY) pri=Math.max(pri,ax[i]-original.u[i]);
            double bound=yi>0?original.u[i]:original.l[i];
            if(yi!=0 && Math.abs(bound)<OSQP.OSQP_INFTY) {
                support+=yi*bound; comp=Math.max(comp,Math.abs(yi*(bound-ax[i])));
            } else if(yi!=0) cone=Math.max(cone,Math.abs(yi));
        }
        double dual=0;
        for(int i=0;i<original.n;i++) dual=Math.max(dual,Math.abs(px[i]+original.q[i]+aty[i]));
        boolean finite=Double.isFinite(norm(x)) && Double.isFinite(norm(y));
        double independentObj=0.5*xpx+qx+original.offset;
        double reportedObj=(solver.getObjectiveValue()-original.offset)/costScale+original.offset;
        String export = System.getProperty("audit.export");
        if (export != null) exportData(original, export);
        System.out.println("file,n,m,offset,cost_scale,status,status_polish,iterations,rho_updates,setup_ms,solve_ms,objective,independent_objective,primal_violation,stationarity,complementarity,dual_cone_error,gap_expression,finite,x_hash,y_hash,linear_calls,linear_backward_error,linear_corrections,linear_scaled_residual");
        System.out.printf("%s,%d,%d,%.17g,%.17g,%s,%d,%d,%d,%.6f,%.6f,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%s,%d,%d,%.0f,%.17g,%.0f,%.17g%n",
            args[0],original.n,original.m,original.offset,costScale,solver.getInfo().status,solver.getInfo().status_polish,solver.getInfo().iter,
            solver.getInfo().rho_updates,setup,solve,reportedObj,independentObj,pri,dual,comp,cone,xpx+qx+support,finite,
            java.util.Arrays.hashCode(x),java.util.Arrays.hashCode(y),
            diagnostic(solver,"linearCalls"),diagnostic(solver,"maxBackwardError"),diagnostic(solver,"linearCorrections"),diagnostic(solver,"maxScaledResidual"));
    }

    static double diagnostic(OSQP solver, String field) {
        try {
            java.lang.reflect.Field f=LinSys.class.getDeclaredField(field); f.setAccessible(true);
            return ((Number)f.get(solver.work.linsys_solver)).doubleValue();
        } catch (NoSuchFieldException ex) { return Double.NaN; }
        catch (IllegalAccessException ex) { throw new IllegalStateException(ex); }
    }

    /** Exact original matrices for native solvers, before JOSQP scaling. Big-endian doubles/ints. */
    static void exportData(OSQP.Data d, String path) {
        try (java.io.DataOutputStream out = new java.io.DataOutputStream(new java.io.FileOutputStream(path))) {
            out.writeInt(0x4a515031); out.writeInt(d.n); out.writeInt(d.m); out.writeDouble(d.offset);
            for (CSCMatrix a : new CSCMatrix[]{d.P,d.A}) {
                int nz=a.Ap[a.n]; out.writeInt(nz);
                for (int i=0;i<=a.n;i++) out.writeInt(a.Ap[i]);
                for (int i=0;i<nz;i++) out.writeInt(a.Ai[i]);
                for (int i=0;i<nz;i++) out.writeDouble(a.Ax[i]);
            }
            for(double[] v : new double[][]{d.q,d.l,d.u}) for(double z:v) out.writeDouble(z);
        } catch (java.io.IOException ex) { throw new java.io.UncheckedIOException(ex); }
    }
}
