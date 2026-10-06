package com.quantego.josqp;

import java.util.Locale;
import java.util.logging.Level;

/** Repeated numeric factorization of an unchanged real KKT matrix; not a full-solver benchmark. */
public final class FactorReplayBenchmark {
    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.ROOT); OSQP.LOG.setLevel(Level.OFF);
        OSQP.Settings settings=new OSQP.Settings();settings.verbose=false;
        OSQP solver=new OSQP(Parser.readQmps(args[0]).getData(),settings);
        LinSys lin=solver.work.linsys_solver;QDLDL q=lin.qdldl;
        long deadline=System.nanoTime()+2_000_000_000L;
        int warmups=0;
        do { q.refactor(lin.kkt);warmups++; } while(System.nanoTime()<deadline || warmups<3);
        int samples=3;
        double[] times=new double[samples];
        for(int i=0;i<samples;i++) { long t=System.nanoTime();q.refactor(lin.kkt);times[i]=(System.nanoTime()-t)/1e6; }
        java.security.MessageDigest sha=java.security.MessageDigest.getInstance("SHA-256");
        java.nio.ByteBuffer buffer=java.nio.ByteBuffer.allocate(8);
        for(double[] array:new double[][]{q._ldl.L.Ax,q._ldl.D,q._ldl.Dinv}) for(double v:array) {
            buffer.clear();buffer.putLong(Double.doubleToLongBits(v));sha.update(buffer.array());
        }
        StringBuilder hash=new StringBuilder();for(byte b:sha.digest()) hash.append(String.format("%02x",b&255));
        System.out.println("file,n,factor_nnz,warmups,sample,factor_ms,factor_sha256");
        for(int i=0;i<samples;i++) System.out.printf("%s,%d,%d,%d,%d,%.6f,%s%n",args[0],lin.kkt.n,q._ldl.L.Ap[lin.kkt.n],warmups,i,times[i],hash);
    }
}
