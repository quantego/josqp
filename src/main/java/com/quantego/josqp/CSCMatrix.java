package com.quantego.josqp;

public class CSCMatrix {
	public final int n;
	public final int m;
	public final int nz;
	public final int nzmax;
	public final int[] Ap;
    public final int[] Ai;
    public final double[] Ax;
	public CSCMatrix(int m, int n, int nzmax, int[] ap, int[] ai, double[] ax) {
		this.n = n;
		this.m = m;
		this.nz = ax.length;
		this.nzmax = nzmax;
		this.Ap = ap;
		this.Ai = ai;
		this.Ax = ax;
	}

}