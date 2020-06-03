package com.quantego.josqp;

public class CSC {
	
	final int nzmax; ///< maximum number of entries
	final int m;     ///< number of rows
	final int n;     ///< number of columns
	final int[] p;     ///< column pointers (size n+1); col indices (size nzmax) start from 0 when using triplet format (direct KKT matrix formation)
	final int[] i;     ///< row indices, size nzmax starting from 0
	final double[] x;     ///< numerical values, size nzmax
	final int nz;    ///< number of entries in triplet matrix, -1 for csc
	
	public CSC(int nzmax, int m, int n, int[] p, int[] i, double[] x, int nz) {
		this.nzmax = nzmax;
		this.m = m;
		this.n = n;
		this.p = p;
		this.i = i;
		this.x = x;
		this.nz = nz;
	}

}
