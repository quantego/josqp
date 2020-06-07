package com.quantego.josqp;

public class CSCMatrix {
	public final int n;
	public final int m;
	public int nz;
	public int nzmax;
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
	
	public CSCMatrix(int m, int n, int nzmax, boolean values, boolean triplet) {
		this.n = n;
		this.m = m;
		this.nz = triplet ? 0 : -1;
		this.nzmax = Math.max(nzmax,  1);
		this.Ap = new int[triplet ? nzmax : n+1];
		this.Ai = new int[nzmax];
		this.Ax = values ? new double[nzmax] : null;
	}
	
	public static CSCMatrix triplet_to_csc(CSCMatrix T, int[] TtoC) {
		int m = T.m, n=T.n, nz=T.nz, p, k; 
		int[] Ti = T.Ai, Tj=T.Ap;
		double[] Tx = T.Ax;
		CSCMatrix C = new CSCMatrix(m,n,nz,T.Ax!=null,false);
		int[] Cp = C.Ap;
		int[] Ci = C.Ai;
		double[] Cx = C.Ax;
		int[] w = new int[n];

		  for (k = 0; k < nz; k++) w[Tj[k]]++;  /* column counts */
		  csc_cumsum(Cp, w, n);                 /* column pointers */

		  for (k = 0; k < nz; k++) {
		    Ci[p = w[Tj[k]]++] = Ti[k];         /* A(i,j) is the pth entry in C */

		    if (Cx!=null) {
		      Cx[p] = Tx[k];

		      if (TtoC != null) TtoC[k] = p;  // Assign vector of indices
		    }
		  }
		  return C;     /* success; free w and return C */
	}
	
	public static CSCMatrix triplet_to_csr(CSCMatrix T, int[] TtoC) {
		int m = T.m, n=T.n, nz=T.nz, p, k; 
		int[] Ti = T.Ai, Tj=T.Ap;
		double[] Tx = T.Ax;
		CSCMatrix C = new CSCMatrix(m,n,nz,T.Ax!=null,false);
		int[] Cp = C.Ap;
		int[] Cj = C.Ai;
		double[] Cx = C.Ax;
		int[] w = new int[n];

		for (k = 0; k < nz; k++) w[Ti[k]]++;  /* row counts */
		  csc_cumsum(Cp, w, m);                 /* row pointers */

		  for (k = 0; k < nz; k++) {
		    Cj[p = w[Ti[k]]++] = Tj[k];         /* A(i,j) is the pth entry in C */

		    if (Cx!=null) {
		      Cx[p] = Tx[k];

		      if (TtoC != null) TtoC[k] = p;  // Assign vector of indices
		    }
		  }
		  return C;     /* success; free w and return C */
	}
	
	public static int csc_cumsum(int[] p, int[] c, int n) {
		  int nz = 0;
		  for (int i = 0; i < n; i++) {
		    p[i] = nz;
		    nz  += c[i];
		    c[i] = p[i];
		  }
		  p[n] = nz;
		  return nz;
	}
	
	public static int[] csc_pinv(int[] p, int n) {
		int[] pinv = new int[n];
		for (int k=0; k<n; k++)
			pinv[p[k]] = k;
		return pinv;
	}
	
	public static CSCMatrix csc_symperm(CSCMatrix A, int[] pinv, int[] AtoC, boolean values) {
		  int i, j, p, q, i2, j2, n;
		  int[] Ap, Ai, Cp, Ci, w;
		  double[] Cx, Ax;
		  CSCMatrix  C;

		  n  = A.n;
		  Ap = A.Ap;
		  Ai = A.Ai;
		  Ax = A.Ax;
		  C  = new CSCMatrix(n, n, Ap[n], values && (Ax != null),false);                                /* alloc result*/
		  w = new int[n];      /* get workspace */


		  Cp = C.Ap;
		  Ci = C.Ai;
		  Cx = C.Ax;

		  for (j = 0; j < n; j++)    /* count entries in each column of C */
		  {
		    j2 = pinv!=null ? pinv[j] : j; /* column j of A is column j2 of C */

		    for (p = Ap[j]; p < Ap[j + 1]; p++) {
		      i = Ai[p];

		      if (i > j) continue;     /* skip lower triangular part of A */
		      i2 =  pinv!=null ? pinv[i] : i; /* row i of A is row i2 of C */
		      w[Math.max(i2, j2)]++;      /* column count of C */
		    }
		  }
		  csc_cumsum(Cp, w, n);        /* compute column pointers of C */

		  for (j = 0; j < n; j++) {
		    j2 = pinv!=null ? pinv[j] : j;   /* column j of A is column j2 of C */

		    for (p = Ap[j]; p < Ap[j + 1]; p++) {
		      i = Ai[p];

		      if (i > j) continue;                             /* skip lower triangular
		                                                          part of A*/
		      i2                         = pinv!=null ? pinv[i] : i; /* row i of A is row i2
		                                                          of C */
		      Ci[q = w[Math.max(i2, j2)]++] = Math.min(i2, j2);

		      if (Cx!=null) Cx[q] = Ax[p];

		      if (AtoC!=null) { // If vector AtoC passed, store values of the mappings
		        AtoC[p] = q;
		      }
		    }
		  }
		  return C; /* success; free workspace, return C */
		}
	
	public static CSCMatrix copy_csc_mat(CSCMatrix A) {
		CSCMatrix B = new CSCMatrix(A.m, A.n, A.Ap[A.n], true, false);
		  LinAlg.prea_int_vec_copy(A.Ap, B.Ap);
		  LinAlg.prea_int_vec_copy(A.Ai, B.Ai);
		  LinAlg.prea_vec_copy(A.Ax, B.Ax);

		  return B;
		}
	
	public static void prea_copy_csc_mat(CSCMatrix A, CSCMatrix B) {
		LinAlg.prea_int_vec_copy(A.Ap, B.Ap);
		LinAlg.prea_int_vec_copy(A.Ai, B.Ai);
		LinAlg.prea_vec_copy(A.Ax, B.Ax);
		B.nzmax = A.nzmax;
		}
	
	CSCMatrix csc_to_triu(CSCMatrix M) {
		CSCMatrix  M_trip;    // Matrix in triplet format
		CSCMatrix  M_triu;    // Resulting upper triangular matrix
		  int nnzorigM;  // Number of nonzeros from original matrix M
		  int nnzmaxM;   // Estimated maximum number of elements of upper triangular M
		  int n;         // Dimension of M
		  int ptr, i, j; // Counters for (i,j) and index in M
		  int z_M = 0;   // Counter for elements in M_trip


		  // Check if matrix is square
		  if (M.m != M.n) {
			  throw new IllegalStateException("Matrix is not square.");
		  }
		  n = M.n;

		  // Get number of nonzeros full M
		  nnzorigM = M.Ap[n];

		  // Estimate nnzmaxM
		  // Number of nonzero elements in original M + diagonal part.
		  // . Full matrix M as input: estimate is half the number of total elements +
		  // diagonal = .5 * (nnzorigM + n)
		  // . Upper triangular matrix M as input: estimate is the number of total
		  // elements + diagonal = nnzorigM + n
		  // The maximum between the two is nnzorigM + n
		  nnzmaxM = nnzorigM + n;

		  // OLD
		  // nnzmaxM = n*(n+1)/2;  // Full upper triangular matrix (This version
		  // allocates too much memory!)
		  // nnzmaxM = .5 * (nnzorigM + n);  // half of the total elements + diagonal

		  // Allocate M_trip
		  M_trip = new CSCMatrix(n, n, nnzmaxM, true, true); // Triplet format


		  // Fill M_trip with only elements in M which are in the upper triangular
		  for (j = 0; j < n; j++) { // Cycle over columns
		    for (ptr = M.Ap[j]; ptr < M.Ap[j + 1]; ptr++) {
		      // Get row index
		      i = M.Ai[ptr];

		      // Assign element only if in the upper triangular
		      if (i <= j) {
		        // c_print("\nM(%i, %i) = %.4f", M.i[ptr], j, M.x[ptr]);

		        M_trip.Ai[z_M] = i;
		        M_trip.Ap[z_M] = j;
		        M_trip.Ax[z_M] = M.Ax[ptr];

		        // Increase counter for the number of elements
		        z_M++;
		      }
		    }
		  }

		  // Set number of nonzeros
		  M_trip.nz = z_M;

		  // Convert triplet matrix to csc format
		  M_triu = triplet_to_csc(M_trip, null);

		  // Assign number of nonzeros of full matrix to triu M
		  M_triu.nzmax = nnzmaxM;


		  // Return matrix in triplet form
		  return M_triu;
		}
}