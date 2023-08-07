package com.quantego.josqp;

import java.lang.reflect.Array;
import java.util.Arrays;

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
		  C.nz = Cx.length;
		  return C;     /* success; free w and return C */
	}
	
	public static CSCMatrix triplet_to_csc(int m, int n, int nz, int[] Ti, int[] Tj, double[] Tx, int[] TtoC) {
		int p, k, i, iChk;
		CSCMatrix C = new CSCMatrix(m,n,nz,Tx!=null,false);
		int[] Cp = C.Ap;
		int[] Ci = C.Ai;
		double[] Cx = C.Ax;
		int[] w = new int[n];
		int iooRcvCap = 1, iooRcvSize = 0;
		int[] iooRcv = new int[iooRcvCap];

		// determining the column pointers
		for (k = 0; k < nz; k++) w[Tj[k]]++;  /* column counts */
		csc_cumsum(Cp, w, n);                 /* column pointers */

		// determining the row indices
		for (k = 0; k < nz; k++) {
			if (k == 0 || Tj[k] != Tj[k-1]) { // if colum number is changed (assuming the RCV data is already sorted based on the column number
				if (iooRcvCap < Cp[Tj[k] + 1] - Cp[Tj[k]]) {
					iooRcvCap = 2 * (Cp[Tj[k] + 1] - Cp[Tj[k]]);
					iooRcv = new int[iooRcvCap];
				} else {
					Arrays.fill(iooRcv, 0);
				}
				iooRcvSize = Cp[Tj[k] + 1] - Cp[Tj[k]];
				// sorting based on the row number (for each column)
				for (i = 1; i < iooRcvSize; i++) {
					iChk = i - 1;
					while (iChk != -1 && Ti[k + iooRcv[iChk]] > Ti[k + i]) {
						iooRcv[iChk + 1] = iooRcv[iChk];
						iChk--;
					}
					iooRcv[iChk + 1] = i;
				}
				for (i = 0; i < iooRcvSize; i++) {
					Ci[k + i] = Ti[k + iooRcv[i]];
					if (Cx != null) {
						Cx[k + i] = Tx[k + iooRcv[i]];
					}
				}
			}
			//Ci[p = w[Tj[k]]++] = Ti[k];         /* A(i,j) is the pth entry in C */
		  //if (Cx!=null) {
		  //  Cx[p] = Tx[k];
			//	if (TtoC != null) TtoC[k] = p;  // Assign vector of indices
			//}
		}
		C.nz = Cx.length;
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
		  C.nz = Cx.length;
		  return C; /* success; free workspace, return C */
		}
	
	public static CSCMatrix copy_csc_mat(CSCMatrix A) {
		CSCMatrix B = new CSCMatrix(A.m, A.n, A.Ap[A.n], true, false);
		  LinAlg.prea_int_vec_copy(A.Ap, B.Ap, A.n+1);
		  LinAlg.prea_int_vec_copy(A.Ai, B.Ai,A.Ap[A.n]);
		  LinAlg.prea_vec_copy(A.Ax, B.Ax, A.Ap[A.n]);
		  B.nz = A.nz;
		  return B;
		}
	
	public static void prea_copy_csc_mat(CSCMatrix A, CSCMatrix B) {
		LinAlg.prea_int_vec_copy(A.Ap, B.Ap, A.n + 1);
		LinAlg.prea_int_vec_copy(A.Ai, B.Ai, A.Ap[A.n]);
		LinAlg.prea_vec_copy(A.Ax, B.Ax, A.Ap[A.n]);
		B.nzmax = A.nzmax;
		}
	
	public static CSCMatrix csc_to_triu(CSCMatrix M) {
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
	
	public String toString() {
		double[][] mat = new double[m][n];
		for (int j = 0; j < this.n; j++) {
		    for (int p = Ap[j]; p < Ap[j + 1]; p++) {
		      int i = Ai[p];
		      mat[i][j] = Ax[p];
		    }
		}
		StringBuilder sb = new StringBuilder();
		for (int i=0;i<m;i++)
			sb.append(Arrays.toString(mat[i])).append("\n");
		return sb.toString();
	}
}