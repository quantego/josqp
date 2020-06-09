package com.quantego.josqp;

public class KKT {
	


	public static CSCMatrix form_KKT(CSCMatrix P,
			CSCMatrix A,
	              int       format,
	              double     param1,
	              double[]    param2,
	              int[]      PtoKKT,
	              int[]     AtoKKT,
	              int[][]     Pdiag_idx,
	              int[]      Pdiag_n,
	              int[]      param2toKKT) {
	  int  nKKT, nnzKKTmax; // Size, number of nonzeros and max number of nonzeros
	                          // in KKT matrix
	  CSCMatrix   KKT_trip, KKT;  // KKT matrix in triplet format and CSC format
	  int  ptr, i, j;       // Counters for elements (i,j) and index pointer
	  int  zKKT = 0;        // Counter for total number of elements in P and in
	                          // KKT
	  int[] KKT_TtoC;        // Pointer to vector mapping from KKT in triplet form
	                          // to CSC

	  // Get matrix dimensions
	  nKKT = P.m + A.m;

	  // Get maximum number of nonzero elements (only upper triangular part)
	  nnzKKTmax = P.Ap[P.n] + // Number of elements in P
	              P.m +       // Number of elements in param1 * I
	              A.Ap[A.n] + // Number of nonzeros in A
	              A.m;        // Number of elements in - diag(param2)

	  // Preallocate KKT matrix in triplet format
	  KKT_trip = new CSCMatrix(nKKT, nKKT, nnzKKTmax, true, true);


	  // Allocate vector of indices on the diagonal. Worst case it has m elements
	  if (Pdiag_idx != null) {
	    Pdiag_idx[0] = new int[P.m];
	    Pdiag_n  [0]   = 0; // Set 0 diagonal elements to start
	  }

	  // Allocate Triplet matrices
	  // P + param1 I
	  for (j = 0; j < P.n; j++) { // cycle over columns
	    // No elements in column j => add diagonal element param1
	    if (P.Ap[j] == P.Ap[j + 1]) {
	      KKT_trip.Ai[zKKT] = j;
	      KKT_trip.Ap[zKKT] = j;
	      KKT_trip.Ax[zKKT] = param1;
	      zKKT++;
	    }

	    for (ptr = P.Ap[j]; ptr < P.Ap[j + 1]; ptr++) { // cycle over rows
	      // Get current row
	      i = P.Ai[ptr];

	      // Add element of P
	      KKT_trip.Ai[zKKT] = i;
	      KKT_trip.Ap[zKKT] = j;
	      KKT_trip.Ax[zKKT] = P.Ax[ptr];

	      if (PtoKKT != null) PtoKKT[ptr] = zKKT;  // Update index from P to
	                                                    // KKTtrip

	      if (i == j) {                                 // P has a diagonal element,
	                                                    // add param1
	        KKT_trip.Ax[zKKT] += param1;

	        // If index vector pointer supplied . Store the index
	        if (Pdiag_idx != null) {
	          Pdiag_idx[0][Pdiag_n[0]] = ptr;
	          Pdiag_n[0]++;
	        }
	      }
	      zKKT++;

	      // Add diagonal param1 in case
	      if ((i < j) &&                  // Diagonal element not reached
	          (ptr + 1 == P.Ap[j + 1])) { // last element of column j
	        // Add diagonal element param1
	        KKT_trip.Ai[zKKT] = j;
	        KKT_trip.Ap[zKKT] = j;
	        KKT_trip.Ax[zKKT] = param1;
	        zKKT++;
	      }
	    }
	  }

	  if (Pdiag_idx != null) {
	    // Realloc Pdiag_idx so that it contains exactly *Pdiag_n diagonal elements
		  Pdiag_idx[0] = new int[Pdiag_n[0]];
	  }


	  // A' at top right
	  for (j = 0; j < A.n; j++) {                      // Cycle over columns of A
	    for (ptr = A.Ap[j]; ptr < A.Ap[j + 1]; ptr++) {
	      KKT_trip.Ap[zKKT] = P.m + A.Ai[ptr];         // Assign column index from
	                                                    // row index of A
	      KKT_trip.Ai[zKKT] = j;                        // Assign row index from
	                                                    // column index of A
	      KKT_trip.Ax[zKKT] = A.Ax[ptr];                // Assign A value element

	      if (AtoKKT != null) AtoKKT[ptr] = zKKT;  // Update index from A to
	                                                    // KKTtrip
	      zKKT++;
	    }
	  }

	  // - diag(param2) at bottom right
	  for (j = 0; j < A.m; j++) {
	    KKT_trip.Ai[zKKT] = j + P.n;
	    KKT_trip.Ap[zKKT] = j + P.n;
	    KKT_trip.Ax[zKKT] = -param2[j];

	    if (param2toKKT != null) param2toKKT[j] = zKKT;  // Update index from
	                                                          // param2 to KKTtrip
	    zKKT++;
	  }

	  // Allocate number of nonzeros
	  KKT_trip.nz = zKKT;

	  // Convert triplet matrix to csc format
	  if (PtoKKT==null && AtoKKT==null && param2toKKT==null) {
	    // If no index vectors passed, do not store KKT mapping from Trip to CSC/CSR
	    if (format == 0) KKT = CSCMatrix.triplet_to_csc(KKT_trip, null);
	    else KKT = CSCMatrix.triplet_to_csr(KKT_trip, null);
	  }
	  else {
	    // Allocate vector of indices from triplet to csc
	    KKT_TtoC = new int[zKKT];


	    // Store KKT mapping from Trip to CSC/CSR
	    if (format == 0)
	      KKT = CSCMatrix.triplet_to_csc(KKT_trip, KKT_TtoC);
	    else
	      KKT = CSCMatrix.triplet_to_csr(KKT_trip, KKT_TtoC);

	    // Update vectors of indices from P, A, param2 to KKT (now in CSC format)
	    if (PtoKKT != null) {
	      for (i = 0; i < P.Ap[P.n]; i++) {
	        PtoKKT[i] = KKT_TtoC[PtoKKT[i]];
	      }
	    }

	    if (AtoKKT != null) {
	      for (i = 0; i < A.Ap[A.n]; i++) {
	        AtoKKT[i] = KKT_TtoC[AtoKKT[i]];
	      }
	    }

	    if (param2toKKT != null) {
	      for (i = 0; i < A.m; i++) {
	        param2toKKT[i] = KKT_TtoC[param2toKKT[i]];
	      }
	    }


	  }


	  return KKT;
	}


	public static void update_KKT_P(CSCMatrix KKT,
			CSCMatrix P,
	        int[] PtoKKT,
	        final double param1,
	        final  int[] Pdiag_idx,
	        final int Pdiag_n) {
	  int i, j; // Iterations

	  // Update elements of KKT using P
	  for (i = 0; i < P.Ap[P.n]; i++) {
	    KKT.Ax[PtoKKT[i]] = P.Ax[i];
	  }

	  // Update diagonal elements of KKT by adding sigma
	  for (i = 0; i < Pdiag_n; i++) {
	    j                  = Pdiag_idx[i]; // Extract index of the element on the
	                                       // diagonal
	    KKT.Ax[PtoKKT[j]] += param1;
	  }
	}

	public static void update_KKT_A(CSCMatrix KKT, final CSCMatrix A, final int[] AtoKKT) {
	  int i; // Iterations

	  // Update elements of KKT using A
	  for (i = 0; i < A.Ap[A.n]; i++) {
	    KKT.Ax[AtoKKT[i]] = A.Ax[i];
	  }
	}

	public static void update_KKT_param2(CSCMatrix KKT, final double[] param2,
	                       final int[] param2toKKT, final int m) {
	  int i; // Iterations

	  // Update elements of KKT using param2
	  for (i = 0; i < m; i++) {
	    KKT.Ax[param2toKKT[i]] = -param2[i];
	  }
	}



}
