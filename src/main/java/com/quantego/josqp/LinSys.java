package com.quantego.josqp;

public class LinSys  {
	


	public static CSCMatrix permute_KKT(CSCMatrix KKT, LinSys p, int Pnz, int Anz, int m, int[] PtoKKT, int[] AtoKKT, int[] rhotoKKT){
	    AMD.Status amd_status;
	    int[] Pinv;
	    CSCMatrix KKT_temp;
	    int[] KtoPKPt;
	    int i; // Indexing


	    // Compute permutation matrix P using AMD

	    amd_status = AMD.amd_order(KKT.n, KKT.Ap, KKT.Ai, p.P);
	    if (amd_status == AMD.Status.AMD_INVALID) {
	        // Free Amd info and return an error
	        throw new IllegalStateException("AMD Error");
	    }


	    // Inverse of the permutation vector
	    Pinv = CSCMatrix.csc_pinv(p.P, KKT.n);

	    // Permute KKT matrix
	    if (PtoKKT==null && AtoKKT==null && rhotoKKT==null){  // No vectors to be stored
	        // Assign values of mapping
	        KKT_temp = CSCMatrix.csc_symperm(KKT, Pinv, null, true);
	    }
	    else {
	        // Allocate vector of mappings from unpermuted to permuted
	        KtoPKPt = new int[KKT.Ap[KKT.n]];
	        KKT_temp = CSCMatrix.csc_symperm(KKT, Pinv, KtoPKPt, true);

	        // Update vectors PtoKKT, AtoKKT and rhotoKKT
	        if (PtoKKT!=null){
	            for (i = 0; i < Pnz; i++){
	                PtoKKT[i] = KtoPKPt[PtoKKT[i]];
	            }
	        }
	        if (AtoKKT!=null){
	            for (i = 0; i < Anz; i++){
	                AtoKKT[i] = KtoPKPt[AtoKKT[i]];
	            }
	        }
	        if (rhotoKKT!=null){
	            for (i = 0; i < m; i++){
	                rhotoKKT[i] = KtoPKPt[rhotoKKT[i]];
	            }
	        }

	    }

	    // Cleanup
	    // Free previous KKT matrix and assign pointer to new one
	   return KKT_temp;
	}
	
	QDLDL qdldl;
	int n;
	int m;
	double sigma;
	boolean polish;
	CSCMatrix L;
	double[] rho_inv_vec;
	CSCMatrix kkt;
	public int[] P;
	double[] sol;
	double[] bp;
	int[] PtoKKT;
	int[] AtoKKT;
	int[][] Pdiag_idx = new int[1][];
	int[] Pdiag_n = new int[1];
	int[] rhotoKKT;

	// Initialize LDL Factorization structure
	public LinSys(CSCMatrix P, CSCMatrix A, double sigma, double[] rho_vec, boolean polish){

	    // Define Variables
	    CSCMatrix KKT_temp;     // Temporary KKT pointer
	    int i;            // Loop counter
	    int n_plus_m;     // Define n_plus_m dimension

	    // Allocate private structure to store KKT factorization
	    qdldl = new QDLDL() ;

	    // Size of KKT
	    this.n = P.n;
	    this.m = A.m;
	    n_plus_m = this.n + this.m;

	    // Sigma parameter
	    this.sigma = sigma;

	    // Polishing flag
	    this.polish = polish;



	    // Sparse matrix L (lower triangular)
	    // NB: We don not allocate L completely (CSC elements)
	    //      L will be allocated during the factorization depending on the
	    //      resulting number of elements.
//	    this.L = new CSCMatrix(n_plus_m,n_plus_m,-1);

	    // Diagonal matrix stored as a vector D
//	    s.Dinv = (QDLDL_float *)c_malloc(sizeof(QDLDL_float) * n_plus_m);
//	    s.D    = (QDLDL_float *)c_malloc(sizeof(QDLDL_float) * n_plus_m);

	    // Permutation vector P
	    this.P    = new int[n_plus_m];

	    // Working vector
	    this.bp   = new double[n_plus_m];

	    // Solution vector
	    this.sol  = new double[n_plus_m];

	    // Parameter vector
	    this.rho_inv_vec = new double[this.m];

	    // Elimination tree workspace
//	    s.etree = (QDLDL_int *)c_malloc(n_plus_m * sizeof(QDLDL_int));
//	    s.Lnz   = (QDLDL_int *)c_malloc(n_plus_m * sizeof(QDLDL_int));

	    // Preallocate L matrix (Lx and Li are sparsity dependent)
//	    s.L.p = (int *)c_malloc((n_plus_m+1) * sizeof(QDLDL_int));

	    // Lx and Li are sparsity dependent, so set them to
	    // null initially so we don't try to free them prematurely
//	    s.L.i = OSQP_NULL;
//	    s.L.x = OSQP_NULL;

	    // Preallocate workspace
//	    s.iwork = (QDLDL_int *)c_malloc(sizeof(QDLDL_int)*(3*n_plus_m));
//	    s.bwork = (QDLDL_bool *)c_malloc(sizeof(QDLDL_bool)*n_plus_m);
//	    s.fwork = (QDLDL_float *)c_malloc(sizeof(QDLDL_float)*n_plus_m);

	    // Form and permute KKT matrix
	    if (polish){ // Called from polish()
	        // Use s.rho_inv_vec for storing param2 = vec(delta)
	        for (i = 0; i < A.m; i++){
	            this.rho_inv_vec[i] = sigma;
	        }

	        KKT_temp = KKT.form_KKT(P, A, 0, sigma, this.rho_inv_vec, null, null, null, null, null);

	        // Permute matrix
	        if (KKT_temp!=null)
	            KKT_temp = permute_KKT(KKT_temp, this, 0, 0, 0, null, null, null);
	    }
	    else { // Called from ADMM algorithm

	        // Allocate vectors of indices
	        this.PtoKKT = new int[P.Ap[P.n]];
	        this.AtoKKT = new int[A.Ap[A.n]];
	        this.rhotoKKT = new int[A.m];

	        // Use p.rho_inv_vec for storing param2 = rho_inv_vec
	        for (i = 0; i < A.m; i++){
	            this.rho_inv_vec[i] = 1. / rho_vec[i];
	        }

	        KKT_temp = KKT.form_KKT(P, A, 0, sigma, this.rho_inv_vec,
	                            this.PtoKKT, this.AtoKKT,
	                            this.Pdiag_idx, this.Pdiag_n, this.rhotoKKT);

	        // Permute matrix
	        if (KKT_temp!=null)
	        	KKT_temp = permute_KKT(KKT_temp, this, P.Ap[P.n], A.Ap[A.n], A.m, this.PtoKKT, this.AtoKKT, this.rhotoKKT);
	    }

	    // Check if matrix has been created

	    qdldl.factor(KKT_temp);
	    // Factorize the KKT matrix


	    if (!polish){ // If KKT passed, assign it to KKT_temp
	        this.kkt = KKT_temp;
	    }


	    // No error
	}



	// Permute x = P*b using P
	public static void permute_x(int n, double[] x, double[] b, int[] P) {
	    int j;
	    for (j = 0 ; j < n ; j++) x[j] = b[P[j]];
	}

	// Permute x = P'*b using P
	public static void permutet_x(int n, double[] x, double[] b, int[] P) {
	    int j;
	    for (j = 0 ; j < n ; j++) x[P[j]] = b[j];
	}


	void LDLSolve(double[] x, double[] b, int[] P, double[] bp) {
	    /* solves P'LDL'P x = b for x */
	    permute_x(bp.length, bp, b, P);
	    bp = qdldl.solve(bp);
	    permutet_x(bp.length, x, bp, P);

	}


	public void solve(double[] b) {
	    int j;

//	    if (this.polish) {
//	        /* stores solution to the KKT system in b */
//	        LDLSolve(b, b, this.P, this.bp);
//	    } else {
	        /* stores solution to the KKT system in s.sol */
	        LDLSolve(this.sol, b, this.P, this.bp);

	        /* copy x_tilde from s.sol */
	        for (j = 0 ; j < this.n ; j++) {
	            b[j] = this.sol[j];
	        }

	        /* compute z_tilde from b and s.sol */
	        for (j = 0 ; j < this.m ; j++) {
	            b[j + this.n] += this.rho_inv_vec[j] * this.sol[j + this.n];
	        }
//	    }

	}


	// Update private structure with new P and A
	void update_solver_matrices(CSCMatrix P, CSCMatrix A) {

	    // Update KKT matrix with new P
	    KKT.update_KKT_P(this.kkt, P, this.PtoKKT, this.sigma, this.Pdiag_idx[0], this.Pdiag_n[0]);

	    // Update KKT matrix with new A
	    KKT.update_KKT_A(this.kkt, A, this.AtoKKT);

	    qdldl.factor(this.kkt);

	}


	boolean update_rho_vec(final double[] rho_vec){
	    int i;

	    // Update internal rho_inv_vec
	    for (i = 0; i < this.m; i++){
	        this.rho_inv_vec[i] = 1. / rho_vec[i];
	    }

	    // Update KKT matrix with new rho_vec
	    KKT.update_KKT_param2(this.kkt, this.rho_inv_vec, this.rhotoKKT, this.m);
	    return qdldl.factor(this.kkt);
	}
	


}
