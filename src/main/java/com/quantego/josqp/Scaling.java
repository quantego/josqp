package com.quantego.josqp;

public class Scaling {


	// Set values lower than threshold SCALING_REG to 1
	public static void limit_scaling(double[] D) {
	  int i;

	  for (i = 0; i < D.length; i++) {
	    D[i] = D[i] < OSQP.MIN_SCALING ? 1.0 : D[i];
	    D[i] = D[i] > OSQP.MAX_SCALING ? OSQP.MAX_SCALING : D[i];
	  }
	}
	
	public static double limit_scaling(double D) {
		D = D < OSQP.MIN_SCALING ? 1.0 : D;
		D = D > OSQP.MAX_SCALING ? OSQP.MAX_SCALING : D;
		return D;
	}
	

	/**
	 * Compute infinite norm of the columns of the KKT matrix without forming it
	 *
	 * The norm is stored in the vector v = (D, E)
	 *
	 * @param P        Cost matrix
	 * @param A        Constraints matrix
	 * @param D        Norm of columns related to variables
	 * @param D_temp_A Temporary vector for norm of columns of A
	 * @param E        Norm of columns related to constraints
	 * @param n        Dimension of KKT matrix
	 */
	public static void compute_inf_norm_cols_KKT(CSCMatrix P, CSCMatrix A,
	                               double[] D, double[] D_temp_A,
	                               double[] E, int n) {
	  // First half
	  //  [ P ]
	  //  [ A ]
	  LinAlg.mat_inf_norm_cols_sym_triu(P, D);
	  LinAlg.mat_inf_norm_cols(A, D_temp_A);
	  LinAlg.vec_ew_max_vec(D, D_temp_A, D, n);

	  // Second half
	  //  [ A']
	  //  [ 0 ]
	  LinAlg.mat_inf_norm_rows(A, E);
	}

	public static void scale_data(OSQP.Workspace work) {
	  // Scale KKT matrix
	  //
	  //    [ P   A']
	  //    [ A   0 ]
	  //
	  // with diagonal matrix
	  //
	  //  S = [ D    ]
	  //      [    E ]
	  //

	  int   i;          // Iterations index
	  double c_temp;     // Cost function scaling
	  double inf_norm_q; // Infinity norm of q

	  int n = work.data.n;
	  int m = work.data.m;

	  // Initialize scaling to 1
	  work.scaling.c = 1.0;
	  LinAlg.vec_set_scalar(work.scaling.D,    1., work.data.n);
	  LinAlg.vec_set_scalar(work.scaling.Dinv, 1., work.data.n);
	  LinAlg.vec_set_scalar(work.scaling.E,    1., work.data.m);
	  LinAlg.vec_set_scalar(work.scaling.Einv, 1., work.data.m);


	  for (i = 0; i < work.settings.scaling; i++) {
	    //
	    // First Ruiz step
	    //

	    // Compute norm of KKT columns
	    compute_inf_norm_cols_KKT(work.data.P, work.data.A,
	                              work.D_temp, work.D_temp_A,
	                              work.E_temp, n);

	    // Set to 1 values with 0 norms (avoid crazy scaling)
	    limit_scaling(work.D_temp);
	    limit_scaling(work.E_temp);

	    // Take square root of norms
	    LinAlg.vec_ew_sqrt(work.D_temp, n);
	    LinAlg.vec_ew_sqrt(work.E_temp, m);

	    // Divide scalings D and E by themselves
	    LinAlg.vec_ew_recipr(work.D_temp, work.D_temp, n);
	    LinAlg.vec_ew_recipr(work.E_temp, work.E_temp, m);

	    // Equilibrate matrices P and A and vector q
	    // P <- DPD
	    
	    LinAlg.mat_premult_diag(work.data.P, work.D_temp);
	    LinAlg.mat_postmult_diag(work.data.P, work.D_temp);

	    // A <- EAD
	    LinAlg.mat_premult_diag(work.data.A, work.E_temp);
	    LinAlg.mat_postmult_diag(work.data.A, work.D_temp);

	    // q <- Dq
	    LinAlg.vec_ew_prod(work.D_temp, work.data.q, work.data.q, n);

	    // Update equilibration matrices D and E
	    LinAlg.vec_ew_prod(work.scaling.D, work.D_temp,  work.scaling.D, n);
	    LinAlg.vec_ew_prod(work.scaling.E, work.E_temp,  work.scaling.E, m);

	    //
	    // Cost normalization step
	    //

	    // Compute avg norm of cols of P
	    LinAlg.mat_inf_norm_cols_sym_triu(work.data.P, work.D_temp);
	    c_temp = LinAlg.vec_mean(work.D_temp, n);

	    // Compute inf norm of q
	    inf_norm_q = LinAlg.vec_norm_inf(work.data.q, n);

	    // If norm_q == 0, set it to 1 (ignore it in the scaling)
	    // NB: Using the same function as with vectors here
	    inf_norm_q = limit_scaling(inf_norm_q);

	    // Compute max between avg norm of cols of P and inf norm of q
	    c_temp = Math.max(c_temp, inf_norm_q);

	    // Limit scaling (use same function as with vectors)
	    c_temp = limit_scaling(c_temp);

	    // Invert scaling c = 1 / cost_measure
	    c_temp = 1. / c_temp;

	    // Scale P
	    LinAlg.mat_mult_scalar(work.data.P, c_temp);

	    // Scale q
	    LinAlg.vec_mult_scalar(work.data.q, c_temp, n);

	    // Update cost scaling
	    work.scaling.c *= c_temp;
	  }


	  // Store cinv, Dinv, Einv
	  work.scaling.cinv = 1. / work.scaling.c;
	  LinAlg.vec_ew_recipr(work.scaling.D, work.scaling.Dinv, n);
	  LinAlg.vec_ew_recipr(work.scaling.E, work.scaling.Einv, m);


	  // Scale problem vectors l, u
	  LinAlg. vec_ew_prod(work.scaling.E, work.data.l, work.data.l, m);
	  LinAlg.vec_ew_prod(work.scaling.E, work.data.u, work.data.u, m);

	}


	public static void unscale_data(OSQP.Workspace work) {
	  // Unscale cost
		LinAlg.mat_mult_scalar(work.data.P, work.scaling.cinv);
		LinAlg.mat_premult_diag(work.data.P, work.scaling.Dinv);
		LinAlg.mat_postmult_diag(work.data.P, work.scaling.Dinv);
		LinAlg.vec_mult_scalar(work.data.q, work.scaling.cinv, work.data.n);
		LinAlg.vec_ew_prod(work.scaling.Dinv, work.data.q, work.data.q, work.data.n);

	  // Unscale constraints
		LinAlg.mat_premult_diag(work.data.A, work.scaling.Einv);
		LinAlg.mat_postmult_diag(work.data.A, work.scaling.Dinv);
		LinAlg.vec_ew_prod(work.scaling.Einv, work.data.l, work.data.l, work.data.m);
		LinAlg.vec_ew_prod(work.scaling.Einv, work.data.u, work.data.u, work.data.m);

	}

	public static void unscale_solution(OSQP.Workspace work) {
	  // primal
		LinAlg.vec_ew_prod(work.scaling.D,
	              work.solution.x,
	              work.solution.x, work.data.n);

	  // dual
		LinAlg.vec_ew_prod(work.scaling.E,
	              work.solution.y,
	              work.solution.y, work.data.m);
		LinAlg.vec_mult_scalar(work.solution.y, work.scaling.cinv, work.data.m);

	}


}
