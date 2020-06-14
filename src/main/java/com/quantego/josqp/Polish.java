package com.quantego.josqp;

import java.util.Arrays;

public class Polish {

	/**
	 * Form reduced matrix A that contains only rows that are active at the
	 * solution.
	 * Ared = vstack[Alow, Aupp]
	 * Active constraints are guessed from the primal and dual solution returned by
	 * the ADMM.
	 * @param  work Workspace
	 * @return      Number of rows in Ared, negative if error
	 */
	public static int form_Ared(OSQP.Workspace work) {
	  int j, ptr;
	  int Ared_nnz = 0;

	  // Initialize counters for active constraints
	  work.pol.n_low = 0;
	  work.pol.n_upp = 0;

	  /* Guess which linear constraints are lower-active, upper-active and free
	   *    A_to_Alow[j] = -1    (if j-th row of A is not inserted in Alow)
	   *    A_to_Alow[j] =  i    (if j-th row of A is inserted at i-th row of Alow)
	   * Aupp is formed in the equivalent way.
	   * Ared is formed by stacking vertically Alow and Aupp.
	   */
	  for (j = 0; j < work.data.m; j++) {
	    if (work.z[j] - work.data.l[j] < -work.y[j]) { // lower-active
	      work.pol.Alow_to_A[work.pol.n_low] = j;
	      work.pol.A_to_Alow[j]                = work.pol.n_low++;
	    } else {
	      work.pol.A_to_Alow[j] = -1;
	    }
	  }

	  for (j = 0; j < work.data.m; j++) {
	    if (work.data.u[j] - work.z[j] < work.y[j]) { // upper-active
	      work.pol.Aupp_to_A[work.pol.n_upp] = j;
	      work.pol.A_to_Aupp[j]                = work.pol.n_upp++;
	    } else {
	      work.pol.A_to_Aupp[j] = -1;
	    }
	  }

	  // Check if there are no active constraints
	  if (work.pol.n_low + work.pol.n_upp == 0) {
	    // Form empty Ared
	    work.pol.Ared = new CSCMatrix(0,work.data.n,0,true,false);
	    LinAlg.int_vec_set_scalar(work.pol.Ared.Ap, 0, work.data.n + 1);
	    return 0; // mred = 0
	  }

	  // Count number of elements in Ared
	  for (j = 0; j < work.data.A.Ap[work.data.A.n]; j++) {
	    if ((work.pol.A_to_Alow[work.data.A.Ai[j]] != -1) ||
	        (work.pol.A_to_Aupp[work.data.A.Ai[j]] != -1)) Ared_nnz++;
	  }

	  // Form Ared
	  // Ared = vstack[Alow, Aupp]
	  work.pol.Ared = new CSCMatrix(work.pol.n_low + work.pol.n_upp,
	                                work.data.n, Ared_nnz, true, false);
	  Ared_nnz = 0; // counter

	  for (j = 0; j < work.data.n; j++) { // Cycle over columns of A
	    work.pol.Ared.Ap[j] = Ared_nnz;

	    for (ptr = work.data.A.Ap[j]; ptr < work.data.A.Ap[j + 1]; ptr++) {
	      // Cycle over elements in j-th column
	      if (work.pol.A_to_Alow[work.data.A.Ai[ptr]] != -1) {
	        // Lower-active rows of A
	        work.pol.Ared.Ai[Ared_nnz] =
	          work.pol.A_to_Alow[work.data.A.Ai[ptr]];
	        work.pol.Ared.Ax[Ared_nnz++] = work.data.A.Ax[ptr];
	      } else if (work.pol.A_to_Aupp[work.data.A.Ai[ptr]] != -1) {
	        // Upper-active rows of A
	        work.pol.Ared.Ai[Ared_nnz] = work.pol.A_to_Aupp[work.data.A.Ai[ptr]] 
	                                       + work.pol.n_low;
	        work.pol.Ared.Ax[Ared_nnz++] = work.data.A.Ax[ptr];
	      }
	    }
	  }

	  // Update the last element in Ared.p
	  work.pol.Ared.Ap[work.data.n] = Ared_nnz;

	  // Return number of rows in Ared
	  return work.pol.n_low + work.pol.n_upp;
	}

	/**
	 * Form reduced right-hand side rhs_red = vstack[-q, l_low, u_upp]
	 * @param  work Workspace
	 * @param  rhs  right-hand-side
	 * @return      reduced rhs
	 */
	public static void form_rhs_red(OSQP.Workspace work, double[]rhs) {
	  int j;

	  // Form the rhs of the reduced KKT linear system
	  for (j = 0; j < work.data.n; j++) { // -q
	    rhs[j] = -work.data.q[j];
	  }

	  for (j = 0; j < work.pol.n_low; j++) { // l_low
	    rhs[work.data.n + j] = work.data.l[work.pol.Alow_to_A[j]];
	  }

	  for (j = 0; j < work.pol.n_upp; j++) { // u_upp
	    rhs[work.data.n + work.pol.n_low + j] =
	      work.data.u[work.pol.Aupp_to_A[j]];
	  }
	}

	/**
	 * Perform iterative refinement on the polished solution:
	 *    (repeat)
	 *    1. (K + dK) * dz = b - K*z
	 *    2. z <- z + dz
	 * @param  work Solver workspace
	 * @param  p    Private variable for solving linear system
	 * @param  z    Initial z value
	 * @param  b    RHS of the linear system
	 */
	public static int iterative_refinement(OSQP.Workspace work,
	                                  LinSys  p,
	                                  double[]      z,
	                                  double[]      b) {
	  int i, j, n;
	  double[] rhs;

	  if (work.settings.polish_refine_iter > 0) {

		  // Assign dimension n
		  n = work.data.n + work.pol.Ared.m;

		  // Allocate rhs vector
	      rhs = new double[n];

	      for (i = 0; i < work.settings.polish_refine_iter; i++) {
	        // Form the RHS for the iterative refinement:  b - K*z
	        LinAlg.prea_vec_copy(b, rhs, n);
	
	        // Upper Part: R^{n}
	        // -= Px (upper triang)
	        LinAlg.mat_vec(work.data.P, z, rhs, 0, 0, -1);
	
	        // -= Px (lower triang)
	        LinAlg.mat_tpose_vec(work.data.P, z, rhs, 0, 0, -1, true);
	
	        // -= Ared'*y_red
	        LinAlg.mat_tpose_vec(work.pol.Ared, z, rhs, work.data.n, 0, -1, false);
	
	        // Lower Part: R^{m}
	        LinAlg.mat_vec(work.pol.Ared, z, rhs, 0, work.data.n, -1);
	
	        // Solve linear system. Store solution in rhs
	        p.solve(rhs);
	
	        // Update solution
	        for (j = 0; j < n; j++) {
	          z[j] += rhs[j];
	        }
	      }
	     
		}
	  return 0;
	}

	/**
	 * Compute dual variable y from yred
	 * @param work Workspace
	 * @param yred Dual variables associated to active constraints
	 */
	public static void get_ypol_from_yred(OSQP.Workspace work, double[] yred, int start) {
	  int j;

	  // If there are no active constraints
	  if (work.pol.n_low + work.pol.n_upp == 0) {
		  LinAlg.vec_set_scalar(work.pol.y, 0., work.data.m);
	    return;
	  }

	  // NB: yred = vstack[ylow, yupp]
	  for (j = 0; j < work.data.m; j++) {
	    if (work.pol.A_to_Alow[j] != -1) {
	      // lower-active
	      work.pol.y[j] = yred[work.pol.A_to_Alow[j]+start];
	    } else if (work.pol.A_to_Aupp[j] != -1) {
	      // upper-active
	      work.pol.y[j] = yred[work.pol.A_to_Aupp[j] + work.pol.n_low + start];
	    } else {
	      // inactive
	      work.pol.y[j] = 0.0;
	    }
	  }
	}

	public static int polish(OSQP.Workspace work) {
	  int mred;
	  boolean polish_successful;
	  LinSys plsh;
	  double[] pol_sol; // Polished solution


	  // Form Ared by assuming the active constraints and store in work.pol.Ared
	  mred = form_Ared(work);
	  if (mred < 0) { // work.pol.red = OSQP_NULL
	    // Polishing failed
	    work.info.status_polish = -1;

	    return -1;
	  }

	  // Form and factorize reduced KKT
	  plsh = new LinSys(work.data.P, work.pol.Ared, work.settings.delta, null, true);



	  // Form reduced right-hand side rhs_red
	  double[] rhs_red = new double[work.data.n*mred];

	  form_rhs_red(work, rhs_red);

	  pol_sol = Arrays.copyOf(rhs_red, work.data.n + mred);
	  

	  // Solve the reduced KKT system
	  plsh.solve(pol_sol);

	  // Perform iterative refinement to compensate for the regularization error
	  iterative_refinement(work, plsh, pol_sol, rhs_red);



	  // Store the polished solution (x,z,y)
	  LinAlg.prea_vec_copy(pol_sol, work.pol.x, work.data.n);   // pol.x
	  LinAlg.mat_vec(work.data.A, work.pol.x, work.pol.z, 0, 0, 0); // pol.z
	  get_ypol_from_yred(work, pol_sol, work.data.n);     // pol.y

	  // Ensure (z,y) satisfies normal cone constraint
	  Projection.project_normalcone(work, work.pol.z, work.pol.y);

	  // Compute primal and dual residuals at the polished solution
	  OSQP.update_info(work, 0, true, true);

	  // Check if polish was successful
	  polish_successful = (work.pol.pri_res < work.info.pri_res &&
	                       work.pol.dua_res < work.info.dua_res) || // Residuals
	                                                                    // are
	                                                                    // reduced
	                      (work.pol.pri_res < work.info.pri_res &&
	                       work.info.dua_res < 1e-10) ||              // Dual
	                                                                    // residual
	                                                                    // already
	                                                                    // tiny
	                      (work.pol.dua_res < work.info.dua_res &&
	                       work.info.pri_res < 1e-10);                // Primal
	                                                                    // residual
	                                                                    // already
	                                                                    // tiny

	  if (polish_successful) {
	    // Update solver information
	    work.info.obj_val       = work.pol.obj_val;
	    work.info.pri_res       = work.pol.pri_res;
	    work.info.dua_res       = work.pol.dua_res;
	    work.info.status_polish = 1;

	    // Update (x, z, y) in ADMM iterations
	    // NB: z needed for warm starting
	    LinAlg.prea_vec_copy(work.pol.x, work.x, work.data.n);
	    LinAlg.prea_vec_copy(work.pol.z, work.z, work.data.m);
	    LinAlg.prea_vec_copy(work.pol.y, work.y, work.data.m);

	  } else { // Polishing failed
	    work.info.status_polish = -1;

	    // TODO: Try to find a better solution on the line connecting ADMM
	    //       and polished solution
	  }

	 

	  return 0;
	}

}
