package com.quantego.josqp;

public class OSQP {
	
	static final double RHO = 0.1;
	static final double SIGMA = 1E-06;
	static final int MAX_ITER = 4000;
	static final double EPS_ABS = 1E-3;
	static final double EPS_REL = 1E-3;
	static final double EPS_PRIM_INF = 1E-4;
	static final double EPS_DUAL_INF = 1E-4;
	static final double ALPHA = 1.6;
	static final LinSys.TYPE LINSYS_SOLVER = LinSys.TYPE.QLDL;

	static final double RHO_MIN = 1e-06;
	static final double RHO_MAX = 1e06;
	static final double HO_EQ_OVER_RHO_INEQ = 1e03;
	static final double RHO_TOL = 1e-04; ///< tolerance for detecting if an inequality is set to equality


	static final double DELTA = 1.0e-6;
	static final boolean POLISH = false;
	static final int POLISH_REFINE_ITER = 3;
	static final boolean VERBOSE = true;

	static final int SCALED_TERMINATION = 0;
	static final int CHECK_TERMINATION = 25;
	static final int WARM_START = 1;
	static final int SCALING = 10;

	static final double MIN_SCALING = 1.0e-04; ///< minimum scaling value
	static final double MAX_SCALING = 1.0e+04; ///< maximum scaling value


	static final double OSQP_NULL = 0.0;
	static final double OSQP_NAN = Double.NaN;
	static final double OSQP_INFTY = Double.POSITIVE_INFINITY;

	static final boolean ADAPTIVE_RHO = true;
	static final int ADAPTIVE_RHO_INTERVAL = 0;
	static final double ADAPTIVE_RHO_FRACTION= 0.4;         ///< fraction of setup time after which we update rho
	static final int ADAPTIVE_RHO_MULTIPLE_TERMINATION = 4; ///< multiple of check_termination after which we update rho (if PROFILING disabled)
	static final int ADAPTIVE_RHO_FIXED = 100;             ///< number of iterations after which we update rho if termination_check  and PROFILING are disabled
	static final int ADAPTIVE_RHO_TOLERANCE = 5;          ///< tolerance for adopting new rho; minimum ratio between new rho and the current one

	static final int TIME_LIMIT = 0;  
	
	public static class Scaling {
		double  c;    ///< cost function scaling
		double[] D;    ///< primal variable scaling
		double[] E;    ///< dual variable scaling
		double  cinv; ///< cost function rescaling
		double[] Dinv; ///< primal variable rescaling
		double[] Einv; ///< dual variable rescaling
	}
	
	public enum Error {
		OSQP_DATA_VALIDATION_ERROR,  /* Start errors from 1 */
	    OSQP_SETTINGS_VALIDATION_ERROR,
	    OSQP_LINSYS_SOLVER_LOAD_ERROR,
	    OSQP_LINSYS_SOLVER_INIT_ERROR,
	    OSQP_NONCVX_ERROR,
	    OSQP_MEM_ALLOC_ERROR,
	    OSQP_WORKSPACE_NOT_INIT_ERROR
	}
	
	public enum Status {
		DUAL_INFEASIBLE_INACCURATE,
		PRIMAL_INFEASIBLE_INACCURATE, 
		SOLVED_INACCURATE,
		SOLVED,
		MAX_ITER_REACHED,
		PRIMAL_INFEASIBLE,   /* primal infeasible  */
		DUAL_INFEASIBLE,  /* dual infeasible */
		SIGINT,             /* interrupted by user */
		TIME_LIMIT_REACHED,
		NON_CVX,            /* problem non convex */
		UNSOLVED  
	}
	
	public static class Info {
		
		public int iter;          ///< number of iterations taken
		public Status status;     ///< status string, e.g. 'solved'
		public int status_val;    ///< status as c_int, defined in constants.h

		public int status_polish; ///< polish status: successful (1), unperformed (0), (-1) unsuccessful

		public double obj_val;     ///< primal objective
		public double pri_res;     ///< norm of primal residual
		public double dua_res;    ///< norm of dual residual

		public double setup_time;  ///< time taken for setup phase (seconds)
		public double solve_time;  ///< time taken for solve phase (seconds)
		public double update_time; ///< time taken for update phase (seconds)
		public double polish_time; ///< time taken for polish phase (seconds)
		public double run_time;    ///< total time  (seconds)

		public int rho_updates;  ///< number of rho updates
		public double rho_estimate; ///< best rho estimate so far from residuals
	}
	
	public static class Polish {
		CSC Ared;          ///< active rows of A
		///<    Ared = vstack[Alow, Aupp]
		int n_low;     ///< number of lower-active rows
		int n_upp;     ///< number of upper-active rows
		int[] A_to_Alow; ///< Maps indices in A to indices in Alow
		int[] A_to_Aupp; ///< Maps indices in A to indices in Aupp
		int[] Alow_to_A; ///< Maps indices in Alow to indices in A
		int[] Aupp_to_A; ///< Maps indices in Aupp to indices in A
		double[] x;         ///< optimal x-solution obtained by polish
		double[] z;         ///< optimal z-solution obtained by polish
		double[] y;         ///< optimal y-solution obtained by polish
		double obj_val;   ///< objective value at polished solution
		double pri_res;   ///< primal residual at polished solution
		double dua_res;   ///< dual residual at polished solution
	}
	
	public static class Settings {
		double rho;                    ///< ADMM step rho
		double sigma;                  ///< ADMM step sigma
		int scaling;                ///< heuristic data scaling iterations; if 0, then disabled.

		boolean adaptive_rho;           ///< boolean, is rho step size adaptive?
		int   adaptive_rho_interval;  ///< number of iterations between rho adaptations; if 0, then it is automatic
		double adaptive_rho_tolerance; ///< tolerance X for adapting rho. The new rho has to be X times larger or 1/X times smaller than the current one to trigger a new factorization.
		double adaptive_rho_fraction;  ///< interval for adapting rho (fraction of the setup time)

		int max_iter;      ///< maximum number of iterations
		double eps_abs;       ///< absolute convergence tolerance
		double eps_rel;       ///< relative convergence tolerance
		double eps_prim_inf;  ///< primal infeasibility tolerance
		double eps_dual_inf;  ///< dual infeasibility tolerance
		double  alpha;         ///< relaxation parameter
		LinSys linsys_solver; ///< linear system solver to use

		double delta;                         ///< regularization parameter for polishing
		boolean polish;                        ///< boolean, polish ADMM solution
		int polish_refine_iter;            ///< number of iterative refinement steps in polishing

		boolean verbose;                         ///< boolean, write out progress

		boolean scaled_termination;              ///< boolean, use scaled termination criteria
		int check_termination;               ///< integer, check termination interval; if 0, then termination checking is disabled
		boolean warm_start;                      ///< boolean, warm start

		double time_limit;                    ///< maximum number of seconds allowed to solve the problem; if 0, then disabled
	}
	
	public static class Data {
		int n; ///< number of variables n
		int m; ///< number of constraints m
		CSC     P; ///< the upper triangular part of the quadratic cost matrix P in csc format (size n x n).
		CSC     A; ///< linear constraints matrix A in csc format (size m x n)
		double[] q; ///< dense array for linear part of cost function (size n)
		double[] l; ///< dense array for lower bound (size m)
		double[] u; ///< dense array for upper bound (size m)
	}
	
	public static class Solution {
		double[] x; ///< primal solution
		double[] y; ///< Lagrange multiplier associated to \f$l <= Ax <= u\f$
	}
	
	public static class Workspace {
		/// Problem data to work on (possibly scaled)
		Data data;

		  /// Linear System solver structure
		LinSys linsys_solver;

		  /// Polish structure
		Polish pol;

		  /**
		   * @name Vector used to store a vectorized rho parameter
		   * @{
		   */
		double[] rho_vec;     ///< vector of rho values
		double[] rho_inv_vec; ///< vector of inv rho values

		  /** @} */

		int[] constr_type; ///< Type of constraints: loose (-1), equality (1), inequality (0)

		  /**
		   * @name Iterates
		   * @{
		   */
		double[] x;        ///< Iterate x
		double[] y;        ///< Iterate y
		double[] z;        ///< Iterate z
		double[] xz_tilde; ///< Iterate xz_tilde

		double[] x_prev;   ///< Previous x

	    /**< NB: Used also as workspace vector for dual residual */
		double[] z_prev;   ///< Previous z

	    /**< NB: Used also as workspace vector for primal residual */

	    /**
	     * @name Primal and dual residuals workspace variables
	     *
	     * Needed for residuals computation, tolerances computation,
	     * approximate tolerances computation and adapting rho
	     * @{
	     */
		double[] Ax;  ///< scaled A * x
		double[] Px;  ///< scaled P * x
		double[] Aty; ///< scaled A * x

		/**
		 * @name Primal infeasibility variables
		 * @{
		 */
		double[] delta_y;   ///< difference between consecutive dual iterates
		double[] Atdelta_y; ///< A' * delta_y

		/**
		 * @name Dual infeasibility variables
		 * @{
		 */
		double[] delta_x;  ///< difference between consecutive primal iterates
		double[] Pdelta_x; ///< P * delta_x
		double[] Adelta_x; ///< A * delta_x

		/**
		 * @name Temporary vectors used in scaling
		 * @{
		 */

		double[] D_temp;   ///< temporary primal variable scaling vectors
		double[] D_temp_A; ///< temporary primal variable scaling vectors storing norms of A columns
		double[] E_temp;   ///< temporary constraints scaling vectors storing norms of A' columns


		/** @} */

		Settings settings; ///< problem settings
		Scaling  scaling;  ///< scaling vectors
		Solution solution; ///< problem solution
		Info     info;     ///< solver information


		  /// flag indicating whether the solve function has been run before
		boolean first_run;

		  /// flag indicating whether the update_time should be cleared
		int clear_update_time;

		  /// flag indicating that osqp_update_rho is called from osqp_solve function
		int rho_update_from_solve;

	}
	
	Data data;
	Settings settings;
	Workspace work;
	
	public OSQP(Data data, Settings settings) {
		//TODO
	}
	
	public void setWorkspace(Workspace work) {
		this.work = work;
	}
	
	public Workspace getWorkspace() {
		return work;
	}
	
	public int solve() {
		int exitflag = 0;
		int iter = 0;
		boolean compute_cost_function = work.settings.verbose; // Boolean: compute the cost function in the loop or not
		boolean can_check_termination; // Boolean: check termination or not
		double temp_run_time;       // Temporary variable to store current run time
		if (work.clear_update_time == 1)
		    work.info.update_time = 0.0;
		work.rho_update_from_solve = 1;

		//TODO: timer


		// Initialize variables (cold start or warm start depending on settings)
		if (!work.settings.warm_start) Auxil.cold_start(work);  // If not warm start ->
		                                                      // set x, z, y to zero

		// Main ADMM algorithm
		for (iter = 1; iter <= work.settings.max_iter; iter++) {
		    // Update x_prev, z_prev (preallocated, no malloc)
		    //swap_vectors(&(work.x), &(work.x_prev));
		    //swap_vectors(&(work.z), &(work.z_prev));

		    /* ADMM STEPS */
		    /* Compute \tilde{x}^{k+1}, \tilde{z}^{k+1} */
		    Auxil.update_xz_tilde(work);

		    /* Compute x^{k+1} */
		    Auxil.update_x(work);

		    /* Compute z^{k+1} */
		    Auxil.update_z(work);

		    /* Compute y^{k+1} */
		    Auxil.update_y(work);
		    
		    can_check_termination = work.settings.check_termination>0 &&
                    (iter % work.settings.check_termination == 0);
		    
		    if (can_check_termination) {
		        // Update information and compute also objective value
		        update_info(work, iter, compute_cost_function, 0);

		        // Check algorithm termination
		        if (check_termination(work, 0)) {
		          // Terminate algorithm
		          break;
		        }
		      }
		    // Adapt rho
		    if (work.settings.adaptive_rho &&
		        work.settings.adaptive_rho_interval>0 &&
		        (iter % work.settings.adaptive_rho_interval == 0))  {
		      // Update info with the residuals if it hasn't been done before
			    if (!can_check_termination) {
			        // Information has not been computed before for termination check
			        update_info(work, iter, compute_cost_function, 0);
			      }
	
			      // Actually update rho
			      if (adapt_rho(work)) 
			        throw new IllegalStateException("Failed rho update");
		    }


		  }        // End of ADMM for loop


		  // Update information and check termination condition if it hasn't been done
		  // during last iteration (max_iter reached or check_termination disabled)
		  if (!can_check_termination) {
		    /* Update information */


		    // If no printing is enabled, update info directly
		    update_info(work, iter - 1, compute_cost_function, 0);



		    /* Check whether a termination criterion is triggered */
		    check_termination(work, 0);
		  }

		  // Compute objective value in case it was not
		  // computed during the iterations
		  if (!compute_cost_function && has_solution(work.info)){
		    work.info.obj_val = compute_obj_val(work, work.x);
		  }



		  /* if max iterations reached, change status accordingly */
		  if (work.info.status_val == Status.UNSOLVED) {
		    if (!check_termination(work, 1)) { // Try to check for approximate
		      work.info = Status.MAX_ITER_REACHED;
		    }
		  }



		  /* Update rho estimate */
		  work.info.rho_estimate = compute_rho_estimate(work);




		  // Polish the obtained solution
		  if (work.settings.polish && (work.info.status_val == OSQP_SOLVED))
		    polish(work);





		  // Store solution
		  store_solution(work);




		  return exitflag;
	}
	
	int check_termination(boolean approximate) {
		  double eps_prim, eps_dual, eps_prim_inf, eps_dual_inf;
		  int   exitflag;
		  boolean   prim_res_check, dual_res_check, prim_inf_check, dual_inf_check;
		  double eps_abs, eps_rel;

		  // Initialize variables to 0
		  exitflag       = 0;
		  prim_res_check = false; dual_res_check = false;
		  prim_inf_check = false; dual_inf_check = false;

		  // Initialize tolerances
		  eps_abs      = work.settings.eps_abs;
		  eps_rel      = work.settings.eps_rel;
		  eps_prim_inf = work.settings.eps_prim_inf;
		  eps_dual_inf = work.settings.eps_dual_inf;

		  // If residuals are too large, the problem is probably non convex
		  if ((work.info.pri_res > OSQP_INFTY) ||
		      (work.info.dua_res > OSQP_INFTY)){
		    // Looks like residuals are diverging. Probably the problem is non convex!
		    // Terminate and report it
		    work.info.status = Status.NON_CVX;
		    work.info.obj_val = OSQP_NAN;
		    return 1;
		  }

		  // If approximate solution required, increase tolerances by 10
		  if (approximate) {
		    eps_abs      *= 10;
		    eps_rel      *= 10;
		    eps_prim_inf *= 10;
		    eps_dual_inf *= 10;
		  }

		  // Check residuals
		  if (work.data.m == 0) {
		    prim_res_check = true; // No constraints -> Primal feasibility always satisfied
		  }
		  else {
		    // Compute primal tolerance
		    eps_prim = compute_pri_tol(work, eps_abs, eps_rel);

		    // Primal feasibility check
		    if (work.info.pri_res < eps_prim) {
		      prim_res_check = true;
		    } else {
		      // Primal infeasibility check
		      prim_inf_check = is_primal_infeasible(work, eps_prim_inf);
		    }
		  } // End check if m == 0

		  // Compute dual tolerance
		  eps_dual = compute_dua_tol(work, eps_abs, eps_rel);

		  // Dual feasibility check
		  if (work.info.dua_res < eps_dual) {
		    dual_res_check = true;
		  } else {
		    // Check dual infeasibility
		    dual_inf_check = is_dual_infeasible(work, eps_dual_inf);
		  }

		  // Compare checks to determine solver status
		  if (prim_res_check && dual_res_check) {
		    // Update final information
		    if (approximate) {
		      work.info = Status.SOLVED_INACCURATE;
		    } else {
		      work.info = Status.SOLVED;
		    }
		    exitflag = 1;
		  }
		  else if (prim_inf_check) {
		    // Update final information
		    if (approximate) {
		      work.info = Status.PRIMAL_INFEASIBLE_INACCURATE;
		    } else {
		      work.info = Status.PRIMAL_INFEASIBLE;
		    }

		    if (work.settings.scaling!=0 && !work.settings.scaled_termination) {
		      // Update infeasibility certificate
		      vec_ew_prod(work.scaling.E, work.delta_y, work.delta_y, work.data.m);
		    }
		    work.info.obj_val = OSQP_INFTY;
		    exitflag            = 1;
		  }
		  else if (dual_inf_check) {
		    // Update final information
		    if (approximate) {
		    	work.info = Status.DUAL_INFEASIBLE_INACCURATE;
		    } else {
		    	work.info = Status.DUAL_INFEASIBLE;
		    }

		    if (work.settings.scaling!=0 && !work.settings.scaled_termination) {
		      // Update infeasibility certificate
		      vec_ew_prod(work.scaling.D, work.delta_x, work.delta_x, work.data.n);
		    }
		    work.info.obj_val = -OSQP_INFTY;
		    exitflag            = 1;
		  }

		  return exitflag;
		}

}
