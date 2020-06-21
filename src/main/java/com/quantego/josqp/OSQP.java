package com.quantego.josqp;

import java.util.Arrays;

public class OSQP {
	
	static final double RHO = 0.1;
	static final double SIGMA = 1E-06;
	static final int MAX_ITER = 4000;
	static final double EPS_ABS = 1E-3;
	static final double EPS_REL = 1E-3;
	static final double EPS_PRIM_INF = 1E-4;
	static final double EPS_DUAL_INF = 1E-4;
	static final double ALPHA = 1.6;

	static final double RHO_MIN = 1e-06;
	static final double RHO_MAX = 1e06;
	static final double RHO_EQ_OVER_RHO_INEQ = 1e03;
	static final double RHO_TOL = 1e-04; ///< tolerance for detecting if an inequality is set to equality


	static final double DELTA = 1.0e-6;
	static final boolean POLISH = false;
	static final int POLISH_REFINE_ITER = 3;
	static final boolean VERBOSE = true;

	static final boolean SCALED_TERMINATION = false;
	static final int CHECK_TERMINATION = 25;
	static final boolean WARM_START = true;
	static final int SCALING = 10;

	static final double MIN_SCALING = 1.0e-04; ///< minimum scaling value
	static final double MAX_SCALING = 1.0e+04; ///< maximum scaling value


	static final double OSQP_NULL = 0.0;
	static final double OSQP_NAN = Double.NaN;
	static final double OSQP_INFTY = 1.0e30;

	static final boolean ADAPTIVE_RHO = true;
	static final int ADAPTIVE_RHO_INTERVAL = 0;
	static final double ADAPTIVE_RHO_FRACTION = 0.4;         ///< fraction of setup time after which we update rho
	static final int ADAPTIVE_RHO_MULTIPLE_TERMINATION = 4; ///< multiple of check_termination after which we update rho (if PROFILING disabled)
	static final int ADAPTIVE_RHO_FIXED = 100;             ///< number of iterations after which we update rho if termination_check  and PROFILING are disabled
	static final int ADAPTIVE_RHO_TOLERANCE = 5;          ///< tolerance for adopting new rho; minimum ratio between new rho and the current one

	static final int TIME_LIMIT = 0;  
	
	public static class ScaledProblemData {
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
		TIME_LIMIT_REACHED,
		NON_CVX,            /* problem non convex */
		UNSOLVED,
		ERROR
	}
	
	public static class Info {
		
		public int iter;          ///< number of iterations taken
		public Status status;     ///< status string, e.g. 'solved'
//		public Status status_val;    ///< status as c_int, defined in constants.h

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
	
	public static class PolishedData {
		CSCMatrix Ared;          ///< active rows of A
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
		public double rho = RHO;                    ///< ADMM step rho
		public double sigma = SIGMA;                  ///< ADMM step sigma
		public int scaling = SCALING;                ///< heuristic data scaling iterations; if 0, then disabled.

		public boolean adaptive_rho = ADAPTIVE_RHO;           ///< boolean, is rho step size adaptive?
		public int   adaptive_rho_interval = ADAPTIVE_RHO_INTERVAL;  ///< number of iterations between rho adaptations; if 0, then it is automatic
		public double adaptive_rho_tolerance = ADAPTIVE_RHO_TOLERANCE; ///< tolerance X for adapting rho. The new rho has to be X times larger or 1/X times smaller than the current one to trigger a new factorization.
		public double adaptive_rho_fraction = ADAPTIVE_RHO_FRACTION;  ///< interval for adapting rho (fraction of the setup time)

		public int max_iter = MAX_ITER;      ///< maximum number of iterations
		public double eps_abs = EPS_ABS;       ///< absolute convergence tolerance
		double eps_rel = EPS_REL;       ///< relative convergence tolerance
		public double eps_prim_inf = EPS_PRIM_INF;  ///< primal infeasibility tolerance
		public double eps_dual_inf = EPS_DUAL_INF;  ///< dual infeasibility tolerance
		public double  alpha = ALPHA;         ///< relaxation parameter
//		LinSys linsys_solver; ///< linear system solver to use

		public double delta = DELTA;                         ///< regularization parameter for polishing
		public boolean polish = POLISH;                        ///< boolean, polish ADMM solution
		public int polish_refine_iter = POLISH_REFINE_ITER;            ///< number of iterative refinement steps in polishing

		public boolean verbose = VERBOSE;                         ///< boolean, write out progress

		public boolean scaled_termination = SCALED_TERMINATION;              ///< boolean, use scaled termination criteria
		public int check_termination = CHECK_TERMINATION;               ///< integer, check termination interval; if 0, then termination checking is disabled
		public boolean warm_start = WARM_START;                      ///< boolean, warm start

		public double time_limit = TIME_LIMIT;                    ///< maximum number of seconds allowed to solve the problem; if 0, then disabled
	}
	
	public static class Data {
		//TODO: make this private again once all tests pass
		public final int n; ///< number of variables n
		public final int m; ///< number of constraints m
		public final CSCMatrix     P; ///< the upper triangular part of the quadratic cost matrix P in csc format (size n x n).
		public final CSCMatrix     A; ///< linear constraints matrix A in csc format (size m x n)
		public final double[] q; ///< dense array for linear part of cost function (size n)
		public final double[] l; ///< dense array for lower bound (size m)
		public final double[] u; ///< dense array for upper bound (size m)
		public int getN() {
			return n;
		}
		public int getM() {
			return m;
		}
		public Data(int n, int m, CSCMatrix P, CSCMatrix A, double[] q, double[] l, double[] u) {
			this.n = n;
			this.m = m;
			this.P = CSCMatrix.copy_csc_mat(P);
			this.A = CSCMatrix.copy_csc_mat(A);
			this.q = q.clone();
			this.l = l.clone();
			this.u = u.clone();
		}
	}
	
	public static class Solution {
		double[] x; ///< primal solution
		double[] y; ///< Lagrange multiplier associated to \f$l <= Ax <= u\f$
	}
	
	
	public static class Workspace {
		/// Problem data to work on (possibly scaled)
		public Data data;

		  /// Linear System solver structure
		LinSys linsys_solver;

		  /// Polish structure
		PolishedData pol;

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

		public Settings settings; ///< problem settings
		ScaledProblemData  scaling;  ///< scaling vectors
		Solution solution; ///< problem solution
		public Info     info;     ///< solver information


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
		this.work = new Workspace();
		work.data = data;
		//TODO: C code creates a deep copy of the data, needed?
		work.rho_vec = new double[data.m];
		work.rho_inv_vec = new double[data.m];
		work.constr_type = new int[data.m]; //maybe use byte->char here or an enum
		work.x = new double[data.n];
		work.z = new double[data.m];
		work.xz_tilde = new double[data.n+data.m];
		work.x_prev = new double[data.n];
		work.z_prev = new double[data.m];
		work.y = new double[data.m];
		//cold_start(work); not needed, since vars are zero by default
		// Primal and dual residuals variables
		work.Ax = new double[data.m];
		work.Px = new double[data.n];
		work.Aty = new double[data.n];
		 // Primal infeasibility variables
		work.delta_y = new double[data.m];
		work.Atdelta_y = new double[data.n];
		// Dual infeasibility variables
		work.delta_x = new double[data.n];
		work.Pdelta_x = new double[data.n];
		work.Adelta_x = new double[data.m];
		// Settings
		work.settings = settings;
		// Perform scaling
		if(settings.scaling!=0) {
			// Allocate scaling structure
		    work.scaling = new OSQP.ScaledProblemData();
		    work.scaling.D    = new double[data.n];
		    work.scaling.Dinv = new double[data.n];
		    work.scaling.E    = new double[data.m];
		    work.scaling.Einv = new double[data.m];

		    // Allocate workspace variables used in scaling
		    work.D_temp   = new double[data.n];
		    work.D_temp_A = new double[data.n];
		    work.E_temp   = new double[data.m];

		    // Scale data
		    
		    Scaling.scale_data(work);
		  }
		// Set type of constraints
		set_rho_vec(work);
		// Load linear system solver
		work.linsys_solver = new LinSys(
				work.data.P, work.data.A,
                work.settings.sigma, work.rho_vec, false
                );
		
		// Initialize active constraints structure
		work.pol = new OSQP.PolishedData();
		work.pol.Aupp_to_A = new int[data.m];
		work.pol.Alow_to_A = new int[data.m];
		work.pol.A_to_Alow = new int[data.m];
		work.pol.A_to_Aupp = new int[data.m];
		work.pol.x = new double[data.n];
		work.pol.z = new double[data.m];
		work.pol.y = new double[data.m];
		
		// Allocate solution
		work.solution = new OSQP.Solution();
		work.solution.x = new double[data.n];
		work.solution.y = new double[data.m];
		
		work.info = new OSQP.Info();
		work.info.status = Status.UNSOLVED;
		
		//for profiling
		work.first_run = true;
		work.info.setup_time = System.currentTimeMillis();
		
		work.info.rho_estimate = work.settings.rho;
		
		work.settings.adaptive_rho_interval = Math.max(
		          work.settings.adaptive_rho_interval,
		          work.settings.check_termination);
		
	}
	
	public void setWorkspace(Workspace work) {
		this.work = work;
	}
	
	public Workspace getWorkspace() {
		return work;
	}
	
	public OSQP.Status solve() {
		int iter = 0;
		boolean compute_cost_function = work.settings.verbose; // Boolean: compute the cost function in the loop or not
		boolean can_check_termination = false; // Boolean: check termination or not
		if (work.clear_update_time == 1)
		    work.info.update_time = 0.0;
		work.rho_update_from_solve = 1;



		// Initialize variables (cold start or warm start depending on settings)
		if (!work.settings.warm_start) 
			cold_start(work);  // If not warm start ->
		                                                      // set x, z, y to zero

		// Main ADMM algorithm
		for (iter = 1; iter <= work.settings.max_iter; iter++) {
		    // Update x_prev, z_prev (preallocated, no malloc)
			
			double[] temp;
			temp = work.x_prev;
			work.x_prev = work.x;
			work.x = temp;
			
			temp = work.z_prev;
			work.z_prev = work.z;
			work.z = temp;

		    /* ADMM STEPS */
		    /* Compute \tilde{x}^{k+1}, \tilde{z}^{k+1} */
		    update_xz_tilde(work);
		    /* Compute x^{k+1} */
		    update_x(work);
		    /* Compute z^{k+1} */
		    update_z(work);
		    /* Compute y^{k+1} */
		    update_y(work);
		    can_check_termination = work.settings.check_termination>0 &&
                    (iter % work.settings.check_termination == 0);
		    if (can_check_termination) {
		        // Update information and compute also objective value
		        update_info(work, iter, compute_cost_function, false);
		        // Check algorithm termination
		        if (check_termination(work, false)) {
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
			        update_info(work, iter, compute_cost_function, false);
			      }
			      // Actually update rho
			      if (adapt_rho(work))
			    	  return Status.ERROR;
		    }
		  }        // End of ADMM for loop


		  // Update information and check termination condition if it hasn't been done
		  // during last iteration (max_iter reached or check_termination disabled)
		  if (!can_check_termination) {
		    /* Update information */
		    // If no printing is enabled, update info directly
		    update_info(work, iter - 1, compute_cost_function, false);
		    /* Check whether a termination criterion is triggered */
		    check_termination(work, false);
		  }
		  // Compute objective value in case it was not
		  // computed during the iterations
		  if (!compute_cost_function && has_solution(work.info)){
		    work.info.obj_val = compute_obj_val(work, work.x);
		  }

		  /* if max iterations reached, change status accordingly */
		  if (work.info.status == Status.UNSOLVED) {
		    if (!check_termination(work, true)) { // Try to check for approximate
		      work.info.status = Status.MAX_ITER_REACHED;
		    }
		  }

		  /* Update rho estimate */
		  work.info.rho_estimate = compute_rho_estimate(work);

		  // Polish the obtained solution
		  if (work.settings.polish && (work.info.status == Status.SOLVED))
		    Polish.polish(work);

		  // Store solution
		  store_solution(work);

		  return work.info.status;
	}
	
//	static int check_termination(OSQP.Workspace work, boolean approximate) {
//		  double eps_prim, eps_dual, eps_prim_inf, eps_dual_inf;
//		  int   exitflag;
//		  boolean   prim_res_check, dual_res_check, prim_inf_check, dual_inf_check;
//		  double eps_abs, eps_rel;
//
//		  // Initialize variables to 0
//		  exitflag       = 0;
//		  prim_res_check = false; dual_res_check = false;
//		  prim_inf_check = false; dual_inf_check = false;
//
//		  // Initialize tolerances
//		  eps_abs      = work.settings.eps_abs;
//		  eps_rel      = work.settings.eps_rel;
//		  eps_prim_inf = work.settings.eps_prim_inf;
//		  eps_dual_inf = work.settings.eps_dual_inf;
//
//		  // If residuals are too large, the problem is probably non convex
//		  if ((work.info.pri_res > OSQP_INFTY) ||
//		      (work.info.dua_res > OSQP_INFTY)){
//		    // Looks like residuals are diverging. Probably the problem is non convex!
//		    // Terminate and report it
//		    work.info.status = Status.NON_CVX;
//		    work.info.obj_val = OSQP_NAN;
//		    return 1;
//		  }
//
//		  // If approximate solution required, increase tolerances by 10
//		  if (approximate) {
//		    eps_abs      *= 10;
//		    eps_rel      *= 10;
//		    eps_prim_inf *= 10;
//		    eps_dual_inf *= 10;
//		  }
//
//		  // Check residuals
//		  if (work.data.m == 0) {
//		    prim_res_check = true; // No constraints -> Primal feasibility always satisfied
//		  }
//		  else {
//		    // Compute primal tolerance
//		    eps_prim = compute_pri_tol(work, eps_abs, eps_rel);
//
//		    // Primal feasibility check
//		    if (work.info.pri_res < eps_prim) {
//		      prim_res_check = true;
//		    } else {
//		      // Primal infeasibility check
//		      prim_inf_check = is_primal_infeasible(work, eps_prim_inf);
//		    }
//		  } // End check if m == 0
//
//		  // Compute dual tolerance
//		  eps_dual = compute_dua_tol(work, eps_abs, eps_rel);
//
//		  // Dual feasibility check
//		  if (work.info.dua_res < eps_dual) {
//		    dual_res_check = true;
//		  } else {
//		    // Check dual infeasibility
//		    dual_inf_check = is_dual_infeasible(work, eps_dual_inf);
//		  }
//
//		  // Compare checks to determine solver status
//		  if (prim_res_check && dual_res_check) {
//		    // Update final information
//		    if (approximate) {
//		      work.info.status = Status.SOLVED_INACCURATE;
//		    } else {
//		      work.info.status = Status.SOLVED;
//		    }
//		    exitflag = 1;
//		  }
//		  else if (prim_inf_check) {
//		    // Update final information
//		    if (approximate) {
//		      work.info.status = Status.PRIMAL_INFEASIBLE_INACCURATE;
//		    } else {
//		      work.info.status = Status.PRIMAL_INFEASIBLE;
//		    }
//
//		    if (work.settings.scaling!=0 && !work.settings.scaled_termination) {
//		      // Update infeasibility certificate
//		    	LinAlg.vec_ew_prod(work.scaling.E, work.delta_y, work.delta_y);
//		    }
//		    work.info.obj_val = OSQP_INFTY;
//		    exitflag            = 1;
//		  }
//		  else if (dual_inf_check) {
//		    // Update final information
//		    if (approximate) {
//		    	work.info.status = Status.DUAL_INFEASIBLE_INACCURATE;
//		    } else {
//		    	work.info.status = Status.DUAL_INFEASIBLE;
//		    }
//
//		    if (work.settings.scaling!=0 && !work.settings.scaled_termination) {
//		      // Update infeasibility certificate
//		      LinAlg.vec_ew_prod(work.scaling.D, work.delta_x, work.delta_x);
//		    }
//		    work.info.obj_val = -OSQP_INFTY;
//		    exitflag            = 1;
//		  }
//
//		  return exitflag;
//		}
	
	public static void cold_start(OSQP.Workspace work) {
		  Arrays.fill(work.x, 0);
		  Arrays.fill(work.z, 0);
		  Arrays.fill(work.y, 0);
		}
	
	static double compute_rho_estimate(Workspace work) {
		  double pri_res, dua_res;           // Primal and dual residuals
		  double pri_res_norm, dua_res_norm; // Normalization for the residuals
		  double temp_res_norm;              // Temporary residual norm
		  double rho_estimate;               // Rho estimate value
		  
		// Get problem dimensions
		  int n = work.data.n;
		  int m = work.data.m;


		  // Get primal and dual residuals
		  pri_res = LinAlg.vec_norm_inf(work.z_prev, m);
		  dua_res = LinAlg.vec_norm_inf(work.x_prev, n);

		  // Normalize primal residual
		  pri_res_norm  = LinAlg.vec_norm_inf(work.z, m);           // ||z||
		  temp_res_norm = LinAlg.vec_norm_inf(work.Ax, m);          // ||Ax||
		  pri_res_norm  = Math.max(pri_res_norm, temp_res_norm); // max (||z||,||Ax||)
		  pri_res      /= (pri_res_norm + 1e-10);             // Normalize primal
		                                                      // residual (prevent 0
		                                                      // division)

		  // Normalize dual residual
		  dua_res_norm  = LinAlg.vec_norm_inf(work.data.q,n);     // ||q||
		  temp_res_norm = LinAlg.vec_norm_inf(work.Aty,n);         // ||A' y||
		  dua_res_norm  = Math.max(dua_res_norm, temp_res_norm);
		  temp_res_norm = LinAlg.vec_norm_inf(work.Px,n);          //  ||P x||
		  dua_res_norm  = Math.max(dua_res_norm, temp_res_norm); // max(||q||,||A' y||,||P
		                                                      // x||)
		  dua_res      /= (dua_res_norm + 1e-10);             // Normalize dual residual
		                                                      // (prevent 0 division)


		  // Return rho estimate
		  rho_estimate = work.settings.rho * Math.sqrt(pri_res / (dua_res + 1e-10)); // (prevent
		                                                                            // 0
		                                                                            // division)
		  rho_estimate = Math.min(Math.max(rho_estimate, RHO_MIN), RHO_MAX);              // Constrain
		                                                                            // rho
		                                                                            // values
		  return rho_estimate;
		}
	
	public static boolean osqp_update_rho(Workspace work, double rho_new) {
		  int i;
		  boolean exitflag = false;

		  // Check value of rho
		  if (rho_new <= 0) {
		    throw new IllegalStateException("rho must be positive");
		  }


		  // Update rho in settings
		  work.settings.rho = Math.min(Math.max(rho_new, RHO_MIN), RHO_MAX);

		  // Update rho_vec and rho_inv_vec
		  for (i = 0; i < work.data.m; i++) {
		    if (work.constr_type[i] == 0) {
		      // Inequalities
		      work.rho_vec[i]     = work.settings.rho;
		      work.rho_inv_vec[i] = 1. / work.settings.rho;
		    }
		    else if (work.constr_type[i] == 1) {
		      // Equalities
		      work.rho_vec[i]     = RHO_EQ_OVER_RHO_INEQ * work.settings.rho;
		      work.rho_inv_vec[i] = 1. / work.rho_vec[i];
		    }
		  }

		  // Update rho_vec in KKT matrix
		  exitflag = work.linsys_solver.update_rho_vec(work.rho_vec);
		  return exitflag;

	}
	
	static boolean adapt_rho(Workspace work) {
		  double rho_new = compute_rho_estimate(work);
		  boolean exitflag = false;

		  // Set rho estimate in info
		  work.info.rho_estimate = rho_new;

		  // Check if the new rho is large or small enough and update it in case
		  if ((rho_new > work.settings.rho * work.settings.adaptive_rho_tolerance) ||
		      (rho_new < work.settings.rho /  work.settings.adaptive_rho_tolerance)) {
			exitflag = osqp_update_rho(work, rho_new);
		    work.info.rho_updates ++;
		  }
		  return exitflag;
		}
	
	static void set_rho_vec(Workspace work) {
		  int i;

		  work.settings.rho = Math.min(Math.max(work.settings.rho, RHO_MIN), RHO_MAX);

		  for (i = 0; i < work.data.m; i++) {
		    if ((work.data.l[i] < -OSQP_INFTY * MIN_SCALING) &&
		        (work.data.u[i] > OSQP_INFTY * MIN_SCALING)) {
		      // Loose bounds
		      work.constr_type[i] = -1;
		      work.rho_vec[i]     = RHO_MIN;
		    } else if (work.data.u[i] - work.data.l[i] < RHO_TOL) {
		      // Equality constraints
		      work.constr_type[i] = 1;
		      work.rho_vec[i]     = RHO_EQ_OVER_RHO_INEQ * work.settings.rho;
		    } else {
		      // Inequality constraints
		      work.constr_type[i] = 0;
		      work.rho_vec[i]     = work.settings.rho;
		    }
		    work.rho_inv_vec[i] = 1. / work.rho_vec[i];
		  }
	}
	
	public void update_rho_vec() {
		  int i, constr_type_changed;

		  constr_type_changed = 0;

		  for (i = 0; i < work.data.m; i++) {
		    if ((work.data.l[i] < -OSQP_INFTY * MIN_SCALING) &&
		        (work.data.u[i] > OSQP_INFTY * MIN_SCALING)) {
		      // Loose bounds
		      if (work.constr_type[i] != -1) {
		        work.constr_type[i] = -1;
		        work.rho_vec[i]     = RHO_MIN;
		        work.rho_inv_vec[i] = 1. / RHO_MIN;
		        constr_type_changed  = 1;
		      }
		    } else if (work.data.u[i] - work.data.l[i] < RHO_TOL) {
		      // Equality constraints
		      if (work.constr_type[i] != 1) {
		        work.constr_type[i] = 1;
		        work.rho_vec[i]     = RHO_EQ_OVER_RHO_INEQ * work.settings.rho;
		        work.rho_inv_vec[i] = 1. / work.rho_vec[i];
		        constr_type_changed  = 1;
		      }
		    } else {
		      // Inequality constraints
		      if (work.constr_type[i] != 0) {
		        work.constr_type[i] = 0;
		        work.rho_vec[i]     = work.settings.rho;
		        work.rho_inv_vec[i] = 1. / work.settings.rho;
		        constr_type_changed  = 1;
		      }
		    }
		  }
		// Update rho_vec in KKT matrix if constraints type has changed
		  if (constr_type_changed == 1) {
		    work.linsys_solver.update_rho_vec(work.rho_vec);
		  }

	}
	
	static void compute_rhs(Workspace work) {
		  int i; // Index

		  for (i = 0; i < work.data.n; i++) {
		    // Cycle over part related to x variables
		    work.xz_tilde[i] = work.settings.sigma * work.x_prev[i] -
		                        work.data.q[i];
		  }

		  for (i = 0; i < work.data.m; i++) {
		    // Cycle over dual variable in the first step (nu)
		    work.xz_tilde[i + work.data.n] = work.z_prev[i] - work.rho_inv_vec[i] *
		                                        work.y[i];
		  }
		}
	
	static void update_xz_tilde(Workspace work) {
		  // Compute right-hand side
		  compute_rhs(work);

		  // Solve linear system
		  work.linsys_solver.solve(work.linsys_solver.sol,work.xz_tilde);
		}

	static void update_x(Workspace work) {
		  int i;

		  // update x
		  for (i = 0; i < work.data.n; i++) {
		    work.x[i] = work.settings.alpha * work.xz_tilde[i] +
		                 (1.0 - work.settings.alpha) * work.x_prev[i];
		  }

		  // update delta_x
		  for (i = 0; i < work.data.n; i++) {
		    work.delta_x[i] = work.x[i] - work.x_prev[i];
		  }
		}

		static void update_z(Workspace work) {
		  int i;

		  // update z
		  for (i = 0; i < work.data.m; i++) {
		    work.z[i] = work.settings.alpha * work.xz_tilde[i + work.data.n] +
		                 (1.0 - work.settings.alpha) * work.z_prev[i] +
		                 work.rho_inv_vec[i] * work.y[i];
		  }

		  // project z
		  Projection.project(work, work.z);
		}

		static void update_y(Workspace work) {
		  int i; // Index

		  for (i = 0; i < work.data.m; i++) {
		    work.delta_y[i] = work.rho_vec[i] *
		                       (work.settings.alpha *
		                        work.xz_tilde[i + work.data.n] +
		                        (1.0 - work.settings.alpha) * work.z_prev[i] -
		                        work.z[i]);
		    work.y[i] += work.delta_y[i];
		  }
		}
	
		static double compute_obj_val(Workspace work, double[] x) {
			  double obj_val;

			  obj_val = LinAlg.quad_form(work.data.P, x) +
					  LinAlg.vec_prod(work.data.q, x, x.length);

			  if (work.settings.scaling>0) {
			    obj_val *= work.scaling.cinv;
			  }

			  return obj_val;
		}

		static double compute_pri_res(Workspace work, double[] x, double[] z) {
		  // NB: Use z_prev as working vector
		  // pr = Ax - z

			LinAlg.mat_vec(work.data.A, x, work.Ax, 0,0, 0); // Ax
			LinAlg.vec_add_scaled(work.z_prev, work.Ax, z, work.data.m, -1.0);

		  // If scaling active . rescale residual
		  if (work.settings.scaling>0 && !work.settings.scaled_termination) {
		    return LinAlg.vec_scaled_norm_inf(work.scaling.Einv, work.z_prev, work.data.m);
		  }

		  // Return norm of the residual
		  return LinAlg.vec_norm_inf(work.z_prev, work.data.m);
		}

		static double compute_pri_tol(Workspace work, double eps_abs, double eps_rel) {
		  double max_rel_eps, temp_rel_eps;

		  // max_rel_eps = max(||z||, ||A x||)
		  if (work.settings.scaling>0 && !work.settings.scaled_termination) {
		    // ||Einv * z||
		    max_rel_eps =
		    		LinAlg.vec_scaled_norm_inf(work.scaling.Einv, work.z, work.data.m);

		    // ||Einv * A * x||
		    temp_rel_eps = LinAlg.vec_scaled_norm_inf(work.scaling.Einv,
		                                       work.Ax, work.data.m);

		    // Choose maximum
		    max_rel_eps = Math.max(max_rel_eps, temp_rel_eps);
		  } else { // No unscaling required
		    // ||z||
		    max_rel_eps = LinAlg.vec_norm_inf(work.z, work.data.m);

		    // ||A * x||
		    temp_rel_eps = LinAlg.vec_norm_inf(work.Ax, work.data.m);

		    // Choose maximum
		    max_rel_eps = Math.max(max_rel_eps, temp_rel_eps);
		  }

		  // eps_prim
		  return eps_abs + eps_rel * max_rel_eps;
		}

		static double compute_dua_res(Workspace work, double[] x, double[] y) {
		  // NB: Use x_prev as temporary vector
		  // NB: Only upper triangular part of P is stored.
		  // dr = q + A'*y + P*x

		  // dr = q
			LinAlg.prea_vec_copy(work.data.q, work.x_prev, work.data.n);

		  // P * x (upper triangular part)
			LinAlg.mat_vec(work.data.P, x, work.Px, 0,0, 0);

		  // P' * x (lower triangular part with no diagonal)
			LinAlg.mat_tpose_vec(work.data.P, x, work.Px,0,0, 1, true);

		  // dr += P * x (full P matrix)
			LinAlg.vec_add_scaled(work.x_prev, work.x_prev, work.Px, work.data.n, 1.0);

		  // dr += A' * y
		  if (work.data.m > 0) {
			  LinAlg.mat_tpose_vec(work.data.A, y, work.Aty, 0, 0, 0, false);
			  LinAlg.vec_add_scaled(work.x_prev, work.x_prev, work.Aty, work.data.n, 1);
		  }

		  // If scaling active . rescale residual
		  if (work.settings.scaling>0 && !work.settings.scaled_termination) {
		    return work.scaling.cinv * LinAlg.vec_scaled_norm_inf(work.scaling.Dinv,
		                                                     work.x_prev, work.data.n);
		  }

		  return LinAlg.vec_norm_inf(work.x_prev, work.data.n);
		}

		static double compute_dua_tol(Workspace work, double eps_abs, double eps_rel) {
			double max_rel_eps, temp_rel_eps;

		  // max_rel_eps = max(||q||, ||A' y|, ||P x||)
		  if (work.settings.scaling>0 && !work.settings.scaled_termination) {
		    // || Dinv q||
		    max_rel_eps = LinAlg.vec_scaled_norm_inf(work.scaling.Dinv,
		                                      work.data.q, work.data.n);

		    // || Dinv A' y ||
		    temp_rel_eps = LinAlg.vec_scaled_norm_inf(work.scaling.Dinv,
		                                       work.Aty, work.data.n);
		    max_rel_eps = Math.max(max_rel_eps, temp_rel_eps);

		    // || Dinv P x||
		    temp_rel_eps = LinAlg.vec_scaled_norm_inf(work.scaling.Dinv,
		                                       work.Px, work.data.n);
		    max_rel_eps = Math.max(max_rel_eps, temp_rel_eps);

		    // Multiply by cinv
		    max_rel_eps *= work.scaling.cinv;
		  } else { // No scaling required
		    // ||q||
		    max_rel_eps = LinAlg.vec_norm_inf(work.data.q, work.data.n);

		    // ||A'*y||
		    temp_rel_eps = LinAlg.vec_norm_inf(work.Aty, work.data.n);
		    max_rel_eps  = Math.max(max_rel_eps, temp_rel_eps);

		    // ||P*x||
		    temp_rel_eps = LinAlg.vec_norm_inf(work.Px, work.data.n);
		    max_rel_eps  = Math.max(max_rel_eps, temp_rel_eps);
		  }

		  // eps_dual
		  return eps_abs + eps_rel * max_rel_eps;
		}

		static boolean is_primal_infeasible(Workspace work, double eps_prim_inf) {
		  // This function checks for the primal infeasibility termination criteria.
		  //
		  // 1) A' * delta_y < eps * ||delta_y||
		  //
		  // 2) u'*max(delta_y, 0) + l'*min(delta_y, 0) < -eps * ||delta_y||
		  //

		  int i; // Index for loops
		  double norm_delta_y;
		  double ineq_lhs = 0.0;

		  // Project delta_y onto the polar of the recession cone of [l,u]
		  for (i = 0; i < work.data.m; i++) {
		    if (work.data.u[i] > OSQP_INFTY * MIN_SCALING) {          // Infinite upper bound
		      if (work.data.l[i] < -OSQP_INFTY * MIN_SCALING) {       // Infinite lower bound
		        // Both bounds infinite
		        work.delta_y[i] = 0.0;
		      } else {
		        // Only upper bound infinite
		        work.delta_y[i] = Math.min(work.delta_y[i], 0.0);
		      }
		    } else if (work.data.l[i] < -OSQP_INFTY * MIN_SCALING) {  // Infinite lower bound
		      // Only lower bound infinite
		      work.delta_y[i] = Math.max(work.delta_y[i], 0.0);
		    }
		  }

		  // Compute infinity norm of delta_y (unscale if necessary)
		  if (work.settings.scaling>0 && !work.settings.scaled_termination) {
		    // Use work.Adelta_x as temporary vector
		    LinAlg.vec_ew_prod(work.scaling.E, work.delta_y, work.Adelta_x, work.data.m);
		    norm_delta_y = LinAlg.vec_norm_inf(work.Adelta_x, work.data.m);
		  } else {
		    norm_delta_y = LinAlg.vec_norm_inf(work.delta_y, work.data.m);
		  }

		  if (norm_delta_y > eps_prim_inf) { // ||delta_y|| > 0

		    for (i = 0; i < work.data.m; i++) {
		      ineq_lhs += work.data.u[i] * Math.max(work.delta_y[i], 0) + 
		                  work.data.l[i] * Math.max(work.delta_y[i], 0);
		    }

		    // Check if the condition is satisfied: ineq_lhs < -eps
		    if (ineq_lhs < -eps_prim_inf * norm_delta_y) {
		      // Compute and return ||A'delta_y|| < eps_prim_inf
		    	LinAlg.mat_tpose_vec(work.data.A, work.delta_y, work.Atdelta_y, 0, 0,0, false);

		      // Unscale if necessary
		      if (work.settings.scaling>0 && !work.settings.scaled_termination) {
		        LinAlg.vec_ew_prod(work.scaling.Dinv,
		                    work.Atdelta_y,
		                    work.Atdelta_y, work.data.n);
		      }

		      return LinAlg.vec_norm_inf(work.Atdelta_y, work.data.n) < eps_prim_inf * norm_delta_y;
		    }
		  }

		  // Conditions not satisfied . not primal infeasible
		  return false;
		}

		static boolean is_dual_infeasible(Workspace work, double eps_dual_inf) {
		  // This function checks for the scaled dual infeasibility termination
		  // criteria.
		  //
		  // 1) q * delta_x < - eps * || delta_x ||
		  //
		  // 2) ||P * delta_x || < eps * || delta_x ||
		  //
		  // 3) . (A * delta_x)_i > -eps * || delta_x ||,    l_i != -inf
		  //    . (A * delta_x)_i <  eps * || delta_x ||,    u_i != inf
		  //


		  int   i; // Index for loops
		  double norm_delta_x;
		  double cost_scaling;

		  // Compute norm of delta_x
		  if (work.settings.scaling>0 && !work.settings.scaled_termination) { // Unscale
		                                                                        // if
		                                                                        // necessary
		    norm_delta_x = LinAlg.vec_scaled_norm_inf(work.scaling.D,
		                                       work.delta_x, work.data.n);
		    cost_scaling = work.scaling.c;
		  } else {
		    norm_delta_x = LinAlg.vec_norm_inf(work.delta_x, work.data.n);
		    cost_scaling = 1.0;
		  }

		  // Prevent 0 division || delta_x || > 0
		  if (norm_delta_x > eps_dual_inf) {
		    // Normalize delta_x by its norm

		    /* vec_mult_scalar(work.delta_x, 1./norm_delta_x, work.data.n); */

		    // Check first if q'*delta_x < 0
		    if (LinAlg.vec_prod(work.data.q, work.delta_x, work.data.n) <
		        -cost_scaling * eps_dual_inf * norm_delta_x) {
		      // Compute product P * delta_x (NB: P is store in upper triangular form)
		    LinAlg.mat_vec(work.data.P, work.delta_x, work.Pdelta_x, 0,0,0);
		    LinAlg.mat_tpose_vec(work.data.P, work.delta_x, work.Pdelta_x, 0, 0, 1, true);

		      // Scale if necessary
		      if (work.settings.scaling>0 && !work.settings.scaled_termination) {
		        LinAlg.vec_ew_prod(work.scaling.Dinv,
		                    work.Pdelta_x,
		                    work.Pdelta_x, work.data.n);
		      }

		      // Check if || P * delta_x || = 0
		      if (LinAlg.vec_norm_inf(work.Pdelta_x, work.data.n) <
		          cost_scaling * eps_dual_inf * norm_delta_x) {
		        // Compute A * delta_x
		    	  LinAlg.mat_vec(work.data.A, work.delta_x, work.Adelta_x, 0,0, 0);

		        // Scale if necessary
		        if (work.settings.scaling>0 && !work.settings.scaled_termination) {
		        	LinAlg.vec_ew_prod(work.scaling.Einv,
		                      work.Adelta_x,
		                      work.Adelta_x, work.data.m);
		        }

		        // De Morgan Law Applied to dual infeasibility conditions for A * x
		        // NB: Note that MIN_SCALING is used to adjust the infinity value
		        //     in case the problem is scaled.
		        for (i = 0; i < work.data.m; i++) {
		          if (((work.data.u[i] < OSQP_INFTY * MIN_SCALING) &&
		               (work.Adelta_x[i] >  eps_dual_inf * norm_delta_x)) ||
		              ((work.data.l[i] > -OSQP_INFTY * MIN_SCALING) &&
		               (work.Adelta_x[i] < -eps_dual_inf * norm_delta_x))) {
		            // At least one condition not satisfied . not dual infeasible
		            return false;
		          }
		        }

		        // All conditions passed . dual infeasible
		        return true;
		      }
		    }
		  }

		  // Conditions not satisfied . not dual infeasible
		  return false;
		}

		static boolean has_solution(Info info){

		  return ((info.status != Status.PRIMAL_INFEASIBLE) &&
		      (info.status != Status.PRIMAL_INFEASIBLE_INACCURATE) &&
		      (info.status != Status.DUAL_INFEASIBLE) &&
		      (info.status != Status.DUAL_INFEASIBLE_INACCURATE) &&
		      (info.status != Status.NON_CVX));

		}

		static void store_solution(Workspace work) {

		  if (has_solution(work.info)) {
			  LinAlg.prea_vec_copy(work.x, work.solution.x, work.data.n); // primal
			  LinAlg.prea_vec_copy(work.y, work.solution.y, work.data.m); // dual

		    // Unscale solution if scaling has been performed
		    if (work.settings.scaling>0)
		    	Scaling.unscale_solution(work);
		  } else {
		    // No solution present. Solution is NaN
			  LinAlg.vec_set_scalar(work.solution.x, OSQP_NAN, work.data.n);
			  LinAlg.vec_set_scalar(work.solution.y, OSQP_NAN, work.data.m);


		    // Normalize infeasibility certificates if embedded is off
		    // NB: It requires a division
		    if ((work.info.status == Status.PRIMAL_INFEASIBLE) ||
		        ((work.info.status == Status.PRIMAL_INFEASIBLE_INACCURATE))) {
		      double norm_vec = LinAlg.vec_norm_inf(work.delta_y, work.data.m);
		      LinAlg.vec_mult_scalar(work.delta_y, 1. / norm_vec, work.data.m);
		    }

		    if ((work.info.status == Status.DUAL_INFEASIBLE) ||
		        ((work.info.status == Status.DUAL_INFEASIBLE_INACCURATE))) {
		      double norm_vec = LinAlg.vec_norm_inf(work.delta_x, work.data.n);
		      LinAlg.vec_mult_scalar(work.delta_x, 1. / norm_vec, work.data.n);
		    }


		    // Cold start iterates to 0 for next runs (they cannot start from NaN)
		    cold_start(work);
		  }
		}

		protected static void update_info(Workspace work,
		                 int          iter,
		                 boolean          compute_objective,
		                 boolean          polish) {
		  double[] x, z, y;                   // Allocate pointers to variables


		  if (polish) {
		    x       = work.pol.x;
		    y       = work.pol.y;
		    z       = work.pol.z;

		  } else {
			  x                = work.x;
			  y                = work.y;
			  z                = work.z;
			  work.info.iter = iter; // Update iteration number
		}



		  // Compute the objective if needed
		  if (compute_objective) {
			  if (polish)
				  work.pol.obj_val = compute_obj_val(work, x);
			  else
				  work.info.obj_val = compute_obj_val(work, x);
		  }

		  // Compute primal residual
		  if (work.data.m == 0) {
		    // No constraints . Always primal feasible
			  if (polish)
				  work.pol.pri_res = 0.;
			  else
				  work.info.pri_res = 0.;
				  
		  } else {
			  if (polish)
				  work.pol.pri_res = compute_pri_res(work, x, z);
			  else
				  work.info.pri_res = compute_pri_res(work, x, z);
		  }

		  // Compute dual residual
		  work.info.dua_res = compute_dua_res(work, x, y);

		}


		void reset_info(Info info) {

		  info.status = Status.UNSOLVED; // Problem is unsolved

		  info.rho_updates = 0;              // Rho updates are now 0
		}


		static boolean check_termination(Workspace work, boolean approximate) {
		  double eps_prim, eps_dual, eps_prim_inf, eps_dual_inf;
		  boolean   exitflag = false;
		  boolean   prim_res_check=false, dual_res_check=false, prim_inf_check=false, dual_inf_check=false;
		  double eps_abs, eps_rel;

		  // Initialize variables to 0


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
		    return true;
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
		    prim_res_check = true; // No constraints . Primal feasibility always satisfied
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
		      work.info.status = Status.SOLVED_INACCURATE;
		    } else {
		      work.info.status = Status.SOLVED;
		    }
		    exitflag = true;
		  }
		  else if (prim_inf_check) {
		    // Update final information
		    if (approximate) {
		      work.info.status = Status.PRIMAL_INFEASIBLE_INACCURATE;
		    } else {
		      work.info.status = Status.PRIMAL_INFEASIBLE;
		    }

		    if (work.settings.scaling>0 && !work.settings.scaled_termination) {
		      // Update infeasibility certificate
		      LinAlg.vec_ew_prod(work.scaling.E, work.delta_y, work.delta_y, work.data.m);
		    }
		    work.info.obj_val = OSQP_INFTY;
		    exitflag            = true;
		  }
		  else if (dual_inf_check) {
		    // Update final information
		    if (approximate) {
		      work.info.status = Status.DUAL_INFEASIBLE_INACCURATE;
		    } else {
		      work.info.status = Status.DUAL_INFEASIBLE;
		    }

		    if (work.settings.scaling>0 && !work.settings.scaled_termination) {
		      // Update infeasibility certificate
		      LinAlg.vec_ew_prod(work.scaling.D, work.delta_x, work.delta_x, work.data.n);
		    }
		    work.info.obj_val = -OSQP_INFTY;
		    exitflag            = true;
		  }

		  return exitflag;
		}
		
		public void update_lin_cost(double[] q_new) {

			  // Replace q by the new vector
			LinAlg.prea_vec_copy(q_new, work.data.q, work.data.n);

			  // Scaling
			  if (work.settings.scaling>0) {
			    LinAlg.vec_ew_prod(work.scaling.D, work.data.q, work.data.q, work.data.n);
			    LinAlg.vec_mult_scalar(work.data.q, work.scaling.c, work.data.n);
			  }

			  // Reset solver information
			  reset_info(work.info);
			}
		
		
		public void update_bounds(double[] l_new, double[] u_new) {
			int i;

		
			// Check if lower bound is smaller than upper bound
			for (i = 0; i < work.data.m; i++) {
				if (l_new[i] > u_new[i]) {
					throw new IllegalArgumentException("Lower bound is greater than upper bound in row "+i);
				}
			}
			
			// Replace l and u by the new vectors
			LinAlg.prea_vec_copy(l_new, work.data.l, work.data.m);
			LinAlg.prea_vec_copy(u_new, work.data.u, work.data.m);
			
			// Scaling
			if (work.settings.scaling>0) {
				LinAlg.vec_ew_prod(work.scaling.E, work.data.l, work.data.l, work.data.m);
				LinAlg.vec_ew_prod(work.scaling.E, work.data.u, work.data.u, work.data.m);
			}
			
			// Reset solver information
			reset_info(work.info);

		}
		
		public void warm_start_x(double[] x) {

			  // Update warm_start setting to true
			  if (!work.settings.warm_start) work.settings.warm_start = true;

			  // Copy primal variable into the iterate x
			  LinAlg.prea_vec_copy(x, work.x, work.data.n);

			  // Scale iterate
			  if (work.settings.scaling>0) {
				  LinAlg.vec_ew_prod(work.scaling.Dinv, work.x, work.x, work.data.n);
			  }

			  // Compute Ax = z and store it in z
			  LinAlg.mat_vec(work.data.A, work.x, work.z, 0, 0, 0);

			}

		public void warm_start_y(double[] y) {
			  // Update warm_start setting to true
			  if (!work.settings.warm_start) work.settings.warm_start = true;

			  // Copy primal variable into the iterate y
			  LinAlg.prea_vec_copy(y, work.y, work.data.m);

			  // Scale iterate
			  if (work.settings.scaling>0) {
				  LinAlg.vec_ew_prod(work.scaling.Einv, work.y, work.y, work.data.m);
				  LinAlg.vec_mult_scalar(work.y, work.scaling.c, work.data.m);
			  }

			}
		
		public void  update_P(double[] Px_new,
                int[]   Px_new_idx) {
			

			if (work.settings.scaling>0) {
			// Unscale data
				Scaling.unscale_data(work);
			}
			
			// Update P elements
			if (Px_new_idx!= null) { // Change only Px_new_idx
				for (int i = 0; i < Px_new_idx.length; i++) {
					work.data.P.Ax[Px_new_idx[i]] = Px_new[i];
				}
			}
			else // Change whole P
			{
				int nnzP = work.data.P.Ap[work.data.P.n];
				for (int i = 0; i < nnzP; i++) {
					work.data.P.Ax[i] = Px_new[i];
				}
			}
			
			if (work.settings.scaling>0) {
			// Scale data
				Scaling.scale_data(work);
			}
			
			// Update linear system structure with new data
			work.linsys_solver.update_solver_matrices(work.data.P,work.data.A);
			
			// Reset solver information
			reset_info(work.info);
		}
		
		public void update_A(double[] Ax_new,
                 int[]   Ax_new_idx) {
			int i;        // For indexing
			int nnzA;     // Number of nonzeros in A
			
			
			nnzA = work.data.A.Ap[work.data.A.n];
			
		
			
			if (work.settings.scaling>0) {
				// Unscale data
				Scaling.unscale_data(work);
			}
			
			// Update A elements
			if (Ax_new_idx!=null) { // Change only Ax_new_idx
			for (i = 0; i < Ax_new_idx.length; i++) {
			  work.data.A.Ax[Ax_new_idx[i]] = Ax_new[i];
			}
			}
			else { // Change whole A
			for (i = 0; i < nnzA; i++) {
			  work.data.A.Ax[i] = Ax_new[i];
			}
			}
			
			if (work.settings.scaling>0) {
				// Scale data
				Scaling.scale_data(work);
			}
			
			// Update linear system structure with new data
			work.linsys_solver.update_solver_matrices(work.data.P,
			                                              work.data.A);
			
			// Reset solver information
			reset_info(work.info);

	
			}
		
		public void update_P_A(double[] Px_new,
                int[]   Px_new_idx,
                double[] Ax_new,
                int[]   Ax_new_idx) {
			int i;          // For indexing
			int nnzP, nnzA; // Number of nonzeros in P and A
			
			// Check if workspace has been initialized
			
	
			
			nnzP = work.data.P.Ap[work.data.P.n];
			nnzA = work.data.A.Ap[work.data.A.n];
			
		
			
			if (work.settings.scaling>0) {
			// Unscale data
				Scaling.unscale_data(work);
			}
			
			// Update P elements
			if (Px_new_idx!=null) { // Change only Px_new_idx
				for (i = 0; i < Px_new_idx.length; i++) {
					work.data.P.Ax[Px_new_idx[i]] = Px_new[i];
				}
			}
			else // Change whole P
			{
				for (i = 0; i < nnzP; i++) {
					work.data.P.Ax[i] = Px_new[i];
				}
			}
			
			// Update A elements
			if (Ax_new_idx!=null) { // Change only Ax_new_idx
				for (i = 0; i < Ax_new_idx.length; i++) {
					work.data.A.Ax[Ax_new_idx[i]] = Ax_new[i];
				}
			}
			else { // Change whole A
				for (i = 0; i < nnzA; i++) {
					work.data.A.Ax[i] = Ax_new[i];
				}
			}
			
			if (work.settings.scaling>0) {
			// Scale data
				Scaling.scale_data(work);
			}
			
			// Update linear system structure with new data
			work.linsys_solver.update_solver_matrices(work.data.P,
			                                            work.data.A);
			
			// Reset solver information
			reset_info(work.info);
			
	
			
			}
		
		public double[] getPrimalSolution() {
			return work.solution.x;
		}
		
		public double[] getDualSolution() {
			return work.solution.y;
		}
		
		public double getObjectiveValue() {
			return work.info.obj_val;
		}

}
