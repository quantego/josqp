package com.quantego.josqp;

public class OSQP {
	
	public static class Scaling {
		double  c;    ///< cost function scaling
		double[] D;    ///< primal variable scaling
		double[] E;    ///< dual variable scaling
		double  cinv; ///< cost function rescaling
		double[] Dinv; ///< primal variable rescaling
		double[] Einv; ///< dual variable rescaling
	}
	
	public static class Info {
		public enum STATUS {
			SOLVED, UNBOUNDED, INFEASIBLE, NOT_CONVERGED, ERROR
		}
		
		public int iter;          ///< number of iterations taken
		public STATUS status;     ///< status string, e.g. 'solved'
		public int status_val;    ///< status as c_int, defined in constants.h

		public int status_polish; ///< polish status: successful (1), unperformed (0), (-1) unsuccessful

		public float obj_val;     ///< primal objective
		public float pri_res;     ///< norm of primal residual
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

		int scaled_termination;              ///< boolean, use scaled termination criteria
		int check_termination;               ///< integer, check termination interval; if 0, then termination checking is disabled
		int warm_start;                      ///< boolean, warm start

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

}
