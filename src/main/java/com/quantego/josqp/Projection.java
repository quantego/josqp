package com.quantego.josqp;

public class Projection {
	
	public static void project(OSQP.Workspace work, double[] z) {
		  int i, m;

		  m = work.data.m;

		  for (i = 0; i < m; i++) {
		    z[i] = Math.min(Math.max(z[i],
		                       work.data.l[i]), // Between lower
		                 work.data.u[i]);       // and upper bounds
		  }
		}

	public static void project_normalcone(OSQP.Workspace work, double[] z, double[] y) {
		  int i, m;

		  // NB: Use z_prev as temporary vector

		  m = work.data.m;

		  for (i = 0; i < m; i++) {
		    work.z_prev[i] = z[i] + y[i];
		    z[i]            = Math.min(Math.max(work.z_prev[i], work.data.l[i]),
		                            work.data.u[i]);
		    y[i] = work.z_prev[i] - z[i];
		  }
		}

}
