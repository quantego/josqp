package com.quantego.josqp;

class Projection {
	
	public static void project(OSQP.Workspace work, double[] z){
        NumericKernels.INSTANCE.project(work, z);
    }

	public static void project_normalcone(OSQP.Workspace work, double[] z, double[] y){
        NumericKernels.INSTANCE.project_normalcone(work, z, y);
    }

}
