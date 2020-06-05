package com.quantego.josqp;

public interface LinSys {
	
	public enum TYPE {
		QLDL
	}
	
	public double[] solve(double[] b);

}
