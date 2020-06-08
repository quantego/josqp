package josqp;

import com.quantego.josqp.CSCMatrix;
import com.quantego.josqp.OSQP;

public class Demo {
	
	public static void main(String... args) {
		double[] P_x = { 4.0, 1.0, 2.0 };
		  int  P_nnz  = 3;
		  int[]   P_i = { 0, 0, 1 };
		  int[]   P_p = { 0, 1, 3 };
		  double[] q   = { 1.0, 1.0 };
		  double[] A_x = { 1.0, 1.0, 1.0, 1.0 };
		  int   A_nnz  = 4;
		  int[]   A_i = { 0, 1, 0, 2 };
		  int[]   A_p = { 0, 2, 4 };
		  double[] l   = { 1.0, 0.0, 0.0 };
		  double[] u   = { 1.0, 0.7, 0.7 };
		  int n = 2;
		  int m = 3;
		  
		  OSQP.Settings settings = new OSQP.Settings();
			CSCMatrix A = new CSCMatrix(m,n,A_nnz,A_p,A_i,A_x);
			CSCMatrix P = new CSCMatrix(n,n,P_nnz,P_p,P_i,P_x);
			OSQP.Data data = new OSQP.Data(n,m,P,A,q,l,u);
			OSQP opt = new OSQP(data,settings);
			opt.solve();
			OSQP.Info info = opt.getWorkspace().info;
			System.out.println(info.obj_val+" "+info.status.toString());
	}
	

}
