package josqp;

import java.util.Arrays;

import com.quantego.josqp.AMD;

public class AMDSimpleDemo {
	
	public static void main(String... args) {
		int n=5;
		int[] Ap = {0,2,6,10,12,14};
		int[] Ai = {0,1,0,1,2,4,1,2,3,4,2,3,1,4};
		int[] P = new int[5];
		AMD.amd_order(n,Ap,Ai,P);
		System.out.println(Arrays.toString(P));
	}

}
