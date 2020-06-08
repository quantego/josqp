package com.quantego.josqp;

import java.util.Arrays;

public class QDLDL  {
    final static int QDLDL_UNKNOWN = -1;
    final static boolean QDLDL_USED = true;
    final static boolean QDLDL_UNUSED = false;
    final static int QDLDL_INT_MAX = Integer.MAX_VALUE; //TODO: what is this?
    
    LDL _ldl;
    
    /*public QDLDL(int An, int Am, int[] Ap, int[] Ai, double[] Ax) {
    	Etree e = computeEtree(An,Ap,Ai);
    	CSCMatrix A = new CSCMatrix(An,Am,Ax.length,Ap,Ai,Ax);
    	_ldl = QDLDL.decompose(A, e);
    }*/
    
    public QDLDL() {
    	
    }
    
    public boolean factor(CSCMatrix A) {
    	Etree e = computeEtree(A.n,A.Ap,A.Ai);
    	_ldl = QDLDL.decompose(A, e);
    	return true; //TODO
    }
    
    public double[] solve(double[] b) {
    	return QDLDL.solve(_ldl, b);
    }
    
    public static Etree computeEtree(int n, int[] Ap, int[] Ai) {
    	final int[] work = new int[n];
    	final int[] Lnz = new int[n];
    	final int[] etree = new int[n];
    	
    	int sumLnz = 0;
    	int i=0, j=0, p=0;


        for(i = 0; i < n; i++){
        // zero out Lnz and work.  Set all etree values to unknown
            etree[i] = QDLDL_UNKNOWN;

            //Abort if A doesn't have at least one entry
            //one entry in every column
            if(Ap[i] == Ap[i+1]){
            	throw new IllegalStateException("A doesn't have at least one entry");
            }
        }

        for(j = 0; j < n; j++){
            work[j] = j;
            for(p = Ap[j]; p < Ap[j+1]; p++){
                i = Ai[p];
                if(i > j){
                	throw new IllegalStateException("Entries on lower traingle");
                	}; //abort if entries on lower triangle
                while(work[i] != j){
                    if(etree[i] == QDLDL_UNKNOWN){
                        etree[i] = j;
                    }
                    Lnz[i]++;         //nonzeros in this column
                    work[i] = j;
                    i = etree[i];
                }
            }
        }

        //compute the total nonzeros in L.  This much
        //space is required to store Li and Lx.  Return
        //error code -2 if the nonzero count will overflow
        //its unteger type.
        sumLnz  = 0;
        for(i = 0; i < n; i++){
            if(sumLnz > QDLDL_INT_MAX - Lnz[i]){
                sumLnz = -2;
                break;
            }
            else{
                sumLnz += Lnz[i];
            }
        }

        return new Etree(Lnz, etree, sumLnz);

    }
    
    public static LDL decompose (CSCMatrix A, Etree e){
    	final int[] Lnz = e.Lnz;
    	final int[] etree = e.etree;
    	final int n = A.n;
    	final int m = A.m;
    	final int[] Ap = A.Ap;
    	final int[] Ai = A.Ai;
    	final double[] Ax = A.Ax;
    	final int[] Lp = new int[n+1];
    	final int[] Li = new int[e.sumLnz];
    	final double[] Lx = new double[e.sumLnz]; 
    	final double[] D = new double[n];
    	final double[] Dinv = new double[n];
    	//final boolean[] bwork = new boolean[n];
    	//final int[] iwork = new int[3*n];
    	//final double[] fwork = new double[n];
    	
    	int i,j,k,nnzY, bidx, cidx, nextIdx, nnzE, tmpIdx;
    	final int[] yIdx = new int[n]; 
    	final int[] elimBuffer = new int[n]; 
    	final int[] LNextSpaceInCol = new int[n];
    	final double[] yVals = new double[n];
    	double yVals_cidx;
    	final boolean[] yMarkers = new boolean[n];
    	int positiveValuesInD = 0;
    	
    	Lp[0] = 0; //first column starts at index zero

    	for(i = 0; i < n; i++){
	    	//compute L column indices
	    	Lp[i+1] = Lp[i] + Lnz[i];   //cumsum, total at the end
	
	    	// set all Yidx to be 'unused' initially
	    	//in each column of L, the next available space
	    	//to start is just the first space in the column
	    	yMarkers[i]  = QDLDL_UNUSED;
	    	yVals[i]     = 0.0;
	    	D[i]         = 0.0;
	    	LNextSpaceInCol[i] = Lp[i];
    	}

    	// First element of the diagonal D.
    	D[0]     = Ax[0];
    	if(D[0] == 0.0){
    		throw new IllegalStateException("First diagonal entry is zero.");
    		}
    	if(D[0]  > 0.0){positiveValuesInD++;}
    	Dinv[0] = 1/D[0];

    	//Start from 1 here. The upper LH corner is trivially 0
    	//in L b/c we are only computing the subdiagonal elements
    	for(k = 1; k < n; k++){

    	//NB : For each k, we compute a solution to
    	//y = L(0:(k-1),0:k-1))\b, where b is the kth
    	//column of A that sits above the diagonal.
    	//The solution y is then the kth row of L,
    	//with an implied '1' at the diagonal entry.

    	//number of nonzeros in this row of L
    	nnzY = 0;  //number of elements in this row

    	//This loop determines where nonzeros
    	//will go in the kth row of L, but doesn't
    	//compute the actual values
    	tmpIdx = Ap[k+1];

    	for(i = Ap[k]; i < tmpIdx; i++){

	    	bidx = Ai[i];   // we are working on this element of b
	
	    	//Initialize D[k] as the element of this column
	    	//corresponding to the diagonal place.  Don't use
	    	//this element as part of the elimination step
	    	//that computes the k^th row of L
	    	if(bidx == k){
	    		D[k] = Ax[i];
	    		continue;
	    	}
	
	    	yVals[bidx] = Ax[i];   // initialise y(bidx) = b(bidx)
	
	    	// use the forward elimination tree to figure
	    	// out which elements must be eliminated after
	    	// this element of b
	    	nextIdx = bidx;
	
	    	if(yMarkers[nextIdx] == QDLDL_UNUSED){   //this y term not already visited
	
	    	  yMarkers[nextIdx] = QDLDL_USED;     //I touched this one
	    	  elimBuffer[0]     = nextIdx;  // It goes at the start of the current list
	    	  nnzE              = 1;         //length of unvisited elimination path from here
	
	    	  nextIdx = etree[bidx];
	
	    	  while(nextIdx != QDLDL_UNKNOWN && nextIdx < k){
	    	    if(yMarkers[nextIdx] == QDLDL_USED) break;
	
	    	    yMarkers[nextIdx] = QDLDL_USED;   //I touched this one
	    	    elimBuffer[nnzE] = nextIdx; //It goes in the current list
	    	    nnzE++;                     //the list is one longer than before
	    	    nextIdx = etree[nextIdx];   //one step further along tree
	
	    	  } //end while
	
	    	  // now I put the buffered elimination list into
	    	  // my current ordering in reverse order
	    	  while(nnzE>0){
	    	    yIdx[nnzY++] = elimBuffer[--nnzE];
	    	  } //end while
	    	} //end if
	
	    	} //end for i
	
	    	//This for loop places nonzeros values in the k^th row
	    	for(i = (nnzY-1); i >=0; i--){
	
		    	//which column are we working on?
		    	cidx = yIdx[i];
		
		    	// loop along the elements in this
		    	// column of L and subtract to solve to y
		    	tmpIdx = LNextSpaceInCol[cidx];
		    	yVals_cidx = yVals[cidx];
		    	for(j = Lp[cidx]; j < tmpIdx; j++){
		    	  yVals[Li[j]] -= Lx[j]*yVals_cidx;
		    	}
		
		    	//Now I have the cidx^th element of y = L\b.
		    	//so compute the corresponding element of
		    	//this row of L and put it into the right place
		    	Li[tmpIdx] = k;
		    	Lx[tmpIdx] = yVals_cidx *Dinv[cidx];
		
		    	//D[k] -= yVals[cidx]*yVals[cidx]*Dinv[cidx];
		    	D[k] -= yVals_cidx*Lx[tmpIdx];
		    	LNextSpaceInCol[cidx]++;
		
		    	//reset the yvalues and indices back to zero and QDLDL_UNUSED
		    	//once I'm done with them
		    	yVals[cidx]     = 0.0;
		    	yMarkers[cidx]  = QDLDL_UNUSED;
	
	    	} //end for i
	
	    	//Maintain a count of the positive entries
	    	//in D.  If we hit a zero, we can't factor
	    	//this matrix, so abort
	    	if(D[k] == 0.0){
	    		throw new IllegalStateException("Diagonal entry zero.");
	    	}
	    	if(D[k]  > 0.0){
	    		positiveValuesInD++;
	    	}
	
	    	//compute the inverse of the diagonal
	    	Dinv[k]= 1/D[k];

    	} //end for k
    	return new LDL(new CSCMatrix(n,m,Lx.length,Lp,Li,Lx), D, Dinv, positiveValuesInD);
    }
    
    static void Lsolve(int n, int[] Lp, int[] Li, double[] Lx, double[] x) {
		int i,j;
		for(i = 0; i < n; i++){
		    for(j = Lp[i]; j < Lp[i+1]; j++){
		        x[Li[j]] -= Lx[j]*x[i];
		    }
		}
    }
    
    static void Ltsolve(int n, int[] Lp, int[] Li, double[] Lx, double[] x) {
    	int i,j;
    	  for(i = n-1; i>=0; i--){
    	      for(j = Lp[i]; j < Lp[i+1]; j++){
    	          x[i] -= Lx[j]*x[Li[j]];
    	      }
    	  }
    }
    
    public static double[] solve(LDL ldl, double[] b) {
    	final double x[] = b.clone();
    	final int n = ldl.L.n;
    	final int[] Lp = ldl.L.Ap;
    	final int[] Li = ldl.L.Ai;
    	final double[] Lx = ldl.L.Ax;
    	int i;
    	Lsolve(n,Lp,Li,Lx,x);
    	for(i = 0; i < n; i++) x[i] *= ldl.Dinv[i];
    	Ltsolve(n,Lp,Li,Lx,x);
    	return x;
    }
    
    public static class Etree {
    	public final int[] Lnz;
    	public final int[] etree;
    	public final int sumLnz;
		public Etree(int[] Lnz, int[] etree, int sumLnz) {
			this.Lnz = Lnz;
			this.etree = etree;
			this.sumLnz = sumLnz;
		}
    }
    
    public static class LDL {
    	public final CSCMatrix L;
    	public final double[] D;
    	public final double[] Dinv;
    	public final int positiveValuesInD;
		public LDL(CSCMatrix L, double[] D, double[] Dinv, int positiveValuesInD) {
			this.L = L;
			this.D = D;
			this.Dinv = Dinv;
			this.positiveValuesInD = positiveValuesInD;
		}
    }
    
    public static void main(String... args) {
    	int An   = 10;
    	int[] Ap = {0, 1, 2, 4, 5, 6, 8, 10, 12, 14, 17};
    	int[] Ai = {0, 1, 1, 2, 3, 4, 1, 5, 0, 6, 3, 7, 6, 8, 1, 2, 9};
    	double[] Ax = {1.0, 0.460641, -0.121189, 0.417928, 0.177828, 0.1,
    	                       -0.0290058, -1.0, 0.350321, -0.441092, -0.0845395,
    	                       -0.316228, 0.178663, -0.299077, 0.182452, -1.56506, -0.1};
    	double[] b = {1,2,3,4,5,6,7,8,9,10};
    	QDLDL q = new QDLDL();
    	CSCMatrix A = new CSCMatrix(An, An,Ax.length,Ap,Ai,Ax);
    	q.factor(A);
    	double[] x = q.solve(b);
    	
    	Etree e = computeEtree(An,Ap,Ai);
//    	CSCMatrix A = new CSCMatrix(An,An,Ax.length,Ap,Ai,Ax);
    	LDL ldl = q._ldl;
    	System.out.println();
    	System.out.println(Arrays.toString(Ap));
    	System.out.println(Arrays.toString(Ai));
    	System.out.println(Arrays.toString(Ax));
    	System.out.println();
    	System.out.println(Arrays.toString(e.etree));
    	System.out.println(Arrays.toString(e.Lnz));
    	System.out.println();
    	System.out.println(Arrays.toString(ldl.L.Ap));
    	System.out.println(Arrays.toString(ldl.L.Ai));
    	System.out.println(Arrays.toString(ldl.L.Ax));
    	System.out.println();
    	System.out.println(Arrays.toString(ldl.D));
    	System.out.println(Arrays.toString(ldl.Dinv));
    	System.out.println();
    	System.out.println(Arrays.toString(x));
    	
    }
}
