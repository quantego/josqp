package com.quantego.josqp;

import java.util.Arrays;

/**
 * Java implementation of the QDLDL factorization routine by Bartolomeo Stellato, Goran Banjac, Paul Goulart, Alberto Bemporad, and
 * Stephen Boyd. The implementation closely follows the original C code that can be found at https://github.com/osqp/qdldl.
 *
 *  @author Paul Goulart
 *  @author Bartolomeo Stellato
 *  @author Goran Banjac
 *  @author Nils Loehndorf
 *
 */
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
    
    private Workspace workspace;
    boolean useSchedule = Boolean.getBoolean("josqp.factorSchedule");
    boolean useBlocks = Boolean.getBoolean("josqp.factorBlocks");

    /** A fresh factorization; previously produced factors remain independent. */
    public boolean factor(CSCMatrix A) {
        Workspace next = new Workspace(A, computeEtree(A.n, A.Ap, A.Ai));
        LDL result = decompose(A, next, useSchedule, useBlocks);
        workspace = next;
        _ldl = result;
        return true;
    }

    /** Internal solver path: reuse analysis and arrays when the CSC pattern is unchanged.
     * A failed in-place numeric factorization invalidates the factor until a successful retry.
     * QDLDL instances, like their owning LinSys, are not thread-safe.
     */
    boolean refactor(CSCMatrix A) {
        _ldl = null;
        if (workspace == null || !workspace.matches(A)) return factor(A);
        _ldl = decompose(A, workspace, useSchedule, useBlocks);
        return true;
    }

    public double[] solve(double[] b) {
        if (_ldl == null) throw new IllegalStateException("No valid LDL factorization");
        return QDLDL.solve(_ldl, b);
    }

    private static final class Workspace {
        final Etree analysis;
        final int[] patternAp, patternAi;
        final CSCMatrix L;
        final double[] D, Dinv, yVals;
        final int[] yIdx, elimBuffer, nextSpace;
        final boolean[] yMarkers;
        int[] scheduleColumns, schedulePositions, scheduleRows;
        Blocks[] blocks;
        boolean completed;

        Workspace(CSCMatrix A, Etree analysis) {
            this.analysis = analysis;
            int n = A.n;
            patternAp = Arrays.copyOf(A.Ap, n + 1);
            patternAi = Arrays.copyOf(A.Ai, A.Ap[n]);
            int[] lp = new int[n + 1];
            for (int i = 0; i < n; i++) lp[i + 1] = lp[i] + analysis.Lnz[i];
            L = new CSCMatrix(A.m, n, analysis.sumLnz, lp,
                    new int[analysis.sumLnz], new double[analysis.sumLnz]);
            D = new double[n]; Dinv = new double[n]; yVals = new double[n];
            yIdx = new int[n]; elimBuffer = new int[n]; nextSpace = new int[n];
            yMarkers = new boolean[n];
        }

        boolean matches(CSCMatrix A) {
            if (A.n != L.n || A.m != L.m || A.Ap.length < patternAp.length
                    || A.Ai.length < patternAi.length || A.Ax.length < patternAi.length) return false;
            for (int i = 0; i < patternAp.length; i++) if (A.Ap[i] != patternAp[i]) return false;
            for (int i = 0; i < patternAi.length; i++) if (A.Ai[i] != patternAi[i]) return false;
            return true;
        }
    }

    /** Long contiguous runs only; short sparse columns retain the indexed loop. */
    private static final class Blocks {
        final int[] starts, rows, lengths;
        Blocks(int[] li, int start, int end, int count) {
            starts=new int[count]; rows=new int[count]; lengths=new int[count];
            int b=0;
            for(int i=start;i<end;) {
                int j=i+1; while(j<end && li[j]==li[j-1]+1) j++;
                starts[b]=i; rows[b]=li[i]; lengths[b++]=j-i; i=j;
            }
        }
    }

    private static void buildBlocks(Workspace w) {
        w.blocks=new Blocks[w.L.n];
        for(int c=0;c<w.L.n;c++) {
            int start=w.L.Ap[c], end=w.L.Ap[c+1], count=0;
            long squares=0;
            for(int i=start;i<end;) {
                int j=i+1; while(j<end && w.L.Ai[j]==w.L.Ai[j-1]+1) j++;
                long len=j-i; squares+=len*len; count++; i=j;
            }
            if(end>start && squares/(double)(end-start)>=16)
                w.blocks[c]=new Blocks(w.L.Ai,start,end,count);
        }
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
                throw new IllegalStateException(String.format("There are more than %d non-zero elements, which exceeds the capacity of LDL factorization.",Integer.MAX_VALUE));
            }
            else{
                sumLnz += Lnz[i];
            }
        }

        return new Etree(Lnz, etree, sumLnz);

    }
    
    public static LDL decompose(CSCMatrix A, Etree e) {
        return decompose(A, new Workspace(A, e));
    }

    private static LDL decompose(CSCMatrix A, Workspace w) {
        return decompose(A,w,false,false);
    }

    private static LDL decompose(CSCMatrix A, Workspace w, boolean schedule, boolean blocks) {
        final int[] etree = w.analysis.etree;
        final int n = A.n;
        final int[] Ap = A.Ap, Ai = A.Ai;
        final double[] Ax = A.Ax;
        final int[] Lp = w.L.Ap, Li = w.L.Ai;
        final double[] Lx = w.L.Ax, D = w.D, Dinv = w.Dinv;
        int i,j,k,nnzY, bidx, cidx, nextIdx, nnzE, tmpIdx;
        final int[] yIdx = w.yIdx, elimBuffer = w.elimBuffer, LNextSpaceInCol = w.nextSpace;
        final double[] yVals = w.yVals;
        double yVals_cidx;
        final boolean[] yMarkers = w.yMarkers;
        boolean replay = schedule && w.scheduleRows != null;
        boolean record = schedule && w.completed && !replay;
        int[] recordedColumns = record ? new int[w.analysis.sumLnz] : null;
        int[] recordedPositions = record ? new int[w.analysis.sumLnz] : null;
        int[] recordedRows = record ? new int[n+1] : null;
        int cursor=0;
    	int positiveValuesInD = 0;
    	

    	for(i = 0; i < n; i++){
	
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
            if(record) recordedRows[k]=cursor;

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
            if(replay) continue;
	
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
        int steps = replay ? w.scheduleRows[k+1]-w.scheduleRows[k] : nnzY;
        for(i = 0; i < steps; i++){
	
		    	//which column are we working on?
                cidx = replay ? w.scheduleColumns[w.scheduleRows[k]+i] : yIdx[nnzY-1-i];
		
		    	// loop along the elements in this
		    	// column of L and subtract to solve to y
                tmpIdx = replay ? w.schedulePositions[w.scheduleRows[k]+i] : LNextSpaceInCol[cidx];
                if(record) {
                    if(Li[tmpIdx]!=k) throw new IllegalStateException("LDL schedule disagrees with factor pattern");
                    recordedColumns[cursor]=cidx; recordedPositions[cursor++]=tmpIdx;
                }
		    	yVals_cidx = yVals[cidx];
                Blocks runs=blocks && w.blocks!=null ? w.blocks[cidx] : null;
                if(runs==null) {
                    for(j=Lp[cidx];j<tmpIdx;j++) yVals[Li[j]]-=Lx[j]*yVals_cidx;
                } else {
                    for(int b=0;b<runs.starts.length && runs.starts[b]<tmpIdx;b++) {
                        int len=Math.min(runs.lengths[b],tmpIdx-runs.starts[b]);
                        NumericKernels.INSTANCE.subtractSegment(yVals,runs.rows[b],Lx,runs.starts[b],len,yVals_cidx);
                    }
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
        if(record) {
            recordedRows[n]=cursor;
            if(cursor!=w.analysis.sumLnz) throw new IllegalStateException("Incomplete LDL schedule");
            w.scheduleColumns=recordedColumns; w.schedulePositions=recordedPositions; w.scheduleRows=recordedRows;
        }
        if(blocks && w.blocks==null) buildBlocks(w);
        w.completed=true;
        return new QDLDL.LDL(w.L, D, Dinv, positiveValuesInD);
    }
    
    static void Lsolve(int n, int[] Lp, int[] Li, double[] Lx, double[] x) {
		int i,j;
		for(i = 0; i < n; i++){
		    final double xi = x[i]; // L is strictly lower triangular: this entry cannot change here.
		    for(j = Lp[i]; j < Lp[i+1]; j++){
		        x[Li[j]] -= Lx[j]*xi;
		    }
		}
    }
    
    static void Ltsolve(int n, int[] Lp, int[] Li, double[] Lx, double[] x) {
    	int i,j;
    	  for(i = n-1; i>=0; i--){
              double xi = x[i];
    	      for(j = Lp[i]; j < Lp[i+1]; j++){
		          xi -= Lx[j]*x[Li[j]];
    	      }
              x[i] = xi;
    	  }
    }

    /**
     * Backward solve with the diagonal scaling folded in. L is strictly lower
     * triangular, so every entry Li[j] has already been processed (and scaled)
     * when column i is reached in descending order. This is the same arithmetic
     * as scaling all of x first, with one pass over x less.
     */
    static void LtsolveScaled(int n, int[] Lp, int[] Li, double[] Lx, double[] Dinv, double[] x) {
    	int i,j;
    	for(i = n-1; i>=0; i--){
            double xi = x[i]*Dinv[i];
    		for(j = Lp[i]; j < Lp[i+1]; j++){
    			xi -= Lx[j]*x[Li[j]];
    		}
            x[i] = xi;
    	}
    }
    
    /** Short-column kernels preserve the generic solve's exact update order.
     * The fallback also remains correct if a caller changes a public factor's pattern.
     */
    private static void solveShort(QDLDL.LDL l,double[] x) {
        int[] p=l.L.Ap,r=l.L.Ai;double[] a=l.L.Ax;
        for(int i=0;i<l.L.n;i++) {
            int j=p[i],len=p[i+1]-j;double xi=x[i];
            switch(len) {
                case 0:break;
                case 1:x[r[j]]-=a[j]*xi;break;
                case 2:x[r[j]]-=a[j]*xi;x[r[j+1]]-=a[j+1]*xi;break;
                case 3:x[r[j]]-=a[j]*xi;x[r[j+1]]-=a[j+1]*xi;x[r[j+2]]-=a[j+2]*xi;break;
                case 4:x[r[j]]-=a[j]*xi;x[r[j+1]]-=a[j+1]*xi;x[r[j+2]]-=a[j+2]*xi;x[r[j+3]]-=a[j+3]*xi;break;
                default:for(int end=p[i+1];j<end;j++) x[r[j]]-=a[j]*xi;
            }
        }
        for(int i=l.L.n-1;i>=0;i--) {
            int j=p[i],len=p[i+1]-j;double xi=x[i]*l.Dinv[i];
            switch(len) {
                case 0:break;
                case 1:xi-=a[j]*x[r[j]];break;
                case 2:xi-=a[j]*x[r[j]];xi-=a[j+1]*x[r[j+1]];break;
                case 3:xi-=a[j]*x[r[j]];xi-=a[j+1]*x[r[j+1]];xi-=a[j+2]*x[r[j+2]];break;
                case 4:xi-=a[j]*x[r[j]];xi-=a[j+1]*x[r[j+1]];xi-=a[j+2]*x[r[j+2]];xi-=a[j+3]*x[r[j+3]];break;
                default:for(int end=p[i+1];j<end;j++) xi-=a[j]*x[r[j]];
            }
            x[i]=xi;
        }
    }
    public static double[] solve(LDL ldl, double[] b) {
    	final double x[] = b.clone();
        solveInPlace(ldl, x);
        return x;
    }

    /** Solve using caller-owned workspace; no per-solve array allocation. */
    public static void solveInPlace(LDL ldl, double[] x) {
        if (ldl == null) throw new IllegalStateException("No valid LDL factorization");
    	final int n = ldl.L.n;
        if (x.length != n) {
            throw new IllegalArgumentException("Right-hand side length must match factorization");
        }
    	final int[] Lp = ldl.L.Ap;
    	final int[] Li = ldl.L.Ai;
    	final double[] Lx = ldl.L.Ax;
        if (ldl.shortColumns) { solveShort(ldl,x); return; }
    	Lsolve(n,Lp,Li,Lx,x);
    	LtsolveScaled(n,Lp,Li,Lx,ldl.Dinv,x);
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
        /** Select once per factor, without an extra plan array or per-solve allocation. */
        final boolean shortColumns;
		public LDL(CSCMatrix L, double[] D, double[] Dinv, int positiveValuesInD) {
			this.L = L;
			this.D = D;
			this.Dinv = Dinv;
			this.positiveValuesInD = positiveValuesInD;
            boolean shortPattern=true;
            for(int i=0;i<L.n;i++) if(L.Ap[i+1]-L.Ap[i]>4) {shortPattern=false;break;}
            this.shortColumns=shortPattern;
		}
    }
    
    /*public static void main(String... args) {
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
    	
    }*/
}
