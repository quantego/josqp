package com.quantego.josqp;

import java.util.Random;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.assertArrayEquals;

/** Independent dense products plus exact comparison to the original arithmetic order. */
public final class ShortColumnSolveTest {
    static void check(int n,int width,Random rng) {
        int[] p=new int[n+1];for(int c=0;c<n;c++) p[c+1]=p[c]+Math.min(width,n-c-1);
        int[] rows=new int[p[n]];double[] values=new double[p[n]],d=new double[n],inv=new double[n];
        double[][] lower=new double[n][n];
        for(int c=0;c<n;c++) {
            lower[c][c]=1;d[c]=(c%2==0?1:-1)*(1+rng.nextDouble());inv[c]=1/d[c];
            for(int k=p[c];k<p[c+1];k++) {int r=c+1+k-p[c];rows[k]=r;values[k]=(rng.nextDouble()-.5)*.2;lower[r][c]=values[k];}
        }
        double[][] kkt=new double[n][n];
        for(int i=0;i<n;i++) for(int j=0;j<n;j++) for(int k=0;k<n;k++) kkt[i][j]+=lower[i][k]*d[k]*lower[j][k];
        double[] known=new double[n],rhs=new double[n];for(int i=0;i<n;i++) known[i]=rng.nextGaussian();
        for(int i=0;i<n;i++) for(int j=0;j<n;j++) rhs[i]+=kkt[i][j]*known[j];
        double[] reference=rhs.clone(),actual=rhs.clone();
        for(int c=0;c<n;c++) for(int k=p[c];k<p[c+1];k++) reference[rows[k]]-=values[k]*reference[c];
        for(int c=n-1;c>=0;c--) {double x=reference[c]*inv[c];for(int k=p[c];k<p[c+1];k++)x-=values[k]*reference[rows[k]];reference[c]=x;}
        QDLDL.LDL factor=new QDLDL.LDL(new CSCMatrix(n,n,values.length,p,rows,values),d,inv,(n+1)/2);
        QDLDL.solveInPlace(factor,actual);
        for(int i=0;i<n;i++) {
            if(Double.doubleToRawLongBits(reference[i])!=Double.doubleToRawLongBits(actual[i])) throw new AssertionError("arithmetic order differs");
            if(Math.abs(actual[i]-known[i])>1e-12) throw new AssertionError("known dense solution differs");
            double product=0;for(int j=0;j<n;j++) product+=kkt[i][j]*actual[j];
            if(Math.abs(product-rhs[i])>1e-12) throw new AssertionError("dense KKT residual differs");
        }
    }
    @Test public void shortAndGeneralPatternsPreserveArithmeticAndSolveIndependentDenseSystems() {
        Random r=new Random(82916);int count=0;
        for(int n:new int[]{0,1,2,7,31,73}) for(int width:new int[]{0,1,2,3,4,7}) for(int sample=0;sample<10;sample++) {check(n,width,r);count++;}
        if(count!=360) throw new AssertionError("Missing factor cases");
    }

    @Test public void changedPublicFactorPatternKeepsGenericFallbackCorrect() {
        int n=7;
        CSCMatrix lower=new CSCMatrix(n,n,0,new int[n+1],new int[0],new double[0]);
        double[] d=new double[n],inv=new double[n];
        java.util.Arrays.fill(d,1);java.util.Arrays.fill(inv,1);
        QDLDL.LDL factor=new QDLDL.LDL(lower,d,inv,n);
        // The public CSC content can change after the short-pattern decision was cached.
        int[] p=new int[n+1],rows=new int[n-1];double[] values=new double[n-1];
        for(int i=1;i<=n;i++) p[i]=n-1;
        for(int i=0;i<n-1;i++) {rows[i]=i+1;values[i]=.02*(i+1);}
        lower.replaceContent(n-1,p,rows,values);
        double[] known={1,2,3,4,5,6,7},rhs=new double[n];
        // For L=I+v*e0', K=L*L': b0=x0+v'*x; bi=xi+vi*b0.
        rhs[0]=known[0];
        for(int i=1;i<n;i++) rhs[0]+=values[i-1]*known[i];
        for(int i=1;i<n;i++) rhs[i]=known[i]+values[i-1]*rhs[0];
        QDLDL.solveInPlace(factor,rhs);
        assertArrayEquals(known,rhs,1e-12);
    }
}
