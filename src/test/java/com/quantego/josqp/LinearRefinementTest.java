package com.quantego.josqp;

import org.junit.Test;
import static org.junit.Assert.*;

public class LinearRefinementTest {
    private LinSys system() {
        CSCMatrix p=new CSCMatrix(2,2,3,new int[]{0,1,3},new int[]{0,0,1},new double[]{4,1,3});
        CSCMatrix a=new CSCMatrix(1,2,2,new int[]{0,1,2},new int[]{0,0},new double[]{1,2});
        return new LinSys(p,a,1e-6,new double[]{0.1},false);
    }

    @Test public void diagnosticPreservesSolutionBitsAndMeasuresOriginalKkt() {
        LinSys plain=system(), audit=system(); audit.configureRefinement(0,true);
        double[] rhs={1,2,3}, expected=new double[3], actual=new double[3];
        plain.LDLSolve(expected,rhs,plain.P,plain.bp);
        audit.LDLSolve(actual,rhs,audit.P,audit.bp);
        for(int i=0;i<3;i++) assertEquals(Double.doubleToLongBits(expected[i]),Double.doubleToLongBits(actual[i]));
        assertEquals(1,audit.linearCalls); assertEquals(0,audit.linearCorrections);
        assertTrue(audit.maxBackwardError<1e-14);
        assertEquals(1,(4+1e-6)*actual[0]+actual[1]+actual[2],1e-14);
        assertEquals(2,actual[0]+(3+1e-6)*actual[1]+2*actual[2],1e-14);
        assertEquals(3,actual[0]+2*actual[1]-10*actual[2],1e-14);
    }

    @Test public void residualControlledCorrectionRecoversFromInexactFactors() {
        LinSys exact=system(), inexact=system(); inexact.configureRefinement(6,true);
        // Controlled factor error: corrections must use residuals of K, not of LDL.
        for(int i=0;i<inexact.qdldl._ldl.Dinv.length;i++) inexact.qdldl._ldl.Dinv[i]*=1.001;
        double[] rhs={1,2,3}, expected=new double[3], actual=new double[3];
        exact.LDLSolve(expected,rhs,exact.P,exact.bp);
        inexact.LDLSolve(actual,rhs,inexact.P,inexact.bp);
        assertTrue(inexact.maxBackwardError>1e-5);
        assertTrue(inexact.linearCorrections>0);
        assertArrayEquals(expected,actual,1e-11);
        assertTrue(inexact.backwardError(inexact.bp)<1e-12);
    }

    @Test public void zeroRhsDoesNotDivideByZeroOrRequestCorrections() {
        LinSys system=system(); system.configureRefinement(3,true);
        system.LDLSolve(new double[3],new double[3],system.P,system.bp);
        assertEquals(0,system.maxBackwardError,0); assertEquals(0,system.linearCorrections);
        assertThrows(IllegalArgumentException.class,()->system.configureRefinement(-1,true));
    }
}
