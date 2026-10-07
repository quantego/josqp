package com.quantego.josqp;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

public class FactorScheduleTest {
    private void same(QDLDL.LDL a,QDLDL.LDL b) {
        assertArrayEquals(a.L.Ap,b.L.Ap); assertArrayEquals(a.L.Ai,b.L.Ai);
        for(double[][] pair:new double[][][]{ {a.L.Ax,b.L.Ax},{a.D,b.D},{a.Dinv,b.Dinv} })
            for(int i=0;i<pair[0].length;i++) assertEquals(Double.doubleToLongBits(pair[0][i]),Double.doubleToLongBits(pair[1][i]));
        assertEquals(a.positiveValuesInD,b.positiveValuesInD);
    }

    @Test public void scheduledAndBlockedFactorsMatchIndependentFactorsOnEveryUpdate() {
        int n=80,nz=n*(n+1)/2,k=0;
        int[] ap=new int[n+1],ai=new int[nz]; double[] ax=new double[nz];
        for(int c=0;c<n;c++) {
            ap[c]=k;
            for(int r=0;r<=c;r++) { ai[k]=r; ax[k++]=r==c?n+1:0.01*((r+c)%7-3); }
        }
        ap[n]=k; CSCMatrix a=new CSCMatrix(n,n,nz,ap,ai,ax);
        for(boolean blocks:new boolean[]{false,true}) {
            QDLDL q=new QDLDL();q.useSchedule=true;q.useBlocks=blocks;
            for(int step=0;step<12;step++) {
                for(int c=0;c<n;c++) ax[ap[c+1]-1]=n+1+step*0.1;
                q.refactor(a);same(QDLDL.decompose(a,QDLDL.computeEtree(n,ap,ai)),q._ldl);
            }
        }
    }

    @Test public void scheduleRebuildsForChangedPatternAndRecoversAfterFailedValues() {
        CSCMatrix a=new CSCMatrix(3,3,5,new int[]{0,1,3,5},new int[]{0,0,1,1,2},new double[]{4,1,4,1,4});
        QDLDL q=new QDLDL();q.useSchedule=q.useBlocks=true;
        q.refactor(a);q.refactor(a);q.refactor(a);
        a.Ax[0]=0;
        assertThrows(IllegalStateException.class,()->q.refactor(a));assertNull(q._ldl);
        a.Ax[0]=4;q.refactor(a);same(QDLDL.decompose(a,QDLDL.computeEtree(3,a.Ap,a.Ai)),q._ldl);
        a.Ai[3]=0;q.refactor(a);q.refactor(a);q.refactor(a);
        same(QDLDL.decompose(a,QDLDL.computeEtree(3,a.Ap,a.Ai)),q._ldl);
    }
}
