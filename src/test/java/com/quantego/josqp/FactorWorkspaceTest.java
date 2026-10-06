package com.quantego.josqp;

import java.util.Random;
import org.junit.Test;
import static org.junit.Assert.*;

public class FactorWorkspaceTest {
    private CSCMatrix matrix() {
        return new CSCMatrix(4, 4, 8, new int[]{0,1,3,5,8},
                new int[]{0,0,1,0,2,1,2,3}, new double[]{8,1,9,-2,10,2,1,11});
    }

    private void bits(double[] expected, double[] actual) {
        assertEquals(expected.length, actual.length);
        for (int i = 0; i < expected.length; i++)
            assertEquals("entry " + i, Double.doubleToLongBits(expected[i]), Double.doubleToLongBits(actual[i]));
    }

    private void check(CSCMatrix a, QDLDL q) {
        QDLDL.LDL fresh = QDLDL.decompose(a, QDLDL.computeEtree(a.n, a.Ap, a.Ai));
        assertArrayEquals(fresh.L.Ap, q._ldl.L.Ap);
        assertArrayEquals(fresh.L.Ai, q._ldl.L.Ai);
        bits(fresh.L.Ax, q._ldl.L.Ax);
        bits(fresh.D, q._ldl.D); bits(fresh.Dinv, q._ldl.Dinv);
        assertEquals(fresh.positiveValuesInD, q._ldl.positiveValuesInD);
        double[] rhs = new double[a.n];
        for (int i = 0; i < rhs.length; i++) rhs[i] = i + 1;
        double[] x = q.solve(rhs);
        bits(QDLDL.solve(fresh, rhs), x);
        double[] ax = new double[a.n];
        for (int c = 0; c < a.n; c++) for (int k = a.Ap[c]; k < a.Ap[c+1]; k++) {
            int r = a.Ai[k];
            ax[r] += a.Ax[k] * x[c];
            if (r != c) ax[c] += a.Ax[k] * x[r];
        }
        assertArrayEquals(rhs, ax, 1e-12);
    }

    @Test public void repeatedNumericUpdatesReuseFactorsAndMatchFreshDecomposition() {
        CSCMatrix a = matrix(); QDLDL q = new QDLDL(); q.refactor(a);
        QDLDL.LDL first = q._ldl;
        Random random = new Random(529);
        for (int step = 0; step < 100; step++) {
            for (int c = 0; c < a.n; c++) for (int k = a.Ap[c]; k < a.Ap[c+1]; k++)
                a.Ax[k] = a.Ai[k] == c ? 10 + random.nextDouble() : random.nextDouble() - 0.5;
            q.refactor(a);
            assertSame(first.L, q._ldl.L);
            assertSame(first.D, q._ldl.D); assertSame(first.Dinv, q._ldl.Dinv);
            check(a, q);
        }
        a.Ax[0] = -10;
        q.refactor(a); check(a, q);
        assertEquals(3, q._ldl.positiveValuesInD);
    }

    @Test public void inPlacePatternMutationAndDimensionChangeRebuildAnalysis() {
        CSCMatrix a = matrix(); QDLDL q = new QDLDL(); q.refactor(a);
        CSCMatrix old = q._ldl.L;
        a.Ai[3] = 1; // Same arrays and nnz, different elimination pattern.
        q.refactor(a); assertNotSame(old, q._ldl.L); check(a, q);
        old = q._ldl.L;
        CSCMatrix two = new CSCMatrix(2,2,3,new int[]{0,1,3},new int[]{0,0,1},new double[]{4,1,3});
        q.refactor(two); assertNotSame(old,q._ldl.L); check(two,q);
    }

    @Test public void failedNumericReuseIsInvalidUntilSuccessfulRecovery() {
        CSCMatrix a = matrix(); QDLDL q = new QDLDL(); q.refactor(a);
        a.Ax[0] = 1; a.Ax[1] = 1; a.Ax[2] = 1; // Singular second pivot, after partial writes.
        assertThrows(IllegalStateException.class, () -> q.refactor(a));
        assertNull(q._ldl);
        assertThrows(IllegalStateException.class, () -> q.solve(new double[4]));
        a.Ax = matrix().Ax;
        q.refactor(a); check(a,q);
    }

    @Test public void publicFactorResultsRemainIndependentAndSurviveRejectedReplacement() {
        CSCMatrix a = matrix(); QDLDL q = new QDLDL(); q.factor(a);
        QDLDL.LDL first = q._ldl;
        double[] firstD = first.D.clone();
        a.Ax[0] = 12; q.factor(a);
        assertNotSame(first.L, q._ldl.L); bits(firstD,first.D);
        QDLDL.LDL good = q._ldl;
        a.Ax[0] = 0;
        assertThrows(IllegalStateException.class, () -> q.factor(a));
        assertSame(good,q._ldl);
    }
}
