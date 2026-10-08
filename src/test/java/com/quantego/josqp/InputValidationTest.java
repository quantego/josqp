package com.quantego.josqp;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

class InputValidationTest {
    private static CSCMatrix matrix(double p) {
        return new CSCMatrix(1,1,1,new int[]{0,1},new int[]{0},new double[]{p});
    }
    private static OSQP.Data data(double lower, double upper) {
        return new OSQP.Data(1,1,matrix(8),matrix(1),new double[]{-8},
                new double[]{lower},new double[]{upper},0);
    }
    private static OSQP.Settings settings() {
        OSQP.Settings s = new OSQP.Settings(); s.verbose=false;
        s.eps_abs=s.eps_rel=1e-9; return s;
    }
    @Test void invalidDataIsRejectedBeforeSolving() {
        assertThrows(IllegalArgumentException.class, () -> data(2,1));
        assertThrows(IllegalArgumentException.class, () -> data(Double.NaN,1));
        assertThrows(IllegalArgumentException.class, () -> data(0,Double.NEGATIVE_INFINITY));
        OSQP.Data d = data(0,10); d.q[0]=Double.NaN;
        assertThrows(IllegalArgumentException.class, () -> new OSQP(d,settings()));
        d.q[0]=-8; d.P.Ax[0]=Double.POSITIVE_INFINITY;
        assertThrows(IllegalArgumentException.class, () -> new OSQP(d,settings()));
        assertThrows(IllegalArgumentException.class, () -> new OSQP(null,settings()));
    }
    @Test void malformedCSCAndVectorLengthsAreRejected() {
        OSQP.Data d = data(0,10);
        d.P.Ap[1]=2; assertThrows(IllegalArgumentException.class, () -> new OSQP(d,settings()));
        d.P.Ap[1]=1; d.P.Ai[0]=-1;
        assertThrows(IllegalArgumentException.class, () -> new OSQP(d,settings()));
        d.P.Ai[0]=0; d.A.Ap[0]=1;
        assertThrows(IllegalArgumentException.class, () -> new OSQP(d,settings()));
        assertThrows(IllegalArgumentException.class, () -> new OSQP.Data(1,1,matrix(8),matrix(1),
                new double[0],new double[]{0},new double[]{10},0));
    }
    @Test void nonFiniteSettingsAreRejected() {
        for (double value : new double[]{Double.NaN,Double.POSITIVE_INFINITY,0,-1}) {
            OSQP.Settings s = settings(); s.rho=value;
            assertThrows(IllegalArgumentException.class, () -> new OSQP(data(0,10),s));
        }
        OSQP.Settings s=settings(); s.eps_abs=Double.NaN; OSQP.Settings invalid=s;
        assertThrows(IllegalArgumentException.class, () -> new OSQP(data(0,10),invalid));
        s=settings(); s.device=1; OSQP.Settings device=s;
        assertThrows(IllegalArgumentException.class, () -> new OSQP(data(0,10),device));
    }
    @Test void dataAndSettingsCanBeReusedWithoutBeingScaledOrMutated() {
        OSQP.Data d=data(0,10); OSQP.Settings s=settings();
        OSQP first=new OSQP(d,s); first.solve();
        OSQP second=new OSQP(d,s); second.solve();
        assertArrayEquals(new double[]{8},d.P.Ax,0);
        assertArrayEquals(new double[]{-8},d.q,0);
        assertArrayEquals(first.getPrimalSolution(),second.getPrimalSolution(),0);
        assertEquals(1,second.getPrimalSolution()[0],1e-7);
        s.max_iter=0; d.q[0]=100;
        assertEquals(OSQP.Status.SOLVED,first.solve());
        assertEquals(1,first.getPrimalSolution()[0],1e-7);
    }
    @Test void rejectedCombinedVectorUpdateLeavesTheOldProblemIntact() {
        OSQP q=new OSQP(data(0,10),settings()); q.solve();
        double[] old=q.work.data.q.clone();
        assertThrows(IllegalArgumentException.class, () -> q.update_data_vec(new double[]{-80},
                new double[]{20},new double[]{10}));
        assertArrayEquals(old,q.work.data.q,0);
        assertThrows(IllegalArgumentException.class, () -> q.update_data_vec(new double[]{Double.NaN},null,null));
        assertThrows(IllegalArgumentException.class, () -> q.update_P(new double[]{Double.NaN},null));
        assertThrows(IllegalArgumentException.class, () -> q.warm_start_x(new double[]{Double.NaN}));
        assertEquals(OSQP.Status.SOLVED,q.solve());
        assertEquals(1,q.getPrimalSolution()[0],1e-7);
    }
    @Test void failedFactorUpdateCannotProduceASolvedResultAndCanBeRepaired() {
        OSQP q=new OSQP(data(0,10),settings()); q.solve();
        assertThrows(IllegalStateException.class, () -> q.update_P(new double[]{-8},null));
        assertEquals(OSQP.Status.ERROR,q.getInfo().status);
        assertThrows(IllegalStateException.class,q::solve);
        q.update_P(new double[]{8},null);
        assertEquals(OSQP.Status.SOLVED,q.solve());
        assertEquals(1,q.getPrimalSolution()[0],1e-7);
    }
    @Test void unboundedRowsAndEmptyConstraintMatricesAreAllowed() {
        assertEquals(OSQP.Status.SOLVED,new OSQP(data(Double.NEGATIVE_INFINITY,Double.POSITIVE_INFINITY),settings()).solve());
        OSQP.Data d=new OSQP.Data(1,0,matrix(8),new CSCMatrix(0,1,0,new int[]{0,0},new int[0],new double[0]),
                new double[]{-8},new double[0],new double[0],0);
        OSQP q=new OSQP(d,settings()); assertEquals(OSQP.Status.SOLVED,q.solve());
        assertEquals(1,q.getPrimalSolution()[0],1e-7);
    }
    @Test void automaticRhoIntervalUsesTheOsqpIterationCadence() {
        OSQP.Settings auto=settings(); auto.adaptive_rho_interval=0; auto.check_termination=25;
        auto.rho=1e-6; auto.scaling=0; auto.max_iter=10000;
        OSQP.Settings fixed=auto.copy(); fixed.adaptive_rho_interval=100;
        OSQP q=new OSQP(data(0,0.01),auto), f=new OSQP(data(0,0.01),fixed);
        q.solve(); f.solve();
        assertArrayEquals(f.getPrimalSolution(),q.getPrimalSolution(),0);
        assertEquals(f.getInfo().rho_updates,q.getInfo().rho_updates);
        assertTrue(q.getInfo().rho_updates>0, "The comparison must exercise rho adaptation.");
        assertEquals(0,auto.adaptive_rho_interval);
    }
}
