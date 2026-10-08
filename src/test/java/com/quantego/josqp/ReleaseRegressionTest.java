package com.quantego.josqp;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

/** Regression cases from the 1.0 release review, with independent analytical optima. */
class ReleaseRegressionTest {
    private static OSQP.Settings tight() {
        OSQP.Settings s = new OSQP.Settings();
        s.verbose = false;
        s.eps_abs = s.eps_rel = 1e-9;
        s.max_iter = 20000;
        return s;
    }
    private static Model quadratic() {
        Model m = new Model();
        m.addVar(); m.setQbj(0, 0, 4); m.setObj(0, -8);
        return m;
    }
    private static void optimum(Model m, double x, double objective) {
        assertEquals(OSQP.Status.SOLVED, m.solve(tight()));
        assertEquals(x, m.getSolution(0), 1e-7);
        assertEquals(objective, m.getObj(), 1e-7);
    }
    @Test void addingRowsDoesNotRescaleTheCachedObjective() {
        Model m = quadratic(); optimum(m, 1, -4);
        m.addConstraintRow(new int[]{0}, new double[]{2}, 0, 10);
        optimum(m, 1, -4);
        m.addConstraintRow(new int[]{0}, new double[]{3}, 0, 10);
        optimum(m, 1, -4);
    }
    @Test void changingSparsityPreservesUnchangedMatrices() {
        Model m = quadratic(); m.addVar(); m.setQbj(1, 1, 4); m.setObj(1, -8);
        int row = m.addConstraintRow(new int[]{0}, new double[]{2}, 0, 10);
        optimum(m, 1, -8);
        m.setLhs(row, 1, 2); optimum(m, 1, -8);
        // Change only P's pattern; the new problem has optimum (8/9,8/9).
        m.setQbj(0, 1, 1);
        optimum(m, 8.0/9, -64.0/9);
    }
    @Test void addingVariablesPreservesTheOldObjective() {
        Model m = quadratic(); optimum(m, 1, -4);
        m.addVar(); m.setQbj(1, 1, 1); m.setObj(1, -2);
        optimum(m, 1, -5);
        assertEquals(1, m.getSolution(1), 1e-7);
    }
    @Test void rebuildAfterInfeasibilityDoesNotWarmStartFromNaN() {
        Model m=quadratic(); m.setUb(0,1);
        int row=m.addConstraintRow(new int[]{0},new double[]{1},2,OSQP.OSQP_INFTY);
        assertEquals(OSQP.Status.PRIMAL_INFEASIBLE,m.solve(tight()));
        m.setConstraintBounds(row,0,OSQP.OSQP_INFTY);
        m.addVar(); m.setQbj(1,1,1);
        optimum(m,1,-4);
    }
    @Test void reverseOrderedCrossTermsAccumulate() {
        Model.Builder b = Model.getBuilder();
        Model.Variable x = b.addVariable(), y = b.addVariable();
        b.setObjective().add(1,x,x).add(1,y,y).add(-0.25,x,y).add(-0.75,y,x)
                .add(-1,x).add(-1,y).minimize();
        Model m = b.build(); optimum(m, 1, -1);
        assertEquals(1, m.getSolution(1), 1e-7);
    }
    @Test void offsetChangesNeedNoRefactorisation() {
        Model m = quadratic(); m.setOffset(3); optimum(m, 1, -1);
        OSQP solver = m._solver;
        m.setOffset(7); assertEquals(3, m.getObj(), 1e-7);
        optimum(m, 1, 3); assertSame(solver, m._solver);
    }
    @Test void changingSenseChangesTheProblemAndSignedOffset() {
        Model m = new Model(); m.addVariables(1, new double[]{0}, new double[]{10});
        m.setObj(0, 1); m.setOffset(3); optimum(m, 0, 3);
        m.setSense(true); optimum(m, 10, 13);
        m.setOffset(5); optimum(m, 10, 15);
        m.setSense(false); optimum(m, 0, 5);
    }
    @Test void changingSenseCannotReuseAConvexFactorForANonconvexProblem() {
        Model m = quadratic(); optimum(m, 1, -4);
        m.setSense(true);
        assertThrows(IllegalStateException.class, () -> m.solve(tight()));
    }
    @Test void fixedPatternUpdatesDoNotRetainWriteHistory() {
        CSCMatrixBuilder b = new CSCMatrixBuilder();
        b.set(0,0,1).set(2,0,2).set(1,1,3); b.build(false);
        for (int k = 0; k < 100000; k++) b.set(0,0,k).set(2,0,k+1).set(1,1,k+2);
        assertEquals(3, b.getNumEntries());
        CSCMatrix negative = b.build(true);
        assertArrayEquals(new double[]{-99999,-100000,-100001}, negative.Ax);
        negative.Ax[0] = 42;
        assertEquals(99999, b.build(false).Ax[0]);
    }
    @Test void duplicateTripletsCompactAndRemainEditable() {
        CSCMatrixBuilder b = new CSCMatrixBuilder();
        for (int k = 0; k < 100; k++) b.set(3,1,k).set(0,0,k+1);
        b.build(false); assertEquals(2, b.getNumEntries());
        b.set(3,1,5).set(1,0,7).set(0,0,9);
        assertArrayEquals(new double[]{9,7,5}, b.build(false).Ax);
        assertEquals(3, b.getNumEntries());
    }
    @Test void subsequentSettingsAreValidatedAndApplied() {
        Model m = quadratic(); optimum(m, 1, -4);
        OSQP initial = m._solver;
        OSQP.Settings bad = tight(); bad.max_iter = 0;
        assertThrows(IllegalArgumentException.class, () -> m.solve(bad));
        OSQP.Settings rejected = tight(); rejected.linsys_solver = OSQP.LinSysSolver.OSQP_INDIRECT_SOLVER;
        assertThrows(IllegalArgumentException.class, () -> m.solve(rejected));
        OSQP.Settings one = tight(); one.max_iter = 1; one.warm_start = false;
        assertEquals(OSQP.Status.MAX_ITER_REACHED, m.solve(one));
        assertSame(initial, m._solver);
        assertEquals(1, m._solver.getInfo().iter);
        assertFalse(m._solver.work.settings.warm_start);
        optimum(m, 1, -4);
    }
    @Test void setupSettingsRebuildFromUnscaledData() {
        Model m = quadratic(); optimum(m, 1, -4);
        OSQP initial = m._solver;
        OSQP.Settings s = tight(); s.scaling = 0; s.sigma = 1e-5;
        assertEquals(OSQP.Status.SOLVED, m.solve(s));
        assertNotSame(initial, m._solver);
        assertEquals(1, m.getSolution(0), 1e-7);
        assertEquals(-4, m.getObj(), 1e-7);
    }
    @Test void getParamSettingsApplyToTheDefaultSolve() {
        Model m = quadratic(); m.getParam().verbose = false;
        m.getParam().max_iter = 1; m.getParam().warm_start = false;
        assertEquals(OSQP.Status.MAX_ITER_REACHED, m.solve());
        m.getParam().max_iter = 20000;
        assertEquals(OSQP.Status.SOLVED, m.solve());
    }
    @Test void combinedModelUpdatesMatchAFreshProblem() {
        Model m = quadratic(); int row = m.addConstraintRow(new int[]{0}, new double[]{1}, 0, 2);
        optimum(m, 1, -4); OSQP initial = m._solver;
        for (int k = 0; k < 25; k++) {
            double a = 2 + k/10.0, p = 3 + k/10.0;
            m.setLhs(row, 0, a); m.setQbj(0,0,p);
            double x = Math.min(4/p, 2/a);
            optimum(m, x, p*x*x-8*x);
            assertSame(initial, m._solver);
            assertEquals(1, m._P.getNumEntries());
            assertEquals(2, m._A.getNumEntries());
        }
    }
    @Test void invalidBulkInputDoesNotPartiallyChangeAModel() {
        Model m = quadratic(); optimum(m, 1, -4);
        assertThrows(IllegalArgumentException.class, () -> m.setQuadraticCost(
                new int[]{0,-1},new int[]{0,0},new double[]{50,1}));
        assertThrows(IllegalArgumentException.class, () -> m.setLinearCost(new int[]{0},new double[]{Double.NaN}));
        assertThrows(IllegalArgumentException.class, () -> m.setConstraintBounds(0,Double.NaN,1));
        assertThrows(IllegalArgumentException.class, () -> m.addVariables(2,new double[]{0,2},new double[]{1,1}));
        assertEquals(1, m.getNumVariables()); optimum(m, 1, -4);
    }
    @Test void polishCandidateDoesNotReplacePublicDiagnostics() {
        Model m = quadratic(); optimum(m, 1, -4);
        double before = m._solver.getInfo().rel_kkt_error;
        OSQP.update_info(m._solver.work, 0, true, true);
        assertEquals(before, m._solver.getInfo().rel_kkt_error, 0);
        OSQP.Settings s = tight(); s.polish = true; m.solve(s);
        OSQP.Info info = m._solver.getInfo();
        assertEquals(Math.max(Math.max(info.pri_res,info.dua_res),Math.abs(info.duality_gap)),
                info.rel_kkt_error, 1e-15);
    }
}
