package com.quantego.josqp;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

public class DualityGapTest {
    private OSQP problem(double q, double lower, double upper, int scaling) {
        CSCMatrix one = new CSCMatrix(1,1,1,new int[]{0,1},new int[]{0},new double[]{1});
        OSQP.Settings s = new OSQP.Settings();
        s.verbose = false; s.scaling = scaling; s.eps_abs = s.eps_rel = 1e-8;
        return new OSQP(new OSQP.Data(1,1,one,one,new double[]{q},new double[]{lower},new double[]{upper},9),s);
    }

    @Test public void feasibleStationaryButNonComplementaryPointCannotTerminate() {
        OSQP p = problem(-1000,0,1,0);
        p.work.settings.eps_abs = p.work.settings.eps_rel = 0.001;
        p.work.x[0] = 0; p.work.z[0] = 0; p.work.y[0] = 1000;
        OSQP.update_info(p.work,25,true,false);
        assertEquals(0,p.work.info.pri_res,0); assertEquals(0,p.work.info.dua_res,0);
        assertEquals(1000,p.work.info.duality_gap,0);
        assertFalse(OSQP.check_termination(p.work,false));
        p.work.settings.check_dualgap = false;
        assertTrue(OSQP.check_termination(p.work,false));
    }

    @Test public void optimalPointHasCorrectPrimalDualCostsAndOffsetCancels() {
        OSQP p = problem(-2,0,1,0);
        p.work.x[0]=1; p.work.z[0]=1; p.work.y[0]=1;
        OSQP.update_info(p.work,25,true,false);
        assertEquals(-1.5,p.work.info.obj_val,0);
        assertEquals(-1.5,p.work.info.dual_obj_val,0);
        assertEquals(0,p.work.info.duality_gap,0);
        assertTrue(OSQP.check_termination(p.work,false));
    }

    /** The OSQP 1.0 info fields follow the gap expression of every info update. */
    @Test public void infoTracksPrimDualIntegralAndRelKktError() {
        OSQP optimal = problem(-2,0,1,0);
        optimal.work.x[0]=1; optimal.work.z[0]=1; optimal.work.y[0]=1;
        OSQP.update_info(optimal.work,25,true,false);
        assertEquals(0,optimal.work.info.primdual_int,0);
        assertEquals(0,optimal.work.info.rel_kkt_error,0);

        OSQP stuck = problem(-1000,0,1,0);
        stuck.work.x[0]=0; stuck.work.z[0]=0; stuck.work.y[0]=1000;
        OSQP.update_info(stuck.work,25,true,false);
        assertEquals(1000,stuck.work.info.primdual_int,0);
        assertEquals(1000,stuck.work.info.rel_kkt_error,0);
        OSQP.update_info(stuck.work,50,true,false);
        assertEquals(2000, stuck.work.info.primdual_int, 0, "the integral accumulates over the info updates");
    }

    @Test public void infiniteBoundsAndTinyMultipliersDoNotCreateNaNs() {
        OSQP p = problem(-2,Double.NEGATIVE_INFINITY,Double.POSITIVE_INFINITY,0);
        p.work.x[0]=2; p.work.z[0]=2; p.work.y[0]=1e-20;
        OSQP.update_info(p.work,25,true,false);
        assertEquals(-2,p.work.info.dual_obj_val,0);
        assertEquals(0,p.work.info.duality_gap,0);
    }

    @Test public void scalingAndScaledTerminationReturnSameAnalyticSolution() {
        for (int scaling : new int[]{0,10}) for (boolean scaled : new boolean[]{false,true}) {
            OSQP p = problem(-2,0,1,scaling); p.work.settings.scaled_termination=scaled;
            assertEquals(OSQP.Status.SOLVED,p.solve());
            assertEquals(1,p.getPrimalSolution()[0],1e-7);
            assertEquals(7.5,p.getObjectiveValue(),1e-7);
            assertEquals(p.getInfo().obj_val-p.getInfo().dual_obj_val,p.getInfo().duality_gap,1e-14);
            assertTrue(Math.abs(p.getInfo().duality_gap)<1e-7);
        }
    }

    @Test public void polishedCandidateKeepsAdmmGapAndToleranceSeparate() {
        OSQP p = problem(-2,0,1,0); assertEquals(OSQP.Status.SOLVED,p.solve());
        double gap=p.work.info.duality_gap, scale=p.work.info.gap_scale;
        p.work.pol.x[0]=0; p.work.pol.y[0]=1; p.work.pol.z[0]=0;
        OSQP.update_info(p.work,0,true,true);
        assertEquals(gap,p.work.info.duality_gap,0); assertEquals(scale,p.work.info.gap_scale,0);
        assertEquals(1,p.work.pol.duality_gap,0);
    }

    @Test public void nonFiniteGapCannotProduceSolvedStatus() {
        OSQP p=problem(-2,0,1,0);
        p.work.info.pri_res=0; p.work.info.dua_res=0; p.work.info.duality_gap=Double.NaN;
        assertFalse(OSQP.check_termination(p.work,false));
    }

    @Test public void polishCannotReplaceSolvedPointWithNoncomplementaryCandidate() {
        for (int scaling : new int[]{0,10}) for (boolean scaled : new boolean[]{false,true}) {
            OSQP p=problem(-1000,0,1,scaling);
            p.work.settings.scaled_termination=scaled;
            p.work.settings.eps_abs=p.work.settings.eps_rel=1e-6;
            assertEquals(OSQP.Status.SOLVED,p.solve());
            // x=0 is feasible and stationary with y=1000, but y acts on the inactive upper bound.
            p.work.pol.x[0]=0; p.work.pol.z[0]=0;
            p.work.pol.y[0]=1000*(scaling>0 ? p.work.scaling.c*p.work.scaling.Einv[0] : 1);
            OSQP.update_info(p.work,0,true,true);
            assertEquals(0,p.work.pol.pri_res,1e-12);
            assertEquals(0,p.work.pol.dua_res,1e-9);
            assertFalse(Polish.candidateMeetsTolerances(p.work));
            p.work.settings.check_dualgap=false;
            assertTrue(Polish.candidateMeetsTolerances(p.work));
        }
    }

    @Test public void optimalPolishCandidateRetainsSolutionQuality() {
        OSQP p=problem(-2,0,1,10); assertEquals(OSQP.Status.SOLVED,p.solve());
        p.work.pol.x[0]=p.work.scaling.Dinv[0]; p.work.pol.z[0]=p.work.scaling.E[0];
        p.work.pol.y[0]=p.work.scaling.c*p.work.scaling.Einv[0];
        OSQP.update_info(p.work,0,true,true);
        assertTrue(Polish.candidateMeetsTolerances(p.work));
        p.work.pol.duality_gap=Double.NaN;
        assertFalse(Polish.candidateMeetsTolerances(p.work));
    }

    private OSQP timeLimitedUnconstrained(double q) {
        CSCMatrix p=new CSCMatrix(1,1,1,new int[]{0,1},new int[]{0},new double[]{1});
        CSCMatrix a=new CSCMatrix(0,1,0,new int[]{0,0},new int[0],new double[0]);
        OSQP.Settings s=new OSQP.Settings();s.verbose=false;s.scaling=0;
        s.eps_abs=1e-6;s.eps_rel=0;s.time_limit=Double.MIN_VALUE;
        return new OSQP(new OSQP.Data(1,0,p,a,new double[]{q},new double[0],new double[0],0),s);
    }

    @Test public void timeLimitCanReturnAQualifiedInaccurateSolution() {
        OSQP p=timeLimitedUnconstrained(1e-5);
        assertEquals(OSQP.Status.SOLVED_INACCURATE,p.solve());
        assertEquals(1,p.getInfo().iter);
        assertTrue(p.getInfo().dua_res>1e-6 && p.getInfo().dua_res<1e-5);
        assertTrue(Math.abs(p.getInfo().duality_gap)<1e-5);
    }

    @Test public void timeLimitCanReturnAnExactStationarySolution() {
        OSQP p=timeLimitedUnconstrained(0);
        assertEquals(OSQP.Status.SOLVED,p.solve());
        assertEquals(1,p.getInfo().iter);
        assertEquals(0,p.getInfo().dua_res,0);
    }

    @Test public void timeLimitRetainsItsStatusWhenQualityIsInsufficient() {
        OSQP p=timeLimitedUnconstrained(1);
        assertEquals(OSQP.Status.TIME_LIMIT_REACHED,p.solve());
        assertTrue(p.getInfo().dua_res>1e-5);
    }

    @Test public void relativeOnlyUnconstrainedPolishNeedsNoPrimalTolerance() {
        OSQP p=timeLimitedUnconstrained(-2);
        p.work.settings.eps_abs=0;p.work.settings.eps_rel=1e-8;
        p.work.settings.time_limit=0;
        assertEquals(OSQP.Status.SOLVED,p.solve());
        p.work.pol.x[0]=2;
        OSQP.update_info(p.work,0,true,true);
        assertTrue(Polish.candidateMeetsTolerances(p.work));
    }
}
