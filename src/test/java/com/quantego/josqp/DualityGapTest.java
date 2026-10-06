package com.quantego.josqp;

import org.junit.Test;
import static org.junit.Assert.*;

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
        assertEquals("the integral accumulates over the info updates",2000,stuck.work.info.primdual_int,0);
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
}
