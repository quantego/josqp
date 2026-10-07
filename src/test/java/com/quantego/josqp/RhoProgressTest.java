package com.quantego.josqp;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

public class RhoProgressTest {
    private OSQP problem() {
        CSCMatrix one=new CSCMatrix(1,1,1,new int[]{0,1},new int[]{0},new double[]{1});
        OSQP.Settings s=new OSQP.Settings();s.verbose=false;s.scaling=0;s.adaptive_rho_kkt=true;
        // Exercise the progress gate with the tighter, termination-synchronised cadence.
        s.adaptive_rho_interval=0;
        return new OSQP(new OSQP.Data(1,1,one,one,new double[]{-2},new double[]{0},new double[]{1},0),s);
    }
    @Test public void insufficientProgressSkipsRefactorAndLeavesRhoUnchanged() {
        OSQP p=problem();p.work.last_rho_kkt_error=1;
        p.work.info.pri_res=10;p.work.info.dua_res=1;p.work.info.duality_gap=-20;
        double rho=p.work.settings.rho;
        assertTrue(OSQP.adapt_rho(p.work));assertEquals(0,p.work.info.rho_updates);
        assertEquals(rho,p.work.settings.rho,0);assertEquals(1,p.work.last_rho_kkt_error,0);
    }
    @Test public void enabledProgressGateSolvesAnalyticQpAndResetsOnRepeatedSolve() {
        OSQP p=problem();p.work.settings.eps_abs=p.work.settings.eps_rel=1e-8;
        assertEquals(OSQP.Status.SOLVED,p.solve());assertEquals(1,p.getPrimalSolution()[0],1e-7);
        p.work.last_rho_kkt_error=0;p.work.settings.warm_start=false;
        assertEquals(OSQP.Status.SOLVED,p.solve());assertEquals(1,p.getPrimalSolution()[0],1e-7);
    }
}
