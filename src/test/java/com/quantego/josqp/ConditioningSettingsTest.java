package com.quantego.josqp;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

public class ConditioningSettingsTest {
    private OSQP.Data data() {
        CSCMatrix p=new CSCMatrix(2,2,2,new int[]{0,1,2},new int[]{0,1},new double[]{2,2});
        CSCMatrix a=new CSCMatrix(3,2,4,new int[]{0,2,4},new int[]{0,1,0,2},new double[]{1,1,1,1});
        return new OSQP.Data(2,3,p,a,new double[]{-4,-4},new double[]{3,0,0},
                new double[]{3,OSQP.OSQP_INFTY,OSQP.OSQP_INFTY},7);
    }
    private OSQP.Settings settings() {
        OSQP.Settings s=new OSQP.Settings();s.verbose=false;s.adaptive_rho=false;
        s.eps_abs=s.eps_rel=1e-8;s.polish=true;s.max_iter=100000;
        s.scaling_cost=1000;s.rho_eq_over_rho_ineq=1e6;return s;
    }
    private void check(OSQP p,double objective,double multiplier) {
        assertEquals(OSQP.Status.SOLVED,p.solve());
        assertArrayEquals(new double[]{1.5,1.5},p.getPrimalSolution(),1e-6);
        assertEquals(multiplier,p.getDualSolution()[0],1e-6);
        assertEquals(objective,p.getObjectiveValue(),1e-6);
        assertEquals(p.getInfo().obj_val-p.getInfo().dual_obj_val,p.getInfo().duality_gap,1e-10);
    }
    @Test public void conditioningReturnsOriginalCostsAndDualUnits() {
        OSQP.Data d=data();check(new OSQP(d,settings()),-.5,1);
        assertArrayEquals(new double[]{2,2},d.P.Ax,0);
        assertArrayEquals(new double[]{-4,-4},d.q,0);
    }
    @Test public void conditionedMatrixAndVectorUpdatesMatchFreshSolver() {
        OSQP p=new OSQP(data(),settings());check(p,-.5,1);
        p.update_data_mat(new double[]{4,4},null,null,null);
        p.update_data_vec(new double[]{-8,-8},null,null);check(p,-8,2);
        OSQP.Data fresh=data();fresh.P.Ax[0]=fresh.P.Ax[1]=4;
        fresh.q[0]=fresh.q[1]=-8;check(new OSQP(fresh,settings()),-8,2);
    }
    @Test public void ratioChangesApplyToRuntimeRhoAndBoundsUpdates() {
        OSQP.Settings initial=settings();initial.rho_eq_over_rho_ineq=1000;
        OSQP p=new OSQP(data(),initial);check(p,-.5,1);
        p.update_settings(settings());assertEquals(1e5,p.work.rho_vec[0],0);
        p.update_rho(.2);assertEquals(2e5,p.work.rho_vec[0],0);
        p.update_bounds(new double[]{2,0,0},new double[]{4,OSQP.OSQP_INFTY,OSQP.OSQP_INFTY});
        assertEquals(.2,p.work.rho_vec[0],0);
        p.update_bounds(new double[]{3,0,0},new double[]{3,OSQP.OSQP_INFTY,OSQP.OSQP_INFTY});
        assertEquals(2e5,p.work.rho_vec[0],0);check(p,-.5,1);
    }
    @Test public void invalidConditioningAndSetupChangesAreRejected() {
        OSQP.Settings s=settings();s.scaling=0;
        assertThrows(IllegalArgumentException.class,()->new OSQP(data(),s));
        s.scaling=10;s.scaling_cost=Double.NaN;
        assertThrows(IllegalArgumentException.class,()->new OSQP(data(),s));
        s.scaling_cost=1;s.rho_eq_over_rho_ineq=Double.MAX_VALUE;
        assertThrows(IllegalArgumentException.class,()->new OSQP(data(),s));
        OSQP p=new OSQP(data(),settings());OSQP.Settings changed=settings();changed.scaling_cost=100;
        assertThrows(IllegalArgumentException.class,()->p.update_settings(changed));
        OSQP.Data tooLarge=data();tooLarge.q[0]=Double.MAX_VALUE;
        assertThrows(IllegalArgumentException.class,()->new OSQP(tooLarge,settings()));
        OSQP.Settings tiny=settings();tiny.scaling_cost=Double.MIN_VALUE;
        assertThrows(IllegalArgumentException.class,()->new OSQP(data(),tiny));
        tiny.scaling_cost=1e-200;OSQP.Data underflow=data();underflow.P.Ax[0]=1e-200;
        assertThrows(IllegalArgumentException.class,()->new OSQP(underflow,tiny));
    }
    @Test public void modelRebuildsWhenInitialCostScalingChanges() {
        Model m=Model.getBuilder().addVariables(1,new double[]{0},new double[]{10}).linearObjective(new int[]{0},new double[]{-4})
                .quadraticObjective(new int[]{0},new int[]{0},new double[]{1}).minimize().build();
        OSQP.Settings s=settings();s.scaling_cost=1;
        assertEquals(OSQP.Status.SOLVED,m.solve(s));double expected=m.getObj();
        s.scaling_cost=1000;assertEquals(OSQP.Status.SOLVED,m.solve(s));
        assertEquals(expected,m.getObj(),1e-6);
    }

    @Test public void overflowingUpdatesLeaveConditionedDataAndFactorUsable() {
        OSQP p=new OSQP(data(),settings());check(p,-.5,1);
        assertThrows(IllegalArgumentException.class,()->p.update_P(new double[]{Double.MAX_VALUE,2},null));
        check(p,-.5,1);
        assertThrows(IllegalArgumentException.class,()->p.update_data_vec(new double[]{Double.MAX_VALUE,-4},
                new double[]{4,0,0},new double[]{4,OSQP.OSQP_INFTY,OSQP.OSQP_INFTY}));
        check(p,-.5,1);
    }
}
