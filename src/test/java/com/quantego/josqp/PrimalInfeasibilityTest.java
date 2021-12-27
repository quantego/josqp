package com.quantego.josqp;


import static org.junit.Assert.assertEquals;

import org.junit.Test;


public class PrimalInfeasibilityTest {

    @Test
    public void testPrimalInfeasibleQPSolve() {
        // Populate data
        final OSQP.Data data = PrimalInfeasibilityTestDataGenerator.generateProblem();
        final PrimalInfeasibilityTestSolsData sols_data = PrimalInfeasibilityTestDataGenerator
                .generateData();

        // Define Solver settings as default
        final OSQP.Settings settings = new OSQP.Settings();
        settings.max_iter = 10000;
        settings.alpha = 1.6;
        settings.polish = true;
        settings.scaling = 0;
        settings.verbose = true;
        settings.warm_start = false;

        // Setup workspace
        final OSQP osqp = new OSQP(data, settings);

        // Solve Problem
        osqp.solve();

        // Compare solver statuses
        assertEquals("Primal infeasible QP test solve: Error in solver status!",sols_data.status_test, osqp.work.info.status);
    }

}
