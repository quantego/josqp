package com.quantego.josqp;

import static org.junit.Assert.assertEquals;

import org.junit.Test;


public class UnconstrainedTest {

    @Test
    public void testUnconstrainedSolve() {

        // Populate data
        final OSQP.Data data = UnconstrainedTestDataGenerator.generateProblem();
        final UnconstrainedTestSolsData sols_data = UnconstrainedTestDataGenerator.generateData();

        // Define Solver settings as default
        final OSQP.Settings settings = new OSQP.Settings();
        settings.verbose = true;

        // Setup workspace
        final OSQP osqp = new OSQP(data, settings);

        // Solve Problem first time
        osqp.solve();

        // Compare solver statuses
        assertEquals("Unconstrained test solve: Error in solver status!",sols_data.status_test, osqp.work.info.status
                );

        // Compare primal solutions
        assertEquals("Unconstrained test solve: Error in primal solution!",
                LinAlg.vec_norm_inf_diff(osqp.work.solution.x, sols_data.x_test, data.n),
                OSQPTester.TESTS_TOL);

        // Compare objective values
        assertEquals("Unconstrained test solve: Error in objective value!",sols_data.obj_value_test, osqp.work.info.obj_val, OSQPTester.TESTS_TOL
                );
    }

}
