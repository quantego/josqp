package com.quantego.josqp;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.lessThan;
import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

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
        assertEquals(sols_data.status_test, osqp.work.info.status,
                "Unconstrained test solve: Error in solver status!");

        // Compare primal solutions
        assertThat("Unconstrained test solve: Error in primal solution!",
                LinAlg.vec_norm_inf_diff(osqp.work.solution.x, sols_data.x_test, data.n),
                lessThan(OSQPTester.TESTS_TOL));

        // Compare objective values
        assertEquals(sols_data.obj_value_test, osqp.work.info.obj_val, OSQPTester.TESTS_TOL,
                "Unconstrained test solve: Error in objective value!");
    }

}
