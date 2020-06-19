package com.quantego.josqp;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.in;
import static org.hamcrest.Matchers.lessThan;
import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

public class PrimalDualInfeasibilityTest {

    @Test
    public void testOptimal() {
        // Load problem data
        final PrimalDualInfeasibilityTestSolsData data = PrimalDualInfeasibilityTestDataGenerator
                .generateData();

        // Populate problem data
        final OSQP.Data problem = new OSQP.Data(data.P.n, data.A12.m, data.P, data.A12, data.q,
                data.l, data.u1);

        // Define Solver settings as default
        final OSQP.Settings settings = new OSQP.Settings();
        settings.max_iter = 2000;
        settings.alpha = 1.6;
        settings.polish = true;
        settings.scaling = 0;
        settings.verbose = true;

        // Setup workspace
        final OSQP osqp = new OSQP(problem, settings);

        // Solve Problem
        osqp.solve();

        // Compare solver statuses
        assertEquals(OSQP.Status.SOLVED, osqp.work.info.status,
                "Primal dual infeasibility test 1: Error in solver status!");

        // Compare primal solutions
        assertThat("Primal dual infeasibility test 1: Error in primal solution!",
                LinAlg.vec_norm_inf_diff(osqp.work.solution.x, data.x1, problem.n),
                lessThan(OSQPTester.TESTS_TOL));

        // Compare dual solutions
        assertThat("Primal dual infeasibility test 1: Error in dual solution!",
                LinAlg.vec_norm_inf_diff(osqp.work.solution.y, data.y1, problem.m),
                lessThan(OSQPTester.TESTS_TOL));

        // Compare objective values
        assertEquals(osqp.work.info.obj_val, data.obj_value1, OSQPTester.TESTS_TOL,
                "Primal dual infeasibility test 1: Error in objective value!");
    }

    @Test
    public void testPrimInfeas() {
        // Load problem data
        final PrimalDualInfeasibilityTestSolsData data = PrimalDualInfeasibilityTestDataGenerator
                .generateData();

        // Populate problem data
        final OSQP.Data problem = new OSQP.Data(data.P.n, data.A12.m, data.P, data.A12, data.q,
                data.l, data.u2);

        // Define Solver settings as default
        final OSQP.Settings settings = new OSQP.Settings();
        settings.max_iter = 2000;
        settings.alpha = 1.6;
        settings.polish = false;
        settings.scaling = 0;
        settings.verbose = true;

        // Setup workspace
        final OSQP osqp = new OSQP(problem, settings);

        // Solve Problem
        osqp.solve();

        // Compare solver statuses
        assertEquals(OSQP.Status.PRIMAL_INFEASIBLE, osqp.work.info.status,
                "Primal dual infeasibility test 2: Error in solver status!");
    }

    @Test
    public void testDualInfeas() {
        // Load problem data
        final PrimalDualInfeasibilityTestSolsData data = PrimalDualInfeasibilityTestDataGenerator
                .generateData();

        // Populate problem data
        final OSQP.Data problem = new OSQP.Data(data.P.n, data.A34.m, data.P, data.A34, data.q,
                data.l, data.u3);

        // Define Solver settings as default
        final OSQP.Settings settings = new OSQP.Settings();
        settings.max_iter = 2000;
        settings.alpha = 1.6;
        settings.polish = false;
        settings.scaling = 0;
        settings.verbose = true;

        // Setup workspace
        final OSQP osqp = new OSQP(problem, settings);

        // Solve Problem
        osqp.solve();

        // Compare solver statuses
        assertEquals(OSQP.Status.DUAL_INFEASIBLE, osqp.work.info.status,
                "Primal dual infeasibility test 3: Error in solver status!");
    }

    @Test
    public void testPrimalDualInfeas() {
        // Load problem data
        final PrimalDualInfeasibilityTestSolsData data = PrimalDualInfeasibilityTestDataGenerator
                .generateData();

        // Populate problem data
        final OSQP.Data problem = new OSQP.Data(data.P.n, data.A34.m, data.P, data.A34, data.q,
                data.l, data.u4);

        // Define Solver settings as default
        final OSQP.Settings settings = new OSQP.Settings();
        settings.max_iter = 2000;
        settings.alpha = 1.6;
        settings.polish = false;
        settings.scaling = 0;
        settings.verbose = true;

        // Setup workspace
        final OSQP osqp = new OSQP(problem, settings);

        // Solve Problem
        osqp.solve();

        // Compare solver statuses
        assertThat("Primal dual infeasibility test 4: Error in solver status!",
                osqp.work.info.status, in(new OSQP.Status[] { OSQP.Status.PRIMAL_INFEASIBLE,
                        OSQP.Status.DUAL_INFEASIBLE }));
    }
}
