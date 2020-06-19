package com.quantego.josqp;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.lessThan;
import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import org.junit.jupiter.api.Test;

public class BasicQPTest {

    @Test
    public void testSolve() {

        // Populate data
        final OSQP.Data data = BasicQPTestDataGenerator.generateProblem();
        final BasicQPTestSolsData sols_data = BasicQPTestDataGenerator.generateData();

        // Problem settings
        final OSQP.Settings settings = new OSQP.Settings();
        settings.max_iter = 2000;
        settings.alpha = 1.6;
        settings.polish = true;
        settings.scaling = 0;
        settings.verbose = true;
        settings.warm_start = false;

        // Setup workspace
        OSQP osqp = new OSQP(data, settings);

        // Solve Problem
        osqp.solve();

        // Compare solver statuses
        assertEquals(sols_data.status_test, osqp.work.info.status,
                "Basic QP test solve: Error in solver status!");

        // Compare primal solutions
        assertThat("Basic QP test solve: Error in primal solution!",
                LinAlg.vec_norm_inf_diff(osqp.work.solution.x, sols_data.x_test, data.n),
                lessThan(OSQPTester.TESTS_TOL));

        // Compare dual solutions
        assertThat("Basic QP test solve: Error in dual solution!",
                LinAlg.vec_norm_inf_diff(osqp.work.solution.y, sols_data.y_test, data.m),
                lessThan(OSQPTester.TESTS_TOL));

        // Compare objective values
        assertEquals(sols_data.obj_value_test, osqp.work.info.obj_val, OSQPTester.TESTS_TOL,
                "Basic QP test solve: Error in objective value!");

        // Try to set wrong settings
        assertThrows(IllegalArgumentException.class, () -> {
            OSQP.osqp_update_rho(osqp.work, -0.1);
        }, "Basic QP test solve: Wrong value of rho not caught!");

//      mu_assert("Basic QP test solve: Wrong value of max_iter not caught!",
//            OSQP.osqp_update_max_iter(work, -1) == 1);
//
//      mu_assert("Basic QP test solve: Wrong value of eps_abs not caught!",
//            osqp_update_eps_abs(work, -1.) == 1);
//
//      mu_assert("Basic QP test solve: Wrong value of eps_rel not caught!",
//            osqp_update_eps_rel(work, -1.) == 1);
//
//      mu_assert("Basic QP test solve: Wrong value of eps_prim_inf not caught!",
//            osqp_update_eps_prim_inf(work, -0.1) == 1);
//
//      mu_assert("Basic QP test solve: Wrong value of eps_dual_inf not caught!",
//            osqp_update_eps_dual_inf(work, -0.1) == 1);
//
//      mu_assert("Basic QP test solve: Wrong value of alpha not caught!",
//            osqp_update_alpha(work, 2.0) == 1);
//
//      mu_assert("Basic QP test solve: Wrong value of warm_start not caught!",
//            osqp_update_warm_start(work, -1) == 1);
//
//      mu_assert("Basic QP test solve: Wrong value of scaled_termination not caught!",
//            osqp_update_scaled_termination(work, 2) == 1);
//
//      mu_assert("Basic QP test solve: Wrong value of check_termination not caught!",
//            osqp_update_check_termination(work, -1) == 1);
//
//      mu_assert("Basic QP test solve: Wrong value of delta not caught!",
//            osqp_update_delta(work, 0.) == 1);
//
//      mu_assert("Basic QP test solve: Wrong value of polish not caught!",
//            osqp_update_polish(work, 2) == 1);
//
//      mu_assert("Basic QP test solve: Wrong value of polish_refine_iter not caught!",
//            osqp_update_polish_refine_iter(work, -1) == 1);
//
//      mu_assert("Basic QP test solve: Wrong value of verbose not caught!",
//            osqp_update_verbose(work, 2) == 1);

        /*
         * ============================= SETUP WITH WRONG SETTINGS
         * =============================
         */

        // Setup workspace with empty settings
        assertThrows(IllegalArgumentException.class, () -> {
            final OSQP osqp1 = new OSQP(data, null);
            osqp1.validate(); // avoid code elimination
        }, "Basic QP test solve: Setup should result in error due to empty settings");

        // Setup workspace with a wrong number of scaling iterations
        int tmp_int = settings.scaling;
        settings.scaling = -1;
        assertThrows(IllegalArgumentException.class, () -> {
            final OSQP osqp1 = new OSQP(data, settings);
            osqp1.validate(); // avoid code elimination
        }, "Basic QP test solve: Setup should result in error due to a negative number"
                + " of scaling iterations");
        settings.scaling = tmp_int;

        // Setup workspace with wrong settings.adaptive_rho_interval
        tmp_int = settings.adaptive_rho_interval;
        settings.adaptive_rho_interval = -1;
        assertThrows(IllegalArgumentException.class, () -> {
            final OSQP osqp1 = new OSQP(data, settings);
            osqp1.validate(); // avoid code elimination
        }, "Basic QP test solve: Setup should result in error due to negative"
                + " settings.adaptive_rho_interval");
        settings.adaptive_rho_interval = tmp_int;

        // Setup workspace with wrong settings.adaptive_rho_tolerance
        double tmp_float = settings.adaptive_rho_tolerance;
        settings.adaptive_rho_tolerance = 0.5;
        assertThrows(

                IllegalArgumentException.class, () -> {
                    final OSQP osqp1 = new OSQP(data, settings);
                    osqp1.validate(); // avoid code elimination
                }, "Basic QP test solve: Setup should result in error due to wrong"
                        + " settings.adaptive_rho_tolerance");
        settings.adaptive_rho_tolerance = tmp_float;

        // Setup workspace with wrong settings.polish_refine_iter
        tmp_int = settings.polish_refine_iter;
        settings.polish_refine_iter = -3;
        assertThrows(IllegalArgumentException.class, () -> {
            final OSQP osqp1 = new OSQP(data, settings);
            osqp1.validate(); // avoid code elimination
        }, "Basic QP test solve: Setup should result in error due to negative"
                + " settings.polish_refine_iter");
        settings.polish_refine_iter = tmp_int;

        // Setup workspace with wrong settings.rho
        tmp_float = settings.rho;
        settings.rho = 0.0;
        assertThrows(IllegalArgumentException.class, () -> {
            final OSQP osqp1 = new OSQP(data, settings);
            osqp1.validate(); // avoid code elimination
        }, "Basic QP test solve: Setup should result in error due to non-positive"
                + " settings.rho");
        settings.rho = tmp_float;

        // Setup workspace with wrong settings.sigma
        tmp_float = settings.sigma;
        settings.sigma = -0.1;
        assertThrows(IllegalArgumentException.class, () -> {
            final OSQP osqp1 = new OSQP(data, settings);
            osqp1.validate(); // avoid code elimination
        }, "Basic QP test solve: Setup should result in error due to non-positive"
                + " settings.sigma");
        settings.sigma = tmp_float;

        // Setup workspace with wrong settings.delta
        tmp_float = settings.delta;
        settings.delta = -1.1;
        assertThrows(IllegalArgumentException.class, () -> {
            final OSQP osqp1 = new OSQP(data, settings);
            osqp1.validate(); // avoid code elimination
        }, "Basic QP test solve: Setup should result in error due to non-positive"
                + " settings.delta");
        settings.delta = tmp_float;

        // Setup workspace with wrong settings.max_iter
        tmp_int = settings.max_iter;
        settings.max_iter = 0;
        assertThrows(IllegalArgumentException.class, () -> {
            final OSQP osqp1 = new OSQP(data, settings);
            osqp1.validate(); // avoid code elimination
        }, "Basic QP test solve: Setup should result in error due to non-positive"
                + " settings.max_iter");
        settings.max_iter = tmp_int;

        // Setup workspace with wrong settings.eps_abs
        tmp_float = settings.eps_abs;
        settings.eps_abs = -1.1;
        assertThrows(IllegalArgumentException.class, () -> {
            final OSQP osqp1 = new OSQP(data, settings);
            osqp1.validate(); // avoid code elimination
        }, "Basic QP test solve: Setup should result in error due to negative"
                + " settings.eps_abs");
        settings.eps_abs = tmp_float;

        // Setup workspace with wrong settings.eps_rel
        tmp_float = settings.eps_rel;
        settings.eps_rel = -0.1;
        assertThrows(IllegalArgumentException.class, () -> {
            final OSQP osqp1 = new OSQP(data, null);
            osqp1.validate(); // avoid code elimination
        }, "Basic QP test solve: Setup should result in error due to negative"
                + " settings.eps_rel");
        settings.eps_rel = tmp_float;

        // Setup workspace with wrong settings.eps_prim_inf
        tmp_float = settings.eps_prim_inf;
        settings.eps_prim_inf = -0.1;
        assertThrows(IllegalArgumentException.class, () -> {
            final OSQP osqp1 = new OSQP(data, settings);
            osqp1.validate(); // avoid code elimination
        }, "Basic QP test solve: Setup should result in error due to non-positive"
                + " settings.eps_prim_inf");
        settings.eps_prim_inf = tmp_float;

        // Setup workspace with wrong settings.eps_dual_inf
        tmp_float = settings.eps_dual_inf;
        settings.eps_dual_inf = 0.0;
        assertThrows(IllegalArgumentException.class, () -> {
            final OSQP osqp1 = new OSQP(data, settings);
            osqp1.validate(); // avoid code elimination
        }, "Basic QP test solve: Setup should result in error due to non-positive"
                + " settings.eps_dual_inf");
        settings.eps_dual_inf = tmp_float;

        // Setup workspace with wrong settings.alpha
        tmp_float = settings.alpha;
        settings.alpha = 2.0;
        assertThrows(IllegalArgumentException.class, () -> {
            final OSQP osqp1 = new OSQP(data, settings);
            osqp1.validate(); // avoid code elimination
        }, "Basic QP test solve: Setup should result in error due to wrong" + " settings.alpha");
        settings.alpha = tmp_float;

        // Setup workspace with wrong settings.check_termination
        tmp_int = settings.check_termination;
        settings.check_termination = -1;
        assertThrows(IllegalArgumentException.class, () -> {
            final OSQP osqp1 = new OSQP(data, settings);
            osqp1.validate(); // avoid code elimination
        }, "Basic QP test solve: Setup should result in error due to non-boolean"
                + " settings.check_termination");
        settings.check_termination = tmp_int;

        /*
         * ========================= SETUP WITH WRONG DATA =========================
         */

        // Setup workspace with empty data
        assertThrows(IllegalArgumentException.class, () -> {
            final OSQP osqp1 = new OSQP(null, settings);
            osqp1.validate(); // avoid code elimination
        }, "Basic QP test solve: Setup should result in error due to empty data");

        // Setup workspace with wrong data.m
        final OSQP.Data tmp_data1 = new OSQP.Data(data.m - 1, data.n, data.P, data.A, data.q,
                data.l, data.u);
        assertThrows(IllegalArgumentException.class, () -> {
            final OSQP osqp1 = new OSQP(tmp_data1, settings);
            osqp1.validate(); // avoid code elimination
        }, "Basic QP test solve: Setup should result in error due to wrong data.m");

        // Setup workspace with wrong data.n
        tmp_int = data.n;
        final OSQP.Data tmp_data2 = new OSQP.Data(data.m, data.n + 1, data.P, data.A, data.q,
                data.l, data.u);
        assertThrows(IllegalArgumentException.class, () -> {
            final OSQP osqp1 = new OSQP(tmp_data2, settings);
            osqp1.validate(); // avoid code elimination
        }, "Basic QP test solve: Setup should result in error due to wrong data.n");

        // Setup workspace with zero data.n
        final OSQP.Data tmp_data3 = new OSQP.Data(data.m, 0, data.P, data.A, data.q, data.l,
                data.u);
        assertThrows(IllegalArgumentException.class, () -> {
            final OSQP osqp1 = new OSQP(tmp_data3, settings);
            osqp1.validate(); // avoid code elimination
        }, "Basic QP test solve: Setup should result in error due to zero data.n");

        // Setup workspace with wrong P.m
        CSCMatrix P_tmp = new CSCMatrix(data.n + 1, data.P.n, data.P.nzmax, data.P.Ap, data.P.Ai,
                data.P.Ax);
        final OSQP.Data tmp_data4 = new OSQP.Data(data.m, data.n, P_tmp, data.A, data.q, data.l,
                data.u);
        assertThrows(IllegalArgumentException.class, () -> {
            final OSQP osqp1 = new OSQP(tmp_data4, settings);
            osqp1.validate();
        }, "Basic QP test solve: Setup should result in error due to wrong P.m");

        // Setup workspace with wrong P.n
        P_tmp = new CSCMatrix(data.P.m, data.n + 1, data.P.nzmax, data.P.Ap, data.P.Ai, data.P.Ax);
        final OSQP.Data tmp_data5 = new OSQP.Data(data.m, data.n, P_tmp, data.A, data.q, data.l,
                data.u);
        assertThrows(IllegalArgumentException.class, () -> {
            final OSQP osqp1 = new OSQP(tmp_data5, settings);
            osqp1.validate();
        }, "Basic QP test solve: Setup should result in error due to wrong P.n");

        // Setup workspace with non-upper-triangular P
        final int P_tmp_m = 2;
        final int P_tmp_n = 2;
        final int P_tmp_nzmax = 4;
        final double[] P_tmp_x = new double[] { 4.0, 1.0, 1.0, 2.0 };
        final int[] P_tmp_i = new int[] { 0, 1, 0, 1 };
        final int[] P_tmp_p = new int[] { 0, 2, 4 };
        P_tmp = new CSCMatrix(P_tmp_m, P_tmp_n, P_tmp_nzmax, P_tmp_p, P_tmp_i, P_tmp_x);
        final OSQP.Data tmp_data6 = new OSQP.Data(data.m, data.n, P_tmp, data.A, data.q, data.l,
                data.u);
        assertThrows(IllegalArgumentException.class, () -> {
            final OSQP osqp1 = new OSQP(tmp_data6, settings);
            osqp1.validate();
        }, "Basic QP test solve: Setup should result in error due to non-triu structure of P");

        // Setup workspace with non-consistent bounds
        data.l[0] = data.u[0] + 1.0;
        assertThrows(IllegalArgumentException.class, () -> {
            final OSQP osqp1 = new OSQP(data, settings);
            osqp1.validate();
        }, "Basic QP test solve: Setup should result in error due to non-consistent bounds");
    }

    @Test
    public void testUpdate() {

        // Populate data
        final OSQP.Data data = BasicQPTestDataGenerator.generateProblem();
        final OSQP.Data originalData = new OSQP.Data(data);
        final BasicQPTestSolsData sols_data = BasicQPTestDataGenerator.generateData();

        // Define Solver settings as default
        final OSQP.Settings settings = new OSQP.Settings();
        settings.max_iter = 200;
        settings.alpha = 1.6;
        settings.polish = true;
        settings.scaling = 0;
        settings.verbose = true;
        settings.warm_start = false;

        // Setup workspace
        final OSQP osqp = new OSQP(data, settings);

        // ====================================================================
        // Update data
        // ====================================================================

        // Update linear cost
        osqp.update_lin_cost(sols_data.q_new);
        assertThat("Basic QP test update: Error in updating linear cost!",
                LinAlg.vec_norm_inf_diff(osqp.work.data.q, sols_data.q_new, data.n),
                lessThan(OSQPTester.TESTS_TOL));

        // UPDATE BOUND
        // Try to update with non-consistent values
        assertThrows(Exception.class, () -> osqp.update_bounds(sols_data.u_new, sols_data.l_new),
                "Basic QP test update: Error in bounds update ordering not caught!");

        // Now update with correct values
        assertDoesNotThrow(() -> osqp.update_bounds(sols_data.l_new, sols_data.u_new),
                "Basic QP test update: Error in bounds update ordering!");

        assertThat("Basic QP test update: Error in bounds update, lower bound!",
                LinAlg.vec_norm_inf_diff(osqp.work.data.l, sols_data.l_new, data.m),
                lessThan(OSQPTester.TESTS_TOL));

        assertThat("Basic QP test update: Error in bounds update, upper bound!",
                LinAlg.vec_norm_inf_diff(osqp.work.data.u, sols_data.u_new, data.m),
                lessThan(OSQPTester.TESTS_TOL));

        // Return original values
        osqp.update_bounds(originalData.l, originalData.u);

        // UPDATE LOWER BOUND
        // Try to update with non-consistent values
        assertThrows(Exception.class, () -> osqp.update_lower_bound(sols_data.u_new),
                "Basic QP test update: Error in lower bound update ordering not caught!");

        // Now update with correct values
        assertDoesNotThrow(() -> osqp.update_lower_bound(sols_data.l_new),
                "Basic QP test update: Error in lower bound update ordering!");

        assertThat("Basic QP test update: Error in updating lower bound!",
                LinAlg.vec_norm_inf_diff(osqp.work.data.l, sols_data.l_new, data.m),
                lessThan(OSQPTester.TESTS_TOL));

        // Return original values
        osqp.update_lower_bound(originalData.l);

        // UPDATE UPPER BOUND
        // Try to update with non-consistent values
        assertThrows(Exception.class, () -> osqp.update_upper_bound(sols_data.l_new),
                "Basic QP test update: Error in upper bound update: ordering not caught!");

        // Now update with correct values
        assertDoesNotThrow(() -> osqp.update_upper_bound(sols_data.u_new),
                "Basic QP test update: Error in upper bound update: ordering!");

        assertThat("Basic QP test update: Error in updating upper bound!",
                LinAlg.vec_norm_inf_diff(osqp.work.data.u, sols_data.u_new, data.m),
                lessThan(OSQPTester.TESTS_TOL));
    }

    @Test
    public void testCheckTermination() {
        // Populate data
        final OSQP.Data data = BasicQPTestDataGenerator.generateProblem();
        final BasicQPTestSolsData sols_data = BasicQPTestDataGenerator.generateData();

        // Define Solver settings as default
        final OSQP.Settings settings = new OSQP.Settings();
        settings.max_iter = 200;
        settings.alpha = 1.6;
        settings.polish = false;
        settings.scaling = 0;
        settings.verbose = true;
        settings.check_termination = 0;
        settings.warm_start = false;

        // Setup workspace
        final OSQP osqp = new OSQP(data, settings);

        // Solve Problem
        osqp.solve();

        // Check if iter == max_iter
        assertEquals(osqp.work.settings.max_iter, osqp.work.info.iter,
                "Basic QP test check termination: Error in number of iterations taken!");

        // Compare solver statuses
        assertEquals(sols_data.status_test, osqp.work.info.status,
                "Basic QP test check termination: Error in solver status!");

        // Compare primal solutions
        assertThat("Basic QP test check termination: Error in primal solution!",
                LinAlg.vec_norm_inf_diff(osqp.work.solution.x, sols_data.x_test, data.n),
                lessThan(OSQPTester.TESTS_TOL));

        // Compare dual solutions
        // print_vec(work.solution.y, data.m, "y_sol");
        // print_vec(sols_data.y_test, data.m, "y_test");
        assertThat("Basic QP test check termination: Error in dual solution!",
                LinAlg.vec_norm_inf_diff(osqp.work.solution.y, sols_data.y_test, data.m),
                lessThan(OSQPTester.TESTS_TOL));

        // Compare objective values
        assertEquals(sols_data.obj_value_test, osqp.work.info.obj_val, OSQPTester.TESTS_TOL,
                "Basic QP test check termination: Error in objective value!");
    }

    @Test
    public void testUpdateRho() {

        // Populate data
        final OSQP.Data data = BasicQPTestDataGenerator.generateProblem();
        final BasicQPTestSolsData sols_data = BasicQPTestDataGenerator.generateData();

        // Define Solver settings as default
        double rho = 0.7;
        final OSQP.Settings settings = new OSQP.Settings();
        settings.rho = rho;
        settings.adaptive_rho = false; // Disable adaptive rho for this test
        settings.eps_abs = 5e-05;
        settings.eps_rel = 5e-05;
        settings.check_termination = 1;

        // Setup workspace
        final OSQP osqp1 = new OSQP(data, settings);

        // Solve Problem
        osqp1.solve();

        // Store number of iterations
        final int n_iter_new_solver = osqp1.work.info.iter;

        // Compare solver statuses
        assertEquals(sols_data.status_test, osqp1.work.info.status,
                "Update rho test solve: Error in solver status!");

        // Compare primal solutions
        assertThat("Update rho test solve: Error in primal solution!",
                LinAlg.vec_norm_inf_diff(osqp1.work.solution.x, sols_data.x_test, data.n)
                        / LinAlg.vec_norm_inf(sols_data.x_test, data.n),
                lessThan(OSQPTester.TESTS_TOL));

        // Compare dual solutions
        assertThat("Update rho test solve: Error in dual solution!",
                LinAlg.vec_norm_inf_diff(osqp1.work.solution.y, sols_data.y_test, data.m)
                        / LinAlg.vec_norm_inf(sols_data.y_test, data.m),
                lessThan(OSQPTester.TESTS_TOL));

        // Compare objective values
        assertEquals(sols_data.obj_value_test, osqp1.work.info.obj_val, OSQPTester.TESTS_TOL,
                "Update rho test solve: Error in objective value!");

        // Create new problem with different rho and update it
        final OSQP.Settings settings2 = new OSQP.Settings();
        settings2.rho = 0.1;
        settings2.adaptive_rho = false;
        settings2.check_termination = 1;
        settings2.eps_abs = 5e-05;
        settings2.eps_rel = 5e-05;
        final OSQP.Data data2 = BasicQPTestDataGenerator.generateProblem();

        // Setup workspace
        final OSQP osqp2 = new OSQP(data2, settings2);

        // Update rho
        assertDoesNotThrow(() -> osqp2.update_rho(rho),
                "Basic QP test update rho: Error update rho!");

        // Solve Problem
        osqp2.solve();

        // Compare solver statuses
        assertEquals(sols_data.status_test, osqp2.work.info.status,
                "Basic QP test update rho: Error in solver status!");

        // Compare primal solutions
        assertThat("Basic QP test update rho: Error in primal solution!",
                LinAlg.vec_norm_inf_diff(osqp2.work.solution.x, sols_data.x_test, data.n)
                        / LinAlg.vec_norm_inf(sols_data.x_test, data.n),
                lessThan(OSQPTester.TESTS_TOL));

        // Compare dual solutions
        assertThat("Basic QP test update rho: Error in dual solution!",
                LinAlg.vec_norm_inf_diff(osqp2.work.solution.y, sols_data.y_test, data.m)
                        / LinAlg.vec_norm_inf(sols_data.y_test, data.m),
                lessThan(OSQPTester.TESTS_TOL));

        // Compare objective values
        assertEquals(sols_data.obj_value_test, osqp2.work.info.obj_val, OSQPTester.TESTS_TOL,
                "Basic QP test update rho: Error in objective value!");

        // Assert same number of iterations
        assertEquals(n_iter_new_solver, osqp2.work.info.iter,
                "Basic QP test update rho: Error in number of iterations!");
    }

    @Test
    public void testWarmStart() {
        // Cold started variables
        double[] x0 = new double[] { 0.0, 0.0 };
        double[] y0 = new double[] { 0.0, 0.0, 0.0, 0.0 };

        // Optimal solution
        double[] xopt = new double[] { 0.3, 0.7, };
        double[] yopt = new double[] { -2.9, 0.0, 0.2, 0.0, };

        // Populate data
        final OSQP.Data data = BasicQPTestDataGenerator.generateProblem();
        // final basic_qp_sols_data sols_data = generate_problem_basic_qp_sols_data();

        // Define Solver settings as default
        final OSQP.Settings settings = new OSQP.Settings();
        settings.check_termination = 1;

        // Setup workspace
        final OSQP osqp = new OSQP(data, settings);

        // Solve Problem
        osqp.solve();
        final int iter = osqp.work.info.iter;

        // Cold start and solve again
        osqp.warm_start(x0, y0);
        osqp.solve();

        // Check if the number of iterations is the same
        assertEquals(iter, osqp.work.info.iter, "Basic QP test warm start: Cold start error!");

        // Warm start from the solution and solve again
        osqp.warm_start_x(xopt);
        osqp.warm_start_y(yopt);
        osqp.solve();

        // Check that the number of iterations equals 1
        assertEquals(1, osqp.work.info.iter, "Basic QP test warm start: Warm start error!");
    }

   /* @Test
    public void test_basic_qp_solve_pardiso() {

        int exitFlag;
        OSQP.Settings settings = new OSQP.Settings();
        OSQP.Data data;
        DataGenerator.basic_qp_sols_data qp_sols_data;

        data = DataGenerator.generate_problem_basic_qp();
        qp_sols_data = DataGenerator.generate_problem_basic_qp_sols_data();

        settings.max_iter = 2000;
        settings.alpha = 1.6;
        settings.polish = true;
        settings.scaling = 0;
        settings.verbose = true;
        settings.warm_start = false;

    } */
}
