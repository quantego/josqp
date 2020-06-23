package com.quantego.josqp;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.lessThan;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

public class UpdateMatricesTest {

    @Test
    public void testFormKKT() {
        // Load problem data
        final UpdateMatricesTestSolsData data = UpdateMatricesTestDataGenerator.generateData();

        // Define rho_vec and sigma to form KKT
        final double sigma = data.test_form_KKT_sigma;
        final int m = data.test_form_KKT_A.m;
        final double[] rho_vec = new double[m];
        final double[] rho_inv_vec = new double[m];
        LinAlg.vec_add_scalar(rho_vec, data.test_form_KKT_rho, m);
        LinAlg.vec_ew_recipr(rho_vec, rho_inv_vec, m);

        // Allocate vectors of indices
        final int[] PtoKKT = new int[data.test_form_KKT_Pu.Ap[data.test_form_KKT_Pu.n]];
        final int[] AtoKKT = new int[data.test_form_KKT_A.Ap[data.test_form_KKT_A.n]];

        // Form KKT matrix storing the index vectors
        final int[][] Pdiag_idx = new int[1][];
        final int[] Pdiag_n = new int[1];
        final CSCMatrix kkt = KKT.form_KKT(data.test_form_KKT_Pu, data.test_form_KKT_A, 0, sigma,
                rho_inv_vec, PtoKKT, AtoKKT, Pdiag_idx, Pdiag_n, null);

        // Assert if KKT matrix is the same as predicted one
        assertTrue(OSQPTester.is_eq_csc(kkt, data.test_form_KKT_KKTu, OSQPTester.TESTS_TOL),
                "Update matrices: error in forming KKT matrix!");

        // Update KKT matrix with new P and new A
        KKT.update_KKT_P(kkt, data.test_form_KKT_Pu_new, PtoKKT, sigma, Pdiag_idx[0], Pdiag_n[0]);
        KKT.update_KKT_A(kkt, data.test_form_KKT_A_new, AtoKKT);

        // Assert if KKT matrix is the same as predicted one
        assertTrue(OSQPTester.is_eq_csc(kkt, data.test_form_KKT_KKTu_new, OSQPTester.TESTS_TOL),
                "Update matrices: error in updating KKT matrix!");
    }

    @Test
    public void testUpdate() {

        // Load problem data
        final UpdateMatricesTestSolsData data = UpdateMatricesTestDataGenerator.generateData();

        // Generate first problem data
        final OSQP.Data problem = new OSQP.Data(data.test_solve_Pu.n, data.test_solve_A.m,
                data.test_solve_Pu, data.test_solve_A, data.test_solve_q, data.test_solve_l,
                data.test_solve_u);

        // Define Solver settings as default
        // Problem settings
        final OSQP.Settings settings = new OSQP.Settings();
        settings.max_iter = 1000;
        settings.alpha = 1.6;
        settings.verbose = true;

        // Setup workspace
        OSQP osqp = new OSQP(problem, settings);

        // Solve Problem
        osqp.solve();

        // Compare solver statuses
        assertEquals(data.test_solve_status, osqp.work.info.status,
                "Update matrices: original problem, error in solver status!");

        // Compare primal solutions
        assertThat("Update matrices: original problem, error in primal solution!",
                LinAlg.vec_norm_inf_diff(osqp.work.solution.x, data.test_solve_x, data.n),
                lessThan(OSQPTester.TESTS_TOL));

        // Compare dual solutions
        assertThat("Update matrices: original problem, error in dual solution!",
                LinAlg.vec_norm_inf_diff(osqp.work.solution.y, data.test_solve_y, data.m),
                lessThan(OSQPTester.TESTS_TOL));

        // Update P
        final int nnzP = data.test_solve_Pu.Ap[data.test_solve_Pu.n];
        final int[] Px_new_idx = new int[nnzP];

        // Generate indices going from beginning to end of P
        for (int i = 0; i < nnzP; i++) {
            Px_new_idx[i] = i;
        }

        osqp.update_P(data.test_solve_Pu_new.Ax, Px_new_idx);

        // Solve Problem
        osqp.solve();

        // Compare solver statuses
        assertEquals(data.test_solve_P_new_status, osqp.work.info.status,
                "Update matrices: problem with updating P, error in solver status!");

        // Compare primal solutions
        assertThat("Update matrices: problem with updating P, error in primal solution!",
                LinAlg.vec_norm_inf_diff(osqp.work.solution.x, data.test_solve_P_new_x,
                        data.n) < OSQPTester.TESTS_TOL);

        // Compare dual solutions
        assertThat("Update matrices: problem with updating P, error in dual solution!",
                LinAlg.vec_norm_inf_diff(osqp.work.solution.y, data.test_solve_P_new_y, data.m),
                lessThan(OSQPTester.TESTS_TOL));

        // Cleanup and setup workspace
        osqp = new OSQP(problem, settings);

        // Update P (all indices)
        osqp.update_P(data.test_solve_Pu_new.Ax, null);

        // Solve Problem
        osqp.solve();

        // Compare solver statuses
        assertEquals(data.test_solve_P_new_status, osqp.work.info.status,
                "Update matrices: problem with updating P (all indices), error in solver status!");

        // Compare primal solutions
        assertThat(
                "Update matrices: problem with updating P (all indices), error in primal solution!",
                LinAlg.vec_norm_inf_diff(osqp.work.solution.x, data.test_solve_P_new_x, data.n),
                lessThan(OSQPTester.TESTS_TOL));

        // Compare dual solutions
        assertThat(
                "Update matrices: problem with updating P (all indices), error in dual solution!",
                LinAlg.vec_norm_inf_diff(osqp.work.solution.y, data.test_solve_P_new_y, data.m),
                lessThan(OSQPTester.TESTS_TOL));

        // Cleanup and setup workspace
        osqp = new OSQP(problem, settings);

        // Update A
        final int nnzA = data.test_solve_A.Ap[data.test_solve_A.n];
        final int[] Ax_new_idx = new int[nnzA];

        // Generate indices going from beginning to end of A
        for (int i = 0; i < nnzA; i++) {
            Ax_new_idx[i] = i;
        }

        osqp.update_A(data.test_solve_A_new.Ax, Ax_new_idx);

        // Solve Problem
        osqp.solve();

        // Compare solver statuses
        assertEquals(data.test_solve_A_new_status, osqp.work.info.status,
                "Update matrices: problem with updating A, error in solver status!");

        // Compare primal solutions
        assertThat("Update matrices: problem with updating A, error in primal solution!",
                LinAlg.vec_norm_inf_diff(osqp.work.solution.x, data.test_solve_A_new_x, data.n),
                lessThan(OSQPTester.TESTS_TOL));

        // Compare dual solutions
        assertThat("Update matrices: problem with updating A, error in dual solution!",
                LinAlg.vec_norm_inf_diff(osqp.work.solution.y, data.test_solve_A_new_y, data.m),
                lessThan(OSQPTester.TESTS_TOL));

        // Cleanup and setup workspace
        osqp = new OSQP(problem, settings);

        // Update A (all indices)
        osqp.update_A(data.test_solve_A_new.Ax, null);

        // Solve Problem
        osqp.solve();

        // Compare solver statuses
        assertEquals(data.test_solve_A_new_status, osqp.work.info.status,
                "Update matrices: problem with updating A (all indices), error in solver status!");

        // Compare primal solutions
        assertThat(
                "Update matrices: problem with updating A (all indices), error in primal solution!",
                LinAlg.vec_norm_inf_diff(osqp.work.solution.x, data.test_solve_A_new_x, data.n),
                lessThan(OSQPTester.TESTS_TOL));

        // Compare dual solutions
        assertThat(
                "Update matrices: problem with updating A (all indices), error in dual solution!",
                LinAlg.vec_norm_inf_diff(osqp.work.solution.y, data.test_solve_A_new_y, data.m),
                lessThan(OSQPTester.TESTS_TOL));

        // Cleanup and setup workspace
        osqp = new OSQP(problem, settings);

        // Update P and A
        osqp.update_P_A(data.test_solve_Pu_new.Ax, Px_new_idx, data.test_solve_A_new.Ax,
                Ax_new_idx);

        // Solve Problem
        osqp.solve();

        // Compare solver statuses
        assertEquals(data.test_solve_P_A_new_status, osqp.work.info.status,
                "Update matrices: problem with updating P and A, error in solver status!");

        // Compare primal solutions
        assertThat("Update matrices: problem with updating P and A, error in primal solution!",
                LinAlg.vec_norm_inf_diff(osqp.work.solution.x, data.test_solve_P_A_new_x, data.n),
                lessThan(OSQPTester.TESTS_TOL));

        // Compare dual solutions
        assertThat("Update matrices: problem with updating P and A, error in dual solution!",
                LinAlg.vec_norm_inf_diff(osqp.work.solution.y, data.test_solve_P_A_new_y, data.m),
                lessThan(OSQPTester.TESTS_TOL * OSQPTester.TESTS_TOL));

        // Cleanup and setup workspace
        osqp = new OSQP(problem, settings);

        // Update P and A (all indices)
        osqp.update_P_A(data.test_solve_Pu_new.Ax, null, data.test_solve_A_new.Ax, null);

        // Solve Problem
        osqp.solve();

        // Compare solver statuses
        assertEquals(data.test_solve_P_A_new_status, osqp.work.info.status,
                "Update matrices: problem with updating P and A (all indices), error in solver status!");

        // Compare primal solutions
        assertThat(
                "Update matrices: problem with updating P and A (all indices), error in primal solution!",
                LinAlg.vec_norm_inf_diff(osqp.work.solution.x, data.test_solve_P_A_new_x, data.n),
                lessThan(OSQPTester.TESTS_TOL));

        // Compare dual solutions
        assertThat(
                "Update matrices: problem with updating P and A (all indices), error in dual solution!",
                LinAlg.vec_norm_inf_diff(osqp.work.solution.y, data.test_solve_P_A_new_y, data.m),
                lessThan(OSQPTester.TESTS_TOL * OSQPTester.TESTS_TOL));
    }
}
