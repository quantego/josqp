package com.quantego.josqp;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.lessThan;

import org.junit.jupiter.api.Test;

public class SolveLinSysTest {
    @Test
    public void testSolveKKT() {
        final SolveLinSysTestSolsData data = SolveLinSysTestDataGenerator.generateData();

        // Settings
        final OSQP.Settings settings = new OSQP.Settings();
        settings.rho = data.test_solve_KKT_rho;
        settings.sigma = data.test_solve_KKT_sigma;

        // Set rho_vec
        final int m = data.test_solve_KKT_A.m;
        final double[] rho_vec = new double[m];
        LinAlg.vec_add_scalar(rho_vec, settings.rho, m);

        // Form and factorize KKT matrix
        final LinSys s = new LinSys(data.test_solve_KKT_Pu, data.test_solve_KKT_A, settings.sigma,
                rho_vec, false);

        // Solve KKT x = b via LDL given factorization
        s.solve(data.test_solve_KKT_rhs);

        assertThat("Linear systems solve tests: error in forming and solving KKT system!",
                LinAlg.vec_norm_inf_diff(data.test_solve_KKT_rhs, data.test_solve_KKT_x,
                        data.test_solve_KKT_m + data.test_solve_KKT_n),
                lessThan(OSQPTester.TESTS_TOL));
    }
}
