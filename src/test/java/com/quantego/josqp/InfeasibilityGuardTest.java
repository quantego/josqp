package com.quantego.josqp;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

/**
 * The infeasibility tests reuse the last dual (primal infeasibility) or primal (dual infeasibility)
 * step as the certificate direction. On a feasible problem such a step can satisfy the certificate
 * residual tolerance on a single termination check by numerical accident, and on badly scaled
 * problems the residual can plateau just below the tolerance for many checks. JOSQP therefore
 * requires the certificate to hold over {@link OSQP.Settings#infeasibility_streak} consecutive
 * checks and tightens the effective tolerance through {@link OSQP.Settings#infeasibility_margin}.
 *
 * <p>These tests pin that contract down on the reference infeasibility problems instead of on a
 * particular large case: the upstream setting (streak 1, factor 1) is used to learn the first
 * iteration at which the certificate passes, and the guarded solver is then stopped right after it.
 */
public class InfeasibilityGuardTest {

    private static OSQP solve(CSCMatrix A, double[] u, int streak, Double margin, int maxIter) {
        final PrimalDualInfeasibilityTestSolsData data =
                PrimalDualInfeasibilityTestDataGenerator.generateData();
        final OSQP.Data problem = new OSQP.Data(data.P.n, A.m, data.P, A, data.q, data.l, u, 0);
        final OSQP.Settings settings = new OSQP.Settings();
        settings.max_iter = maxIter;
        settings.alpha = 1.6;
        settings.polish = false;
        settings.scaling = 0;
        settings.verbose = false;
        settings.infeasibility_streak = streak;
        if (margin != null) settings.infeasibility_margin = margin;
        final OSQP osqp = new OSQP(problem, settings);
        osqp.solve();
        return osqp;
    }

    /** Upstream geometry: one check and the OSQP tolerance as-is. */
    private static OSQP solveUpstream(CSCMatrix A, double[] u, int maxIter) {
        return solve(A, u, 1, 1.0, maxIter);
    }

    @Test
    public void defaultsAreConservative() {
        final OSQP.Settings settings = new OSQP.Settings();
        assertEquals(2, settings.infeasibility_streak);
        assertEquals(0.1, settings.infeasibility_margin, 0.0);
    }

    @Test
    public void primalCertificateMustSurviveConsecutiveChecks() {
        final PrimalDualInfeasibilityTestSolsData data =
                PrimalDualInfeasibilityTestDataGenerator.generateData();

        final OSQP upstream = solveUpstream(data.A12, data.u2, 2000);
        assertEquals(OSQP.Status.PRIMAL_INFEASIBLE, upstream.work.info.status);
        final int firstCheck = upstream.work.info.iter;

        // A streak of 1 and any non-positive value are equivalent.
        final OSQP clamped = solve(data.A12, data.u2, 0, 1.0, 2000);
        assertEquals(OSQP.Status.PRIMAL_INFEASIBLE, clamped.work.info.status);
        assertEquals(firstCheck, clamped.work.info.iter);

        // The genuine certificate is still recognised, one check interval later.
        final OSQP sustained = solve(data.A12, data.u2, 2, null, 2000);
        assertEquals(OSQP.Status.PRIMAL_INFEASIBLE, sustained.work.info.status);
        assertTrue(
                sustained.work.info.iter > firstCheck,
                "A sustained certificate must be declared later than a single-check one");

        // Stopping on the first passing check must not produce the accurate verdict. The
        // last-resort approximate check may still report an inaccurate hint.
        final OSQP single = solve(data.A12, data.u2, 2, null, sustained.work.info.iter - 25);
        assertEquals(OSQP.Status.PRIMAL_INFEASIBLE_INACCURATE, single.work.info.status);
    }

    @Test
    public void dualCertificateMustSurviveConsecutiveChecks() {
        final PrimalDualInfeasibilityTestSolsData data =
                PrimalDualInfeasibilityTestDataGenerator.generateData();

        final OSQP upstream = solveUpstream(data.A34, data.u3, 2000);
        assertEquals(OSQP.Status.DUAL_INFEASIBLE, upstream.work.info.status);
        final int firstCheck = upstream.work.info.iter;

        final OSQP sustained = solve(data.A34, data.u3, 2, null, 2000);
        assertEquals(OSQP.Status.DUAL_INFEASIBLE, sustained.work.info.status);
        assertTrue(
                sustained.work.info.iter > firstCheck,
                "A sustained certificate must be declared later than a single-check one");

        final OSQP single = solve(data.A34, data.u3, 2, null, sustained.work.info.iter - 25);
        assertEquals(OSQP.Status.DUAL_INFEASIBLE_INACCURATE, single.work.info.status);
    }

    /** The margin must reject plateaus without hiding genuine certificates. */
    @Test
    public void marginKeepsGenuineCertificates() {
        final PrimalDualInfeasibilityTestSolsData data =
                PrimalDualInfeasibilityTestDataGenerator.generateData();
        for (double margin : new double[] { 1.0, 0.1, 0.01 }) {
            assertEquals(
                OSQP.Status.PRIMAL_INFEASIBLE,
                solve(data.A12, data.u2, 2, margin, 20000).work.info.status,
                "primal certificate lost at margin " + margin);
            assertEquals(
                OSQP.Status.DUAL_INFEASIBLE,
                solve(data.A34, data.u3, 2, margin, 20000).work.info.status,
                "dual certificate lost at margin " + margin);
        }
    }

    @Test
    public void feasibleProblemIsUnaffected() {
        final PrimalDualInfeasibilityTestSolsData data =
                PrimalDualInfeasibilityTestDataGenerator.generateData();

        final OSQP one = solveUpstream(data.A12, data.u1, 2000);
        final OSQP two = solve(data.A12, data.u1, 2, null, 2000);
        assertEquals(OSQP.Status.SOLVED, one.work.info.status);
        assertEquals(OSQP.Status.SOLVED, two.work.info.status);
        assertEquals(one.work.info.obj_val, two.work.info.obj_val, 0.0);
        assertEquals(one.work.info.iter, two.work.info.iter);
    }
}
