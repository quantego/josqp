package com.quantego.josqp;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.util.function.Consumer;

import org.junit.jupiter.api.Test;

/**
 * The consolidated OSQP 1.0 update entry points: {@code osqp_update_data_vec} and
 * {@code osqp_update_data_mat}. They must produce the same solver state as the individual 0.6
 * calls, leave untouched vectors alone and reject malformed index arrays.
 *
 * <p>The raw {@link OSQP.Data} API solves {@code min 0.5*x'Px + q'x}, so the P value array is the
 * upper triangle of the Hessian and its second entry is the cross term. The lower bound of the
 * default data is binding, which is what makes an update of A visible in the optimum.</p>
 */
public class UpdateApiTest {

    private static final double[] BASE_P = { 1.0, -0.5, 1.0 };
    private static final double[] BASE_A = { 1.0, 1.0 };
    private static final double[] BASE_Q = { -1.0, -1.0 };
    private static final double[] NEW_P = { 1.0, 0.5, 2.0 };
    private static final double[] NEW_A = { 2.0, 1.0 };
    private static final double[] NEW_Q = { -3.0, 0.5 };
    private static final double[] NEW_L = { -2.0 };
    private static final double[] INF = { OSQP.OSQP_INFTY };
    private static final double BINDING = 6.0;

    private static OSQP.Settings tight() {
        OSQP.Settings settings = new OSQP.Settings();
        settings.eps_abs = 1e-9;
        settings.eps_rel = 1e-9;
        settings.verbose = false;
        return settings;
    }

    private static OSQP.Data data(double[] P, double lower, double[] A, double[] q) {
        CSCMatrix p = new CSCMatrix(2, 2, 3, new int[] { 0, 1, 3 }, new int[] { 0, 0, 1 }, P.clone());
        CSCMatrix a = new CSCMatrix(1, 2, 2, new int[] { 0, 1, 2 }, new int[] { 0, 0 }, A.clone());
        return new OSQP.Data(2, 1, p, a, q.clone(), new double[] { lower },
                new double[] { OSQP.OSQP_INFTY }, 0);
    }

    private static OSQP.Data start() {
        return data(BASE_P, BINDING, BASE_A, BASE_Q);
    }

    private static double[] solutionOf(OSQP.Data data) {
        OSQP solver = new OSQP(data, tight());
        solver.solve();
        return new double[] { solver.getPrimalSolution()[0], solver.getPrimalSolution()[1],
                solver.getObjectiveValue() };
    }

    private static double[] solutionAfterUpdate(OSQP.Data data, Consumer<OSQP> update) {
        OSQP solver = new OSQP(data, tight());
        solver.solve();
        update.accept(solver);
        solver.solve();
        return new double[] { solver.getPrimalSolution()[0], solver.getPrimalSolution()[1],
                solver.getObjectiveValue() };
    }

    private static void assertSameSolution(String message, double[] expected, double[] actual) {
        assertEquals(expected[0], actual[0], 1e-6, message + " (x)");
        assertEquals(expected[1], actual[1], 1e-6, message + " (y)");
        assertEquals(expected[2], actual[2], 1e-6, message + " (objective)");
    }

    @Test
    public void dataMatMatchesSeparateMatrixUpdates() {
        double[] expected = solutionOf(data(NEW_P, BINDING, NEW_A, BASE_Q));
        assertNotEquals(solutionOf(start())[0], expected[0], 1e-6,
                "the test data must differ from the updated data");

        double[] combined = solutionAfterUpdate(start(),
                solver -> solver.update_data_mat(NEW_P, null, NEW_A, null));
        double[] separately = solutionAfterUpdate(start(), solver -> {
            solver.update_P(NEW_P, null);
            solver.update_A(NEW_A, null);
        });
        double[] legacyName = solutionAfterUpdate(start(),
                solver -> solver.update_P_A(NEW_P, null, NEW_A, null));

        assertSameSolution("update_data_mat", expected, combined);
        assertSameSolution("update_P followed by update_A", expected, separately);
        assertSameSolution("update_P_A", expected, legacyName);
    }

    @Test
    public void dataMatSupportsPartialIndexUpdates() {
        double[] expected = solutionOf(data(new double[] { 1.0, 0.25, 1.0 }, BINDING, NEW_A, BASE_Q));
        assertNotEquals(solutionOf(data(BASE_P, BINDING, NEW_A, BASE_Q))[0], expected[0], 1e-6,
                "the changed cross term must move the optimum");
        double[] actual = solutionAfterUpdate(data(BASE_P, BINDING, NEW_A, BASE_Q),
                solver -> solver.update_data_mat(new double[] { 0.25 }, new int[] { 1 }, null, null));
        assertSameSolution("partial P update", expected, actual);
    }

    @Test
    public void dataMatLeavesTheMissingMatrixAlone() {
        double[] expected = solutionOf(data(BASE_P, BINDING, NEW_A, BASE_Q));
        double[] actual = solutionAfterUpdate(start(),
                solver -> solver.update_data_mat(null, null, NEW_A, null));
        assertSameSolution("only A updated", expected, actual);

        OSQP.Data untouchedStart = start();
        double[] untouched = solutionAfterUpdate(untouchedStart,
                solver -> solver.update_data_mat(null, null, null, null));
        assertSameSolution("no-op update", solutionOf(start()), untouched);
    }

    @Test
    public void dataVecMatchesIndividualVectorUpdates() {
        double[] expected = solutionOf(data(BASE_P, NEW_L[0], BASE_A, NEW_Q));
        double[] combined = solutionAfterUpdate(start(),
                solver -> solver.update_data_vec(NEW_Q, NEW_L, INF));
        double[] separately = solutionAfterUpdate(start(), solver -> {
            solver.update_lin_cost(NEW_Q);
            solver.update_bounds(NEW_L, INF);
        });
        assertNotEquals(solutionOf(start())[0], expected[0], 1e-6,
                "the updated vectors must move the optimum");
        assertSameSolution("update_data_vec", expected, combined);
        assertSameSolution("update_lin_cost and update_bounds", expected, separately);
    }

    @Test
    public void dataVecAllowsOneSidedBoundUpdates() {
        double[] expected = solutionOf(data(BASE_P, NEW_L[0], BASE_A, BASE_Q));
        double[] actual = solutionAfterUpdate(start(),
                solver -> solver.update_data_vec(null, NEW_L, null));
        assertSameSolution("lower bound only", expected, actual);

        OSQP.Data untouchedStart = start();
        double[] noop = solutionAfterUpdate(untouchedStart,
                solver -> solver.update_data_vec(null, null, null));
        assertSameSolution("all vectors null", solutionOf(start()), noop);
    }

    @Test
    public void malformedUpdatesAreRejected() {
        OSQP solver = new OSQP(start(), tight());
        assertThrows(IllegalArgumentException.class,
                () -> solver.update_data_mat(new double[] { 1.0, 2.0 }, new int[] { 0 }, null, null));
        assertThrows(IllegalArgumentException.class,
                () -> solver.update_data_mat(new double[] { 1.0 }, new int[] { 7 }, null, null));
        assertThrows(IllegalArgumentException.class,
                () -> solver.update_data_mat(new double[] { 1.0, 2.0, 3.0, 4.0 }, null, null, null));
        assertThrows(IllegalArgumentException.class,
                () -> solver.update_data_mat(null, new int[] { 0 }, null, null));
        assertThrows(IllegalArgumentException.class,
                () -> solver.update_data_vec(null, new double[] { 1.0, 2.0 }, null));
    }
}
