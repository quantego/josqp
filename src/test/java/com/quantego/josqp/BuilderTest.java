package com.quantego.josqp;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;


/** Covers the algebraic model builder, which assembles CSC data and reports objective values. */
public class BuilderTest {

    private static OSQP.Settings tight() {
        OSQP.Settings settings = new OSQP.Settings();
        settings.eps_abs = 1e-9;
        settings.eps_rel = 1e-9;
        return settings;
    }

    @Test
    public void cscAssemblyIsSortedAndLastWriteWins() {
        CSCMatrix matrix = new CSCMatrixBuilder()
                .set(2, 0, 1.0)
                .set(0, 0, 2.0)
                .set(0, 0, 3.0)
                .set(1, 1, 5.0)
                .build(false);
        assertArrayEquals(new int[]{0, 2, 3}, matrix.Ap);
        assertArrayEquals(new int[]{0, 2, 1}, matrix.Ai);
        assertEquals(3.0, matrix.Ax[0], 0.0);
        assertEquals(1.0, matrix.Ax[1], 0.0);
        assertEquals(5.0, matrix.Ax[2], 0.0);
    }

    @Test
    public void cscAssemblySupportsOutOfOrderColumns() {
        CSCMatrix matrix = new CSCMatrixBuilder()
                .set(1, 2, 4.0)
                .set(0, 0, 1.0)
                .set(2, 2, 6.0)
                .set(1, 1, 2.0)
                .build(false);
        assertArrayEquals(new int[]{0, 1, 2, 4}, matrix.Ap);
        assertArrayEquals(new int[]{0, 1, 1, 2}, matrix.Ai);
        assertEquals(1.0, matrix.Ax[0], 0.0);
        assertEquals(2.0, matrix.Ax[1], 0.0);
        assertEquals(4.0, matrix.Ax[2], 0.0);
        assertEquals(6.0, matrix.Ax[3], 0.0);
    }

    @Test
    public void minimizeWithOffsetReportsOffsetOnce() {
        Model.Builder builder = Model.getBuilder();
        Model.Variable x = builder.addVariable();
        Model.Variable y = builder.addVariable();
        x.free();
        y.free();
        Model.Objective objective = builder.setObjective();
        objective.add(1.0, x, x);
        objective.add(1.0, y, y);
        objective.add(7.0);
        objective.minimize();
        Model model = builder.build();

        assertEquals(OSQP.Status.SOLVED, model.solve(tight()));
        assertEquals(0.0, model.getSolution(0), 1e-8);
        assertEquals(0.0, model.getSolution(1), 1e-8);
        assertEquals(7.0, model.getObj(), 1e-8);
    }

    @Test
    public void maximizeWithOffsetUsesSameSignConvention() {
        Model.Builder builder = Model.getBuilder();
        Model.Variable x = builder.addVariable();
        x.free();
        Model.Objective objective = builder.setObjective();
        objective.add(-1.0, x, x);
        objective.add(7.0);
        objective.maximize();
        Model model = builder.build();

        assertEquals(OSQP.Status.SOLVED, model.solve(tight()));
        assertEquals(0.0, model.getSolution(0), 1e-8);
        assertEquals(7.0, model.getObj(), 1e-8);
    }

    @Test
    public void constraintRhsAndTypeCanBeChangedBetweenSolves() {
        Model.Builder builder = Model.getBuilder();
        Model.Variable x = builder.addVariable();
        Model.Variable y = builder.addVariable();
        x.free();
        y.free();
        Model.Objective objective = builder.setObjective();
        objective.add(1.0, x, x);
        objective.add(1.0, y, y);
        objective.minimize();
        Model.Constraint constraint = builder.addConstraint();
        constraint.add(1.0, x);
        constraint.add(1.0, y);
        constraint.geq(1.0);
        Model model = builder.build();

        assertEquals(OSQP.Status.SOLVED, model.solve(tight()));
        assertEquals(0.5, model.getSolution(0), 1e-6);
        assertEquals(0.5, model.getSolution(1), 1e-6);

        constraint.setRhs(3.0);
        assertEquals(OSQP.Status.SOLVED, model.solve(tight()));
        assertEquals(1.5, model.getSolution(0), 1e-6);
        assertEquals(1.5, model.getSolution(1), 1e-6);

        constraint.setType(Model.Constraint.Type.EQ);
        assertEquals(OSQP.Status.SOLVED, model.solve(tight()));
        assertEquals(1.5, model.getSolution(0), 1e-6);
        assertEquals(1.5, model.getSolution(1), 1e-6);
    }

    @Test
    public void variableBoundsAreHonored() {
        Model.Builder builder = Model.getBuilder();
        Model.Variable x = builder.addVariable();
        x.lb(0.0).ub(1.0);
        Model.Objective objective = builder.setObjective();
        objective.add(1.0, x, x);
        objective.add(-10.0, x);
        objective.minimize();
        Model model = builder.build();

        assertEquals(OSQP.Status.SOLVED, model.solve(tight()));
        assertEquals(1.0, model.getSolution(0), 1e-6);
        assertEquals(-9.0, model.getObj(), 1e-6);
    }

    @Test
    public void largeModelKeepsPrimitiveStorage() {
        int n = 20_000;
        Model.Builder builder = Model.getBuilder();
        Model.Variable[] x = builder.addVariables(n);
        Model.Objective objective = builder.setObjective();
        for (int i = 0; i < n; i++) {
            objective.add(1.0, x[i], x[i]);
            objective.add(1e-3, x[i]);
        }
        objective.minimize();
        Model.Constraint constraint = builder.addConstraint();
        for (int i = 0; i < n; i++) {
            constraint.add(1.0, x[i]);
        }
        constraint.leq(1.0);
        Model model = builder.build();

        model.solve();
        // n structural bound rows plus the single user constraint
        assertEquals(n + 1, model.getNumConstraints());
        assertEquals(n, model.getNumVariables());
    }

    @Test
    public void bulkApiMatchesAlgebraicBuilder() {
        // min (x-1)^2 + (y-2)^2 + 5  subject to  x + y == 2
        Model.Builder algebraic = Model.getBuilder();
        Model.Variable x = algebraic.addVariable();
        Model.Variable y = algebraic.addVariable();
        x.free();
        y.free();
        Model.Objective objective = algebraic.setObjective();
        objective.add(1.0, x, x);
        objective.add(-2.0, x);
        objective.add(1.0, y, y);
        objective.add(-4.0, y);
        objective.add(5.0);
        objective.minimize();
        Model.Constraint constraint = algebraic.addConstraint();
        constraint.add(1.0, x);
        constraint.add(1.0, y);
        constraint.eq(2.0);
        Model reference = algebraic.build();

        Model.Builder bulk = Model.getBuilder();
        double[] free = {Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY};
        double[] unbounded = {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
        bulk.addVariables(2, free, unbounded);
        bulk.quadraticObjective(new int[]{0, 1}, new int[]{0, 1}, new double[]{1.0, 1.0});
        bulk.linearObjective(new int[]{0, 1}, new double[]{-2.0, -4.0});
        bulk.offset(5.0);
        bulk.minimize();
        bulk.constraint(new int[]{0, 1}, new double[]{1.0, 1.0}, 2.0, 2.0);
        Model model = bulk.build();

        assertEquals(OSQP.Status.SOLVED, reference.solve(tight()));
        assertEquals(OSQP.Status.SOLVED, model.solve(tight()));
        assertEquals(reference.getNumVariables(), model.getNumVariables());
        assertEquals(reference.getNumConstraints(), model.getNumConstraints());
        assertEquals(reference.getSolution(0), model.getSolution(0), 1e-9);
        assertEquals(reference.getSolution(1), model.getSolution(1), 1e-9);
        assertEquals(reference.getObj(), model.getObj(), 1e-9);
        assertEquals(0.5, model.getSolution(0), 1e-8);
        assertEquals(1.5, model.getSolution(1), 1e-8);
        // (x-1)^2 + (y-2)^2 at (0.5, 1.5); the offset 5 is the constant hidden in the squares.
        assertEquals(0.5, model.getObj(), 1e-8);
    }

    @Test
    public void bulkBoundsAndOffsetAreHonored() {
        Model.Builder builder = Model.getBuilder();
        builder.addVariables(2, new double[]{0.0, -1.0}, new double[]{1.0, 1.0});
        builder.quadraticObjective(new int[]{0, 1}, new int[]{0, 1}, new double[]{1.0, 1.0});
        builder.offset(7.0);
        builder.minimize();
        Model model = builder.build();

        assertEquals(OSQP.Status.SOLVED, model.solve(tight()));
        assertEquals(0.0, model.getSolution(0), 1e-8);
        assertEquals(0.0, model.getSolution(1), 1e-8);
        assertEquals(7.0, model.getObj(), 1e-8);
    }

    @Test
    public void bulkConstraintBoundsCanBeUpdated() {
        Model.Builder builder = Model.getBuilder();
        builder.addVariables(2, new double[]{Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY},
                new double[]{Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY});
        builder.quadraticObjective(new int[]{0, 1}, new int[]{0, 1}, new double[]{1.0, 1.0});
        builder.minimize();
        int row = builder.build() == null ? -1 : 0;   // builder keeps the model, row added below
        builder.constraint(new int[]{0, 1}, new double[]{1.0, 1.0}, 1.0, 1.0);
        Model model = builder.build();

        assertEquals(OSQP.Status.SOLVED, model.solve(tight()));
        assertEquals(0.5, model.getSolution(0), 1e-6);
        assertEquals(0.5, model.getSolution(1), 1e-6);

        model.setConstraintBounds(model.getNumConstraints() - 1, 3.0, 3.0);
        assertEquals(OSQP.Status.SOLVED, model.solve());
        assertEquals(1.5, model.getSolution(0), 1e-4);
        assertEquals(1.5, model.getSolution(1), 1e-4);
        assertTrue(row == 0);
    }

    @Test
    public void bulkApiValidatesItsArguments() {
        Model.Builder builder = Model.getBuilder();
        builder.addVariables(2, null, null);
        try {
            builder.linearObjective(new int[]{0, 1}, new double[]{1.0});
            fail("length mismatch must be rejected");
        } catch (IllegalArgumentException expected) {
        }
        try {
            builder.quadraticObjective(new int[]{0, 3}, new int[]{0, 0}, new double[]{1.0, 1.0});
            fail("unknown variable must be rejected");
        } catch (IllegalArgumentException expected) {
        }
        try {
            builder.constraint(new int[]{0}, new double[]{1.0}, 2.0, 1.0);
            fail("inverted bounds must be rejected");
        } catch (IllegalArgumentException expected) {
        }
    }

    @Test
    public void bulkApiBuildsLargeModelWithoutRunnerObjects() {
        int n = 100_000;
        int[] all = new int[n];
        double[] linear = new double[n];
        double[] quadratic = new double[n];
        for (int i = 0; i < n; i++) {
            all[i] = i;
            linear[i] = 1e-3;
            quadratic[i] = 1.0;
        }
        Model.Builder builder = Model.getBuilder();
        builder.addVariables(n, null, null);
        builder.linearObjective(all, linear);
        builder.quadraticObjective(all, all, quadratic);
        builder.minimize();
        Model model = builder.build();

        assertEquals(n, model.getNumVariables());
        assertEquals(n, model.getNumConstraints());
    }

    @Test
    public void offDiagonalTermUsesTheMonomialCoefficient() {
        // min x^2 + y^2 - x*y - x - y is minimal at (1, 1) with value -1. The coefficient passed to
        // add() must be the coefficient of x*y, exactly as it is the coefficient of x^2 on the
        // diagonal; the old behaviour doubled the off-diagonal term.
        Model.Builder builder = Model.getBuilder();
        Model.Variable x = builder.addVariable();
        Model.Variable y = builder.addVariable();
        builder.setObjective().add(1, x, x).add(1, y, y).add(-1, x, y).add(-1, x).add(-1, y)
                .minimize();
        Model model = builder.build();

        assertEquals(OSQP.Status.SOLVED, model.solve(tight()));
        assertEquals(1.0, model.getSolution(0), 1e-5);
        assertEquals(1.0, model.getSolution(1), 1e-5);
        assertEquals(-1.0, model.getObj(), 1e-5);
    }

    @Test
    public void offDiagonalTermIsOrderIndependent() {
        Model.Builder builder = Model.getBuilder();
        Model.Variable x = builder.addVariable();
        Model.Variable y = builder.addVariable();
        builder.setObjective().add(1, x, x).add(1, y, y).add(-1, y, x).add(-1, x).add(-1, y)
                .minimize();
        Model model = builder.build();

        assertEquals(OSQP.Status.SOLVED, model.solve(tight()));
        assertEquals(1.0, model.getSolution(0), 1e-5);
        assertEquals(1.0, model.getSolution(1), 1e-5);
    }

    @Test
    public void bulkQuadraticObjectiveUsesTheMonomialCoefficient() {
        Model.Builder builder = Model.getBuilder();
        builder.addVariables(2, null, null);
        builder.linearObjective(new int[] { 0, 1 }, new double[] { -1, -1 });
        builder.quadraticObjective(new int[] { 0, 1, 0 }, new int[] { 0, 1, 1 },
                new double[] { 1, 1, -1 });
        builder.minimize();
        Model model = builder.build();

        assertEquals(OSQP.Status.SOLVED, model.solve(tight()));
        assertEquals(1.0, model.getSolution(0), 1e-5);
        assertEquals(1.0, model.getSolution(1), 1e-5);
    }
}
