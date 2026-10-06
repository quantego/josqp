package com.quantego.josqp;

import org.junit.Test;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;


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
}
