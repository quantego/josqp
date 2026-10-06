package com.quantego.josqp;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

/**
 * Repeated solves after changing matrix values.
 *
 * <p>The update path passed CSC row indices where element positions are expected. That is only
 * harmless while the pattern happens to run along the diagonal; otherwise the values land in the
 * wrong entries (or outside the array) without an error. A changed sparsity pattern cannot be applied
 * in place at all and needs a fresh symbolic factorisation. These tests compare a solver reused
 * across updates against a freshly built model for the same final data.
 */
public class ModelUpdateTest {

    private static OSQP.Settings tight() {
        OSQP.Settings settings = new OSQP.Settings();
        settings.eps_abs = 1e-9;
        settings.eps_rel = 1e-9;
        return settings;
    }

    private static void assertSameSolution(Model fresh, Model reused) {
        for (int i = 0; i < fresh.getNumVariables(); i++)
            assertEquals("variable " + i, fresh.getSolution(i), reused.getSolution(i), 1e-5);
        assertEquals("objective", fresh.getObj(), reused.getObj(), 1e-5);
    }

    @Test
    public void changedConstraintCoefficientMatchesFreshModel() {
        // A holds entries in rows 0 and 2 for x0 and rows 1 and 3 for x1, so the CSC row indices
        // differ from the element positions.
        Model.Builder builder = Model.getBuilder();
        Model.Variable x0 = builder.addVariable().lb(0);
        Model.Variable x1 = builder.addVariable().lb(0);
        builder.setObjective().add(1, x0).add(2, x1).minimize();
        Model.Constraint c0 = builder.addConstraint();
        c0.add(1, x1);
        c0.geq(5);
        Model.Constraint c1 = builder.addConstraint();
        c1.add(1, x0);
        c1.geq(5);
        Model reused = builder.build();
        assertEquals(OSQP.Status.SOLVED, reused.solve(tight()));
        assertEquals(5.0, reused.getSolution(0), 1e-5);
        assertEquals(5.0, reused.getSolution(1), 1e-5);

        c0.setLhs(x1, 2.0);
        assertEquals(OSQP.Status.SOLVED, reused.solve(tight()));

        Model.Builder oracleBuilder = Model.getBuilder();
        Model.Variable y0 = oracleBuilder.addVariable().lb(0);
        Model.Variable y1 = oracleBuilder.addVariable().lb(0);
        oracleBuilder.setObjective().add(1, y0).add(2, y1).minimize();
        Model.Constraint d0 = oracleBuilder.addConstraint();
        d0.add(2, y1);
        d0.geq(5);
        Model.Constraint d1 = oracleBuilder.addConstraint();
        d1.add(1, y0);
        d1.geq(5);
        Model oracle = oracleBuilder.build();
        assertEquals(OSQP.Status.SOLVED, oracle.solve(tight()));

        assertEquals(5.0, reused.getSolution(0), 1e-5);
        assertEquals(2.5, reused.getSolution(1), 1e-5);
        assertSameSolution(oracle, reused);
    }

    @Test
    public void changedCrossTermMatchesFreshModel() {
        Model.Builder builder = Model.getBuilder();
        Model.Variable x = builder.addVariable();
        Model.Variable y = builder.addVariable();
        Model.Objective objective = builder.setObjective();
        objective.add(1, x, x).add(1, y, y).add(-1, x, y).add(-1, x).add(-1, y).minimize();
        Model reused = builder.build();
        assertEquals(OSQP.Status.SOLVED, reused.solve(tight()));
        assertEquals(1.0, reused.getSolution(0), 1e-5);

        objective.set(-0.25, x, y);
        assertEquals(OSQP.Status.SOLVED, reused.solve(tight()));

        Model.Builder oracleBuilder = Model.getBuilder();
        Model.Variable u = oracleBuilder.addVariable();
        Model.Variable v = oracleBuilder.addVariable();
        oracleBuilder.setObjective().add(1, u, u).add(1, v, v).add(-0.25, u, v)
                .add(-1, u).add(-1, v).minimize();
        Model oracle = oracleBuilder.build();
        assertEquals(OSQP.Status.SOLVED, oracle.solve(tight()));

        assertEquals(4.0 / 7.0, reused.getSolution(0), 1e-5);
        assertEquals(4.0 / 7.0, reused.getSolution(1), 1e-5);
        assertSameSolution(oracle, reused);
    }

    @Test
    public void addedNonzeroRebuildsTheSolver() {
        // x0 >= 4, minimising 5*x0 + x1 -> (4, 0). Adding the x1 coefficient turns the row into
        // x0 + x1 >= 4, which moves the optimum to (0, 4): the pattern grows and needs a rebuild.
        Model.Builder builder = Model.getBuilder();
        Model.Variable x0 = builder.addVariable().lb(0);
        Model.Variable x1 = builder.addVariable().lb(0);
        builder.setObjective().add(5, x0).add(1, x1).minimize();
        Model.Constraint c0 = builder.addConstraint();
        c0.add(1, x0);
        c0.geq(4);
        Model reused = builder.build();
        assertEquals(OSQP.Status.SOLVED, reused.solve(tight()));
        assertEquals(4.0, reused.getSolution(0), 1e-5);
        assertEquals(0.0, reused.getSolution(1), 1e-5);

        c0.setLhs(x1, 1.0);
        assertEquals(OSQP.Status.SOLVED, reused.solve(tight()));

        Model.Builder oracleBuilder = Model.getBuilder();
        Model.Variable y0 = oracleBuilder.addVariable().lb(0);
        Model.Variable y1 = oracleBuilder.addVariable().lb(0);
        oracleBuilder.setObjective().add(5, y0).add(1, y1).minimize();
        Model.Constraint d0 = oracleBuilder.addConstraint();
        d0.add(1, y0).add(1, y1);
        d0.geq(4);
        Model oracle = oracleBuilder.build();
        assertEquals(OSQP.Status.SOLVED, oracle.solve(tight()));

        assertEquals(0.0, reused.getSolution(0), 1e-5);
        assertEquals(4.0, reused.getSolution(1), 1e-5);
        assertSameSolution(oracle, reused);
    }
}
