package com.quantego.josqp;

import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

public class RegressionTest {
    @Test public void boundsWithAndWithoutNamesHandleValuelessTypes() throws Exception {
        for (String boundName : new String[]{"", "bounds "}) {
            Path file = Files.createTempFile("josqp-bounds-", ".qps");
            try {
                String text = "NAME bounds\nROWS\n N obj\nCOLUMNS\n"
                        + " x obj 1\n y obj 2\n z obj 3\nBOUNDS\n"
                        + " PL " + boundName + "x\n MI " + boundName + "y\n"
                        + " FR " + boundName + "z\n UP " + boundName + "x 7\nENDATA\n";
                Files.write(file, text.getBytes(StandardCharsets.UTF_8));
                OSQP.Data d = Parser.readQmps(file.toString()).getData();
                assertEquals(3, d.n);
                assertArrayEquals(new double[]{0, -OSQP.OSQP_INFTY, -OSQP.OSQP_INFTY}, d.l, 0);
                assertArrayEquals(new double[]{7, OSQP.OSQP_INFTY, OSQP.OSQP_INFTY}, d.u, 0);
            } finally { Files.deleteIfExists(file); }
        }
    }
    private OSQP solver() {
        OSQP.Settings s = new OSQP.Settings();
        s.scaling = 0;
        s.verbose = false;
        s.polish = false;
        s.adaptive_rho = false;
        CSCMatrix identity = new CSCMatrix(1, 1, 1, new int[]{0, 1}, new int[]{0}, new double[]{1});
        return new OSQP(new OSQP.Data(1, 1, identity, identity, new double[]{-2},
                new double[]{0}, new double[]{1}, 0), s);
    }

    @Test public void inPlaceSolveMatchesPublicCopyAndDoesNotMutateInput() {
        QDLDL q = new QDLDL();
        q.factor(new CSCMatrix(2, 2, 3, new int[]{0, 1, 3}, new int[]{0, 0, 1}, new double[]{4, 1, 3}));
        double[] b = {1, 2};
        double[] copy = q.solve(b);
        assertArrayEquals(new double[]{1, 2}, b, 0);
        double[] inPlace = b.clone();
        QDLDL.solveInPlace(q._ldl, inPlace);
        assertArrayEquals(copy, inPlace, 0);
        assertEquals(1, 4 * inPlace[0] + inPlace[1], 1e-14);
        assertEquals(2, inPlace[0] + 3 * inPlace[1], 1e-14);
        assertThrows(IllegalArgumentException.class, () -> QDLDL.solveInPlace(q._ldl, new double[1]));
    }

    @Test public void productDoesNotSwallowInvalidInput() {
        assertThrows(IndexOutOfBoundsException.class,
                () -> LinAlg.vec_ew_prod(new double[1], new double[1], new double[0], 1));
    }

    @Test public void missingInputFilesFailInsteadOfReturningEmptyProblems() {
        String missing = "target/does-not-exist-regression.mps";
        assertThrows(UncheckedIOException.class, () -> Parser.readQmps(missing));
        assertThrows(UncheckedIOException.class, () -> Utils.readDoubleColumn(missing, ",", 0, 0));
        assertThrows(UncheckedIOException.class, () -> Utils.readIntColumn(missing, ",", 0, 0));
    }

    @Test public void mpsNamesAreLiteralAndRemainDistinct() throws Exception {
        Path file = Files.createTempFile("josqp-names-", ".mps");
        try {
            String text = "NAME names\nROWS\n N obj[0]\n L c.1\n L c!_1\n"
                    + "COLUMNS\n x obj[0] 2 c.1 3\n x c!_1 4\n"
                    + "RHS\n rhs c.1 5 c!_1 6\nENDATA\n";
            Files.write(file, text.getBytes(StandardCharsets.UTF_8));
            OSQP.Data d = Parser.readQmps(file.toString()).getData();
            assertEquals(1, d.n);
            assertEquals(3, d.m); // two constraints and the generated variable bound
            assertEquals(2, d.q[0], 0);
            assertArrayEquals(new double[]{3, 4, 1}, d.A.Ax, 0);
            assertArrayEquals(new double[]{5, 6, OSQP.OSQP_INFTY}, d.u, 0);
        } finally { Files.deleteIfExists(file); }
    }

    @Test public void malformedQuadraticEntryDoesNotSilentlyChangeModel() throws Exception {
        Path file = Files.createTempFile("josqp-quad-", ".qps");
        try {
            String text = "NAME quad\nROWS\n N obj\n L c\nCOLUMNS\n x obj 1 c 1\n"
                    + "RHS\n rhs c 2\nQUADOBJ\n x unknown 1\nENDATA\n";
            Files.write(file, text.getBytes(StandardCharsets.UTF_8));
            assertThrows(IllegalArgumentException.class, () -> Parser.readQmps(file.toString()));
        } finally { Files.deleteIfExists(file); }
    }

    @Test public void malformedInputClosesItsFileDescriptor() throws Exception {
        java.lang.management.OperatingSystemMXBean bean = java.lang.management.ManagementFactory.getOperatingSystemMXBean();
        Class<?> unix;
        try { unix = Class.forName("com.sun.management.UnixOperatingSystemMXBean"); }
        catch (ClassNotFoundException e) { org.junit.jupiter.api.Assumptions.abort("UnixOperatingSystemMXBean unavailable: " + e); return; }
        org.junit.jupiter.api.Assumptions.assumeTrue(unix.isInstance(bean));
        java.lang.reflect.Method count = unix.getMethod("getOpenFileDescriptorCount");
        Path file = Files.createTempFile("josqp-descriptor-", ".mps");
        try {
            Files.write(file, "NOT_A_SECTION\n".getBytes(StandardCharsets.UTF_8));
            long before = (Long) count.invoke(bean);
            for (int i = 0; i < 100; i++)
                assertThrows(IllegalStateException.class, () -> Parser.readQmps(file.toString()));
            assertTrue((Long) count.invoke(bean) <= before + 2, "File descriptors leaked");
        } finally { Files.deleteIfExists(file); }
    }

    @Test public void quietTimeLimitWorksWithoutTerminationChecks() {
        OSQP q = solver();
        q.work.settings.time_limit = 1e-12;
        q.work.settings.check_termination = 0;
        q.work.settings.max_iter = 100000;
        assertEquals(OSQP.Status.TIME_LIMIT_REACHED, q.solve());
        assertEquals(1, q.work.info.iter);
        assertTrue(q.work.info.solve_time > 0);
        assertTrue(q.work.info.run_time >= q.work.info.solve_time);
    }

    @Test public void zeroTimeLimitIsDisabled() {
        OSQP q = solver();
        q.work.settings.time_limit = 0;
        assertEquals(OSQP.Status.SOLVED, q.solve());
    }

    @Test public void repeatedSolveDoesNotRetainOldSolvedStatus() {
        OSQP q = solver();
        assertEquals(OSQP.Status.SOLVED, q.solve());
        q.work.settings.warm_start = false;
        q.work.settings.max_iter = 1;
        q.work.settings.eps_abs = 1e-15;
        q.work.settings.eps_rel = 1e-15;
        assertEquals(OSQP.Status.MAX_ITER_REACHED, q.solve());
    }

    @Test public void infiniteBoundsDoNotDisablePrimalInfeasibilityDetection() {
        OSQP.Settings s = new OSQP.Settings();
        s.scaling = 0; s.verbose = false; s.adaptive_rho = false;
        CSCMatrix p = new CSCMatrix(1, 1, 1, new int[]{0, 1}, new int[]{0}, new double[]{1});
        CSCMatrix a = new CSCMatrix(2, 1, 2, new int[]{0, 2}, new int[]{0, 1}, new double[]{1, 1});
        OSQP q = new OSQP(new OSQP.Data(1, 2, p, a, new double[]{0},
                new double[]{Double.NEGATIVE_INFINITY, 0}, new double[]{-1, Double.POSITIVE_INFINITY}, 0), s);
        assertEquals(OSQP.Status.PRIMAL_INFEASIBLE, q.solve());
    }

    @Test public void updatingBothBoundsUpdatesConstraintPenalty() {
        OSQP q = solver();
        double oldRho = q.work.rho_vec[0];
        q.update_bounds(new double[]{0.5}, new double[]{0.5});
        assertEquals(OSQP.RHO_EQ_OVER_RHO_INEQ * oldRho, q.work.rho_vec[0], 0);
        assertEquals(1 / q.work.rho_vec[0], q.work.linsys_solver.rho_inv_vec[0], 0);
        q.update_bounds(new double[]{0}, new double[]{1});
        assertEquals(oldRho, q.work.rho_vec[0], 0);
        assertEquals(OSQP.Status.SOLVED, q.solve());
    }

    @Test public void rejectedBoundsLeaveSolverUsable() {
        for (int scaling : new int[]{0, 10}) {
            OSQP q = solver();
            if (scaling > 0) {
                OSQP.Settings s = new OSQP.Settings();
                s.scaling = scaling; s.verbose = false;
                q = new OSQP(new OSQP.Data(q.work.data), s);
            }
            final OSQP solver = q;
            double[] lower = q.work.data.l.clone(), upper = q.work.data.u.clone();
            assertThrows(IllegalStateException.class, () -> solver.update_lower_bound(new double[]{2}));
            assertArrayEquals(lower, q.work.data.l, 0);
            assertThrows(IllegalStateException.class, () -> solver.update_upper_bound(new double[]{-1}));
            assertArrayEquals(upper, q.work.data.u, 0);
            assertThrows(IllegalArgumentException.class, () -> solver.update_bounds(new double[0], new double[]{1}));
            assertArrayEquals(lower, q.work.data.l, 0);
            assertArrayEquals(upper, q.work.data.u, 0);
            assertEquals(OSQP.Status.SOLVED, q.solve());
        }
    }
}
