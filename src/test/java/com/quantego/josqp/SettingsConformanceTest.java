package com.quantego.josqp;

import org.junit.jupiter.api.Test;

import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.LinkedHashMap;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Keeps {@link OSQP.Settings} compatible with the public {@code OSQPSettings} structure of
 * OSQP 1.0.0 (include/public/osqp_api_types.h, tag v1.0.0).
 *
 * <p>Every setting of the C API must exist here, either under its own name, under a documented
 * alias, or - for the four OSQP 1.0 names that jOSQP keeps as accessor methods - as a bridge to the
 * 0.6 field name. {@link #INERT} lists the settings that exist for API parity but cannot change the
 * behaviour of a pure Java solver; settings that would silently produce a different result than
 * requested are rejected instead ({@link #rejectedSettingsFailLoudly()}).</p>
 */
public class SettingsConformanceTest {

    /** Field names of OSQPSettings in OSQP 1.0.0, in the order of the header. */
    private static final String[] OSQP_1_0_SETTINGS = {
            "device", "linsys_solver",
            "allocate_solution", "verbose", "profiler_level", "warm_starting", "scaling", "polishing",
            "rho", "rho_is_vec", "sigma", "alpha",
            "cg_max_iter", "cg_tol_reduction", "cg_tol_fraction", "cg_precond",
            "adaptive_rho", "adaptive_rho_interval", "adaptive_rho_fraction", "adaptive_rho_tolerance",
            "max_iter", "eps_abs", "eps_rel", "eps_prim_inf", "eps_dual_inf",
            "scaled_termination", "check_termination", "check_dualgap", "time_limit",
            "delta", "polish_refine_iter"
    };

    /** OSQP 1.0 names that jOSQP exposes under the older 0.6 field name. */
    private static final Map<String, String> ALIASES = new LinkedHashMap<>();

    static {
        ALIASES.put("warm_starting", "warm_start");
        ALIASES.put("polishing", "polish");
    }

    /** OSQP 1.0 settings that exist for parity but cannot affect a solve here. */
    private static final Map<String, String> INERT = new LinkedHashMap<>();

    static {
        INERT.put("device", "pure Java solver without device selection");
        INERT.put("allocate_solution", "the Java API always allocates the solution vectors");
        INERT.put("profiler_level", "no profiler annotations in the Java port");
        INERT.put("cg_max_iter", "conjugate-gradient solver not implemented");
        INERT.put("cg_tol_reduction", "conjugate-gradient solver not implemented");
        INERT.put("cg_tol_fraction", "conjugate-gradient solver not implemented");
        INERT.put("cg_precond", "conjugate-gradient solver not implemented");
    }

    private static boolean hasField(String name) {
        for (Field field : OSQP.Settings.class.getFields()) {
            if (field.getName().equals(name))
                return true;
        }
        return false;
    }

    private static boolean hasMethod(String name, Class<?>... parameters) {
        try {
            Method method = OSQP.Settings.class.getMethod(name, parameters);
            return !java.lang.reflect.Modifier.isStatic(method.getModifiers());
        } catch (NoSuchMethodException ex) {
            return false;
        }
    }

    @Test
    public void everyOsqpSettingExistsUnderItsNameOrAlias() {
        for (String name : OSQP_1_0_SETTINGS) {
            String alias = ALIASES.get(name);
            assertTrue(hasField(name) || (alias != null && hasField(alias)), "OSQP 1.0 setting '" + name + "' is missing");
        }
    }

    /** The renamed entries are reachable under both names. */
    @Test
    public void renamedSettingsProvideAccessorBridges() {
        assertTrue(hasMethod("isPolishing"));
        assertTrue(hasMethod("setPolishing", boolean.class));
        assertTrue(hasMethod("isWarmStarting"));
        assertTrue(hasMethod("setWarmStarting", boolean.class));

        OSQP.Settings settings = new OSQP.Settings();
        settings.setPolishing(true);
        assertTrue(settings.polish, "setPolishing must drive the polish field");
        assertTrue(settings.isPolishing());
        settings.polish = false;
        assertFalse(settings.isPolishing());

        settings.setWarmStarting(false);
        assertFalse(settings.warm_start, "setWarmStarting must drive the warm_start field");
        assertFalse(settings.isWarmStarting());
        settings.warm_start = true;
        assertTrue(settings.isWarmStarting());
    }

    @Test
    public void settingsThatCannotBeHonouredAreRejected() {
        Model.Builder builder = Model.getBuilder();
        builder.addVariable();
        Model model = builder.build();

        OSQP.Settings indirect = new OSQP.Settings();
        indirect.linsys_solver = OSQP.LinSysSolver.OSQP_INDIRECT_SOLVER;
        try {
            model.solve(indirect);
            fail("the conjugate-gradient solver must be rejected");
        } catch (IllegalArgumentException expected) {
        }

        OSQP.Settings scalarRho = new OSQP.Settings();
        scalarRho.rho_is_vec = false;
        try {
            model.solve(scalarRho);
            fail("scalar rho must be rejected");
        } catch (IllegalArgumentException expected) {
        }
    }

    @Test
    public void inertSettingsArePresentWithTheOsqpDefaults() {
        OSQP.Settings settings = new OSQP.Settings();
        for (String name : INERT.keySet())
            assertTrue(hasField(name), "'" + name + "' is listed as inert but does not exist");
        assertTrue(settings.device.isEmpty(), "device must default to the CPU device");
        assertTrue(settings.allocate_solution, "solutions must be allocated");
        assertTrue(settings.profiler_level == 0, "no profiler level by default");
        assertTrue(settings.cg_max_iter == OSQP.CG_MAX_ITER, "the CG defaults must match OSQP 1.0");
        assertTrue(settings.cg_tol_reduction == OSQP.CG_TOL_REDUCTION);
        assertTrue(settings.cg_tol_fraction == OSQP.CG_TOL_FRACTION);
        assertTrue(settings.cg_precond == OSQP.Preconditioner.OSQP_DIAGONAL_PRECONDITIONER);
    }

    @Test
    public void supportedSettingsUseTheOsqpNameOrDefaultValue() {
        OSQP.Settings settings = new OSQP.Settings();
        assertTrue(settings.check_dualgap, "check_dualgap must default to true like OSQP 1.0 doubles");
        assertFalse(settings.polish, "OSQP 1.0 defaults polishing to off");
        assertTrue(
                settings.linsys_solver == OSQP.LinSysSolver.OSQP_DIRECT_SOLVER,
                "linsys_solver must default to the direct QDLDL solver");
        assertTrue(settings.rho_is_vec, "rho_is_vec must default to true (the builtin OSQP 1.0 default)");
        assertFalse(settings.check_linear_residual, "the optional KKT diagnostics must stay off by default");
        assertFalse(settings.adaptive_rho_kkt, "the optional rho progress gate must stay off by default");
    }

    /** Defaults of the OSQP 1.0 double/CPU build, pinned so a change is always deliberate. */
    @Test
    public void defaultsMatchOsqp100() {
        OSQP.Settings settings = new OSQP.Settings();
        assertEquals(50, settings.adaptive_rho_interval);
        assertEquals(1e3, OSQP.RHO_EQ_OVER_RHO_INEQ, 0.0);
        assertEquals(1e10, settings.time_limit, 0.0);
        assertTrue(settings.adaptive_rho);
        assertEquals(5.0, settings.adaptive_rho_tolerance, 0.0);
        assertTrue(settings.warm_start);
        assertEquals(25, settings.check_termination);
        assertEquals(4000, settings.max_iter);
        assertEquals(1e-6, settings.delta, 0.0);
        assertEquals(3, settings.polish_refine_iter);
        assertEquals(1e-3, settings.eps_abs, 0.0);
        assertEquals(1e-4, settings.eps_prim_inf, 0.0);
        assertEquals(10, settings.scaling);
    }
}
