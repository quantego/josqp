package com.quantego.josqp;

import org.junit.Test;

import java.lang.reflect.Field;
import java.util.LinkedHashMap;
import java.util.Map;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

/**
 * Keeps {@link OSQP.Settings} compatible with the public {@code OSQPSettings} structure of
 * OSQP 1.0.0 (include/public/osqp_api_types.h, tag v1.0.0).
 *
 * <p>Every setting of the C API must either exist here (under the same name or a documented
 * alias) or be listed in {@link #UNSUPPORTED} with the reason. That makes the gap explicit and
 * forces a decision whenever upstream adds a parameter.</p>
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

    /** OSQP 1.0 names that JOSQP deliberately does not implement, with the reason. */
    private static final Map<String, String> UNSUPPORTED = new LinkedHashMap<>();

    static {
        UNSUPPORTED.put("device", "pure Java solver without GPU backend");
        UNSUPPORTED.put("linsys_solver", "only the direct QDLDL solver is implemented");
        UNSUPPORTED.put("allocate_solution", "solutions are always allocated by the Java API");
        UNSUPPORTED.put("profiler_level", "no profiler annotations in the Java port");
        UNSUPPORTED.put("rho_is_vec", "rho is always the internal scalar step size here");
        UNSUPPORTED.put("cg_max_iter", "conjugate-gradient solver not implemented");
        UNSUPPORTED.put("cg_tol_reduction", "conjugate-gradient solver not implemented");
        UNSUPPORTED.put("cg_tol_fraction", "conjugate-gradient solver not implemented");
        UNSUPPORTED.put("cg_precond", "conjugate-gradient solver not implemented");
    }

    /** OSQP 1.0 names that JOSQP exposes under the older 0.6 name. */
    private static final Map<String, String> ALIASES = new LinkedHashMap<>();

    static {
        ALIASES.put("warm_starting", "warm_start");
        ALIASES.put("polishing", "polish");
    }

    private static boolean hasField(String name) {
        for (Field field : OSQP.Settings.class.getFields()) {
            if (field.getName().equals(name))
                return true;
        }
        return false;
    }

    @Test
    public void everyOsqpSettingIsSupportedOrDocumented() {
        for (String name : OSQP_1_0_SETTINGS) {
            String alias = ALIASES.get(name);
            boolean present = hasField(name) || (alias != null && hasField(alias));
            assertTrue("OSQP 1.0 setting '" + name + "' is neither implemented nor listed as unsupported",
                    present || UNSUPPORTED.containsKey(name));
        }
    }

    @Test
    public void supportedSettingsUseTheOsqpNameOrDefaultValue() throws Exception {
        OSQP.Settings settings = new OSQP.Settings();
        assertTrue("check_dualgap must default to true like OSQP 1.0 doubles",
                OSQP.Settings.class.getField("check_dualgap").getBoolean(settings));
        assertFalse("OSQP 1.0 defaults polishing to off",
                OSQP.Settings.class.getField("polish").getBoolean(settings));
        assertFalse("the optional KKT diagnostics must stay off by default",
                OSQP.Settings.class.getField("check_linear_residual").getBoolean(settings));
        assertFalse("the optional rho progress gate must stay off by default",
                OSQP.Settings.class.getField("adaptive_rho_kkt").getBoolean(settings));
    }

    /** Defaults of the OSQP 1.0 double/CPU build, pinned so a change is always deliberate. */
    @Test
    public void defaultsMatchOsqp100() {
        OSQP.Settings settings = new OSQP.Settings();
        assertTrue("adaptive_rho_interval must be 50 as in OSQP 1.0",
                settings.adaptive_rho_interval == 50);
        assertTrue("rho_eq_over_rho_ineq must be 1e3 as in OSQP 1.0",
                OSQP.RHO_EQ_OVER_RHO_INEQ == 1e3);
        assertTrue("time_limit must default to 1e10 as in OSQP 1.0",
                settings.time_limit == 1e10);
        assertTrue("adaptive_rho must default to true", settings.adaptive_rho);
        assertTrue("adaptive_rho_tolerance must default to 5.0",
                settings.adaptive_rho_tolerance == 5.0);
        assertTrue("warm_start must default to true", settings.warm_start);
        assertTrue("check_termination must default to 25", settings.check_termination == 25);
        assertTrue("max_iter must default to 4000", settings.max_iter == 4000);
        assertTrue("delta must default to 1e-6", settings.delta == 1e-6);
        assertTrue("polish_refine_iter must default to 3", settings.polish_refine_iter == 3);
        assertTrue("eps_abs must default to 1e-3", settings.eps_abs == 1e-3);
        assertTrue("eps_prim_inf must default to 1e-4", settings.eps_prim_inf == 1e-4);
        assertTrue("scaling must default to 10", settings.scaling == 10);
    }

    @Test
    public void documentedUnsupportedSettingsAreAbsent() {
        for (String name : UNSUPPORTED.keySet())
            assertFalse("'" + name + "' is implemented now; remove it from UNSUPPORTED", hasField(name));
    }
}
