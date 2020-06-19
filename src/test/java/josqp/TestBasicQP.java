package josqp;

import com.quantego.josqp.CSCMatrix;
import com.quantego.josqp.DataGenerator;
import com.quantego.josqp.DataGenerator.basic_qp_sols_data;
import com.quantego.josqp.LinAlg;
import com.quantego.josqp.OSQP;
import com.quantego.josqp.OSQP.Data;
import com.quantego.josqp.OSQP.Settings;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class TestBasicQP {
    final static double TESTS_TOL = 0.0001;
    final static int OSQP_TIME_LIMIT_REACHED = -6;

    @Test
    public void testBasicQPSolve() {

        OSQP osqp;
        int exitFlag, tmpInt;
        float tmpFloat;
        CSCMatrix tmpMtrix, pTmp;

        Settings settings = new Settings();

        Data data = null;
        basic_qp_sols_data sols_data;

        data = DataGenerator.generate_problem_basic_qp();
        sols_data = DataGenerator.generate_problem_basic_qp_sols_data();

        settings.max_iter = 2000;
        settings.alpha = 1.6;
        settings.polish = true;
        settings.scaling = 0;
        settings.verbose = true;
        settings.warm_start = false;
        osqp = new OSQP(data, settings);

        exitFlag =0;// osqp_setup(&work, data, settings);

        assertEquals(0, exitFlag, "Basic QP test solve: Setup error!");

        osqp.solve();

        assertEquals(osqp.getWorkspace().info.status, sols_data.getStatus_test(), "Basic QP test solve: Error in solver status!");

        assertTrue((LinAlg.vec_norm_inf_diff(osqp.getPrimalSolution(), sols_data.getX_test(), data.getN()) < TESTS_TOL), "Basic QP test solve: Error in primal solution!");

        assertTrue((LinAlg.vec_norm_inf_diff(osqp.getDualSolution(), sols_data.getY_test(), data.getM()) < TESTS_TOL), "Basic QP test solve: Error in dual solution!");

        assertTrue((Math.abs(osqp.getWorkspace().info.obj_val - sols_data.getObj_value_test()) < 1 ), "Basic QP test solve: Error in objective value!");

        assertEquals(OSQP.osqp_update_rho(osqp.getWorkspace(), -0.1),"Basic QP test solve: Wrong value of rho not caught!");

        //assertEquals(osqp.osqp_update_max_iter(osqp.getWorkspace(), -1),1,"Basic QP test solve: Wrong value of max_iter not caught!");

        //assertEquals(osqp.osqp_update_eps_abs(osqp.getWorkspace(),-1.),1,"Basic QP test solve: Wrong value of eps_abs not caught!");

        //assertEquals(osqp.osqp_update_eps_rel(osqp.getWorkspace(),-1.),1,"Basic QP test solve: Wrong value of eps_rel not caught!");

        //assertEquals(osqp.osqp_update_eps_prim_inf(osqp.getWorkspace(),-0.1),1,"Basic QP test solve: Wrong value of eps_prim_inf not caught!");

        //assertEquals(osqp.osqp_update_eps_dual_inf(osqp.getWorkspace(),-0.1),1,"Basic QP test solve: Wrong value of eps_dual_inf not caught!");

        //assertEquals(osqp.osqp_update_alpha(osqp.getWorkspace(),2.0),1,"Basic QP test solve: Wrong value of alpha not caught!");

        //assertEquals(osqp.osqp_update_warm_start(osqp.getWorkspace(), -1), 1, "Basic QP test solve: Wrong value of warm_start not caught!");

        //assertEquals(osqp.osqp_update_scaled_termination(osqp.getWorkspace(), 2), 1, "Basic QP test solve: Wrong value of scaled_termination not caught!");

        //assertEquals(osqp.osqp_update_check_termination(osqp.getWorkspace(), -1), 1, "Basic QP test solve: Wrong value of check_termination not caught!");

        //assertEquals(osqp.osqp_update_delta(osqp.getWorkspace(), 0.), 1, "Basic QP test solve: Wrong value of delta not caught!");

        //assertEquals(osqp.osqp_update_polish(osqp.getWorkspace(), 2), 1, "Basic QP test solve: Wrong value of polish not caught!");

        //assertEquals(osqp.osqp_update_polish_refine_iter(osqp.getWorkspace(), -1), 1, "Basic QP test solve: Wrong value of polish_refine_iter not caught!");

        //assertEquals(osqp.osqp_update_verbose(osqp.getWorkspace(), 2), 1, "Basic QP test solve: Wrong value of verbose not caught!");

        //osqp_cleanup(osqp.getWorkspace());

    }

    @Test
    public void test_basic_qp_solve_pardiso() {

        int exitFlag;
        Settings settings = new Settings();
        Data data;
        basic_qp_sols_data qp_sols_data;

        data = DataGenerator.generate_problem_basic_qp();
        qp_sols_data = DataGenerator.generate_problem_basic_qp_sols_data();

        settings.max_iter = 2000;
        settings.alpha = 1.6;
        settings.polish = true;
        settings.scaling = 0;
        settings.verbose = true;
        settings.warm_start = false;

    }


    @Test
    public void test_basic_qp_warm_start() {
        int iter;
        double[] x0 = { 0.0f, 0.0f, };
        double[] y0 = { 0.0f, 0.0f, 0.0f, 0.0f, };

        double[] xopt = { 0.3f, 0.7f, };
        double[] yopt = { -2.9f, 0.0f, 0.2f, 0.0f, };

        Settings settings = new Settings();
        Data data = DataGenerator.generate_problem_basic_qp();
        basic_qp_sols_data sols_data = DataGenerator.generate_problem_basic_qp_sols_data();
        settings.check_termination = 1;
        OSQP osqp = new OSQP(data, settings);
        osqp.solve();
        iter = osqp.getWorkspace().info.iter;
        osqp.warm_start_x(x0);
        osqp.warm_start_y(y0);
        osqp.solve();

        assertEquals(osqp.getWorkspace().info.iter,iter,"Basic QP test warm start: Cold start error!");

        osqp.warm_start_x(xopt);
        osqp.warm_start_y(yopt);
        osqp.solve();

        assertEquals(osqp.getWorkspace().info.iter,1,"Basic QP test warm start: Warm start error!");
    }

    @Test
    public void test_basic_qp_check_termination() {
        int exitFlag;
        Settings settings = new Settings();
        Data data = DataGenerator.generate_problem_basic_qp();
        basic_qp_sols_data sols_data = DataGenerator.generate_problem_basic_qp_sols_data();

        settings.max_iter = 200;
        settings.alpha = 1.6;
        settings.polish = false;
        settings.scaling = 0;
        settings.verbose = true;
        settings.check_termination = 0;
        settings.warm_start = false;
        OSQP osqp = new OSQP(data, settings);

        exitFlag =0;// osqp_setup(&work, data, settings);
        assertEquals(exitFlag,0,"Basic QP test solve: Setup error!");
        osqp.solve();

        assertEquals( osqp.getWorkspace().info.iter, settings.max_iter, "Basic QP test check termination: Error in number of iterations taken!");

        assertEquals(osqp.getWorkspace().info.status, sols_data.getStatus_test(), "Basic QP test check termination: Error in solver status!");

        assertTrue(LinAlg.vec_norm_inf_diff(osqp.getPrimalSolution(), sols_data.getX_test(), data.getN()) < TESTS_TOL , "Basic QP test check termination: Error in primal solution!");

        assertTrue(LinAlg.vec_norm_inf_diff(osqp.getDualSolution(), sols_data.getY_test(), data.getM()) < TESTS_TOL , "Basic QP test check termination: Error in dual solution!");

        assertTrue(Math.abs(osqp.getWorkspace().info.obj_val - sols_data.getObj_value_test()) < TESTS_TOL , "Basic QP test check termination: Error in objective value!");

    }

    @Test
    public void test_basic_qp_time_limit() {
        int exitFlag;
        Settings settings = new Settings();
        Data data = DataGenerator.generate_problem_basic_qp();
        basic_qp_sols_data sols_data = DataGenerator.generate_problem_basic_qp_sols_data();

        settings.rho = 20;
        settings.adaptive_rho = false;

        assertEquals(settings.time_limit, 0, "Basic QP test time limit: Default not correct");

        OSQP osqp = new OSQP(data, settings);

        exitFlag =0;// osqp_setup(&work, data, settings);
        assertEquals(exitFlag,0,"Basic QP test time limit: Setup error!");
        osqp.solve();

        assertEquals(osqp.getWorkspace().info.status, sols_data.getStatus_test(), "Basic QP test time limit: Error in no time limit solver status!");

        //osqp.osqp_update_time_limit(1e-5);
        //osqp.osqp_update_eps_rel(1e-09);
        //osqp.osqp_update_eps_abs(1e-09);

        //osqp.osqp_update_time_limit(1e-7);
        //osqp.osqp_update_eps_rel(1e-12);
        //osqp.osqp_update_eps_abs(1e-12);

        //osqp.osqp_update_max_iter((int)2e9);
        //osqp.osqp_update_check_termination(0);

        OSQP.cold_start(osqp.getWorkspace());
        osqp.solve();


        assertEquals(osqp.getWorkspace().info.status, OSQP_TIME_LIMIT_REACHED, "Basic QP test time limit: Error in timed out solver status!");
    }

    public void test_basic_qp_update_rho() {
    }





}