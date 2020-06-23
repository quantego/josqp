package com.quantego.josqp;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import com.quantego.josqp.NonCvxTestDataGenerator;
import com.quantego.josqp.OSQP;
import com.quantego.josqp.OSQP.Settings;

public class NonCvxTest {


    public void test_non_cvx_solve()
    {
        int exitflag = 0;

        // Problem settings
        OSQP.Settings settings = new Settings();

        // Structures
        //OSQP.Workspace work; // Workspace
        OSQP.Data data;      // Data
        NonCvxTestDataGenerator.non_cvx_sols_data sols_data;


        // Populate data
        data = NonCvxTestDataGenerator.generate_problem_non_cvx();
        sols_data = NonCvxTestDataGenerator.generate_problem_non_cvx_sols_data();


        // Define Solver settings as default
        //osqp.osqp_set_default_settings(settings);
        settings.verbose = true;
        settings.sigma = 1e-6;
        OSQP osqp = new OSQP(data,settings);

        // Setup should fail due to (P + sigma I) having a negative eigenvalue
        //assertTrue(exitflag == OSQP.Error.valueOf("OSQP_NONCVX_ERROR").ordinal(),"Non Convex test solve: Setup should have failed!");



        // Update Solver settings
        settings.sigma = sols_data.sigma_new;

        // Setup workspace again
        osqp = new OSQP(data,settings);
        //exitflag = osqp_setup(work, data, settings);

        // Setup should work this time because (P + sigma I) is positive definite
        assertTrue(exitflag == 0,"Non Convex test solve: Setup error!");

        // Solve Problem first time
        OSQP.Status status = osqp.solve();

        //boolean exitone  = osqp.ososqp_update_rho(work,compute_rho_estimate(work));

        // boolean exit = osqp_solve(work);

        // Compare solver statuses'
  /*assertTrue(status.ordinal() == OSQP.Status.valueOf("ERROR").ordinal(),"Non Convex test solve: Error in solver status!"
          );

  // Compare objective values
/* assertTrue(osqp.getWorkspace().info.obj_val == OSQP.Error.valueOf("OSQP_NAN").ordinal(),"Non Convex test solve: Error in objective value!"
   );
*/

        return ;
    }


}