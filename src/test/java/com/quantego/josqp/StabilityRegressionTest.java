package com.quantego.josqp;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.Test;
import static org.junit.Assert.*;

public class StabilityRegressionTest {
    private OSQP.Data read(String rows, String rhs, String ranges) throws Exception {
        Path file=Files.createTempFile("josqp-ranges-", ".qps");
        try {
            Files.write(file,("NAME range\nROWS\n N obj\n"+rows+"COLUMNS\n x obj 1\nRHS\n"+rhs
                    +"RANGES\n"+ranges+"ENDATA\n").getBytes(StandardCharsets.UTF_8));
            return Parser.readQmps(file.toString()).getData();
        } finally { Files.deleteIfExists(file); }
    }
    @Test public void bothRangesOnOneLineAreApplied() throws Exception {
        OSQP.Data d=read(" G first\n G second\n", " rhs first 5 second 4\n", " range first 7 second 7\n");
        assertArrayEquals(new double[]{5,4,0},d.l,0);
        assertArrayEquals(new double[]{12,11,OSQP.OSQP_INFTY},d.u,0);
    }
    @Test public void signedRangesRespectOriginalRowSense() throws Exception {
        OSQP.Data d=read(" L less\n G greater\n E pos\n E neg\n",
                " rhs less 10 greater 10\n rhs pos 10 neg 10\n",
                " range less -3 greater -3\n range pos 3 neg -3\n");
        assertArrayEquals(new double[]{7,10,10,7,0},d.l,0);
        assertArrayEquals(new double[]{10,13,13,10,OSQP.OSQP_INFTY},d.u,0);
    }
    @Test public void interleavedTripletsAreSortedWithCorrectMapping() {
        int[] map=new int[4];
        CSCMatrix c=CSCMatrix.triplet_to_csc(3,2,4,new int[]{2,2,0,1},new int[]{1,0,1,0},
                new double[]{20,21,10,11},map);
        assertArrayEquals(new int[]{0,2,4},c.Ap);
        assertArrayEquals(new int[]{1,2,0,2},c.Ai);
        assertArrayEquals(new double[]{11,21,10,20},c.Ax,0);
        assertArrayEquals(new int[]{3,1,2,0},map);
    }
    @Test public void polishResidualEvaluationPreservesAdmmCandidate() {
        OSQP.Settings s=new OSQP.Settings(); s.scaling=0; s.verbose=false; s.adaptive_rho=false;
        CSCMatrix one=new CSCMatrix(1,1,1,new int[]{0,1},new int[]{0},new double[]{1});
        OSQP solver=new OSQP(new OSQP.Data(1,1,one,one,new double[]{-2},new double[]{0},new double[]{1},0),s);
        assertEquals(OSQP.Status.SOLVED,solver.solve());
        double admmDual=solver.getInfo().dua_res;
        // A deliberately worse polish candidate has exact stationarity residual |0-2+0| = 2.
        solver.work.pol.x[0]=0; solver.work.pol.y[0]=0; solver.work.pol.z[0]=0;
        OSQP.update_info(solver.work,0,true,true);
        assertEquals(admmDual,solver.getInfo().dua_res,0);
        assertEquals(2,solver.work.pol.dua_res,0);
    }

    @Test public void negativeCurvatureIsRejectedDuringSetup() {
        CSCMatrix p=new CSCMatrix(1,1,1,new int[]{0,1},new int[]{0},new double[]{-1});
        CSCMatrix a=new CSCMatrix(0,1,0,new int[]{0,0},new int[]{},new double[]{});
        OSQP.Settings settings=new OSQP.Settings(); settings.scaling=0; settings.verbose=false;
        assertThrows(IllegalStateException.class, () -> new OSQP(new OSQP.Data(1,0,p,a,new double[]{0},
                new double[]{},new double[]{},0),settings));
    }
    @Test public void invalidMatrixUpdateBlocksSolveUntilSuccessfulRefactor() {
        CSCMatrix p=new CSCMatrix(1,1,1,new int[]{0,1},new int[]{0},new double[]{1});
        CSCMatrix a=new CSCMatrix(1,1,1,new int[]{0,1},new int[]{0},new double[]{1});
        LinSys system=new LinSys(p,a,1e-6,new double[]{1},false);
        p.Ax[0]=-3;
        assertThrows(IllegalStateException.class, () -> system.update_solver_matrices(p,a));
        assertNull(system.qdldl._ldl);
        assertThrows(IllegalStateException.class, () -> system.solve(new double[]{1,0}));
        p.Ax[0]=1;system.update_solver_matrices(p,a);
        assertEquals(1,system.qdldl._ldl.positiveValuesInD);
        double[] rhs={1,0};system.solve(rhs);
        assertEquals(1/(2+1e-6),rhs[0],1e-14);
    }

}
