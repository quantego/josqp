package com.quantego.josqp;

public class OSQPTester {

    public static final double TESTS_TOL = 1e-4; //TODO: check why 3e-4 fixes .testUpdateRho()

    public static double[] csc_to_dns(CSCMatrix M) {
        int i, j = 0; // Predefine row index and column index
        int idx;

        // Initialize matrix of zeros
        double[] A = new double[M.m * M.n];

        // Allocate elements
        for (idx = 0; idx < M.Ap[M.n]; idx++) {
            // Get row index i (starting from 1)
            i = M.Ai[idx];

            // Get column index j (increase if necessary) (starting from 1)
            while (M.Ap[j + 1] <= idx) {
                j++;
            }

            // Assign values to A
            A[j * (M.m) + i] = M.Ax[idx];
        }
        return A;
    }

    public static boolean is_eq_csc(CSCMatrix A, CSCMatrix B, double tol) {
        int j, i;

        // If number of columns does not coincide, they are not equal.
        if (A.n != B.n)
            return false;

        for (j = 0; j < A.n; j++) { // Cycle over columns j
            // if column pointer does not coincide, they are not equal
            if (A.Ap[j] != B.Ap[j])
                return false;

            for (i = A.Ap[j]; i < A.Ap[j + 1]; i++) { // Cycle rows i in column j
                if ((A.Ai[i] != B.Ai[i]) || // Different row indices
                        (Math.abs(A.Ax[i] - B.Ax[i]) > tol)) {
                    return false;
                }
            }
        }
        return true;
    }

}
