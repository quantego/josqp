public class TriBench {
    static int n, nz;
    static int[] Lp, Li;
    static double[] Lx;

    /** Current JOSQP shape. */
    static void solveA(int[] Lp, int[] Li, double[] Lx, double[] x) {
        for (int i = 0; i < n; i++) {
            final double xi = x[i];
            for (int j = Lp[i]; j < Lp[i + 1]; j++) x[Li[j]] -= Lx[j] * xi;
        }
        for (int i = n - 1; i >= 0; i--) {
            double xi = x[i] * 0.999;
            for (int j = Lp[i]; j < Lp[i + 1]; j++) xi -= Lx[j] * x[Li[j]];
            x[i] = xi;
        }
    }

    /** Clamped bounds: gives C2 the range facts it needs to drop the Li/Lx checks. */
    static void solveB(int[] Lp, int[] Li, double[] Lx, double[] x) {
        final int len = Li.length;
        for (int i = 0; i < n; i++) {
            final double xi = x[i];
            int start = Lp[i], end = Lp[i + 1];
            if (start < 0) start = 0;
            if (end > len) end = len;
            for (int j = start; j < end; j++) x[Li[j]] -= Lx[j] * xi;
        }
        for (int i = n - 1; i >= 0; i--) {
            double xi = x[i] * 0.999;
            int start = Lp[i], end = Lp[i + 1];
            if (start < 0) start = 0;
            if (end > len) end = len;
            for (int j = start; j < end; j++) xi -= Lx[j] * x[Li[j]];
            x[i] = xi;
        }
    }

    public static void main(String[] args) {
        n = 20000;
        final int percol = 10;
        nz = n * percol;
        Lp = new int[n + 1];
        Li = new int[nz];
        Lx = new double[nz];
        long s = 12345;
        for (int i = 0; i < n; i++) {
            Lp[i] = i * percol;
            for (int k = 0; k < percol; k++) {
                s = s * 6364136223846793005L + 1442695040888963407L;
                int span = n - i - 1;
                Li[i * percol + k] = span > 0 ? i + 1 + (int) ((s >>> 33) % span) : i;
                Lx[i * percol + k] = 1e-4 + ((s >>> 20) % 1000) * 1e-7;
            }
        }
        Lp[n] = nz;
        double[] x = new double[n];
        for (int mode = 0; mode < 2; mode++) {
            double bestA = 1e30, bestB = 1e30, checksum = 0;
            for (int rep = 0; rep < 400; rep++) {
                for (int i = 0; i < n; i++) x[i] = 1.0 + 1e-9 * i;
                long t = System.nanoTime();
                solveA(Lp, Li, Lx, x);
                double dtA = (System.nanoTime() - t) / 1e6;
                if (rep >= 150 && dtA < bestA) bestA = dtA;
                checksum += x[n / 2];

                for (int i = 0; i < n; i++) x[i] = 1.0 + 1e-9 * i;
                t = System.nanoTime();
                solveB(Lp, Li, Lx, x);
                double dtB = (System.nanoTime() - t) / 1e6;
                if (rep >= 150 && dtB < bestB) bestB = dtB;
                checksum += x[n / 2];
            }
            System.out.printf("Java  A %.3f ms (%.2f ns/nnz)   B %.3f ms (%.2f ns/nnz)   B/A %.3f  checksum %.6g%n",
                    bestA, bestA * 1e6 / (2.0 * nz), bestB, bestB * 1e6 / (2.0 * nz), bestB / bestA, checksum);
        }
    }
}
