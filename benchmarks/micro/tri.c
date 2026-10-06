/* Same sparse triangular solve as JOSQP's QDLDL.Lsolve/Ltsolve, in C, as the reference. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int n, nz;
static int *Lp, *Li;
static double *Lx;

static double now(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

int main(int argc, char **argv) {
    n = 20000;
    int percol = 10;
    nz = n * percol;
    Lp = malloc((n + 1) * sizeof(int));
    Li = malloc(nz * sizeof(int));
    Lx = malloc(nz * sizeof(double));
    unsigned long long s = 12345;
    for (int i = 0; i < n; i++) {
        Lp[i] = i * percol;
        for (int k = 0; k < percol; k++) {
            s = s * 6364136223846793005ULL + 1442695040888963407ULL;
            int span = n - i - 1;
            Li[i * percol + k] = span > 0 ? i + 1 + (int)((s >> 33) % span) : i;
            Lx[i * percol + k] = 1e-4 + (double)((s >> 20) % 1000) * 1e-7;
        }
    }
    Lp[n] = nz;
    double *x = malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) x[i] = 1.0 + 1e-9 * i;

    double best = 1e30, checksum = 0;
    for (int rep = 0; rep < 300; rep++) {
        for (int i = 0; i < n; i++) x[i] = 1.0 + 1e-9 * i;
        double t = now();
        for (int i = 0; i < n; i++) {
            double xi = x[i];
            for (int j = Lp[i]; j < Lp[i + 1]; j++) x[Li[j]] -= Lx[j] * xi;
        }
        for (int i = n - 1; i >= 0; i--) {
            double xi = x[i] * 0.999;
            for (int j = Lp[i]; j < Lp[i + 1]; j++) xi -= Lx[j] * x[Li[j]];
            x[i] = xi;
        }
        double dt = now() - t;
        if (rep >= 100 && dt < best) best = dt;
        checksum += x[n / 2];
    }
    printf("C     best %.3f ms  per-nnz %.2f ns  checksum %.6g\n",
           best * 1e3, best * 1e9 / (2.0 * nz), checksum);
    return 0;
}
