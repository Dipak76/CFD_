#include <stdio.h>
#include <math.h>

int main() {
    int M = 31, N = 21;

    /* Defining the grid size */
    double L = 6.0, H = 4.0;
    double dx = L / (M - 1);
    double dy = H / (N - 1);

    /* Defining matrices to store grid point values */
    double psi_k[M][N], psi_k1[M][N];

    /* Defining some points for the function */
    double ai, bi, ci, aj, bj, cj, ap, ae, aw, an, as;
    ae = 1.0 / (pow(dx, 2.0));
    aw = 1.0 / (pow(dx, 2.0));
    an = 1.0 / (pow(dy, 2.0));
    as = 1.0 / (pow(dy, 2.0));
    ap = ae + aw + an + as;
    ai = ap;
    bi = -ae;
    ci = -aw;
    aj = ap;
    bj = -an;
    cj = -as;

    /* Defining boundary conditions */
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            psi_k[i][j] = 0.0; /* Assigning zero values to all the grid points including interior points, Top boundary, Left boundary, and Right boundary up to i=5 */
        }
    }
    for (int i = 6; i < M; i++) {
        psi_k[i][0] = 100.0; /* Bottom boundary */
    }

    /* Alternate Direction Implicit (ADI) method */
    double d1[M - 2], p1[M - 2], q1[M - 2];
    double d2[M - 2], p2[M - 2], q2[M - 2];

    int it = 0;
    double temp, error = 1.0;
    FILE *myfile1 = fopen("Error5.dat", "w"); /* Writing error in the file */
    printf("Iteration\tError\n");
    do {
        error = 0.0;
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                psi_k1[i][j] = psi_k[i][j];
            }
        }
        for (int j = 1; j < N - 1; j++) {
            d1[0] = an * (psi_k[0][j - 1] + psi_k[0][j + 1]);
            p1[0] = -bi / ai;
            q1[0] = d1[0] / ai;
            for (int i = 1; i < M - 1; i++) { /* Finding coefficients */
                d1[i] = an * (psi_k[i][j - 1] + psi_k[i][j + 1]);
                p1[i] = -bi / (ai + ci * p1[i - 1]);
                q1[i] = (d1[i] - ci * q1[i - 1]) / (ai + ci * p1[i - 1]);
            }
            for (int i = M - 2; i > 0; i--) {
                psi_k[i][j] = p1[i] * psi_k[i + 1][j] + q1[i];
            }
        }
        for (int i = 1; i < M - 1; i++) {
            d2[0] = ae * (psi_k[i - 1][0] + psi_k[i + 1][0]);
            p2[0] = -bj / aj;
            q2[0] = d2[0] / aj;
            for (int j = 1; j < N - 1; j++) { /* Finding coefficients */
                d2[j] = ae * (psi_k[i - 1][j] + psi_k[i + 1][j]);
                p2[j] = -bj / (aj + cj * p2[j - 1]);
                q2[j] = (d2[j] - cj * q2[j - 1]) / (aj + cj * p2[j - 1]);
            }
            for (int j = N - 2; j > 0; j--) {
                psi_k[i][j] = p2[j] * psi_k[i][j + 1] + q2[j];
            }
        }
        for (int j = 0; j < N; j++) {
            psi_k[M - 1][j] = psi_k[M - 2][j];
        }
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                error += pow((psi_k[i][j] - psi_k1[i][j]), 2);
            }
        }
        error = sqrt(error / (double)((M - 2) * (N - 2)));
        fprintf(myfile1, "%d\t%.6f\n", it, log10(error));
        printf("%d\t%.6f\n", it, error);
        it++;
    } while (error > 1e-6);
    fclose(myfile1);

    FILE *myfile2 = fopen("Stream5.plt", "w");
    fprintf(myfile2, "VARIABLES = \"X\", \"Y\", \"PHI\"\n");
    fprintf(myfile2, "ZONE T = \"BLOCK1\", I = 31, J = 21, F = POINT\n\n");
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            fprintf(myfile2, "%.6f\t%.6f\t%.6f\n", i * dx, j * dy, psi_k[i][j]);
        }
    }
    fclose(myfile2);

    return 0;
}
