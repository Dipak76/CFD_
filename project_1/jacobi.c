#include <stdio.h>
#include <math.h>

int main() {
    int M = 31, N = 21;

    /* Defining the grid size */
    double L = 6.0, H = 4.0;
    double dx = L / (M - 1);
    double dy = H / (N - 1);

    /* Defining matrices to store grid point values */
    double psiOld[M][N], psiNew[M][N];

    /* Defining some points for the function */
    double beta = dx / dy;
    double B = pow(beta, 2);
    double A = fabs(1 / (2 * (1 + B)));

    /* Defining boundary conditions */
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            psiNew[i][j] = 0.0; /* Assigning zero values to all the grid points including interior points, Top boundary, Left boundary and Right boundary up to i=5 */
        }
    }
    for (int i = 6; i < M; i++) {
        psiNew[i][0] = 100.0; /* Bottom boundary */
    }

    /* Jacobi iterative method */
    int it = 0;
    double error = 1.0;
    FILE *myfile1 = fopen("Error1.dat", "w"); /* Writing error in the file */
    printf("Iteration\tError\n");
    do {
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                psiOld[i][j] = psiNew[i][j];  /* Replacing old function values with new function values */
            }
        }
        for (int i = 1; i < (M - 1); i++) {   /* Iterating for interior points only */
            for (int j = 1; j < (N - 1); j++) {
                psiNew[i][j] = A * (B * psiOld[i][j - 1] + psiOld[i - 1][j] + psiOld[i + 1][j] + B * psiOld[i][j + 1]);
            }
        }
        for (int j = 0; j < N; j++) {   /* For homogeneous Neumann boundary on right wall */
            psiNew[M - 1][j] = psiNew[M - 2][j];
        }
        error = 0.0;
        for (int i = 1; i < M - 2; i++) {
            for (int j = 1; j < N - 2; j++) {
                error = error + pow((psiNew[i][j] - psiOld[i][j]), 2.0);
            }
        }
        error = sqrt(error / (double)((M - 2) * (N - 2)));
        fprintf(myfile1, "%d\t%.6f\n", it, log10(error));
        printf("%d\t%.6f\n", it, error);
        it++;
    } while (error > 1e-6);
    fclose(myfile1);

    FILE *myfile2 = fopen("Stream1.plt", "w");
    fprintf(myfile2, "VARIABLES = \"X\", \"Y\", \"PHI\"\n");
    fprintf(myfile2, "ZONE T = \"BLOCK1\", I = 31, J = 21, F = POINT\n\n");
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            fprintf(myfile2, "%.6f\t%.6f\t%.6f\n", i * dx, j * dy, psiNew[i][j]);
        }
    }
    fclose(myfile2);

    return 0;
}
