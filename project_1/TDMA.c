#include <stdio.h>
#include <math.h>

int main() {
    int M = 31, N = 21;

    /* Defining the grid size */
    double L = 6.0, H = 4.0;
    double dx = L / (M - 1);
    double dy = H / (N - 1);

    /* Defining matrices to store grid point values */
    double psi[M][N];

    /* Defining some points for the function */
    double beta = dx / dy;
    double B = pow(beta, 2);
    double A = fabs(1 / (2 * (1 + B)));

    /* Defining boundary conditions */
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            psi[i][j] = 0.0; /* Assigning zero values to all the grid points including interior points, Top boundary, Left boundary, and Right boundary up to i=5 */
        }
    }
    for (int i = 6; i < M; i++) {
        psi[i][0] = 100.0; /* Bottom boundary */
    }

    /* Tri Diagonal Matrix Algorithm (TDMA) method */
    double d[M - 2], p[M - 2], q[M - 2], ai = (-2 * (1 + B)), bi = 1, ci = 1;

    int it = 0;
    double temp, error = 1.0;
    FILE *myfile1 = fopen("Error4.dat", "w"); /* Writing error in the file */
    printf("Iteration\tError\n");
    do {
        error = 0.0;
        for (int j = 1; j < N - 1; j++) {
            d[0] = -B * (psi[0][j - 1] + psi[0][j + 1]);
            p[0] = -bi / ai;
            q[0] = d[0] / ai;
            for (int i = 1; i < M - 1; i++) { /* Finding coefficients */
                d[i] = -B * (psi[i][j - 1] + psi[i][j + 1]);
                p[i] = -bi / (ai + ci * p[i - 1]);
                q[i] = (d[i] - ci * q[i - 1]) / (ai + ci * p[i - 1]);
            }
            for (int i = M - 2; i > 0; i--) {
                temp = psi[i][j];
                psi[i][j] = p[i] * psi[i + 1][j] + q[i];
                error += pow((psi[i][j] - temp), 2.0);
            }
        }
        for (int j = 0; j < N; j++) { /* Neumann boundary condition */
            psi[M - 1][j] = psi[M - 2][j];
        }
        error = sqrt(error / (double)((M - 2) * (N - 2)));
        fprintf(myfile1, "%d\t%.6f\n", it, log10(error));
        printf("%d\t%.6f\n", it, error);
        it++;
    } while (error > 1e-6);
    fclose(myfile1);

    double x[M], y[N];
    x[0] = y[0] = 0.0;
    for (int i = 0; i < M - 1; i++) {
        x[i + 1] = x[i] + dx;
    }
    for (int i = 0; i < N - 1; i++) {
        y[i + 1] = y[i] + dy;
    }
    FILE *myfile2 = fopen("Stream4.plt", "w");
    fprintf(myfile2, "VARIABLES = \"X\", \"Y\", \"PHI\"\n");
    fprintf(myfile2, "ZONE T = \"BLOCK1\", I = 31, J = 21, F = POINT\n\n");
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            fprintf(myfile2, "%.6f\t%.6f\t%.6f\n", x[i], y[j], psi[i][j]);
        }
    }
    fclose(myfile2);

    return 0;
}
