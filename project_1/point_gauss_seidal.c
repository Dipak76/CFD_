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
    double A = fabs(1 / (2 * (1 + pow(beta, 2))));
    double B = pow(beta, 2);

    /* Defining boundary conditions */
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            psi[i][j] = 0.0; /* Assigning zero values to all the grid points including interior points, Top boundary, Left boundary and Right boundary up to i=5 */
        }
    }
    for (int i = 6; i < M; i++) {
        psi[i][0] = 100.0; /* Bottom boundary */
    }

    /* Point Gauss-Seidel iterative method */
    int it = 0;
    double temp, error = 1.0;
    FILE *myfile1 = fopen("Error2.dat", "w"); /* Writing error in the file */
    printf("Iteration\tError\n");
    do {
        error = 0.0;
        for (int i = 1; i < (M - 1); i++) {   /* Iterating for interior points only */
            for (int j = 1; j < (N - 1); j++) {
                temp = psi[i][j];
                psi[i][j] = A * (B * psi[i][j - 1] + psi[i - 1][j] + psi[i + 1][j] + B * psi[i][j + 1]);
                error += pow((psi[i][j] - temp), 2.0);
            }
        }
        for (int j = 0; j < N; j++) {   /* For homogeneous Neumann boundary on right wall */
            psi[M - 1][j] = psi[M - 2][j];
        }
        error = sqrt(error / (double)((M - 2) * (N - 2)));
        fprintf(myfile1, "%d\t%.6f\n", it, log(error));
        printf("%d\t%.6f\n", it, error);
        it++;
    } while (error > 1e-6);
    fclose(myfile1);

    FILE *myfile2 = fopen("Stream2.plt", "w");
    fprintf(myfile2, "VARIABLES = \"X\", \"Y\", \"PHI\"\n");
    fprintf(myfile2, "ZONE T = \"BLOCK1\", I = 31, J = 21, F = POINT\n\n");
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            fprintf(myfile2, "%.6f\t%.6f\t%.6f\n", i * dx, j * dy, psi[i][j]);
        }
    }
    fclose(myfile2);

    return 0;
}
