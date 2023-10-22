// laplacian_test.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <stdio.h>
#include <stdlib.h>
#define EPS 10e-6
#define RESOLUTION 250
#define LENGTH 1

void abs(double* a) {
    if (*a < 0) {
        *a = -*a;
    }
}

int main() {

    double dx = 1.0 * LENGTH / RESOLUTION;

    double b[RESOLUTION][RESOLUTION];
    double x[RESOLUTION + 2][RESOLUTION + 2];
    double delta[RESOLUTION][RESOLUTION];

    bool finish;

    double omega;


    // Init
    int iter = 0;

    for (int i = 0; i < RESOLUTION; i++) {
        for (int j = 0; j < RESOLUTION; j++) {
            b[i][j] = 5 * dx;
        }
    }


    // Solver
    for (int i = 0; i < RESOLUTION + 2; i++) {
        for (int j = 0; j < RESOLUTION + 2; j++) {
            x[i][j] = 0;
        }
    }

    omega = 1.85;
    finish = 0;

    while (!finish) {
        finish = 1;

        for (int i = 0; i < RESOLUTION; i++) {
            for (int j = 0; j < RESOLUTION; j++) {
                delta[i][j] = x[i + 1][j + 1];
                x[i + 1][j + 1] = (1 - omega) * x[i + 1][j + 1] + omega * (b[i][j] + x[i][j + 1] + x[i + 2][j + 1] + x[i + 1][j] + x[i + 1][j + 2]) / 4;
                delta[i][j] -= x[i + 1][j + 1];
                abs(&delta[i][j]);
                if (delta[i][j] > EPS) {
                    finish = 0;
                }
            }
        }
        iter++;
    }

    for (int i = 0; i < RESOLUTION + 2; i++) {
        for (int j = 0; j < RESOLUTION + 2; j++) {
            printf("%2.1f ", x[i][j]);
        }
        printf("\n");
    }

    printf("\niter = %d\n", iter);
    return 0;
}
