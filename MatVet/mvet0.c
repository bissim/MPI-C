#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "matrixutils.h"

// const int NANOSECONDS_PER_SECOND = 1E9;
#define NANOSECONDS_PER_SECOND 1E9;

int main(int argc, char **argv) {
    int i, j;
    int n;
    double *A, *v, *w;
	double elapsedtime;
	struct timespec start, stop;
    int flag = 0;
    const int MS_IN_S = 1000;

    // take matrix dimension
    // printf("Insert n: ");
    // scanf("%d", &n);
    // printf("\n");
    if (argc < 3) {
		printf("Numero di parametri insufficiente!\n");
        fflush(stdout);
		exit(EXIT_FAILURE);
	}
	else {
		n = atoi(argv[1]);
        flag = atoi(argv[2]);
	}
    printf("Dimensione della matrice: %d (%d elementi)\n", n, n*n);
    fflush(stdout);

    // Alloco spazio di memoria
    A = malloc(n * n * sizeof *A);
    v = malloc(n * sizeof *v);

    // printf("\nGenerating matrix and vector...\n");
    for (i = 0; i < n; i++) {
        v[i] = i;
        for (j = 0; j < n; j++) {
            if (j == 0) {
                A[i*n+j] = 1.0/(i + 1) - 1;
            }
            else {
                A[i*n+j] = 1.0/(i+1) - pow(1.0/2.0, j); 
            }
        }
    }

    if (flag) {
        printf("Matrix A:\n");
        print_matrix_array(A, n, n);
        printf("\n");
        printf("Vector v:\n");
        print_matrix_array(v, n, 1);
        printf("\n");
        fflush(stdout);
    }

    w = malloc(n * sizeof *w);

    // calculate matvet product
	clock_gettime(CLOCK_REALTIME, &start);
    prod_mat_vett(w, A, n, n, v);
	clock_gettime(CLOCK_REALTIME, &stop);

    if (flag) {
        printf("Result vector w:\n");
        print_matrix_array(w, n, 1);
    }

	elapsedtime = (stop.tv_sec - start.tv_sec) + (stop.tv_nsec - start.tv_nsec) / NANOSECONDS_PER_SECOND;
	printf("Elapsed time: %f ms.\n", elapsedtime * MS_IN_S);
    fflush(stdout);

    return EXIT_SUCCESS;
}
