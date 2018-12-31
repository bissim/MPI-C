#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "matrixutils.h"

/**
 * Funzione che esegue la stampa di una matrice di interi
 */
void print_matrix_array(double *array, int rows, int columns) {
    int i;

    if (fmax((float) rows, (float) columns) > 50 || rows * columns > 100) {
        printf("\tToo large to represent (%d elements)!\n", rows * columns);
        return;
    }

    for (i = 0; i < rows * columns; ++i) {
        printf("%6.3f\t", array[i]);
        if (i % columns == columns - 1) {
            printf("\n");
        }
    }
}

/**
 * Funzione che esegue il prodotto matrice vettore
 */
void prod_mat_vett(double result[], double *a, int rows, int cols, double v[]) {
    int i, j;
    
    for (i = 0; i < rows; i++) {
        result[i] = 0;
        for (j = 0; j < cols; j++) {
            result[i] += a[i*cols+j] * v[j];
        }
    }
}

/**
 * Funzione che esegue la trasposizione di una matrice di interi
 */
void transpose_matrix_array(double *array, double *transpose, int rows, int columns) { // TODO restore
    int i, j;
    int source_position, destination_position;
    float source;

    for (i = 0; i < rows; ++i) {
        for (j = 0; j < columns; ++j) {
            source_position = columns * i + j;
            destination_position = rows * j + i;
            source = *(array + source_position);
            *(transpose + destination_position) = source;
        }
    }
}
