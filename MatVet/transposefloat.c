/**
 * Transpose test
 */
#include <stdio.h>
#include <stdlib.h>

void print_matrix_array(float *array, int rows, int columns);
float *transpose_matrix_array(float *array, int rows, int columns);
void transpose_matrix_array_opt(float *array, int rows, int columns);
void swapfloat(float *a, float *b);

int main() {
    float *a, *transpose;
    int r, c, i, j;

    printf("Enter rows of matrix: ");
    fflush(stdout);
    scanf("%d", &r);
    printf("Enter columns of matrix: ");
    fflush(stdout);
    scanf("%d", &c);

    // allocate space for matrix a
    a = (float *) malloc(r * c * sizeof(float));

    // Storing elements of the matrix
    printf("\nGenerating elements of matrix...\n");
    for (i=0; i<r*c; ++i) {
        a[i] = (float) i;
    }

    // Displaying the matrix a
    printf("\nEntered Matrix: \n");
    print_matrix_array(a, r, c);

    // Finding the transpose of matrix a
    transpose = transpose_matrix_array(a, r, c);

    // Displaying the transpose of matrix a
    printf("\nTranspose of Matrix:\n");
    print_matrix_array(transpose, c, r);

    return 0;
}

void print_matrix_array(float *array, int rows, int columns) {
    int i, j;

    for (i = 0; i < rows * columns; ++i) {
        printf("%5.2f\t", array[i]);
        if (i % columns == columns - 1) {
            printf("\n\n");
        }
    }
}

float *transpose_matrix_array(float *array, int rows, int columns) {
    int i, j;
    int source, source_position, destination_position;
    float *transpose;

    // allocate space for transpose matrix
    transpose = (float *) calloc(columns * rows, sizeof(float));

    for (i = 0; i < rows; ++i) {
        // printf("ROW %d\n", i);
        for (j = 0; j < columns; ++j) {
            // printf("Element %d of row %d\n", j, i);
            source_position = columns * i + j;
            destination_position = rows * j + i;
            source = *(array + source_position);
            // printf("Transposing element %f (%d%d) from array cell %d to array cell %d\n", source, i, j, source_position, destination_position);
            *(transpose + destination_position) = source;
            // printf("Cell %d (%d%d) value of transposed matrix is %d\n", destination_position, i, j, transpose[destination_position]);
        }
        // printf("\n");
    }

    return transpose;
}

void transpose_matrix_array_opt(float *array, int rows, int columns) {
    int i, j;
    int source, source_position, destination_position;

    for (i = 0; i < rows; ++i) {
        // printf("ROW %d\n", i);
        for (j = 0; j < columns; ++j) {
            // printf("Element %d of row %d\n", j, i);
            source_position = columns * i + j;
            destination_position = rows * j + i;
            // printf("Transposing element %f (%d%d) from array cell %d to array cell %d\n", array[source_position], i, j, source_position, destination_position);
            swapfloat(array+source_position, array+destination_position);
            // printf("Cell %d (%d%d) value of transposed matrix is %d\n", destination_position, i, j, transpose[destination_position]);
        }
        // printf("\n");
    }

    return;
}

void swapfloat(float *a, float *b) {
    float temp;

    // printf("A is %f, B is %f\n", *a, *b);
    temp = *a;
    *a = *b;
    *b = temp;
    // printf("Now A is %f, B is %f\n\n", *a, *b);

    return;
}
