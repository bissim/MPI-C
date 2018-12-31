#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mpi.h"
#include "matrixutils.h"

int main (int argc, char **argv) {
    int nproc;              // Numero di processi totale
    int me;                 // Il mio id
    int n;                  // Dimensione della matrice
    int local_n;            // Dimensione dei dati locali
    int i, j;               // Iteratori vari

    // Variabili di lavoro
    double *A, *At, *v , *local_v, *localA, *localAt, *local_w, *w;
    int num;

    // initialise pointers
    A = At = v = local_v = localA = local_w = w = NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);

    printf("[P%d] MPI process started!\nRunning program %s\n", me, argv[0]);

    // master generates data
    if (me == 0) {
        if (argc < 2) {
            printf("[P%d] Usage: <programname> <rows> [<columns>]\n", me);
            printf("[P%d] Insert n: ", me); 
            fflush(stdout);
            scanf("%d", &n);
            printf("\n");
        }
        else {
            n = atoi(argv[1]);
            printf("[P%d] n: %d\n", me, n);
        }

        // allocate memory
        A = malloc(sizeof *A * n * n);
        v = malloc(sizeof *v * n);
        w = malloc(sizeof *w * n);

        // create transposed matrix and vector
        for (i = 0; i < n; i++) {
            v[i] = i;
            for (j = 0; j < n; j++) {
                if (j == 0) {
                    A[j*n+i] = 1.0/(i + 1) - 1; // inverted index for transposed matrix
                }
                else {
                    A[j*n+i] = 1.0/(i + 1) - pow(1.0/2.0, j); 
                }
            }
        }

        printf("[P%d] Transposed matrix A:\n", me);
        print_matrix_array(A, n, n);
        printf("\n\n");
        printf("[P%d] Vector v:\n", me);
        print_matrix_array(v, n, 1);
        printf("\n\n");
        fflush(stdout);
    }

    // data part for every process
    local_n = n / nproc;

    // master sends n and local_n
    // to other processes
    if (me == 0){
        printf("[P%d] Sending n = %d, local_n = %d\n", me, n, local_n);
        fflush(stdout);
    }
    MPI_Bcast(
        &n, 1, MPI_INT,
        0, MPI_COMM_WORLD
    );
    MPI_Bcast(
        &local_n, 1, MPI_INT,
        0, MPI_COMM_WORLD
    );
    if (me != 0) {
        printf("[P%d] n: %d, local_n: %d\n", me, n, local_n);
    }

    // everyone allocates local_v
    local_v = malloc(sizeof *local_v * local_n);

    if (me == 0) {
        printf("[P%d] Scattering vector into parts of %d elements...\n", me, local_n);
        fflush(stdout);
    }
    MPI_Scatter(
        v, local_n, MPI_DOUBLE,
        local_v, local_n, MPI_DOUBLE,
        0, MPI_COMM_WORLD
    );
    if (me == 0 ) {
        printf("[P%d] Vector scattered!\n", me);
        fflush(stdout);
    }
    else {
        printf("[P%d] Vector received!\n", me);
        fflush(stdout);
    }
    printf("[P%d] local_v:\n", me);
    print_matrix_array(local_v, local_n, 1);
    printf("\n");
    fflush(stdout);

    // generating transposed matrix before sending
    // if (me == 0) {
    //     At = malloc(sizeof *At * n * n);
    //     transpose_matrix_array(A, At, n, n);
    //     printf("[P%d] Matrix transposed!\n", me);
    //     print_matrix_array(At, n, n);
    //     printf("\n");
    //     fflush(stdout);
    // }

    // everyone allocates localA and local_w
    num = n * local_n;
    localAt = malloc(sizeof *localAt * num);
    MPI_Scatter(
        // At, num, MPI_DOUBLE,
        A, num, MPI_DOUBLE,
        localAt, num, MPI_DOUBLE,
        0, MPI_COMM_WORLD
    );

    if (me == 0) {
        printf("[P%d] Matrix scattered!\n", me);
    }
    else {
        printf("[P%d] Matrix received!\n", me);
    }
    printf("[P%d] localAt:\n", me);
    print_matrix_array(localAt, local_n, n);
    printf("\n");
    fflush(stdout);

    // transpose received matrix
    // localA = malloc(sizeof *localA * num);
    // transpose_matrix_array(localAt, localA, local_n, n);
    // printf("[P%d] localA:\n", me);
    // print_matrix_array(localA, n, local_n);
    // printf("\n");
    // fflush(stdout);

    // perform matvect product
    local_w = malloc(sizeof *local_w * n);
    // prod_mat_vett(local_w, localA, n, local_n, local_v);
    // prod_mat_vett(local_w, localAt, n, local_n, local_v);
    for (i = 0; i < n; i++) {
        local_w[i] = 0;
        for (j = 0; j < n; j++) {
            local_w[i] += localAt[j*n+i] * local_v[j];
        }
    }

    // print local results
    printf("[P%d] Partial result:\n", me);
    print_matrix_array(local_w, n, 1);
    printf("\n");
    fflush(stdout);

    // master collects local results
    MPI_Reduce(
        local_w, w, n, MPI_DOUBLE, MPI_SUM,
        0, MPI_COMM_WORLD
    );
    if (me == 0) {
        printf("[P%d] Results collected!\nw:\n", me);
        print_matrix_array(w, n, 1);
    }
    else {
        printf("[P%d] Results sent!\n", me);
    }
    fflush(stdout);

    MPI_Finalize();
    return EXIT_SUCCESS;
}
