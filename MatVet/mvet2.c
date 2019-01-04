#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mpi.h"
#include "mpiutils.h"
#include "matrixutils.h"

extern int master;

int main(int argc, char **argv) {
    char *progname;
    int nproc;              // Numero di processi totale
    int me;                 // Il mio id
    int n;                  // Dimensione della matrice
    int rows, cols;
    int local_r, local_c;   // Dimensione dei dati locali
    int i, j;               // Iteratori vari 
    int err = MPI_SUCCESS;  // variabile di errore MPI

    // Variabili di lavoro
    double *A = NULL, *At = NULL, *v = NULL, *local_v = NULL, *localA = NULL, *localAt = NULL, *local_w = NULL, *w = NULL;
    int num;

    // Enable MPI
    err = MPI_Init(&argc, &argv);
    checkMPIerror(&me, &err);

    // find out number of total processes
    err = MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    checkMPIerror(&me, &err);

    // Identify every process
    err = MPI_Comm_rank(MPI_COMM_WORLD, &me);
    checkMPIerror(&me, &err);

    MPI_Printf(me, "MPI process started!\n");
    MPI_Printf(me, "Running program ");
    progname = argv[0];
    printf("%s\n\n", progname);

    // Master generates data
    master = 0;
    if (me == master) {
        if (argc < 3) { // inserted rows only, square matrix
            rows = atoi(argv[1]);
            cols = rows;
            n = rows * cols;
        }
        else if (argc < 2) { // no number inserted
            MPI_Printf(me, "Usage: <programname> <rows> [<columns>]\n");
            MPI_Printf(me, "Insert n: "); 
            fflush(stdout);
            scanf("%d", &rows);
            printf("\n");
            rows = cols;
            n = rows * cols;
        }
        else { // non square matrix
            rows = atoi(argv[1]);
            cols = atoi(argv[2]);
            n = rows * cols;
        }

        MPI_Printf(me, "rows: ");
        printf("%d\n", rows);
        MPI_Printf(me, "columns: ");
        printf("%d\n", cols);
        MPI_Printf(me, "n: ");
        printf("%d\n", n);
        
        // Allocate memory
        A = malloc(sizeof *A * n);
        v = malloc(sizeof *v * cols);
        w = malloc(sizeof *w * rows); 

        // Create matrix A and vector v
        for (i = 0; i < rows; i++) {
            for (j = 0; j < cols; j++) {
                v[j] = j;
                if (j == 0) {
                    A[i*rows+j] = 1.0/(i + 1) - 1;
                }
                else {
                    A[i*rows+j] = 1.0/(i + 1) - pow(1.0/2.0, j); 
                }
            }
        }

        MPI_Printf(me, "Matrix A:\n\n");
        print_matrix_array(A, rows, cols);
        printf("\n");
        MPI_Printf(me, "Vector v:\n\n");
        print_matrix_array(v, cols, 1);
        printf("\n");
        fflush(stdout);
    }

    // Data part for every process to deal with
    local_r = rows / nproc;
    local_c = cols / nproc;

    // master sends n and local_n
    // to other processes
    if (me == master) {
        printf("[P%d] Sending rows = %d, columns = %d, local_r = %d, local_c = %d\n", me, rows, cols, local_r, local_c);
        fflush(stdout);
    }
    err = MPI_Bcast(
        &rows, 1, MPI_INT,
        master, MPI_COMM_WORLD
    );
    checkMPIerror(&me, &err);
    err = MPI_Bcast(
        &cols, 1, MPI_INT,
        master, MPI_COMM_WORLD
    );
    checkMPIerror(&me, &err);
    err = MPI_Bcast(
        &local_r, 1, MPI_INT,
        master, MPI_COMM_WORLD
    );
    checkMPIerror(&me, &err);
    err = MPI_Bcast(
        &local_c, 1, MPI_INT,
        master, MPI_COMM_WORLD
    );
    checkMPIerror(&me, &err);

    MPI_Printf(me, "rows: ");
    printf("%d\n", rows);
    MPI_Printf(me, "columns: ");
    printf("%d\n", cols);

    // everyone allocates local_v
    local_v = malloc(sizeof *local_v * local_c);

    if (me == master) {
        printf("Scattering vector into parts of %d elements...\n", local_c);
        fflush(stdout);
    }
    err = MPI_Scatter(
        v, local_c, MPI_DOUBLE,
        local_v, local_c, MPI_DOUBLE,
        master, MPI_COMM_WORLD
    );
    checkMPIerror(&me, &err);

    // free unneeded memory for vector
    // MPI_Free(me, v);

    if (me == master) {
        MPI_Printf(me, "Vector scattered!\n");
        fflush(stdout);
    }
    else {
        MPI_Printf(me, "Vector received!\n");
        fflush(stdout);
    }
    MPI_Printf(me, "local_v:\n");
    print_matrix_array(local_v, local_c, 1);
    printf("\n");
    fflush(stdout);

    // il master crea la trasposta della matrice
    // per poterla inviare per colonne
    if (me == master) {
        At = transpose_matrix_array(A, rows, cols);
        printf("Matrix transposed!\n\n");
        print_matrix_array(At, cols, rows);
        printf("\n");
        fflush(stdout);
    }

    // free unneeded memory for initial matrix
    // MPI_Free(me, A);

    // tutti allocano A locale e il vettore dei risultati
    num = local_c * rows;
    // localAt = malloc(sizeof *localAt * num);
    localAt = malloc(sizeof *localA * num);
    local_w = malloc(sizeof *local_w * local_r);

    // Adesso il master invia a tutti un pezzo della matrice
    err = MPI_Scatter(
        At, num, MPI_DOUBLE,
        // localAt, num, MPI_DOUBLE,
        localAt, num, MPI_DOUBLE,
        master, MPI_COMM_WORLD
    );
    checkMPIerror(&me, &err);

    // free unneeded memory for transposed matrix
    // MPI_Free(me, At);

    if (me == master) {
        MPI_Printf(me, "Matrix scattered!\n");
        fflush(stdout);
    }
    else {
        MPI_Printf(me, "Matrix received!\n");
        fflush(stdout);
    }

    // Print received transposed matrix
    MPI_Printf(me, "localAt:\n");
    print_matrix_array(localAt, local_r, cols);
    printf("\n");
    fflush(stdout);

    // synchronize to start taking time
    // err = MPI_Barrier(MPI_COMM_WORLD);
    // checkMPIerror(&me, &err);

    // la matrice locale ricevuta è in righe
    // va ritrasposta per averla in colonne
    localA = transpose_matrix_array(localAt, local_r, cols);

    // free unneeded memory for original local matrix
    // MPI_Free(me, localAt);

    // Scriviamo la matrice locale ricevuta
    MPI_Printf(me, "localA:\n\n");
    print_matrix_array(localA, local_r, cols);
    printf("\n");
    fflush(stdout);

    // Effettuiamo i calcoli
    prod_mat_vett(local_w, localA, cols, local_r, local_v);

    // free unneeded memory for transposed local matrix
    // MPI_Free(me, localA);

    // free unneeded memory for local vector
    // MPI_Free(me, local_v);

    // mostriamo i risultati parziali
    // err = MPI_Barrier(MPI_COMM_WORLD);
    // checkMPIerror(&me, &err);
    MPI_Printf(me, "Partial result:\n\n");
    print_matrix_array(local_w, cols, 1);
    printf("\n");

    // il master raccoglie i risultati parziali
    err = MPI_Reduce(
        local_w, w, cols, MPI_DOUBLE, MPI_SUM,
        master, MPI_COMM_WORLD
    );
    checkMPIerror(&me, &err);
    if (me == master) {
        MPI_Printf(me, "Results gathered!\n");
    }
    else {
        MPI_Printf(me, "Result sent!\n");
    }

    // free unneeded memory for local result
    // MPI_Free(me, local_w);

    // il master stampa la soluzione
    // err = MPI_Barrier(MPI_COMM_WORLD);
    // checkMPIerror(&me, &err);
    if (me == master) {
        MPI_Printf(me, "Result is w:\n\n"); 
        print_matrix_array(w, cols, 1);
    }
    fflush(stdout);

    // free unneeded memory for result
    // MPI_Free(me, w);

    err = MPI_Finalize();
    checkMPIerror(&me, &err);

    return EXIT_SUCCESS;
}
