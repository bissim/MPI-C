#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mpi.h"

void print_matrix_array(float *, int, int);
void prod_mat_vett(float *, float *, int, int, float*);
float *transpose_matrix_array(float *, int, int);
void transpose_matrix_array_2(float *, float *, int, int);
void checkMPIerror(int *, int *);

int main(int argc, char **argv) {
    // dichiarazione variabili
    int nproc, me;
    int n, numcells, sub_n, sub_nsq;
    int i, j;
    float *A, *subA, *subAt, *localA, *local_At, *v, *local_v, *w, *local_w, *row_w;
    // variabili di lavoro MPI
    MPI_Comm comm_grid, comm_rows, comm_cols;
    int err, dim, reorder, coord_rows, coord_cols, me_grid, me_row, me_col;
    int *ndim, *period, *coord, *belongs;
    double T_inizio, T_fine, T_max;
    int flag = 0;

    // costanti
    const int root = 0;
    const int MS_IN_S = 1000;
    const MPI_Comm comm = MPI_COMM_WORLD;

    // inizializzazione puntatori
    A = subA = subAt = localA = local_At = v = local_v = w = local_w = row_w = NULL;
    ndim = period = coord = belongs = NULL;

    // inizializzazione MPI
    err = MPI_Init(&argc, &argv);
    checkMPIerror(&me, &err);
    err = MPI_Comm_size(comm, &nproc);
    checkMPIerror(&me, &err);
    err = MPI_Comm_rank(comm, &me);
    checkMPIerror(&me, &err);

    if (me == 0) {
        if (argc < 3) {
            printf("Usage: <programName> <rows> <printDebug>\n");
            fflush(stdout);

            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
        else {
            n = atoi(argv[1]);
            flag = atoi(argv[2]);
            if (flag) {
                printf("Debug print enabled.\n");
                fflush(stdout);
            }
        }
    }

    flag = 0; // TODO for debug purposes, set to 1
    if (flag) {
        printf("[P%d] process started.\n\n", me);
        fflush(stdout);
    }

    // creazione matrice e vettore
    if (me == root) {
        // Porzione di dati da processare
        sub_n = n / sqrt(nproc); // if nproc is n, there are sqrt(n) procs for row and col
        // sub_n = n / nproc;
        
        // Alloco spazio di memoria
        A = malloc(sizeof *A * n * n);
        v = malloc(sizeof *v * n);

        // Create matrix A and vector v
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
            printf("\n\n");

            printf("\nVector v:\n");
            print_matrix_array(v, n, 1);
            printf("\n\n");

            fflush(stdout);
        }
    }

    // 1. partizionamento per righe della matrice e distribuzione del vettore
    dim = 2; // topology dimensions
    ndim = calloc(dim, sizeof *ndim);
    ndim[0] = ndim[1] = 2; // TODO modify for non-square matrices
    period = calloc(dim, sizeof *period);
    period[0] = period[1] = 0; // no periodicity over rows nor columns
    reorder = 0; // disable MPI topology optimizations
    if (flag && me == root) {
        printf("Creating topology...");
    }

    // 1.1 creazione topologia cartesiana a griglia non toroidale
    err = MPI_Cart_create(comm, dim, ndim, period, reorder, &comm_grid); // create communication group
    err = MPI_Comm_rank(comm_grid, &me_grid); // process rank into topology
    coord = malloc(sizeof *coord * dim);
    err = MPI_Cart_coords(comm_grid, me_grid, dim, coord); // process coordinates
    checkMPIerror(&me, &err);
    if (flag) {
        printf("[P%d] grid rank %d, coordinates (%d, %d)\n", me, me_grid, coord[0], coord[1]);
        printf("\n");
    }
    // fflush(stdout);
    // if (me == root) {
    //     free(ndim);
    //     free(period);
    // } 

    // invia n e sub_n ai processi
    err = MPI_Bcast(&n, 1, MPI_INT, root, comm_grid);
    checkMPIerror(&me, &err);
    err = MPI_Bcast(&sub_n, 1, MPI_INT, root, comm_grid);
    checkMPIerror(&me, &err);
    if (flag) {
        printf("[P%d] I received n = %d, sub_n = %d\n\n", me, n, sub_n);
        fflush(stdout);
    }

    // 1.2 creazione sottogruppi per righe
    belongs = malloc(sizeof *belongs * dim);
    belongs[0] = 0;
    belongs[1] = 1; // fisso l'indice di riga
    err = MPI_Cart_sub(comm_grid, belongs, &comm_rows);
    err = MPI_Comm_rank(comm_rows, &me_row);
    err = MPI_Cart_coords(comm_rows, me_row, 1, &coord_rows);
    checkMPIerror(&me, &err);
    if (flag) {
        printf("[P%d] (%d, %d) row id %d, row coordinate (%d)\n\n", me_grid, coord[0], coord[1], me_row, coord_rows);
        fflush(stdout);
    }

    // 1.3 creazione sottogruppi per colonne
    belongs[0] = 1; // fisso l'indice di colonna
    belongs[1] = 0;
    err = MPI_Cart_sub(comm_grid, belongs, &comm_cols);
    err = MPI_Comm_rank(comm_cols, &me_col);
    err = MPI_Cart_coords(comm_cols, me_col, 1, &coord_cols);
    checkMPIerror(&me, &err);
    if (flag) {
        printf("[P%d] (%d, %d) column id %d, column coordinate (%d)\n\n", me_grid, coord[0], coord[1], me_col, coord_cols);
        fflush(stdout);
    }
    // if (me == root) {
    //     free(belongs);
    // }

    // sincronizza i processi per far assumere a tutti
    // le proprie coordinate prima di proseguire
    err = MPI_Barrier(comm);
    checkMPIerror(&me, &err);

    // 1.4 scatter matrice tra i processi della prima colonna
    if (flag) {
        printf("[P%d] sub_n = %d\n", me, sub_n);
        fflush(stdout);
    }
    numcells = sub_n * n;
    if (coord[1] == 0) { // processi della prima colonna, hanno l'indice di riga a 0
        subA = calloc(numcells, sizeof *subA);
        if (flag) {
            printf("[P%d] (%d, %d) memory space allocated for subA\n\n", me, coord[0], coord[1]);
            fflush(stdout);
        }
        err = MPI_Scatter(
            A, numcells, MPI_FLOAT, // send numcells cells of A to every P
            subA, numcells, MPI_FLOAT,
            root, comm_cols // every P in the column of root P
        );
        checkMPIerror(&me, &err);
        // verifica la ricezione della sottomatrice
        if (flag) {
            printf("[P%d] (%d, %d) subA:\n", me, coord[0], coord[1]);
            print_matrix_array(subA, sub_n, n);
            printf("\n");
            fflush(stdout);
        }
    }
    // if (me == root) {
    //     free(A);
    // }

    // 1.5 scatter vettore tra i processi della prima riga
    local_v = calloc(sub_n, sizeof *local_v);
    if (coord[0] == 0) {
        err = MPI_Scatter(
            v, sub_n, MPI_FLOAT,
            local_v, sub_n, MPI_FLOAT,
            root, comm_rows // distribuisci lungo le righe
        );
        checkMPIerror(&me, &err);
        if (flag) {
            printf("[P%d] (%d, %d) v scattered to first row processes!\n\n", me, coord[0], coord[1]);
            fflush(stdout);
        }
    }
    // if (me == root) {
    //     free(v);
    // }
    // verifica la ricezione del vettore
    if (flag && coord[0] == 0) {
        printf("[P%d] (%d, %d) local_v:\n", me, coord[0], coord[1]);
        print_matrix_array(local_v, sub_n, 1);
        printf("\n");
        fflush(stdout);
    }

    // 2. partizionamento per colonne delle sottomatrici e per righe del vettore
    // trasponi le sottomatrici
    if (coord[1] == 0) {
        subAt = calloc(numcells, sizeof *subAt); // sub_n * n
        if (flag) {
            printf("[P%d] (%d, %d) memory space allocated for subAt\n\n", me, coord[0], coord[1]);
            fflush(stdout);
        }
    }
    sub_nsq = sub_n * sub_n;
    local_At = calloc(sub_nsq, sizeof *local_At);
    if (flag) {
        printf("[P%d] (%d, %d) memory space allocated for local_At\n\n", me, coord[0], coord[1]);
        fflush(stdout);
    }
    if (coord[1] == 0) { // prima colonna
        transpose_matrix_array_2(subA, subAt, sub_n, n);
        if (flag) {
            printf("[P%d] (%d, %d) subAt:\n", me, coord[0], coord[1]);
            print_matrix_array(subAt, n, sub_n);
            printf("\n");
            fflush(stdout);
        }
    }
    err = MPI_Scatter(
        subAt, sub_nsq, MPI_FLOAT,
        local_At, sub_nsq, MPI_FLOAT,
        0, comm_rows // distribuisci lungo le righe dal primo elemento della riga
    );
    checkMPIerror(&me, &err);
    // if (me == root) {
    //     free(subA);
    //     free(subAt);
    // }
    // verifica ricezione blocchi di matrice
    // locale da tutti i processi
    if (flag) {
        printf("[P%d] (%d, %d) local_At:\n", me, coord[0], coord[1]);
        print_matrix_array(local_At, sub_n, sub_n);
        printf("\n");
        fflush(stdout);
    }
    // ritrasponi i blocchi di matrice
    localA = malloc(sizeof *localA * sub_nsq);
    transpose_matrix_array_2(local_At, localA, sub_n, sub_n);
    if (flag) {
        printf("[P%d](%d, %d) localA:\n", me, coord[0], coord[1]);
        print_matrix_array(localA, sub_n, sub_n);
        printf("\n");
        fflush(stdout);
    }
    // if (me == root){
    //     free(local_At);
    // }
    // broadcast vettore
    err = MPI_Bcast(
        local_v, sub_n, MPI_FLOAT,
        0, comm_cols
    );
    checkMPIerror(&me, &err);
    // verifica ricezione vettore
    if (flag) {
        printf("[P%d](%d, %d) local_v:\n", me, coord[0], coord[1]);
        print_matrix_array(local_v, sub_n, 1);
        printf("\n");
        fflush(stdout);
    }

    // 3. calcolo risultati parziali
    local_w = malloc(sizeof *local_w * sub_n);
    T_inizio = MPI_Wtime(); // inizio del cronometro per il calcolo del tempo di inizio
    prod_mat_vett(local_w, localA, sub_n, sub_n, local_v);
    // verifica calcoli parziali
    if (flag) {
        printf("[P%d](%d, %d) local_w:\n", me, coord[0], coord[1]);
        print_matrix_array(local_w, sub_n, 1);
        printf("\n");
        fflush(stdout);
    }
    // if (me == root) {
    //     free(localA);
    //     free(local_v);
    // }

    // 4. reduce risultati di ogni riga
    // if (coord[1] == 0) {
        row_w = malloc(sizeof *row_w * sub_n);
    // }
    err = MPI_Reduce(
        local_w, row_w, sub_n, MPI_FLOAT,
        MPI_SUM, 0, comm_rows
    );
    checkMPIerror(&me, &err);
    // if (me == root) {
    //     free(local_w);
    // }
    if (flag && coord[1] == 0) {
        printf("[P%d](%d, %d) row_w:\n", me, coord[0], coord[1]);
        print_matrix_array(row_w, sub_n, 1);
        printf("\n");
        fflush(stdout);
    }

    // 5. gather dei sottovettori parziali
    // if (me == root) {
        w = malloc(sizeof *w * n);
    // }
    err = MPI_Gather( // TODO controllare
        row_w, sub_n, MPI_FLOAT,
        w, sub_n, MPI_FLOAT,
        root, comm_cols
    );
    checkMPIerror(&me, &err);

    T_fine = MPI_Wtime() - T_inizio; // calcolo del tempo di fine

    if (flag && me == root) {
        // free(row_w);
        // stampa w
        printf("w:\n");
        print_matrix_array(w, n, 1);
        printf("\n");
        fflush(stdout);
        // free(w);

        printf("Press any key to quit...\n");
        getchar();
    }

    /* calcolo del tempo totale di esecuzione*/
    MPI_Reduce(
        &T_fine, &T_max, 1, MPI_DOUBLE, MPI_MAX,
        root, MPI_COMM_WORLD
    );

    if (me == root) {
        printf("\n");
        printf("Tempo calcolo locale: %.3lf ms\n", T_fine * MS_IN_S);
        printf("MPI_Reduce max time: %.3f ms\n", T_max * MS_IN_S);
    }

    // terminazione MPI
    err = MPI_Finalize();
    checkMPIerror(&me, &err);

    return EXIT_SUCCESS;
}

/**
 * Funzione che esegue la stampa di una matrice di interi
 */
void print_matrix_array(float *array, int rows, int columns) {
    int i;

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
void prod_mat_vett(float result[], float *a, int rows, int cols, float v[]) {
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
float *transpose_matrix_array(float *array, int rows, int columns) {
    int i, j;
    int source_position, destination_position;
    float *transpose, source;

    // print input matrix
    // printf("Print input matrix to transpose:\n");
    // print_matrix_array(array, rows, columns);

    // allocate space for transpose matrix
    transpose = (float *) calloc(columns * rows, sizeof(float));

    for (i = 0; i < rows; ++i) {
        // printf("ROW %d\n", i);
        for (j = 0; j < columns; ++j) {
            // printf("Element %d of row %d\n", j, i);
            source_position = columns * i + j;
            destination_position = rows * j + i;
            source = *(array + source_position);
            // printf("Transposing element %6.3f (%d%d) from array cell %d to array cell %d\n", source, i, j, source_position, destination_position);
            *(transpose + destination_position) = source;
            // printf("Cell %d (%d%d) value of transposed matrix is %f\n", destination_position, i, j, transpose[destination_position]);
        }
        // printf("\n");
    }

    // print transposed matrix
    // printf("Print transposed matrix into function:\n");
    // print_matrix_array(transpose, columns, rows);

    return transpose;
}

/**
 * Funzione che esegue la trasposizione di una matrice di interi
 */
void transpose_matrix_array_2(float *array, float *transpose, int rows, int columns) {
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

    return;
}

void checkMPIerror(int *process, int *error) {
    if (*error != MPI_SUCCESS) {
        fprintf(stderr, "[P%d] MPI call failed!\nError %d.\n", *process, *error);
        fflush(stderr);
    }

    return;
}
