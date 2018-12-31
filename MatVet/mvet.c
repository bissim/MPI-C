#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mpi.h"

/**
 * Funzione che esegue la stampa di una matrice di interi
 */
void print_matrix_array(double *array, int rows, int columns) {
    int i, j;

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
void prod_mat_vett(double w[], double *a, int ROWS, int COLS, double v[]) {
    int i, j;
    
    for (i = 0; i < ROWS; i++) {
        w[i] = 0;
        for (j = 0; j < COLS; j++) {
            w[i] += a[i*COLS+j]* v[j];
        }
    }
}

int main(int argc, char **argv) {
    int nproc;              // Numero di processi totale
    int me;                 // Il mio id
    int n;                  // Dimensione della matrice
    int local_n;            // Dimensione dei dati locali
    int i,j;                // Iteratori vari 

    // Variabili di lavoro
    double *A, *v, *localA,*local_w, *w;

    /*Attiva MPI*/
    MPI_Init(&argc, &argv);
    /*Trova il numero totale dei processi*/
    MPI_Comm_size (MPI_COMM_WORLD, &nproc);
    /*Da ad ogni processo il proprio numero identificativo*/
    MPI_Comm_rank (MPI_COMM_WORLD, &me);

    // Se sono a radice
    if(me == 0) {
        printf("inserire n = "); 
        fflush(stdout);
        scanf("%d",&n);
        printf("\n");

        // Porzione di dati da processare
        local_n = n/nproc;  
        
        // Alloco spazio di memoria
        A = malloc(n * n * sizeof(double));
        v = malloc(n * sizeof(double));
        w = malloc(n * sizeof(double)); 
        
        printf("A = \n"); 
        for (i = 0; i < n; i++) {
            v[i] = i;
            for (j = 0; j < n; j++) {
                if (j == 0) {
                    A[i*n+j] = 1.0/(i + 1) - 1;
                }
                else {
                    A[i*n+j] = 1.0/(i+1) - pow(1.0/2.0, j); 
                }
                printf("%6.3f ", A[i*n+j]);
            }
            printf("\n");
        }
        printf("\n");

        printf("v = \n");
        print_matrix_array(v, n, 1);
        printf("\n");
    }


    // Spedisco n e local_v
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);            
    MPI_Bcast(&local_n, 1, MPI_INT, 0, MPI_COMM_WORLD);            

    // Se sono un figlio alloco v
    if (me != 0) {
        v = malloc(n * sizeof(double));
    }

    MPI_Bcast(&v[0], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);            

    // tutti allocano A locale e il vettore dei risultati
    localA = malloc(local_n * n * sizeof(double));
    local_w = malloc(local_n * sizeof(double));

    // Adesso 0 invia a tutti un pezzo della matrice
    int num = local_n * n;
    MPI_Scatter(
        &A[0], num, MPI_DOUBLE,
        &localA[0], num, MPI_DOUBLE,
        0, MPI_COMM_WORLD
    );

    // Scriviamo la matrice locale ricevuta
    printf("[P%d] localA:\n", me);
    print_matrix_array(localA, local_n, n);
    printf("\n");

    // Effettuiamo i calcoli
    prod_mat_vett(local_w, localA, local_n, n, v);
        
    // 0 raccoglie i risultati parziali
    MPI_Gather(&local_w[0], local_n, MPI_DOUBLE, &w[0], local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // 0 stampa la soluzione
    if (me == 0) {
        printf("w =\n"); 
        print_matrix_array(w, n, 1);
        printf("\n");
    }
    fflush(stdout);

    MPI_Finalize(); /* Disattiva MPI */

    return 0;
}
