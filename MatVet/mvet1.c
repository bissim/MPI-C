#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mpi.h"
#include "mpiutils.h"
#include "matrixutils.h"

extern int master;

int main(int argc, char **argv) {
    int nproc;                      // Numero di processi totale
    int me;                         // Il mio id
    int n;                          // Dimensione della matrice
    int local_n;                    // Dimensione dei dati locali
    int i,j;                        // Iteratori vari
	double T_inizio, T_fine, T_max; // tempi di esecuzione

    // Variabili di lavoro
    double *A, *v, *localA,*local_w, *w;

    /*Attiva MPI*/
    MPI_Init(&argc, &argv);
    /*Trova il numero totale dei processi*/
    MPI_Comm_size (MPI_COMM_WORLD, &nproc);
    /*Da ad ogni processo il proprio numero identificativo*/
    MPI_Comm_rank (MPI_COMM_WORLD, &me);

    MPI_Printf(me, "MPI process started!\n\n");

    // Se sono a radice
    master = 0;
    if (me == master) {
        MPI_Printf(me, "inserire n = "); 
        fflush(stdout);
        scanf("%d",&n);
        printf("\n");

        // Porzione di dati da processare
        local_n = n / nproc;  

        // Alloco spazio di memoria
        A = malloc(n * n * sizeof(double));
        v = malloc(n * sizeof(double));
        w = malloc(n * sizeof(double)); 
 
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
        printf("\n");

        MPI_Printf(me, "A = \n");
        print_matrix_array(A, n, n);
        MPI_Printf(me, "v =\n");
        print_matrix_array(v, n, 1);
        printf("\n");
    }

    // Spedisco n e local_v
    MPI_Bcast(&n, 1, MPI_INT, master, MPI_COMM_WORLD);            
    MPI_Bcast(&local_n, 1, MPI_INT, master, MPI_COMM_WORLD);            

    // Se sono un figlio alloco v
    if (me != master) {
        v = (double *) malloc(n * sizeof(double));
    }

    MPI_Bcast(&v[0], n, MPI_DOUBLE, master, MPI_COMM_WORLD);           

    // tutti allocano A locale e il vettore dei risultati
    localA = (double *) malloc(local_n * n * sizeof(double));
    local_w = (double *) malloc(local_n * sizeof(double));

    // Adesso 0 invia a tutti un pezzo della matrice
    int num = local_n * n;
    MPI_Scatter(
        A, num, MPI_DOUBLE,
        localA, num, MPI_DOUBLE,
        master, MPI_COMM_WORLD
    );

    // synchronize processes at
    // local array receival
    MPI_Barrier(MPI_COMM_WORLD);

    // start computation time
	T_inizio = MPI_Wtime();

    // Scriviamo la matrice locale ricevuta
    MPI_Printf(me, "localA:\n");
    print_matrix_array(localA, local_n, n);
    printf("\n");

    // Effettuiamo i calcoli
    prod_mat_vett(local_w, localA, local_n, n, v);
        
    // 0 raccoglie i risultati parziali
    MPI_Gather(&local_w[0], local_n, MPI_DOUBLE, &w[0], local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD); // sync
	T_fine = MPI_Wtime() - T_inizio; // end time

    // calculate total execution time
	MPI_Reduce(
        &T_fine, &T_max, 1, MPI_DOUBLE,
        MPI_MAX, master, MPI_COMM_WORLD
    );

    // master stampa la soluzione
    if (me == master) {
        MPI_Printf(me, "w =\n"); 
        print_matrix_array(w, n, 1);
        printf("\n-----\n");
		printf("\nProcessori impegnati: %d\n", nproc);
		printf("\nTempo calcolo locale: %lf\n", T_fine);
		printf("\nMPI_Reduce max time: %f\n", T_max);
    }
    fflush(stdout);

    MPI_Finalize(); /* Disattiva MPI */

    return 0;
}
