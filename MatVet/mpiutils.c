#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"
#include "mpiutils.h"

/**
 * Funzione che stampa un messaggio con il
 * rank del processo come prefisso
 */
void MPI_Printf(int process, char *message) {
    printf("[P%d] %s", process, message);
}

/**
 * Funzione per liberare memoria
 * dal processo master
 */
void MPI_Free(int rank, void *memory) {
    if (rank == master && memory != NULL) {
        free(memory);
        memory = NULL;
    }
    else if (rank == master && memory == NULL) {
        MPI_Printf(rank, "Trying to free null memory!\n");
    }
}

/**
 * Funzione per controllare l'esito
 * di una chiamata a funzione MPI
 */
void checkMPIerror(int *process, int *error) {
    if (*error != MPI_SUCCESS) {
        fprintf(stderr, "[P%d] MPI call failed!\nError %d.\n", *process, *error);
        fflush(stderr);

        exit(EXIT_FAILURE);
    }
    // else {
    //     MPI_Printf(*process, "MPI call successful!\n");
    // }
}
