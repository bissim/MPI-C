#include <stdio.h>  
#include <stdlib.h>

#include "mpi.h"

/*
 * 
 * 
 * 
 */
int main(int argc, char **argv) {
	int menum, nproc, menum_grid, row, col;
	int dim, *ndim, reorder, *period;
	int coordinate[2];
	int value;
	MPI_Comm comm_grid; /* definizione di tipo communicator */  

	MPI_Init(&argc, &argv); 
	MPI_Comm_rank(MPI_COMM_WORLD, &menum);  
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	/* Numero di righe della griglia di processo */  
	if (menum == 0)	{	
		printf("Numero di righe della griglia ");
		fflush(stdout);
		scanf("%d", &row);
	}

	/* Spedizione di row da parte di 0 a tutti i processori */  
	MPI_Bcast(&row, 1, MPI_INT, 0, MPI_COMM_WORLD);
	col = nproc / row;	 /* Numero di colonne della griglia */  
	dim = 2;     	 /* Numero di dimensioni della griglia */

	/* vettore contenente le lunghezze di ciascuna dimensione*/  
	ndim = (int*) calloc(dim, sizeof(int));
	ndim[0] = row;
	ndim[1] = col;

	/* vettore contenente la periodicità delle dimensioni */  
	period = (int*) calloc(dim, sizeof(int));
	period[0] = period[1] = 0; 
	reorder = 0;

	/* Definizione della griglia bidimensionale */
	MPI_Cart_create(MPI_COMM_WORLD, dim, ndim, period, reorder, &comm_grid);
	MPI_Comm_rank(comm_grid, &menum_grid);

	/* Definizione delle coordinate di ciascun processo  nella griglia bidimensionale */
	MPI_Cart_coords(comm_grid, menum, dim, coordinate);

	/* Stampa delle coordinate */
//	printf("Processore %d coordinate nella griglia  (%d,%d) \n", menum, *coordinate, *(coordinate + 1));

	/* generazione valori per processore */
	if (coordinate[0] == coordinate[1]) {
		value = coordinate[0] * coordinate[1];
	}
	else {
		value = coordinate[0] + coordinate[1];
	}

	/* Print generated value */
	printf("Process (%d, %d): value generated is %d.\n", coordinate[0], coordinate[1], value);

	MPI_Finalize();

	return 0;
}
