#ifndef MPIUTILS_H_
#define MPIUTILS_H_

int master;

void MPI_Printf(int, char *);
void MPI_Free(int, void *);
void checkMPIerror(int *, int *);

#endif // MPIUTILS_H_
