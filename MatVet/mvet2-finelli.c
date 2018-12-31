#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h> 


/**
* Funzione che esegue il prodotto matrice vettore
*/
void prod_mat_vett(double w[], double *a, int ROWS, int COLS, double v[])
{
    int i, j;
    
    for(i=0;i<ROWS;i++)
    {
        w[i]=0;
        for(j=0;j<COLS;j++)
        { 
            w[i] += a[j*ROWS+i]* v[j];
        } 
    } 
}


int main(int argc, char **argv) {

int nproc;                      // Numero di processi totale
int me;                         // Il mio id
int n;                          // Dimensione della matrice
int local_n;                    // Dimensione dei dati locali
int i,j;                        // Iteratori vari 
double T_inizio,T_fine,T_max;   //Variabili per il calcolo del tempo di esecuzione

// Variabili di lavoro
double *A, *v, *localA, *local_v, *local_w, *w, *T;


/*Attiva MPI*/
MPI_Init(&argc, &argv);
/*Trova il numero totale dei processi*/
MPI_Comm_size (MPI_COMM_WORLD, &nproc);
/*Da ad ogni processo il proprio numero identificativo*/
MPI_Comm_rank (MPI_COMM_WORLD, &me);

// Se sono a radice
if(me == 0)
{
    printf("inserire la dimensione della matrice quadrata = \n"); 
    scanf("%d",&n); 
    // Porzione di dati da processare
    local_n = n/nproc;  
    
    // Alloco spazio di memoria
    A = malloc(n * n * sizeof(double));
    T = malloc(n * n * sizeof(double));
    v = malloc(sizeof(double)*n);
    w = malloc(sizeof(double)*n); 
    
    for (i=0;i<n;i++)
    {
        v[i]=i;  
        for(j=0;j<n;j++)
        {
            if (j==0)
                A[i*n+j]= 1.0/(i+1)-1;
            else
                A[i*n+j]= 1.0/(i+1)-pow(1.0/2.0,j);
        }
    }

    //Effettuiamo la trasposta della matrice
    for (i=0; i<n; i++)
        for( j=0; j<n; j++ )
            T[i*n+j] = A[j*n+i];

    if(n <= 20) {
        // Scriviamo la matrice A
        printf("A = \n"); 
        for (i=0;i<n;i++)
        {
            for(j=0;j<n;j++)
            {   
                if(A[i*n+j] >= 0){
                    printf(" %.4f  ", A[i*n+j] );
                } else {
                    printf("%.4f  ", A[i*n+j] );
                }
            }
            printf("\n");
        }

        printf("\n"); 

        printf("v = \n"); 
        for (i=0;i<n;i++)
        {   
            printf("%f\n", v[i]);
        }
        printf("\n");

        printf("T = \n"); 
        for (i=0;i<n;i++)
        {
            for(j=0;j<n;j++)
            {   
                if(T[i*n+j] >= 0){
                    printf(" %.4f  ", T[i*n+j] );
                } else {
                    printf("%.4f  ", T[i*n+j] );
                }
            }
            printf("\n");
        }

        printf("\n"); 
    }
}     


// Spedisco n e local_v
MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);            
MPI_Bcast(&local_n,1,MPI_INT,0,MPI_COMM_WORLD);

// tutti allocano A locale e il vettore dei risultati
localA = malloc(n * local_n * sizeof(double));
local_v = malloc(local_n * sizeof(double));
local_w = malloc(n * sizeof(double));

//Invio del vettore v(x) ai vari processi
MPI_Scatter(&v[0],local_n,MPI_DOUBLE,&local_v[0],local_n,MPI_DOUBLE,0,MPI_COMM_WORLD); 

// Adesso 0 invia a tutti un pezzo della matrice
int num = local_n*n;
MPI_Scatter(
    &T[0], num, MPI_DOUBLE,
    &localA[0], num, MPI_DOUBLE,
    0, MPI_COMM_WORLD);

if(n <= 20) {
    // Scriviamo la matrice locale ricevuta
    printf("localA %d = \n", me); 
    for(i = 0; i < local_n; i++)
    {
        for(j = 0; j < n; j++) {

            if(localA[i*n+j] >= 0){
                printf(" %.4lf\t", localA[i*n+j]);
            } else {
                printf("%.4lf\t", localA[i*n+j]);
            }
        }
        printf("\n");
    }
}

T_inizio=MPI_Wtime(); //inizio del cronometro per il calcolo del tempo di inizio

// Effettuiamo i calcoli
prod_mat_vett(local_w,localA,n,local_n,local_v);
    
// 0 raccoglie i risultati parziali
MPI_Reduce(&local_w[0],&w[0],n,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

T_fine=MPI_Wtime()-T_inizio; // calcolo del tempo di fine
 
/* calcolo del tempo totale di esecuzione*/
MPI_Reduce(&T_fine,&T_max,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

// 0 stampa la soluzione
if(me==0)
{ 
    if(n <= 400) {
        printf("\n\nw = \n"); 
        for(i = 0; i < n; i++)
            printf("%.4f ", w[i]);
        printf("\n");
    }
    printf("\nTempo calcolo locale: %lf\n", T_fine);
    printf("\nMPI_Reduce max time: %f\n",T_max);
}

MPI_Finalize (); /* Disattiva MPI */
return 0;  
}
