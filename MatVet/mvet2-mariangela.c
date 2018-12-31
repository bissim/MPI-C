#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h> 


/**
* Funzione che esegue il prodotto matrice vettore
*/
void prod_mat_vett(double w[], double *a, int rows, int cols, double v[])
{
    int i, j;
    
    for(i=0;i<rows;i++)
    {
        w[i]=0;
        for(j=0;j<cols;j++)
        { 
            w[i] += a[i*cols+j]* v[j];
        } 
    }    
}

void print_matrix_array(double *array, int rows, int columns) {
    int i;

    for (i = 0; i < rows * columns; ++i) {
        printf("%5.3f \t", array[i]);
        if (i % columns == columns - 1) {
            printf("\n\n");
        }
    }
}

double *transpose_matrix_array(double *array, int rows, int columns) {
    int i, j;
    int source_position, destination_position;
    double source, *transpose;

    // allocate space for transpose matrix
    transpose = (double *) calloc(columns * rows, sizeof(double));

    for (i = 0; i < rows; ++i) {
        // printf("ROW %d\n", i);
        for (j = 0; j < columns; ++j) {
            // printf("Element %d of row %d\n", j, i);
            //source_position = i * columns + j;
            //destination_position = i*columns + j;
            source_position = columns * i + j;
            destination_position = rows * j + i;
            source = *(array + source_position);
            // printf("Transposing element %d (%d%d) from array cell %d to array cell %d\n", source, i, j, source_position, destination_position);
            *(transpose + destination_position) = source;
            // printf("Cell %d (%d%d) value of transposed matrix is %d\n", destination_position, i, j, transpose[destination_position]);
        }
        // printf("\n");
    }

    return transpose;
}

int main(int argc, char **argv) {

int nproc;              // Numero di processi totale
int me;                 // Il mio id
int n;                  // Dimensione della matrice
int local_n;            // Dimensione dei dati locali
int i,j;                    // Iteratori vari 

// Variabili di lavoro
double *A, *v, *localA,*local_w, *w, *local_v, *At, *local_At;
//Variabili per il tempo
double T_inizio,T_fine,T_max;


/*Attiva MPI*/
MPI_Init(&argc, &argv);
/*Trova il numero totale dei processi*/
MPI_Comm_size (MPI_COMM_WORLD, &nproc);
/*Da ad ogni processo il proprio numero identificativo*/
MPI_Comm_rank (MPI_COMM_WORLD, &me);

printf("Processo %d \n",me);

// Master
if(me == 0){
    printf("Inserire n: \n"); 
    scanf("%d",&n); 
    // Porzione di dati da processare
    local_n = n/nproc;  
    
    // Allocazione spazio di memoria
    A = malloc(n * n * sizeof(double));
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
    printf("\n"); 
    //stampa la matrice solo nel caso sia al massimo 10*10
    if(n<10){
        printf("Matrice A: \n");
        print_matrix_array(A,n,n);
    	printf("\n"); 

    	printf("v[] \n"); 
    	for (i=0;i<n;i++){   
        	printf("%f\n", v[i]);
    	}
    	printf("\n");
    } 
}     


// Spedisco n e local_v
MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);            
MPI_Bcast(&local_n,1,MPI_INT,0,MPI_COMM_WORLD);            

//tutti i processi allocano local_v 
local_v = malloc(sizeof(double)*local_n); //local_n è il numero di elementi del vettore v 
											//che è uguale al numero di colonne che inviamo
    
//scatter vettore per righe
MPI_Scatter(&v[0],local_n, MPI_DOUBLE, local_v, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

// tutti allocano A locale e il vettore dei risultati
localA = malloc(local_n * n * sizeof(double));
local_w = malloc(local_n * sizeof(double));

if(me==0){
	//At= malloc(n * n * sizeof(double));// lo spazio è allocato già dalla funzione di trasposizione
	At = transpose_matrix_array(A,n,n);
	if(n<10){
		printf("A traspose  \n"); 
		print_matrix_array(At,n,n);
	}	
	printf("\n");
}

local_At = malloc(n * local_n * sizeof(double));
// Adesso 0 invia a tutti un pezzo della matrice
int num = local_n*n;
MPI_Scatter(
    &At[0], num, MPI_DOUBLE,
    &local_At[0], num, MPI_DOUBLE,
    0, MPI_COMM_WORLD);

// Scriviamo la matrice locale ricevuta se sono al massimo 30 elementi
if(num<30){
	printf("local_At %d: \n", me); 
	for(i = 0; i < local_n; i++){
    	for(j = 0; j < n; j++)
        printf("%5.3lf\t", local_At[i*n+j]);
    printf("\n");
	}
}
printf("\n");

/* sincronizzazione dei processori del contesto MPI_COMM_WORLD*/
	
MPI_Barrier(MPI_COMM_WORLD);
 
T_inizio=MPI_Wtime(); //inizio del cronometro per il calcolo del tempo di inizio

//printf("\n calcoliamo prod_mat_vett [%d] \n",me);
// Effettuiamo la trasposta prima di fare i calcoli
localA = transpose_matrix_array(local_At,local_n,n);
//Calcoliamo il prodotto mat-vet
prod_mat_vett(local_w,localA,n,local_n,local_v);
//printf("\n prodotto calcolato [%d] \n",me);

// 0 raccoglie i risultati parziali
MPI_Reduce(&local_w[0],&w[0],n,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

MPI_Barrier(MPI_COMM_WORLD); // sincronizzazione
T_fine=MPI_Wtime()-T_inizio; // calcolo del tempo di fine
 
/* calcolo del tempo totale di esecuzione*/
MPI_Reduce(&T_fine,&T_max,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

/*stampa a video dei risultati finali*/
if(me==0){	
	printf("\n Il risultato è w [] \n"); 
    for(i = 0; i < n; i++)
        printf("%5.3f ", w[i]);
    printf("\n");
	printf("\nProcessori impegnati: %d\n", nproc);
	printf("\nTempo calcolo locale: %lf\n", T_fine);
	printf("\nMPI_Reduce max time: %f\n",T_max);
}// end if
 
MPI_Finalize (); /* Disattiva MPI */
return 0;  
}
