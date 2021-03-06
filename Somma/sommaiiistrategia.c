/*----------------------------------------------------------------------------
  Inclusione dei file di intestazione usati
----------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
/*--------------------------------------------------------------------------
  Inclusione del file che contiene le definizioni necessarie al preprocessore
  per l'utilizzo di MPI.
-----------------------------------------------------------------------------*/
#include "mpi.h"

int main(int argc, char **argv) {
	/*dichiarazioni variabili*/
    int menum, nproc, tag;
	int n, nloc, i, resto, nlocgen;
	int ind, p, r, sendTo, recvBy;
	int *dist, passi = 0;
	float *vett, *vett_loc;
	float sommaloc = 0.0, tmp;
	double T_inizio, T_fine, T_max;

	MPI_Status info;

	/*Inizializzazione dell'ambiente di calcolo MPI*/
	MPI_Init(&argc, &argv);
	/*assegnazione IdProcessore a menum*/
	MPI_Comm_rank(MPI_COMM_WORLD, &menum);
	/*assegna numero processori a nproc*/
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	/* lettura e inserimento dati*/
	if (menum == 0) {
		if (argc < 2) {
			printf("[P%d] Numero di parametri insufficiente!", menum);
			exit(EXIT_FAILURE);
		}
		else {
			n = atoi(argv[1]);
		}
		printf("[P%d] Addizione di %d numeri reali.\n", menum, n);

	// 	printf("Inserire il numero di elementi da sommare: ");
	// 	fflush(stdout);
	// 	scanf("%d",&n);
		
    	vett = (float*) calloc(n, sizeof(float));
	}

	/*invio del valore di n a tutti i processori appartenenti a MPI_COMM_WORLD*/
	MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);

	/*numero di addendi da assegnare a ciascun processore*/
	nlocgen = n / nproc; // divisione intera
	
    resto = n % nproc; // resto della divisione

	/* Se resto è non nullo, i primi resto processi ricevono un addento in più */
	if (menum < resto) {
		nloc = nlocgen + 1;
	}
	else {
		nloc = nlocgen;
	}

    /*allocazione di memoria del vettore per le somme parziali */
	vett_loc = (float *) calloc(nloc, sizeof(float));

	// P0 genera e assegna gli elementi da sommare
    
	if (menum == 0) {
        /*Inizializza la generazione random degli addendi utilizzando l'ora attuale del sistema*/                
        srand((unsigned int) time(0)); 
		
        for (i = 0; i < n; i++) {
			/*creazione del vettore contenente numeri casuali */
			*(vett + i) = ((float)rand() / (float)(RAND_MAX / 50)) + 1;
		}
		
   		// Stampa del vettore che contiene i dati da sommare, se sono meno di 100 
		// if (n < 100) {
		// 	for (i = 0; i < n; i++) {
		// 		printf("\nElemento %d del vettore: %5.2f\n", i, *(vett + i));
		// 	}
        // }

		// assegnazione dei primi addendi a P0
        for (i = 0; i < nloc; i++) {
			*(vett_loc + i) = *(vett + i);
		}
  
  		// ind è il numero di addendi già assegnati     
		ind = nloc;
        
		/* P0 assegna i restanti addendi agli altri processori */
		for (i = 1; i < nproc; i++) {
			tag = i; /* tag del messaggio uguale all'id del processo che riceve*/
			/*SE ci sono addendi in sovrannumero da ripartire tra i processori*/
            if (i < resto) {
				/*il processore P0 gli invia il corrispondete vettore locale considerando un addendo in piu'*/
				MPI_Send(vett+ind, nloc, MPI_INT, i, tag, MPI_COMM_WORLD);
				ind = ind + nloc;
			}
			else {
				/*il processore P0 gli invia il corrispondete vettore locale*/
				MPI_Send(vett + ind, nlocgen, MPI_INT, i, tag, MPI_COMM_WORLD);
				ind = ind + nlocgen;
			}// end if
		}//end for
	}

    /*SE non siamo il processore P0 riceviamo i dati trasmessi dal processore P0*/
    else {
		// tag è uguale numero di processore che riceve
		tag = menum;

		/*fase di ricezione*/
		MPI_Recv(vett_loc, nloc, MPI_INT, 0, tag, MPI_COMM_WORLD, &info);
	}// end if

	/* sincronizzazione dei processori del contesto MPI_COMM_WORLD*/	
	MPI_Barrier(MPI_COMM_WORLD);

	T_inizio = MPI_Wtime(); //inizio del cronometro per il calcolo del tempo di inizio

	for (i = 0; i < nloc; i++) {
		/* ogni processore effettua la somma parziale */
		sommaloc = sommaloc + *(vett_loc + i);
	}

	// calcolo di p=log_2 (nproc)
	p = nproc;

	while (p != 1) {
		/* shift di un bit a destra */
		p = p >> 1;
		passi++;
	}

	/* creazione del vettore potenze, che contiene le potenze di 2*/
	dist = (int*) calloc(passi + 1, sizeof(int));
		
	for (i = 0; i <= passi; i++) {
		dist[i] = p << i;
	}

	/* calcolo delle somme parziali e combinazione dei risultati parziali */
	for (i = 0; i < passi; i++) {
		// ... calcolo identificativo del processore
		r = menum % dist[i + 1];

		// 
		if (r < dist[i]) {
			// calcolo dell'id del processore a cui spedire la somma locale
			sendTo = menum + dist[i];
			tag = sendTo;
			MPI_Send(&sommaloc, 1, MPI_INT, sendTo, tag, MPI_COMM_WORLD);
			/* ricevo dal processore a cui ho inviato */
			MPI_Recv(&tmp, 1, MPI_INT, sendTo, tag, MPI_COMM_WORLD, &info);
			/*calcolo della somma parziale al passo i*/
			sommaloc = sommaloc + tmp;
		}
		else { //
			recvBy = menum - dist[i];
			tag = menum;
			MPI_Recv(&tmp,1,MPI_INT,recvBy,tag,MPI_COMM_WORLD,&info);
			/* Invio al processore da cui ho ricevuto */
			MPI_Send(&sommaloc,1,MPI_INT,recvBy,tag,MPI_COMM_WORLD);
			/*calcolo della somma parziale al passo i*/
			sommaloc = sommaloc + tmp;
		} //end
	} // end for

	MPI_Barrier(MPI_COMM_WORLD); // sincronizzazione
	T_fine=MPI_Wtime()-T_inizio; // calcolo del tempo di fine

	/* calcolo del tempo totale di esecuzione*/
	MPI_Reduce(&T_fine,&T_max,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

	if (menum == 0) {
		printf("[P%d] La somma e': %5.2f\n", menum, sommaloc);
	}

	/*stampa a video dei risultati finali*/
	if (menum==0) {
		// printf("\nProcessori impegnati: %d\n", nproc);
		printf("Tempo calcolo locale: %f\n", T_fine);
		printf("MPI_Reduce max time: %f\n",T_max);
	} // end if

	/*routine chiusura ambiente MPI*/
	MPI_Finalize();
} // fine programma
