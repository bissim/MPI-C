/*
 * somma.c
 *
 * Created on: 26 set 2018
 * Author: Simone Bisogno (05225/560)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

#define NANOSECONDS_PER_SECOND 1E9
#define SEED 1

int main(int argc, char *argv[]) {
	int N = 0, sum, lower = 1, upper = 50;
	int *nums;
	float elapsedtime;
	struct timespec start, stop;

//	printf("%d arguments provided.\n", argc);
//    for (int i = 0; i < argc; ++i) {
//        printf("argv[%d]: %s\n", i, argv[i]);
//    }

//	N = atoi(argv[1]);
//	if (N == 0) {
		printf("Inserire il numero di elementi da addizionare: ");
		fflush(stdout);
		scanf("%d", &N);
//	}
	printf("Saranno generati %d numeri da %1d a %2d.\n", N, lower, upper);

	srand((unsigned int) time(0));

	nums = (int*) calloc(N, sizeof(int));

	// initialize array
	for (int i = 0; i < N; i++) {
		nums[i] = lower + rand() / (RAND_MAX / (upper - lower + 1) + 1);
	}

	//	printf("{");
	clock_gettime(CLOCK_REALTIME, &start);
	for (int i = 0; i < N; i++) {
//		if (i == N - 1) {
//			printf("%d", randnum);
//		}
//		else {
//			printf("%d, ", randnum);
//		}
		sum += nums[i];
	}
//	printf("}\n");
	clock_gettime(CLOCK_REALTIME, &stop);

	printf("\nLa somma dei numeri generati vale %d.\n", sum);
	elapsedtime = (stop.tv_sec - start.tv_sec) + (stop.tv_nsec - start.tv_nsec) / NANOSECONDS_PER_SECOND;
	printf("Tempo impiegato: %f secondi.\n", elapsedtime);

	return 0;
}
