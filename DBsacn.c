#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define ull unsigned long long

ull Rand(ull N){
	return (ull)((double)rand() / (double)RAND_MAX * (double)N);
}


void gama_remove(ull *gama, ull *gama_length, ull index){
	(*gama_length)--;
	for(ull i = index; i < *gama_length; i++){
		gama[i] = gama[i+1];
	}
}

ull gama_index(ull *gama, ull target, ull low, ull high){
	ull mid = (low + high) / 2;
	if (target == gama[mid])
		return mid;
	if (target > gama[mid])
		return gama_index(gama, target, mid+1, high);
	return gama_index(gama, target, low, mid-1);
}

void findNeighbor(ull *NeighborPts, ull *NeighborPts_length, ull j, double *X, double eps, ull N){

	double tmp;


	for(ull p = 0; p < N; p++){
		tmp = sqrt(pow((X[p*2+0] - X[j*2+0]), 2) + pow((X[p*2+1] - X[j*2+1]), 2));
		if(tmp <= eps){
			NeighborPts[ (*NeighborPts_length)++ ] = p;
		}
	}
}

void dbscan(ull* cluster, double *X, double eps, int min_pts, ull N){
	ull k = 0, i, a, b;

	ull *NeighborPts, NeighborPts_length; 
	ull *Ner_NeighborPts, Ner_NeighborPts_length; 
	ull *gama, gama_length = N; 
	bool *fil;

	NeighborPts = (ull*)malloc(N*sizeof(ull)); 
	Ner_NeighborPts = (ull*)malloc(N*sizeof(ull));
	gama = (ull*)malloc(N*sizeof(ull)); 
	fil = (bool*)malloc(N*sizeof(bool));

	memset(fil, 0, N*1*sizeof(bool));
	for(i = 0; i < N; i++){
		gama[i] = i;
	}

	while(gama_length > 0){
		ull index = Rand(gama_length);
		ull j = gama[ index ];
		gama_remove(gama, &gama_length, index);
		fil[j] = 1;
		NeighborPts_length = 0;
		findNeighbor(NeighborPts, &NeighborPts_length, j, X, eps, N);
		if(NeighborPts_length < min_pts){
			cluster[j] = 0;
		}
		else{
			k += 1;
			cluster[j] = k;
			for(i = 0; i < NeighborPts_length; i++){
				if(fil[ NeighborPts[i] ] == 0){
					gama_remove(gama, &gama_length, gama_index(gama, NeighborPts[i], 0, gama_length-1));
					fil[ NeighborPts[i] ] = 1;
					Ner_NeighborPts_length = 0;
					findNeighbor(Ner_NeighborPts, &Ner_NeighborPts_length, NeighborPts[i], X, eps, N);
					if(Ner_NeighborPts_length >= min_pts){						
						bool *exist = (bool*)malloc(N*sizeof(bool));
						memset(exist, 0, N*1*sizeof(bool));

						for(a = 0; a < NeighborPts_length; a++){
							exist[ NeighborPts[a] ] = 1;
						}

						for(a = 0; a < Ner_NeighborPts_length; a++){
							if(exist[ Ner_NeighborPts[a] ] == 0){
								NeighborPts[NeighborPts_length++] = Ner_NeighborPts[a];
							}
						}
					}
					if(cluster[NeighborPts[i]] == 0){
						cluster[NeighborPts[i]] = k;
					}
				}
			}
		}
	}
	free(NeighborPts); free(Ner_NeighborPts); free(gama); free(fil);
}

int main(int argc, char **argv){
	srand(time(NULL));
    
	double eps = argc > 2 ? atof(argv[1]) : 0.15;
	int min_pts = argc > 2 ? atoi(argv[2]) :  10;

	printf("%lf %d\n", eps, min_pts);

	clock_t start, end;
	double *X;
	ull *C, N;
	int i;

	start = clock();

	
	char *dataname_in = argc > 3 ? argv[3] : "DBSCAN_data.txt";
	FILE *fp_in = fopen(dataname_in, "r");
	fscanf(fp_in, "%llu", &N);
	X = (double*)malloc(N*2*sizeof(double));
	double tmp1,tmp2;
	for(i=0; i<N; i++){
		fscanf(fp_in,  "%lf", &tmp1); X[i*2+0] = tmp1;
		fscanf(fp_in,  "%lf", &tmp2); X[i*2+1] = tmp2;
	}
	fclose(fp_in);
	printf("\nread end\n\n");

	C = (ull*)malloc(N*1*sizeof(ull));
	memset(C, 0, N*1*sizeof(ull));

	dbscan(C, X, eps, min_pts, N);

	end = clock();
    

	printf("write start\n");
	char *dataname_out = argc > 4 ? argv[4] : "DBSCAN_answer.txt";
	FILE *fp_out = fopen(dataname_out, "w");
	for(i=0; i<N; i++)
		fprintf(fp_out, "%lf %lf %llu\n", X[i*2+0], X[i*2+1], C[i]);
	fclose(fp_out);
    
	printf("\nDBSCAN: %lf ms\n", (double)(end-start)/CLOCKS_PER_SEC);
	free(X); free(C);

	return 0;
}
