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

// 模擬 python 的 remove
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

void findNeighbor(ull *NeighborPts, ull *NeighborPts_num, ull j, double eps, ull N, bool *arr){
    ull p;
    for(p=0; p<N; p++){
        if(arr[j*N+p] == 1){
            NeighborPts[*NeighborPts_num] = p;
            *NeighborPts_num = *NeighborPts_num + 1;
        }
    }
}

void dbscan(ull* cluster, double eps, int min_pts, ull N, bool *arr){
    ull k = 0, i, a;

    ull *NeighborPts, NeighborPts_length; // 會存入點的編號，使用 ull
	ull *Ner_NeighborPts, Ner_NeighborPts_length; // 會存入點的編號，使用 ull
	ull *gama, gama_length = N; // 會存入點的編號，使用 ull
	bool *fil; // 判斷點有沒有訪問，值只會是 0(未訪問)，1(已訪問)，所以使用 bool 縮減空間
	bool *exist = (bool*)malloc(N*sizeof(bool));

	NeighborPts = (ull*)malloc(N*sizeof(ull)); // array, 某點領域內的對象
	Ner_NeighborPts = (ull*)malloc(N*sizeof(ull));
	gama = (ull*)malloc(N*sizeof(ull)); // 初始時將所有點標記爲未訪問
	fil = (bool*)malloc(N*sizeof(bool)); // 初始時已訪問對象列表爲空

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

        findNeighbor(NeighborPts, &NeighborPts_length, j, eps, N, arr);

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
                    findNeighbor(Ner_NeighborPts, &Ner_NeighborPts_length, NeighborPts[i], eps, N, arr);
                    if(Ner_NeighborPts_length >= min_pts){						
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
    free(NeighborPts); free(Ner_NeighborPts); free(gama); free(fil); free(exist);
}

__global__ void Euclidean_distance(bool *d_arr, double *d_X, double eps, ull N){
    ull ix = threadIdx.x + blockIdx.x * blockDim.x;
    ull iy = threadIdx.y + blockIdx.y * blockDim.y;

    //N*N 的二維陣列，iy -> j：檢測點，ix -> p：計算與檢測點距離的點
    //陣列內容存放該點(p) 與檢測點(j) 的距離是否小於 eps
    if(ix<N && iy<N){
                //sqrt(pow((X[p*2+0] - X[j*2+0]), 2) + pow((X[p*2+1] - X[j*2+1]), 2))
        if(eps >= sqrt(pow((d_X[ix*2+0] - d_X[iy*2+0]), 2) + pow((d_X[ix*2+1] - d_X[iy*2+1]), 2)) )
            d_arr[iy*N+ix] = 1;
        else
            d_arr[iy*N+ix] = 0;
        
    }
}

int main(int argc, char **argv){
    srand(time(NULL));
    
    double eps = argc > 2 ? atof(argv[1]) : 0.15;
    int min_pts = argc > 2 ? atoi(argv[2]) :  10;

    printf("%lf %d\n", eps, min_pts);

    clock_t start, end;
    double *X, *d_X;
    ull *C, N;
    int i;
    bool *arr, *d_arr;

    start = clock();

    //讀檔
    const char *dataname_in = argc > 3 ? argv[3] : "DBSCAN_data.txt";
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

    arr = (bool*)malloc(N*N*sizeof(bool));

    cudaMalloc((void**) &d_arr, N*N*sizeof(bool));

    cudaMalloc((void**) &d_X, N*2*sizeof(double));
    cudaMemcpy(d_X, X, N*2*sizeof(double), cudaMemcpyHostToDevice);

    dim3 block(32,32);
    dim3 grid((N+block.x-1)/block.x, (N+block.y-1)/block.y);

    Euclidean_distance<<<grid, block>>>(d_arr, d_X, eps, N);

    cudaMemcpy(arr, d_arr, N*N*sizeof(bool), cudaMemcpyDeviceToHost);
    printf("ED end\n\n");

    dbscan(C, eps, min_pts, N, arr);

    end = clock();
    
    /*
    for(i=0; i<N; i++)
        printf("%llu ", C[i]);
    */

    //寫檔
    printf("write end\n");
    const char *dataname_out = argc > 4 ? argv[4] : "DBSCAN_answer.txt";
    FILE *fp_out = fopen(dataname_out, "w");
    for(i=0; i<N; i++)
        fprintf(fp_out, "%lf %lf %llu\n", X[i*2+0], X[i*2+1], C[i]);
    fclose(fp_out);
    
    printf("\nDBSCAN: %lf ms\n", (double)(end-start)/CLOCKS_PER_SEC);
    free(X); free(C); free(arr);
    cudaFree(d_X); cudaFree(d_arr);

    return 0;
}
