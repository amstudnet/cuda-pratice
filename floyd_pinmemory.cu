#include<stdio.h>
#include<string.h>
#include <cuda_runtime.h>
#include <sys/time.h>
#define n 10005
#define inf 10000
#define BLOCK_SIZE 16

#define min(a,b) ((a)<(b))?a:b
#define max(a,b) ((a)>(b))?a:b
__global__ void floyd_warshell_gpu(int *d_e,int k){
 int i=blockIdx.x * blockDim.x + threadIdx.x;
 int j=blockIdx.y * blockDim.y + threadIdx.y;
 // 1d 1d
 
 if(i<n && j<n){
 	d_e[i*n+j]=min(d_e[i*n+k]+d_e[k*n+j],d_e[i*n+j]);
  // d_e[i][j]=min(d_e[i][k]+d_e[k][j],d_e[i][j])
 }

}
double cpuSecond() {
	struct timeval tp;
	gettimeofday(&tp,NULL);
	return ((double)tp.tv_sec + (double)tp.tv_usec*1.e-6);
}
int main(int argc,char ** argv){
  FILE *fp;
  double istart=0;
  double ielaps=0;
  int *h_e;
  //cudaMallocHost((void**)&h_aPinned, nbytes);
  cudaMallocHost((void**)&h_e, (n+1)*(n+1)*sizeof(int));
   int *d_e;
  
  cudaMalloc((void**)&d_e,(n+1)*(n+1)*sizeof(int));
  //h_e = (int*)malloc((n+1)*(n+1)*sizeof(int));
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
        if(i==j)h_e[i*n+j]=0;
        else h_e[i*n+j]=10000;
    }
  }

  fp = fopen("floyd_warshell_testdata-10000.txt", "r");
for (int i=0; i<10000; i++)
 {
       int tmp1,tmp2,tmp3;
        fscanf(fp, "%d", &tmp1);
        fscanf(fp, "%d", &tmp2);
        fscanf(fp, "%d", &tmp3);
        h_e[tmp1*n+tmp2] = tmp3;
        h_e[tmp2*n+tmp1] = tmp3;
 }
 
 
  /*h_e[0*n+0]=0;
  h_e[0*n+1]=5;
  h_e[0*n+2]=10000;
  h_e[0*n+3]=10;

  h_e[1*n+0]=10000;
  h_e[1*n+1]=0;
  h_e[1*n+2]=3;
  h_e[1*n+3]=10000;

  h_e[2*n+0]=10000;
  h_e[2*n+1]=10000;
  h_e[2*n+2]=0;
  h_e[2*n+3]=1;

  h_e[3*n+0]=10000;
  h_e[3*n+1]=10000;
  h_e[3*n+2]=10000;
  h_e[3*n+3]=0;*/
  dim3 dimblock(BLOCK_SIZE,BLOCK_SIZE);
  dim3 dimGrid((n+2+BLOCK_SIZE-1)/BLOCK_SIZE, (n+2+BLOCK_SIZE-1)/BLOCK_SIZE);
  istart=cpuSecond(); 
  for(int k=0;k<n;k++){
     cudaMemcpy(d_e,h_e,(n+1)*(n+1)*sizeof(int),cudaMemcpyHostToDevice);
  	floyd_warshell_gpu<<<dimGrid,dimblock>>>(d_e,k);
  	cudaMemcpy(h_e, d_e, (n+1)*(n+1)*sizeof(int), cudaMemcpyDeviceToHost);
  	
  }

  /*for(int i=0;i<n;i++){

  	for(int j=0;j<n;j++){
  		if(h_e[i*n+j]==10000){
  			 printf("%7s", "INF");
  		}
  		else{
  			printf("%7d",h_e[i*n+j]);
  		}
  		
  	}
  	printf("\n");
  }*/
    ielaps=cpuSecond()-istart;
  printf("10000pinned_memory_floyd_warshell elapsed %7.2f ms\n",ielaps*1000);
  printf("%d\n",h_e[20*n+21] );
  //free(h_e);
  cudaFreeHost(h_e);
  cudaFree(d_e);
  return 0;
}
