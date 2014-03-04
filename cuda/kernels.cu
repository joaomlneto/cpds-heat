#include <math.h>
#include <float.h>
#include <cuda.h>

// First solution with global memory
__global__ void gpu_Heat (float *u, float *utmp, float *residual,int N) {

	// TODO: kernel computation
	int sizey = N;
	int j = blockIdx.x * blockDim.x + threadIdx.x;
	int i = blockIdx.y * blockDim.y + threadIdx.y;
	float diff=0.0;
 	if( i < N-1 && j < N-1 && i > 0 && j > 0) {
  		utmp[i*sizey+j]= 0.25 *
								 (u[ i*sizey     + (j-1) ]+  // left
					           u[ i*sizey     + (j+1) ]+  // right
				              u[ (i-1)*sizey + j     ]+  // top
				              u[ (i+1)*sizey + j     ]); // bottom
      diff = utmp[i*sizey+j] - u[i*sizey + j];
      residual[i*sizey+j] = diff * diff;
	}
}

// Shared memory residual calculation
// Reduction code from CUDA Slides - Mark Harris

__global__ void gpu_HeatReduction (float *res, float *result) {

	extern __shared__ float sdata[];
	unsigned int tid = threadIdx.x;
	unsigned int index= blockIdx.x*blockDim.x+ threadIdx.x;

	sdata[tid] = res[index];
	__syncthreads();

	
	// Reduce the shared table to compute the residual

	for(unsigned int s=blockDim.x/2; s>0; s>>=1) {
		if (tid < s) {
		sdata[tid] += sdata[tid + s];
		}
		__syncthreads();
	}
	if (tid == 0)
	{
		int blockIndex = blockIdx.x;

		result[blockIndex] = sdata[tid];



	}

}
