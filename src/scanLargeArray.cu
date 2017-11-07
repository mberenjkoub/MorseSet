/*
 * Copyright 1993-2006 NVIDIA Corporation.  All rights reserved.
 *
 * NOTICE TO USER:   
 *
 * This source code is subject to NVIDIA ownership rights under U.S. and 
 * international Copyright laws.  
 *
 * NVIDIA MAKES NO REPRESENTATION ABOUT THE SUITABILITY OF THIS SOURCE 
 * CODE FOR ANY PURPOSE.  IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR 
 * IMPLIED WARRANTY OF ANY KIND.  NVIDIA DISCLAIMS ALL WARRANTIES WITH 
 * REGARD TO THIS SOURCE CODE, INCLUDING ALL IMPLIED WARRANTIES OF 
 * MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS FOR A PARTICULAR PURPOSE.   
 * IN NO EVENT SHALL NVIDIA BE LIABLE FOR ANY SPECIAL, INDIRECT, INCIDENTAL, 
 * OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS 
 * OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE 
 * OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE 
 * OR PERFORMANCE OF THIS SOURCE CODE.  
 *
 * U.S. Government End Users.  This source code is a "commercial item" as 
 * that term is defined at 48 C.F.R. 2.101 (OCT 1995), consisting  of 
 * "commercial computer software" and "commercial computer software 
 * documentation" as such terms are used in 48 C.F.R. 12.212 (SEPT 1995) 
 * and is provided to the U.S. Government only as a commercial end item.  
 * Consistent with 48 C.F.R.12.212 and 48 C.F.R. 227.7202-1 through 
 * 227.7202-4 (JUNE 1995), all U.S. Government End Users acquire the 
 * source code with only those rights set forth herein.
 */

#ifdef _WIN32
#  define NOMINMAX 
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime.h>

#include <helper_cuda.h>
#include <helper_functions.h>
#include <helper_timer.h>
#include "stdint.h"
// includes, kernels
//#include <scan.cu>  // defines prescanArray()

////////////////////////////////////////////////////////////////////////////////


// 16 banks on G80
#define NUM_BANKS 16
#define LOG_NUM_BANKS 4

#ifdef ZERO_BANK_CONFLICTS
#define CONFLICT_FREE_OFFSET(index) ((index) >> LOG_NUM_BANKS + (index) >> (2*LOG_NUM_BANKS))
#else
#define CONFLICT_FREE_OFFSET(index) ((index) >> LOG_NUM_BANKS)
#endif

///////////////////////////////////////////////////////////////////////////////
// Work-efficient compute implementation of scan, one thread per 2 elements
// Work-efficient: O(log(n)) steps, and O(n) adds.
// Also shared storage efficient: Uses n + n/NUM_BANKS shared memory -- no ping-ponging
// Also avoids most bank conflicts using single-element offsets every NUM_BANKS elements.
//
// In addition, If ZERO_BANK_CONFLICTS is defined, uses 
//     n + n/NUM_BANKS + n/(NUM_BANKS*NUM_BANKS) 
// shared memory. If ZERO_BANK_CONFLICTS is defined, avoids ALL bank conflicts using 
// single-element offsets every NUM_BANKS elements, plus additional single-element offsets 
// after every NUM_BANKS^2 elements.
//
// Uses a balanced tree type algorithm.  See Blelloch, 1990 "Prefix Sums 
// and Their Applications", or Prins and Chatterjee PRAM course notes:
// http://www.cs.unc.edu/~prins/Classes/203/Handouts/pram.pdf
// 
// This work-efficient version is based on the algorithm presented in Guy Blelloch's
// excellent paper "Prefix sums and their applications".
// http://www-2.cs.cmu.edu/afs/cs.cmu.edu/project/scandal/public/papers/CMU-CS-90-190.html
//
// Pro: Work Efficient, very few bank conflicts (or zero if ZERO_BANK_CONFLICTS is defined)
// Con: More instructions to compute bank-conflict-free shared memory addressing,
// and slightly more shared memory storage used.
//

template <bool isNP2>
__device__ void loadSharedChunkFromMem(uint32_t *s_data,
	const uint32_t *g_idata,
	int n, int baseIndex,
	int& ai, int& bi,
	int& mem_ai, int& mem_bi,
	int& bankOffsetA, int& bankOffsetB)
{
	int thid = threadIdx.x;
	mem_ai = baseIndex + threadIdx.x;
	mem_bi = mem_ai + blockDim.x;

	ai = thid;
	bi = thid + blockDim.x;

	// compute spacing to avoid bank conflicts
	bankOffsetA = CONFLICT_FREE_OFFSET(ai);
	bankOffsetB = CONFLICT_FREE_OFFSET(bi);

	// Cache the computational window in shared memory
	// pad values beyond n with zeros
	s_data[ai + bankOffsetA] = g_idata[mem_ai];

	if (isNP2) // compile-time decision
	{
		s_data[bi + bankOffsetB] = (bi < n) ? g_idata[mem_bi] : 0;
	}
	else
	{
		s_data[bi + bankOffsetB] = g_idata[mem_bi];
	}
}

template <bool isNP2>
__device__ void storeSharedChunkToMem(uint32_t* g_odata,
	const uint32_t* s_data,
	int n,
	int ai, int bi,
	int mem_ai, int mem_bi,
	int bankOffsetA, int bankOffsetB)
{
	__syncthreads();

	// write results to global memory
	g_odata[mem_ai] = s_data[ai + bankOffsetA];
	if (isNP2) // compile-time decision
	{
		if (bi < n)
			g_odata[mem_bi] = s_data[bi + bankOffsetB];
	}
	else
	{
		g_odata[mem_bi] = s_data[bi + bankOffsetB];
	}
}

template <bool storeSum>
__device__ void clearLastElement(uint32_t* s_data,
	uint32_t *g_blockSums,
	int blockIndex)
{
	if (threadIdx.x == 0)
	{
		int index = (blockDim.x << 1) - 1;
		index += CONFLICT_FREE_OFFSET(index);

		if (storeSum) // compile-time decision
		{
			// write this block's total sum to the corresponding index in the blockSums array
			g_blockSums[blockIndex] = s_data[index];
		}

		// zero the last element in the scan so it will propagate back to the front
		s_data[index] = 0;
	}
}



__device__ unsigned int buildSum(uint32_t *s_data)
{
	unsigned int thid = threadIdx.x;
	unsigned int stride = 1;

	// build the sum in place up the tree
	for (int d = blockDim.x; d > 0; d >>= 1)
	{
		__syncthreads();

		if (thid < d)
		{
			int i = __mul24(__mul24(2, stride), thid);
			int ai = i + stride - 1;
			int bi = ai + stride;

			ai += CONFLICT_FREE_OFFSET(ai);
			bi += CONFLICT_FREE_OFFSET(bi);

			s_data[bi] += s_data[ai];
		}

		stride *= 2;
	}

	return stride;
}

__device__ void scanRootToLeaves(uint32_t *s_data, unsigned int stride)
{
	unsigned int thid = threadIdx.x;

	// traverse down the tree building the scan in place
	for (int d = 1; d <= blockDim.x; d *= 2)
	{
		stride >>= 1;

		__syncthreads();

		if (thid < d)
		{
			int i = __mul24(__mul24(2, stride), thid);
			int ai = i + stride - 1;
			int bi = ai + stride;

			ai += CONFLICT_FREE_OFFSET(ai);
			bi += CONFLICT_FREE_OFFSET(bi);

			uint32_t t = s_data[ai];
			s_data[ai] = s_data[bi];
			s_data[bi] += t;
		}
	}
}

template <bool storeSum>
__device__ void prescanBlock(uint32_t *data, int blockIndex, uint32_t *blockSums)
{
	int stride = buildSum(data);               // build the sum in place up the tree
	clearLastElement<storeSum>(data, blockSums,
		(blockIndex == 0) ? blockIdx.x : blockIndex);
	scanRootToLeaves(data, stride);            // traverse down tree to build the scan 
}

template <bool storeSum, bool isNP2>
__global__ void prescan(uint32_t *g_odata,
	const uint32_t *g_idata,
	uint32_t *g_blockSums,
	int n,
	int blockIndex,
	int baseIndex)
{
	int ai, bi, mem_ai, mem_bi, bankOffsetA, bankOffsetB;
	extern __shared__ uint32_t s_data[];

	// load data into shared memory
	loadSharedChunkFromMem<isNP2>(s_data, g_idata, n,
		(baseIndex == 0) ?
		__mul24(blockIdx.x, (blockDim.x << 1)) : baseIndex,
		ai, bi, mem_ai, mem_bi,
		bankOffsetA, bankOffsetB);
	// scan the data in each block
	prescanBlock<storeSum>(s_data, blockIndex, g_blockSums);
	// write results to device memory
	storeSharedChunkToMem<isNP2>(g_odata, s_data, n,
		ai, bi, mem_ai, mem_bi,
		bankOffsetA, bankOffsetB);
}


__global__ void uniformAdd(uint32_t *g_data,
	uint32_t *uniforms,
	int n,
	int blockOffset,
	int baseIndex)
{
	__shared__ uint32_t uni;
	if (threadIdx.x == 0)
		uni = uniforms[blockIdx.x + blockOffset];

	unsigned int address = __mul24(blockIdx.x, (blockDim.x << 1)) + baseIndex + threadIdx.x;

	__syncthreads();

	// note two adds per thread
	g_data[address] += uni;
	g_data[address + blockDim.x] += (threadIdx.x + blockDim.x < n) * uni;
}



//======================================================================================
// includes, kernels



inline bool
isPowerOfTwo(int n)
{
	return ((n&(n - 1)) == 0);
}

inline int
floorPow2(int n)
{
#ifdef WIN32
	// method 2
	return 1 << (int)logb((float)n);
#else
	// method 1
	// float nf = (float)n;
	// return 1 << (((*(int*)&nf) >> 23) - 127); 
	int exp;
	frexp((float)n, &exp);
	return 1 << (exp - 1);
#endif
}

#define BLOCK_SIZE 256

uint32_t** g_scanBlockSums;
unsigned int g_numEltsAllocated = 0;
unsigned int g_numLevelsAllocated = 0;

void preallocBlockSums(unsigned int maxNumElements)
{
	assert(g_numEltsAllocated == 0); // shouldn't be called 

	g_numEltsAllocated = maxNumElements;

	unsigned int blockSize = BLOCK_SIZE; // max size of the thread blocks
	unsigned int numElts = maxNumElements;

	int level = 0;

	do
	{
		unsigned int numBlocks =
			max(1, (int)ceil((float)numElts / (2.f * blockSize)));
		if (numBlocks > 1)
		{
			level++;
		}
		numElts = numBlocks;
	} while (numElts > 1);

	g_scanBlockSums = (uint32_t**)malloc(level * sizeof(uint32_t*));
	g_numLevelsAllocated = level;

	numElts = maxNumElements;
	level = 0;

	do
	{
		unsigned int numBlocks =
			max(1, (int)ceil((float)numElts / (2.f * blockSize)));
		if (numBlocks > 1)
		{
			(cudaMalloc((void**)&g_scanBlockSums[level++],
				numBlocks * sizeof(uint32_t)));
		}
		numElts = numBlocks;
	} while (numElts > 1);

	("preallocBlockSums");
}

void deallocBlockSums()
{
	for (int i = 0; i < g_numLevelsAllocated; i++)
	{
		cudaFree(g_scanBlockSums[i]);
	}

	("deallocBlockSums");

	free((void**)g_scanBlockSums);

	g_scanBlockSums = 0;
	g_numEltsAllocated = 0;
	g_numLevelsAllocated = 0;
}


void prescanArrayRecursive(uint32_t *outArray,
	const uint32_t *inArray,
	int numElements,
	int level)
{
	unsigned int blockSize = BLOCK_SIZE; // max size of the thread blocks
	unsigned int numBlocks =
		max(1, (int)ceil((float)numElements / (2.f * blockSize)));
	unsigned int numThreads;

	if (numBlocks > 1)
		numThreads = blockSize;
	else if (isPowerOfTwo(numElements))
		numThreads = numElements / 2;
	else
		numThreads = floorPow2(numElements);

	unsigned int numEltsPerBlock = numThreads * 2;

	// if this is a non-power-of-2 array, the last block will be non-full
	// compute the smallest power of 2 able to compute its scan.
	unsigned int numEltsLastBlock =
		numElements - (numBlocks - 1) * numEltsPerBlock;
	unsigned int numThreadsLastBlock = max(1, numEltsLastBlock / 2);
	unsigned int np2LastBlock = 0;
	unsigned int sharedMemLastBlock = 0;

	if (numEltsLastBlock != numEltsPerBlock)
	{
		np2LastBlock = 1;

		if (!isPowerOfTwo(numEltsLastBlock))
			numThreadsLastBlock = floorPow2(numEltsLastBlock);

		unsigned int extraSpace = (2 * numThreadsLastBlock) / NUM_BANKS;
		sharedMemLastBlock =
			sizeof(uint32_t)* (2 * numThreadsLastBlock + extraSpace);
	}

	// padding space is used to avoid shared memory bank conflicts
	unsigned int extraSpace = numEltsPerBlock / NUM_BANKS;
	unsigned int sharedMemSize =
		sizeof(uint32_t)* (numEltsPerBlock + extraSpace);

#ifdef DEBUG
	if (numBlocks > 1)
	{
		assert(g_numEltsAllocated >= numElements);
	}
#endif

	// setup execution parameters
	// if NP2, we process the last block separately
	dim3  grid(max(1, numBlocks - np2LastBlock), 1, 1);
	dim3  threads(numThreads, 1, 1);

	// make sure there are no CUDA errors before we start
	("prescanArrayRecursive before kernels");

	// execute the scan
	if (numBlocks > 1)
	{
		prescan<true, false> << < grid, threads, sharedMemSize >> >(outArray,
			inArray,
			g_scanBlockSums[level],
			numThreads * 2, 0, 0);
		getLastCudaError("prescanWithBlockSums");
		if (np2LastBlock)
		{
			prescan<true, true> << < 1, numThreadsLastBlock, sharedMemLastBlock >> >
				(outArray, inArray, g_scanBlockSums[level], numEltsLastBlock,
				numBlocks - 1, numElements - numEltsLastBlock);
			getLastCudaError("prescanNP2WithBlockSums");
		}

		// After scanning all the sub-blocks, we are mostly done.  But now we 
		// need to take all of the last values of the sub-blocks and scan those.  
		// This will give us a new value that must be sdded to each block to 
		// get the final results.
		// recursive (CPU) call
		prescanArrayRecursive(g_scanBlockSums[level],
			g_scanBlockSums[level],
			numBlocks,
			level + 1);

		uniformAdd << < grid, threads >> >(outArray,
			g_scanBlockSums[level],
			numElements - numEltsLastBlock,
			0, 0);
		getLastCudaError("uniformAdd");
		if (np2LastBlock)
		{
			uniformAdd << < 1, numThreadsLastBlock >> >(outArray,
				g_scanBlockSums[level],
				numEltsLastBlock,
				numBlocks - 1,
				numElements - numEltsLastBlock);
			getLastCudaError("uniformAdd");
		}
	}
	else if (isPowerOfTwo(numElements))
	{
		prescan<false, false> << < grid, threads, sharedMemSize >> >(outArray, inArray,
			0, numThreads * 2, 0, 0);
		getLastCudaError("prescan");
	}
	else
	{
		prescan<false, true> << < grid, threads, sharedMemSize >> >(outArray, inArray,
			0, numElements, 0, 0);
		getLastCudaError("prescanNP2");
	}
}

void prescanArray(uint32_t *outArray, uint32_t *inArray, int numElements)
{
	prescanArrayRecursive(outArray, inArray, numElements, 0);
}

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
//int 
//main( int argc, char** argv) 
//{
//    runTest( argc, argv);
//    CUT_EXIT(argc, argv);
//}

////////////////////////////////////////////////////////////////////////////////
//! Run a scan test for CUDA
////////////////////////////////////////////////////////////////////////////////
void
runTest(uint32_t* d_In, uint32_t* d_Out,uint32_t* h_Out, int num_elements)
{
   


    int num_test_iterations = 100;
  //  int num_elements = 1000000; // can support large, non-power-of-2 arrays!


	unsigned int mem_size = sizeof(uint32_t)* num_elements;
    
    unsigned int timerGPU, timerCPU;
  /*  CUT_SAFE_CALL(cutCreateTimer(&timerCPU));
    CUT_SAFE_CALL(cutCreateTimer(&timerGPU));*/

    // allocate host memory to store the input data
    //uint32_t* h_data = (uint32_t*) malloc( mem_size);
    //  
    //// initialize the input data on the host
    //for( unsigned int i = 0; i < num_elements; ++i) 
    //{
    //    h_data[i] = 1.0f;//(int)(10 * rand()/32768.f);
    //}

    // compute reference solution
	uint32_t* reference = (uint32_t*)malloc(mem_size);
  //  cutStartTimer(timerCPU);
  //  for (int i = 0; i < num_test_iterations; i++)
  //  {
  ////      computeGold( reference, h_data, num_elements);
  //  }
  //  cutStopTimer(timerCPU);

    // allocate device memory input and output arrays
    uint32_t* d_idata = NULL;
	uint32_t* d_odata = NULL;

   /* ( cudaMalloc( (void**) &d_idata, mem_size));
    ( cudaMalloc( (void**) &d_odata, mem_size));*/
    
    // copy host memory to device input array
  //  ( cudaMemcpy( d_idata, h_data, mem_size, cudaMemcpyHostToDevice) );
    // initialize all the other device arrays to be safe
  //  ( cudaMemcpy( d_odata, h_data, mem_size, cudaMemcpyHostToDevice) );

    //printf("Running parallel prefix sum (prescan) of %d elements\n", num_elements);
    //printf("This version is work efficient (O(n) adds)\n");
    //printf("and has very few shared memory bank conflicts\n\n");

    preallocBlockSums(num_elements);

    // run once to remove startup overhead
    prescanArray(d_Out, d_In, num_elements);


    // Run the prescan
//    cutStartTimer(timerGPU);
    for (int i = 0; i < num_test_iterations; i++)
    {
        //printf("prescanArray\n");
        prescanArray(d_Out, d_In, num_elements);
    }
   // cutStopTimer(timerGPU);

    deallocBlockSums();    

    // copy result from device to host
	(cudaMemcpy(h_Out, d_Out, sizeof(uint32_t)* num_elements,
                               cudaMemcpyDeviceToHost));

    //// If this is a regression test write the results to a file
    //if( cutCheckCmdLineFlag( argc, (const char**) argv, "regression")) 
    //{
    //    // write file for regression test 
    //    cutWriteFilef( "./data/result.dat", h_data, num_elements, 0.0);
    //}
    //else 
    //{
    //    // custom output handling when no regression test running
    //    // in this case check if the result is equivalent to the expected soluion
    //    unsigned int result_regtest = cutComparef( reference, h_data, num_elements);
    //    printf( "Test %s\n", (1 == result_regtest) ? "PASSED" : "FAILED");
    //    printf( "Average GPU execution time: %f ms\n", cutGetTimerValue(timerGPU) / num_test_iterations);
    //    printf( "CPU execution time:         %f ms\n", cutGetTimerValue(timerCPU) / num_test_iterations);
    //}

    // cleanup memory
 //   cutDeleteTimer(timerCPU);
 //   cutDeleteTimer(timerGPU);
    //free( h_data);
    //free( reference);
    //cudaFree( d_odata);
    //cudaFree( d_idata);
}
