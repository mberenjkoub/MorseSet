

#include <cstdio>
#include <cuda_runtime.h>

#include <helper_cuda.h>
#include <helper_functions.h>
#include <helper_timer.h>
#include<stack>
#include<set>
#include<map>
#include<queue>
#include<math.h>

#include "SCC.h"
#include "scc_kernel.h"
//#include "graph_generator.h"
//#include "parallel_fwd.h"
//#include "hash_table.h"

#ifdef _DEBUG
void bbin_printf(uint32_t elem, int N = 32, int end = 0)
{
	for ( int i = N - 1; i >= end; i-- )
		printf("%d", (bool)(elem & ((uint32_t)1 << i)));
}
#endif

bool _DeviceSet;
//
__global__ void checkTerminateAndSetOldMap(COL_vertex * m, const uint32_t num_rows, uint32_t * propagate)
{
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x + 1;
	if (row >= num_rows)
		return;
	COL_vertex vertex = m[row];
	if (vertex.isPropagate())
		return;

	*propagate = 1;
	vertex.setOldMap(vertex.getMap());
	vertex.setMap(0);

	m[row] = vertex;
}
//
//
__global__ void COL_FWD(const Edge * Fc, const uint32_t * Fr, COL_vertex * m, const uint32_t num_rows, uint32_t * propagate)
{
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x + 1; // + 1 we dont work with undefined vertex
	if (row >= num_rows)
		return;
	COL_vertex vertex = m[row];
	if (vertex.isPropagate() || !vertex.isBackwardVisited())
		return;

	uint32_t row_begin = Fr[row];
	uint32_t row_end = Fr[row + 1];

	int prop = 0;
	for (uint32_t column = row_begin; column < row_end; column++) {
		uint32_t index = Fc[column].getValue();
		COL_vertex p_vertex = m[index];

		if (p_vertex.isBackwardVisited() || p_vertex.getMap() != vertex.getMap())
			continue;
		p_vertex.setBackwardVisitedBit();
		m[index] = p_vertex;
		prop = 1;
	}
	if (prop)
		*propagate = prop;
	vertex.setPropagateBit();
	m[row] = vertex;
}
//
//
__global__ void OBFcomputeSCCs(const OBF_vertex * m, uint32_t * Fr, const uint32_t num_rows)
{
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
	if (row >= num_rows || row == 0)
		return;

	OBF_vertex vertex = m[row];

	Fr[vertex.getOldRange() - 1] = 1;
}

 void OBFcomputeSCCs_c(const OBF_vertex * m,ASF_vertex*a,Dimension*d,  uint32_t * Fr, const uint32_t num_rows)
{
	 uint32_t oldrange[100];
	 uint32_t oldrangeCounter[100];

	//uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
	 int cc = 1;
	 uint32_t range = 0;
	 for (uint32_t row = 0; row < num_rows; row++)
	 {
		 if (row >= num_rows || row == 0)
			 continue;
		 int  i = 1;
		 OBF_vertex vertex = m[row];


		 if (vertex.getOldRange() != row)
		 {
			 for ( i = 1; i < cc; i++)
			 {

				 if (vertex.getOldRange() == oldrange[i])
				 {
					 a[row].setInSCC();
					 range = a[row].getOldRange();
			/*		 a[row].Fr[0] = row + 1;
					 a[row].Fr[1] = row + d->x + 1;
					 a[row].Fr[2] = row + d->x;*/

					 oldrangeCounter[i]++;


					 break;
					 printf(" %d,%d ", row, vertex.getOldRange());

				 }
			 }
			 if (i == cc)
			 {
				 oldrange[cc] = vertex.getOldRange();
				 cc++;
			 }




		 }
		 if(vertex.getOldRange() > 0)
		 Fr[vertex.getOldRange() - 1] = 1;
	 }


}
//
//
__global__ void COLcomputeSCCs(const COL_vertex * m, uint32_t * Fr, const uint32_t num_rows, const ordering o)
{
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x + 1;
	if (row >= num_rows)
		return;

	COL_vertex vertex = m[row];
	uint32_t temp = MAX_VERTEX;
	temp -= vertex.getMap() + 1;
	if (temp < num_rows)
		Fr[temp] = 1;

	//	Fr[ ( o == MAX ) ? (vertex.getMap() - 1) : ((MAX_VERTEX - vertex.getMap()) - 1) ] = 1;
	//	Fr[ 0 ] = 1;
	//	Fr[ (MAX_VERTEX - vertex.getMap()) - 1 ] = 1;
}
//
//
__global__ void setBWDseed(COL_vertex * m, const uint32_t num_rows, const ordering o)
{
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x + 1;
	if (row >= num_rows)
		return;
	COL_vertex vertex = m[row];
	if (vertex.isPropagate())
		return;

	if (o == MAX) {
		if (row == vertex.getMap())
			vertex.setBackwardVisitedBit();
	}
	else {
		if ((MAX_VERTEX - row) == vertex.getMap())
			vertex.setBackwardVisitedBit();
	}
	m[row] = vertex;
}
//
//
__global__ void p_COL_MAP(const Edge * Fc, const uint32_t * Fr, COL_vertex * m, const uint32_t num_rows, const ordering o,
	uint32_t * propagate)
{
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x + 1; // + 1 we dont work with undefined vertex
	if (row >= num_rows)
		return;
	COL_vertex vertex = m[row];
	if (vertex.isPropagate())
		return;

	uint32_t row_begin = Fr[row];
	uint32_t row_end = Fr[row + 1];
	uint32_t max_candidate;

	if (o == MAX)
		max_candidate = row;
	else
		max_candidate = MAX_VERTEX - row;

	for (uint32_t column = row_begin; column < row_end; column++) {
		uint32_t index = Fc[column].getValue();
		COL_vertex p_vertex = m[index];
		if (!p_vertex.isPropagate() && p_vertex.getOldMap() == vertex.getOldMap())
			max_candidate = max(max_candidate, max(vertex.getMap(), p_vertex.getMap()));
	}

	if (max_candidate > vertex.getMap()) {
		vertex.setMap(max_candidate);
		*propagate = 1;
	}
	m[row] = vertex;
}
//
__global__ void OBFtoCOLKernel(const OBF_vertex * m, COL_vertex * cm, const uint32_t num_rows)
{
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
	if (row >= num_rows)
		return;

	cm[row] = m[row];
}
//
__global__ void SetCOL(OBF_vertex * m, const uint32_t * OldRange, const uint32_t num_rows)
{
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
	if (row >= num_rows)
		return;

	OBF_vertex vertex = m[row];
	if (vertex.getOldRange() == *OldRange && !vertex.isInSCC()) {
		vertex.setInCOL();
		m[row] = vertex;
	}
}

//
//
// void OBFCOLcompute_c(const OBF_vertex * m, uint32_t * out_field, const uint32_t * OldRange, const uint32_t num_rows)
//{
//	 for (uint32_t row = 0; row <= num_rows; row++)
//	 {
//		 //uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
//		 if (row >= num_rows)
//			 return;
//
//		 OBF_vertex vertex = m[row];
//		 if (vertex.getOldRange() == *OldRange && !vertex.isInSCC())
//			 out_field[row - 1] = 1;
//	 }
//}
//
//
__global__ void OBFCOLcompute(const OBF_vertex * m, uint32_t * out_field, const uint32_t * OldRange, const uint32_t num_rows)
{
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
	if (row >= num_rows)
		return;
	if (row > 100)
		return;
	OBF_vertex vertex = m[row];
	//printf("old range = %d, %d", vertex.getOldRange(), OldRange);

	if (vertex.getOldRange() == *OldRange && !vertex.isInSCC())
		out_field[row - 1] = 1;
}
//
//
__global__ void FWDOWCTY_pivot2(OBF_vertex * m, const uint32_t * FWDOldRange, const uint32_t OWCTYOldRange, const uint32_t * FWDpivot,
	const uint32_t * OWCTYpivot, const uint32_t num_rows)
{
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
	if (row >= num_rows)
		return;
	OBF_vertex vertex = m[row];
	if (vertex.isInFWD() && vertex.getOldRange() == *FWDOldRange && *FWDpivot == row) {
		vertex.setFWDVisited();
		vertex.setRange(row);
		vertex.setDone2();
		m[row] = vertex;
	}
	else {
		if (vertex.isInOWCTY() && vertex.getOldRange() == OWCTYOldRange && vertex.isReached()) {
			if (row == *OWCTYpivot)
				vertex.setDone2();
			vertex.setRange(*OWCTYpivot);
			m[row] = vertex;
		}
	}
}
//
//
// void cuReduce1_c(uint32_t *g_idata, uint32_t * g_odata)
//{
//	 uint32_t sdata[48];
//	volatile uint32_t * stest = sdata;
//	uint32_t tid = threadIdx.x;
//	stest[tid] = g_idata[tid] + g_idata[tid + 32];
//	//	stest[ tid ] = g_idata[ tid ];
//	// 	sdata[ tid ] += sdata[ tid + 16 ];
//	// 	sdata[ tid ] += sdata[ tid + 8 ];
//	// 	sdata[ tid ] += sdata[ tid + 4 ];
//	// 	sdata[ tid ] += sdata[ tid + 2 ];
//	// 	sdata[ tid ] += sdata[ tid + 1 ];
//
//	// 	if ( tid == 0 )
//	// 		*g_odata = sdata[ 0 ];
//
//	stest[tid] += stest[tid + 16];
//	stest[tid] += stest[tid + 8];
//	stest[tid] += stest[tid + 4];
//	stest[tid] += stest[tid + 2];
//	stest[tid] += stest[tid + 1];
//
//	if (tid == 0)
//		*g_odata = sdata[0];
//}
//
//
__global__ void cuReduce1(uint32_t *g_idata, uint32_t * g_odata)
{
	__shared__ uint32_t sdata[48];
	volatile uint32_t * stest = sdata;
	uint32_t tid = threadIdx.x;
	stest[tid] = g_idata[tid] + g_idata[tid + 32];
	//	stest[ tid ] = g_idata[ tid ];
	// 	sdata[ tid ] += sdata[ tid + 16 ];
	// 	sdata[ tid ] += sdata[ tid + 8 ];
	// 	sdata[ tid ] += sdata[ tid + 4 ];
	// 	sdata[ tid ] += sdata[ tid + 2 ];
	// 	sdata[ tid ] += sdata[ tid + 1 ];
//
	// 	if ( tid == 0 )
	// 		*g_odata = sdata[ 0 ];

	stest[tid] += stest[tid + 16];
	stest[tid] += stest[tid + 8];
	stest[tid] += stest[tid + 4];
	stest[tid] += stest[tid + 2];
	stest[tid] += stest[tid + 1];

	if (tid == 0)
		*g_odata = sdata[0];
}


__global__ void cuReduce(const uint32_t *g_idata, uint32_t *g_odata, const uint32_t n)
{
	__shared__ uint32_t sdata[blockSize];

	uint32_t tid = threadIdx.x;
	if (tid > 128)
	printf("threadidx = %d \t", threadIdx.x);



	uint32_t i = blockIdx.x * (blockSize * 2) + threadIdx.x;
	//printf("i = %d \n", i);
	uint32_t gridSize = blockSize * 2 * gridDim.x;
	//printf("grid size = %d \n", gridSize);

	sdata[tid] = 0;

	while (i < n) {
		sdata[tid] += g_idata[i];
		if (i + blockSize < n)
			sdata[tid] += g_idata[i + blockSize];
		else
			break;
		i += gridSize;
	}
	__syncthreads();

	if (tid < 64)
		sdata[tid] += sdata[tid + 64];
	__syncthreads();

	// 	if ( tid < 32 ) {
	// 		sdata[ tid ] += sdata[ tid + 32 ];
	// 		sdata[ tid ] += sdata[ tid + 16 ];
	// 		sdata[ tid ] += sdata[ tid + 8 ];
	// 		sdata[ tid ] += sdata[ tid + 4 ];
	// 		sdata[ tid ] += sdata[ tid + 2 ];
	// 		sdata[ tid ] += sdata[ tid + 1 ];
	// 	}
	//
	// 	if ( tid == 0 )
	// 		g_odata[ blockIdx.x ] = sdata[ 0 ];
	volatile uint32_t * stest = sdata;
	if (tid < 32) {
		stest[tid] += stest[tid + 32];
		stest[tid] += stest[tid + 16];
		stest[tid] += stest[tid + 8];
		stest[tid] += stest[tid + 4];
		stest[tid] += stest[tid + 2];
		stest[tid] += stest[tid + 1];
	}

	if (tid == 0)
		g_odata[blockIdx.x] = sdata[0];
}


//
//
__global__ void BWDSynchKernel_FWDOWCTYpivot1(OBF_vertex * m, const uint32_t OldRange, uint32_t * NewRange, uint32_t * FWDpivot,
	uint32_t * OWCTYpivot, const uint32_t num_rows)
{
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
	if (row >= num_rows)
		return;
	OBF_vertex vertex = m[row];
	if (!vertex.isInBWD() || vertex.getOldRange() != OldRange)
		return;

	if (row == vertex.getRange())
		*NewRange = row;
	if (vertex.isBWDVisited()) {
		vertex.setOldRange(vertex.getRange());
		vertex.setRange(0);
		vertex.setInFWD();
		if (!*FWDpivot)
			*FWDpivot = row;
	}
	else {
		if (vertex.isBWDPropagate()) {
			*OWCTYpivot = row;
			vertex.setInOWCTY();
			vertex.setReached();
		}
		else
			vertex.setInOWCTY();
	}
	m[row] = vertex;
}
//
//
//
__global__ void BWDTerminateKernel_OWCTY1(OBF_vertex * m, const Edge * Fc, const uint32_t * Fr, const uint32_t OldRange,
	const uint32_t num_rows)
{
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
	if (row >= num_rows)
		return;
	OBF_vertex vertex = m[row];
	if (!vertex.isInBWD() || vertex.getOldRange() != OldRange)
		return;

	if (vertex.isBWDVisited()) {
		uint32_t row_begin = Fr[row];
		uint32_t row_end = Fr[row + 1];

		for (uint32_t column = row_begin; column < row_end; column++) {
			uint32_t index = Fc[column].getValue();
			OBF_vertex p_vertex = m[index];

			if (p_vertex.getOldRange() != OldRange || !p_vertex.isInBWD() || p_vertex.isBWDVisited())
				continue;
//
			p_vertex.setBWDPropagate();
			m[index] = p_vertex;
		}

		if (m[OldRange].getRange() == vertex.getRange())
			m[row].setInSCC();
	}
}
//
//
__global__ void CheckTerminateKernel(const OBF_vertex * m, uint32_t * terminate, const uint32_t num_rows)
{
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
	if (row >= num_rows)
		return;
	OBF_vertex vertex = m[row];
	if (!vertex.isInSCC())
		*terminate = 0;
}
//
//
__global__ void BWD_pivot2(OBF_vertex * m, const uint32_t OldRange, uint32_t * pivot, const uint32_t num_rows)
{
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
	if (row >= num_rows)
		return;
	OBF_vertex vertex = m[row];
	if (!vertex.isInBWD() || vertex.getOldRange() != OldRange)
		return;
//
	if (vertex.isBWDVisited())
		vertex.setRange(*pivot);

	if (row == *pivot) {
		vertex.setDone2();
	}
	m[row] = vertex;
}
//
//
__global__ void OWCTYSynchKernel_BWDpivot1(OBF_vertex * m, const Edge * Bc, const uint32_t * Br, const uint32_t OldRange, uint32_t * pivot,
	const uint32_t num_rows)
{
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
	if (row >= num_rows)
		return;
	OBF_vertex vertex = m[row];
	if (!vertex.isInOWCTY() || vertex.getOldRange() != OldRange)
		return;

	if (vertex.isElim()) {
		vertex.setOldRange(row);
		vertex.setInSCC();
		//testing	vertex.setRange( 0 );
	}
	else {
		if (vertex.isReached()) {
			vertex.setInBWD();
			vertex.setBWDVisited();
			*pivot = row;
		}
		else
			vertex.setInBWD();
	}
	m[row] = vertex;
}
//
//
//
__global__ void FWD_pivot2(OBF_vertex * m, const uint32_t OldRange, uint32_t * pivot, const uint32_t num_rows)
{
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
	if (row >= num_rows)
		return;
	OBF_vertex vertex = m[row];
	if (!vertex.isInFWD() || vertex.getOldRange() != OldRange)
		return;
//
	if (vertex.isFWDVisited()) {
		vertex.setOldRange(vertex.getRange());
		OBF_vertex old_pivot = m[vertex.getRange()];
		if (old_pivot.isInBWD())
			vertex.setInBWD();
		else
			vertex.setInOWCTY();
	}
	else {
		if (*pivot == row) {
			vertex.setFWDVisited();
			vertex.setRange(row);
			vertex.setDone2();
		}
	}
	m[row] = vertex;
}

// void FWDSynchKernel_FWDpivot1_c(OBF_vertex * m, const uint32_t OldRange, const Edge * Bc, const uint32_t * Br, uint32_t * pivot,
//	const uint32_t num_rows, uint32_t * COL_pivot = NULL)
//{
//	//uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
//	 for (uint32_t row = 0; row <= num_rows; row++)
//	 {
//		 if (row >= num_rows)
//			 return;
//		 OBF_vertex vertex = m[row];
//		 if (!vertex.isInFWD() || vertex.getOldRange() != OldRange)
//			 return;
//
//		 if (vertex.isFWDVisited()) {
//			 if (row == vertex.getRange()) {
//				 if (COL_pivot)
//					 *COL_pivot = row;
//
//				 uint32_t row_begin = Br[row];
//				 uint32_t row_end = Br[row + 1];
//
//				 bool skip = false;
//				 for (uint32_t column = row_begin; column < row_end; column++) {
//					 uint32_t index = Bc[column].getValue();
//					 OBF_vertex p_vertex = m[index];
//
//					 if (p_vertex.isInFWD() && p_vertex.getOldRange() == vertex.getOldRange() && p_vertex.isFWDVisited()) {
//						 skip = true;
//						 break;
//					 }
//				 }
//				 if (skip) {
//					 vertex.setInBWD();
//					 vertex.setBWDVisited();
//				 }
//				 else {
//					 vertex.setInOWCTY();
//					 vertex.setReached();
//				 }
//				 vertex.setOldRange(vertex.getRange());
//				 vertex.setDone2();
//
//				 m[row] = vertex;
//			 }
//		 }
//		 else {
//			 if (!*pivot)
//				 *pivot = row;
//		 }
//	 }
//}
//
//
__global__ void FWDSynchKernel_FWDpivot1(OBF_vertex * m, const uint32_t OldRange, const Edge * Bc, const uint32_t * Br, uint32_t * pivot,
	const uint32_t num_rows, uint32_t * COL_pivot = NULL)
{
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
	if (row >= num_rows)
		return;
	OBF_vertex vertex = m[row];
	if (!vertex.isInFWD() || vertex.getOldRange() != OldRange)
		return;

	if (vertex.isFWDVisited()) {
		if (row == vertex.getRange()) {
			if (COL_pivot)
				*COL_pivot = row;

			uint32_t row_begin = Br[row];
			uint32_t row_end = Br[row + 1];

			bool skip = false;
			for (uint32_t column = row_begin; column < row_end; column++) {
				uint32_t index = Bc[column].getValue();
				OBF_vertex p_vertex = m[index];

				if (p_vertex.isInFWD() && p_vertex.getOldRange() == vertex.getOldRange() && p_vertex.isFWDVisited()) {
					skip = true;
					break;
				}
			}
			if (skip) {
				vertex.setInBWD();
				vertex.setBWDVisited();
			}
			else {
				vertex.setInOWCTY();
				vertex.setReached();
			}
			vertex.setOldRange(vertex.getRange());
			vertex.setDone2();

			m[row] = vertex;
		}
	}
	else {
		if (!*pivot)
			*pivot = row;
	}
}


//
//
// void OBFKernel_c(const Edge * Fc, const uint32_t * Fr, const Edge * Bc, const uint32_t * Br, OBF_vertex * m, uint32_t * Temps,
//	const uint32_t num_rows)
//{
//	//uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
//	 for (uint32_t row = 0; row <= num_rows; row++)
//	 {
//
//
//		 if (row >= num_rows)
//			 return;
//		 OBF_vertex vertex = m[row];
//
//		 int phase = vertex.getPhase();
//		 if (row == 32734)
//			 printf("%d", vertex.getRange());
//		 if (phase != 3 && vertex.getRange() == row) {
//			 if (!vertex.isDone1() && !vertex.isDone2())
//				 Temps[phase] = vertex.getOldRange() + (phase == 0);
//			 else {
//				 if (vertex.isDone1()) {
//					 vertex.unsetDone1();
//					 vertex.setDone2();
//					 m[row] = vertex;
//				 }
//				 else
//					 m[row].unsetDone2();
//			 }
//		 }
//
//		 uint32_t row_begin;
//		 uint32_t row_end;
//		 bool prop;
//		 switch (phase) {
//		 case 0: //FWD
//			 if (!vertex.isFWDVisited() || vertex.isFWDPropagate())
//				 continue;
//
//			 row_begin = Fr[row];
//			 row_end = Fr[row + 1];
//
//			 prop = false;
//			 for (uint32_t column = row_begin; column < row_end; column++) {
//				 uint32_t index = Fc[column].getValue();
//				 OBF_vertex p_vertex = m[index];
//
//				 if (!p_vertex.isInFWD() || p_vertex.getOldRange() != vertex.getOldRange() || p_vertex.isFWDVisited())
//					 continue;
//
//				 p_vertex.setFWDVisited();
//				 p_vertex.setRange(vertex.getRange());
//				 m[index] = p_vertex;
//				 prop = true;
//			 }
//			 vertex.setFWDPropagate();
//			 m[row] = vertex;
//			 if (prop)
//				 m[vertex.getRange()].setDone1();
//
//			 break;
//		 case 1: //OWCTY
//			 if (!vertex.isReached() || vertex.isElim())
//				 continue;
//
//			 row_begin = Br[row];
//			 row_end = Br[row + 1];
//
//			 prop = true;
//			 for (uint32_t column = row_begin; column < row_end; column++) {
//				 uint32_t index = Bc[column].getValue();
//				 OBF_vertex p_vertex = m[index];
//
//				 if (p_vertex.isInOWCTY() && p_vertex.getOldRange() == vertex.getOldRange() && !p_vertex.isElim())// && index != row )
//					 prop = false;
//			 }
//			 if (prop) {
//				 row_begin = Fr[row];
//				 row_end = Fr[row + 1];
//
//				 for (uint32_t column = row_begin; column < row_end; column++) {
//					 uint32_t index = Fc[column].getValue();
//					 OBF_vertex p_vertex = m[index];
//
//					 if (!p_vertex.isInOWCTY() || p_vertex.getOldRange() != vertex.getOldRange() || p_vertex.isReached())
//						 continue;
//
//					 p_vertex.setReached();
//					 p_vertex.setRange(vertex.getRange());
//					 m[index] = p_vertex;
//				 }
//				 vertex.setElim();
//				 m[row] = vertex;
//				 m[vertex.getRange()].setDone1();
//			 }
//
//			 break;
//		 case 2: //BWD
//			 if (!vertex.isBWDVisited() || vertex.isBWDPropagate())
//				 continue;
//
//			 row_begin = Br[row];
//			 row_end = Br[row + 1];
//
//			 prop = false;
//			 for (uint32_t column = row_begin; column < row_end; column++) {
//				 uint32_t index = Bc[column].getValue();
//				 OBF_vertex p_vertex = m[index];
//
//				 if (!p_vertex.isInBWD() || p_vertex.getOldRange() != vertex.getOldRange() || p_vertex.isBWDVisited())
//					 continue;
//
//				 p_vertex.setBWDVisited();
//				 p_vertex.setRange(vertex.getRange());
//				 m[index] = p_vertex;
//				 prop = true;
//			 }
//			 vertex.setBWDPropagate();
//			 m[row] = vertex;
//			 if (prop)
//				 m[vertex.getRange()].setDone1();
//
//			 break;
//		 }
//	 }
//}
//
//
__global__ void OBFKernel(const Edge * Fc, const uint32_t * Fr, const Edge * Bc, const uint32_t * Br, OBF_vertex * m, uint32_t * Temps,
	const uint32_t num_rows)
{
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
	if (row >= num_rows)
		return;
	OBF_vertex vertex = m[row];

	int phase = vertex.getPhase();

	if (phase != 3 && vertex.getRange() == row) {
		if (!vertex.isDone1() && !vertex.isDone2())
			Temps[phase] = vertex.getOldRange() + (phase == 0);
		else {
			if (vertex.isDone1()) {
				vertex.unsetDone1();
				vertex.setDone2();
				m[row] = vertex;
			}
			else
				m[row].unsetDone2();
		}
	}
//
	uint32_t row_begin;
	uint32_t row_end;
	bool prop;
	switch (phase) {
	case 0: //FWD
		if (!vertex.isFWDVisited() || vertex.isFWDPropagate())
			return;

		row_begin = Fr[row];
		row_end = Fr[row + 1];
//
		prop = false;
		for (uint32_t column = row_begin; column < row_end; column++) {
			uint32_t index = Fc[column].getValue();
			OBF_vertex p_vertex = m[index];

			if (!p_vertex.isInFWD() || p_vertex.getOldRange() != vertex.getOldRange() || p_vertex.isFWDVisited())
				continue;

			p_vertex.setFWDVisited();
			p_vertex.setRange(vertex.getRange());
			m[index] = p_vertex;
			prop = true;
		}
		vertex.setFWDPropagate();
		m[row] = vertex;
		if (prop)
			m[vertex.getRange()].setDone1();

		break;
	case 1: //OWCTY
		if (!vertex.isReached() || vertex.isElim())
			return;

		row_begin = Br[row];
		row_end = Br[row + 1];
//
		prop = true;
		for (uint32_t column = row_begin; column < row_end; column++) {
			uint32_t index = Bc[column].getValue();
			OBF_vertex p_vertex = m[index];

			if (p_vertex.isInOWCTY() && p_vertex.getOldRange() == vertex.getOldRange() && !p_vertex.isElim())// && index != row )
				prop = false;
		}
		if (prop) {
			row_begin = Fr[row];
			row_end = Fr[row + 1];
//
			for (uint32_t column = row_begin; column < row_end; column++) {
				uint32_t index = Fc[column].getValue();
				OBF_vertex p_vertex = m[index];

				if (!p_vertex.isInOWCTY() || p_vertex.getOldRange() != vertex.getOldRange() || p_vertex.isReached())
					continue;

				p_vertex.setReached();
				p_vertex.setRange(vertex.getRange());
				m[index] = p_vertex;
			}
			vertex.setElim();
			m[row] = vertex;
			m[vertex.getRange()].setDone1();
		}
//
		break;
	case 2: //BWD
		if (!vertex.isBWDVisited() || vertex.isBWDPropagate())
			return;

		row_begin = Br[row];
		row_end = Br[row + 1];

		prop = false;
		for (uint32_t column = row_begin; column < row_end; column++) {
			uint32_t index = Bc[column].getValue();
			OBF_vertex p_vertex = m[index];

			if (!p_vertex.isInBWD() || p_vertex.getOldRange() != vertex.getOldRange() || p_vertex.isBWDVisited())
				continue;

			p_vertex.setBWDVisited();
			p_vertex.setRange(vertex.getRange());
			m[index] = p_vertex;
			prop = true;
		}
		vertex.setBWDPropagate();
		m[row] = vertex;
		if (prop)
			m[vertex.getRange()].setDone1();

		break;
	}
}
//
//
// void _Trimming_c(OBF_vertex * m, const Edge * Bc, const uint32_t * Br, const uint32_t num_rows,
//	uint32_t * terminate)
//{
//	//uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
//	 for (uint32_t row = 0; row <= num_rows; row++)
//	 {
//
//		 if (row >= num_rows)
//			 return;
//		 OBF_vertex vertex = m[row];
//
//		 if (vertex.isInSCC())
//		 {
//			 continue;
//
//		 }
//
//
//		 uint32_t row_begin = Br[row];
//		 uint32_t row_end = Br[row + 1];
//
//		 bool eliminate = true;
//		 for (uint32_t column = row_begin; column < row_end; column++) {
//			 //printf("col = %d ,old range = %d \n", column, vertex.getOldRange());
//
//			 uint32_t index = Bc[column].getValue();
//			 OBF_vertex p_vertex = m[index];
//
//			 if (!p_vertex.isInSCC())
//			 {
//				 eliminate = false;
//			 }
//
//		 }
//		 /*	if ( !eliminate ) {
//		 eliminate = true;
//		 row_begin = Fr[ row ];
//		 row_end = Fr[ row + 1 ];
//
//		 for ( uint32_t column = row_begin; column < row_end; column++ ) {
//		 uint32_t index = Fc[ column ].getValue();
//		 OBF_vertex p_vertex = m[ index ];
//
//		 if ( !p_vertex.isInSCC() )
//		 eliminate = false;
//		 }
//		 }*/
//		 if (eliminate) {
//
//			 vertex.setOldRange(row);
//			 vertex.setInSCC();
//			 m[row] = vertex;
//			 *terminate = 1;
//		 }
//	 }
//	return;
//}
//
//


__global__ void _Trimming(OBF_vertex * m, const Edge * Bc, const uint32_t * Br, const uint32_t num_rows,
	uint32_t * terminate)
{
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
	if (row >= num_rows)
		return;
	OBF_vertex vertex = m[row];

	if (vertex.isInSCC())
	{
		return;

	}


	uint32_t row_begin = Br[row];
	uint32_t row_end = Br[row + 1];

	bool eliminate = true;
	for (uint32_t column = row_begin; column < row_end; column++) {
		//printf("col = %d ,old range = %d \n", column, vertex.getOldRange());

		uint32_t index = Bc[column].getValue();
		OBF_vertex p_vertex = m[index];

		if (!p_vertex.isInSCC())
		{
			eliminate = false;
		}

	}
	/*	if ( !eliminate ) {
	eliminate = true;
	row_begin = Fr[ row ];
	row_end = Fr[ row + 1 ];

	for ( uint32_t column = row_begin; column < row_end; column++ ) {
	uint32_t index = Fc[ column ].getValue();
	OBF_vertex p_vertex = m[ index ];

	if ( !p_vertex.isInSCC() )
	eliminate = false;
	}
	}*/
	if (eliminate) {

		vertex.setOldRange(row);
		vertex.setInSCC();
		m[row] = vertex;
		*terminate = 1;
	}
	return;
}



//
// void CheckTerminateSetKernel_c(const OBF_vertex * m, uint32_t * terminate, const uint32_t num_rows)
//{
//	//uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
//	 for (uint32_t row = 0; row <= num_rows; row++)
//	 {
//		 if (row >= num_rows)
//			 return;
//		 OBF_vertex vertex = m[row];
//		 if (!vertex.isInSCC())
//			 *terminate = row;
//		 else
//			 printf(" %d ", row);
//	 }
//}
//
__global__ void CheckTerminateSetKernel(const OBF_vertex * m, uint32_t * terminate, const uint32_t num_rows)
{
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
	if (row >= num_rows)
		return;
	OBF_vertex vertex = m[row];
	if (!vertex.isInSCC())
		*terminate = row;
}




pair <uint32_t, float> OBF_Decomposition(uint32_t CSize, uint32_t RSize, Edge *Fc, uint32_t *Fr, Edge * Bc, uint32_t * Br, OBF_vertex * m,ASF_vertex* o_m,Dimension*d,
	bool COL_enforced, uint32_t COL_limit, bool trimm, bool quiet)
{
	//-----------GPU initialization---------------------------->
	uint32_t  * d_Fr, *d_Br;
	Edge * d_Fc, *d_Bc;
	OBF_vertex * d_m;
	float temp;

	uint32_t * d_pivot;
	uint32_t * d_Temps;//[0]-d_FWD, [1]-d_OWCTY, [2]-d_BWD
	uint32_t Temps[Temp_count];
	uint32_t terminate = 1;
	int interruptions = 0;

#ifdef _DEBUG
	int FWD_ints = 0;
	int OWCTY_ints = 0;
	int BWD_ints = 0;
	StopWatchInterface* KernelTime = 0;
	StopWatchInterface* IntTime = 0;
	(sdkCreateTimer(&KernelTime));
	(sdkCreateTimer(&IntTime));
#endif

	if (!_DeviceSet) {
		_DeviceSet = true;
		checkCudaErrors(cudaSetDevice(2));
	}

	if (COL_enforced)
		COL_limit = RSize / 100;

	cudaError_t e1, e2, e3, e4, e5, e6, e7;
	checkCudaErrors(e1 = cudaMalloc((void**)&d_Fc, CSize * sizeof(Edge)));
	checkCudaErrors(e2 = cudaMalloc((void**)&d_Fr, RSize * sizeof(uint32_t)));
	checkCudaErrors(e3 = cudaMalloc((void**)&d_Bc, CSize * sizeof(Edge)));
	checkCudaErrors(e4 = cudaMalloc((void**)&d_Br, RSize * sizeof(uint32_t)));
	checkCudaErrors(e5 = cudaMalloc((void**)&d_m, (RSize - 1) * sizeof(OBF_vertex)));
	checkCudaErrors(e6 = cudaMalloc((void**)&d_Temps, Temp_count * sizeof(uint32_t)));
	checkCudaErrors(e7 = cudaMalloc((void**)&d_pivot, sizeof(uint32_t)));

	if (e1 == cudaErrorMemoryAllocation || e2 == cudaErrorMemoryAllocation ||
		e3 == cudaErrorMemoryAllocation || e4 == cudaErrorMemoryAllocation ||
		e5 == cudaErrorMemoryAllocation || e6 == cudaErrorMemoryAllocation ||
		e7 == cudaErrorMemoryAllocation) {
		throw "Error: Not enough memory on GPU\n";
	}

	//col
	//unsigned int COLTime = 0;
	StopWatchInterface* COLTime = 0;
	(sdkCreateTimer(&COLTime));
	uint32_t * d_temp_COL;
	uint32_t * d_temp_COL2;
	uint32_t * d_COL_OldRange;
	bool COL_used = false;
	COL_vertex * d_cm;

	checkCudaErrors(e1 = cudaMalloc((void**)&d_temp_COL, (RSize - 2) * sizeof(uint32_t)));
	checkCudaErrors(e2 = cudaMalloc((void**)&d_temp_COL2, 64 * sizeof(uint32_t)));
	checkCudaErrors(e3 = cudaMalloc((void**)&d_COL_OldRange, sizeof(uint32_t)));
	if (e1 == cudaErrorMemoryAllocation || e2 == cudaErrorMemoryAllocation || e3 == cudaErrorMemoryAllocation) {
		throw "Error: Not enough memory on GPU\n";
	}

	checkCudaErrors(cudaMemcpy(d_Fc, Fc, CSize * sizeof(Edge), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_Fr, Fr, RSize * sizeof(uint32_t), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_Bc, Bc, CSize * sizeof(Edge), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_Br, Br, RSize * sizeof(uint32_t), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_m, m, (RSize - 1) * sizeof(OBF_vertex), cudaMemcpyHostToDevice));

	//unsigned int SCCTime = 0;
	StopWatchInterface* SCCTime = 0;
	(sdkCreateTimer(&SCCTime));
	(sdkStartTimer(&SCCTime));

	dim3 grid(((RSize + 510) / 512), 1, 1);  // (RSize-2) valid vertecis the first one is undefined
	dim3 threads(512, 1, 1);
	dim3 grid1(64, 1, 1);
	dim3 threads1(blockSize, 1, 1);
	dim3 grid2(1, 1, 1);
	dim3 threads2(32, 1, 1);
	//<----------GPU initialization-----------------------------
	//if (!quiet)
	{
		// 		printf("Computing OBF decomposition \n");
		// 		printf("Size: %u, blockdim.x: %u, threaddim.x: %u\n", RSize - 1, grid.x, threads.x);
		printf("Vertices: %u\n", RSize - 2);
		printf("Edges: %u\n", CSize);
	}
	//-----------Main algorithm-------------------------------->
	OBF_vertex initial;
	initial.setInSCC();
	checkCudaErrors(cudaMemcpy(&(d_m[0]), &initial, sizeof(OBF_vertex), cudaMemcpyHostToDevice));
	//-----------Trimming-------------------------------------->
	if (trimm) {
		uint32_t tits = 1;
		do {
			if (!(tits++ & 0x7FF)) {
				temp = sdkGetTimerValue(&SCCTime);
				if (temp > Time_limit)
					break;
			}
			checkCudaErrors(cudaMemset(d_pivot, 0, sizeof(uint32_t)));
			_Trimming << <grid, threads >> >(d_m, d_Bc, d_Br, RSize - 1, d_pivot);
			checkCudaErrors(cudaMemcpy(&terminate, d_pivot, sizeof(uint32_t), cudaMemcpyDeviceToHost));
		} while (terminate);
		checkCudaErrors(cudaMemset(d_pivot, 0, sizeof(uint32_t)));
		CheckTerminateSetKernel << <grid, threads >> >(d_m, d_pivot, RSize - 1);
		checkCudaErrors(cudaMemcpy(&terminate, d_pivot, sizeof(uint32_t), cudaMemcpyDeviceToHost));
	}
	//<----------Trimming---------------------------------------

	//-----------Initial FWD pivot setup----------------------->
	if (terminate) {
		initial.setInFWD();
		initial.setRange(terminate);
		initial.setFWDVisited();
		initial.setDone2();
		checkCudaErrors(cudaMemcpy(&(d_m[terminate]), &initial, sizeof(OBF_vertex), cudaMemcpyHostToDevice));
		terminate = 0;
	}
	else
		terminate = 1;
	//<----------Initial FWD pivot setup------------------------

#ifdef _DEBUG
	long int its = 0;
#endif

	int i;
	while (!terminate) {

#ifdef _DEBUG
		its++;
		//printf("its: %d\n", its);
		cudaThreadSynchronize();
		(sdkStartTimer(&KernelTime));
#endif


		checkCudaErrors(cudaMemset(d_Temps, 0, Temp_count * sizeof(uint32_t)));
		OBFKernel << <grid, threads >> >(d_Fc, d_Fr, d_Bc, d_Br, d_m, d_Temps, RSize - 1);
		checkCudaErrors(cudaMemcpy(Temps, d_Temps, Temp_count * sizeof(uint32_t), cudaMemcpyDeviceToHost));

#ifdef _DEBUG
		cudaThreadSynchronize();
		(sdkStopTimer(&KernelTime));
		(sdkStartTimer(&IntTime));
#endif

		//-----------Interrupt handling---------------------------->
		for (i = 0; i < Temp_count; i++)
		if (Temps[i])
			break;

		if (i == Temp_count) {
			continue;
		}

		interruptions++;
		//printf("inters: %d\n", i);

		//time_limit
		if (!(interruptions & 0x7FF)) {
			temp = sdkGetTimerValue(&SCCTime);
			if (temp > Time_limit)
				break;
		}

		checkCudaErrors(cudaMemset(d_pivot, 0, sizeof(uint32_t)));
		switch (i) {
		case 0:

#ifdef _DEBUG
			FWD_ints++;
#endif

			Temps[0]--;

			if (COL_enforced)
				FWDSynchKernel_FWDpivot1 << <grid, threads >> >(d_m, Temps[0], d_Bc, d_Br, d_pivot, RSize - 1, d_COL_OldRange);
			else
				FWDSynchKernel_FWDpivot1 << <grid, threads >> >(d_m, Temps[0], d_Bc, d_Br, d_pivot, RSize - 1);
			FWD_pivot2 << <grid, threads >> >(d_m, Temps[0], d_pivot, RSize - 1);
			break;
		case 1:

#ifdef _DEBUG
			OWCTY_ints++;
#endif
			OWCTYSynchKernel_BWDpivot1 << <grid, threads >> >(d_m, d_Bc, d_Br, Temps[1], d_pivot, RSize - 1);
			BWD_pivot2 << <grid, threads >> >(d_m, Temps[1], d_pivot, RSize - 1);
			if (COL_enforced)
				checkCudaErrors(cudaMemcpy(d_COL_OldRange, &(d_Temps[1]), sizeof(uint32_t), cudaMemcpyDeviceToDevice));
			break;
		case 2:

#ifdef _DEBUG
			BWD_ints++;
#endif
			BWDTerminateKernel_OWCTY1 << <grid, threads >> >(d_m, d_Fc, d_Fr, Temps[2], RSize - 1);
			BWDSynchKernel_FWDOWCTYpivot1 << <grid, threads >> >(d_m, Temps[2], &(d_Temps[2]), d_pivot, &(d_Temps[0]), RSize - 1);
			FWDOWCTY_pivot2 << <grid, threads >> >(d_m, &(d_Temps[2]), Temps[2], d_pivot, &(d_Temps[0]), RSize - 1);
			if (COL_enforced)
				checkCudaErrors(cudaMemcpy(d_COL_OldRange, &(d_Temps[2]), sizeof(uint32_t), cudaMemcpyDeviceToDevice));
			break;
		}
		if (COL_enforced) {

#ifdef _DEBUG
			cudaThreadSynchronize();
			(sdkStartTimer(&COLTime));
			(sdkStopTimer(&IntTime));
#endif

			checkCudaErrors(cudaMemset(d_temp_COL, 0, (RSize - 2) * sizeof(uint32_t)));
			OBFCOLcompute << <grid, threads >> >(d_m, d_temp_COL, d_COL_OldRange, RSize - 1);
			cuReduce << <grid1, threads1 >> >(d_temp_COL, d_temp_COL2, RSize - 2);
			cuReduce1 << <grid2, threads2 >> >(d_temp_COL2, d_pivot);
			checkCudaErrors(cudaMemcpy(&terminate, d_pivot, sizeof(uint32_t), cudaMemcpyDeviceToHost));
			if (terminate < COL_limit && terminate != 0) {
				COL_used = true;
				SetCOL << <grid, threads >> >(d_m, d_COL_OldRange, RSize - 1);
			}

#ifdef _DEBUG
			cudaThreadSynchronize();
			(sdkStopTimer(&COLTime));
			(sdkStartTimer(&IntTime));
#endif

		}
		//<----------Interrupt handling-----------------------------

		//-----------Termination detection------------------------->
		checkCudaErrors(cudaMemset(d_pivot, 1, sizeof(uint32_t)));
		CheckTerminateKernel << <grid, threads >> >(d_m, d_pivot, RSize - 1);
		checkCudaErrors(cudaMemcpy(&terminate, d_pivot, sizeof(uint32_t), cudaMemcpyDeviceToHost));
		//<----------Termination detection--------------------------

#ifdef _DEBUG
		cudaThreadSynchronize();
		(sdkStopTimer(&IntTime));
#endif

	}
	//<----------Main algorithm---------------------------------

	//-----------COL finish------------------------------------>
	if (COL_enforced && COL_used) {

#ifdef _DEBUG
		cudaThreadSynchronize();
		(sdkStartTimer(&COLTime));
#endif

		checkCudaErrors(cudaMalloc((void**)&d_cm, (RSize - 1) * sizeof(COL_vertex)));
		OBFtoCOLKernel << <grid, threads >> >(d_m, d_cm, RSize - 1);

		int c_its = 0;
		do {
			c_its++;

			//time_limit
		/*	if (!(c_its & 0x7)) {
				temp = sdkGetTimerValue(&SCCTime);
				if (temp > Time_limit)
					break;
			}*/

			do {
				checkCudaErrors(cudaMemset(d_pivot, 0, sizeof(uint32_t)));
				p_COL_MAP << <grid, threads >> >(d_Fc, d_Fr, d_cm, RSize - 1, MIN, d_pivot);
				checkCudaErrors(cudaMemcpy(&terminate, d_pivot, sizeof(uint32_t), cudaMemcpyDeviceToHost));
			} while (terminate != 0);

			setBWDseed << <grid, threads >> >(d_cm, RSize - 1, MIN);
			do {
				checkCudaErrors(cudaMemset(d_pivot, 0, sizeof(uint32_t)));
				COL_FWD << <grid, threads >> >(d_Fc, d_Fr, d_cm, RSize - 1, d_pivot);
				checkCudaErrors(cudaMemcpy(&terminate, d_pivot, sizeof(uint32_t), cudaMemcpyDeviceToHost));
			} while (terminate != 0);

			checkCudaErrors(cudaMemset(d_pivot, 0, sizeof(uint32_t)));
			checkTerminateAndSetOldMap << <grid, threads >> >(d_cm, RSize - 1, d_pivot);
			checkCudaErrors(cudaMemcpy(&terminate, d_pivot, sizeof(uint32_t), cudaMemcpyDeviceToHost));
		} while (terminate != 0);

#ifdef _DEBUG
		cudaThreadSynchronize();
		(sdkStopTimer(&COLTime));
#endif

	}
	//<----------COL finish-------------------------------------

	//-----------Scc extraction-------------------------------->


	uint32_t* Fr_F = new uint32_t[RSize];
	memset(Fr_F, 0, RSize*sizeof(uint32_t));
	checkCudaErrors(cudaMemset(d_Fr, 0, (RSize - 2) * sizeof(uint32_t)));
	if (COL_enforced && COL_used)
		COLcomputeSCCs << <grid, threads >> >(d_cm, d_Fr, RSize - 1, MIN);
	else
	{
		checkCudaErrors(cudaMemcpy(Fr, d_Fr, RSize * sizeof(uint32_t), cudaMemcpyDeviceToHost));
		checkCudaErrors(cudaMemcpy(m, d_m, (RSize - 1) * sizeof(OBF_vertex), cudaMemcpyDeviceToHost));

		//OBFcomputeSCCs << <grid, threads >> >(d_m, d_Fr, RSize - 1);
		OBFcomputeSCCs_c(m,o_m,d, Fr, RSize - 1);

	}
	//checkCudaErrors(cudaMemcpy(Fr, d_Fr, RSize * sizeof(uint32_t), cudaMemcpyDeviceToHost));

	int cc = 0;
	for (int jj = 0; jj < RSize-2; jj++)
	{
		if (Fr[jj] == 0)
			cc++;// printf("%d ", jj);
	}
	/*cuReduce << <grid1, threads1 >> >(d_Fr, d_Br, RSize - 2);
	checkCudaErrors(cudaMemcpy(Br, d_Br, RSize * sizeof(uint32_t), cudaMemcpyDeviceToHost));

	cuReduce1 << <grid2, threads2 >> >(d_Br, d_pivot);
	checkCudaErrors(cudaMemcpy(&terminate, d_pivot, sizeof(uint32_t), cudaMemcpyDeviceToHost));*/
	//<----------Scc extraction---------------------------------
	(sdkStopTimer(&SCCTime));
	float f;
	//if (!quiet)
		//		printf("%u components found.\n", terminate);
		printf("Components: %u\n", cc);

#ifdef _DEBUG
	printf("Kernel iterations = %d\n", its);
	printf("Interruptions =  %d, FWD interruptions = %d, OWCTY interruptions = %d, BWD interruptions = %d \n", interruptions, FWD_ints,
		OWCTY_ints, BWD_ints);
	f = sdkGetTimerValue(&KernelTime);
	printf("Uninterrupted time = %f ms\n", f);
	f = sdkGetTimerValue(&IntTime);
	printf("Interruption time = %f ms\n", f);
	if (COL_enforced) {
		f = sdkGetTimerValue(&COLTime);
		printf("Colouring time = %f ms\n", f);
	}
#endif

	f = sdkGetTimerValue(&SCCTime);
	int min = (int)(f / 60000.0f);
	int sec = (int)(f / 1000.0f) % 60;
	//if (!quiet)
		//		printf("CUDA SCC decomposition time: %d minutes %d seconds (%f ms).\n", min, sec, f);
		printf("Time: %f ms\n", f);

	//CUT_CHECK_ERROR("Kernel execution failed");

	checkCudaErrors(cudaFree(d_Fc));
	checkCudaErrors(cudaFree(d_Fr));
	checkCudaErrors(cudaFree(d_Bc));
	checkCudaErrors(cudaFree(d_Br));
	checkCudaErrors(cudaFree(d_m));
	checkCudaErrors(cudaFree(d_Temps));
	checkCudaErrors(cudaFree(d_pivot));
	if (COL_used)
		checkCudaErrors(cudaFree(d_cm));
	checkCudaErrors(cudaFree(d_temp_COL));
	checkCudaErrors(cudaFree(d_temp_COL2));
	checkCudaErrors(cudaFree(d_COL_OldRange));

	(sdkDeleteTimer(&SCCTime));

#ifdef _DEBUG
	(sdkDeleteTimer(&KernelTime));
	(sdkDeleteTimer(&IntTime));
	if (COL_enforced)
		(sdkDeleteTimer(&COLTime));
#endif

	return make_pair(terminate, f);
}




//
//pair <uint32_t, float> OBF_Decomposition(uint32_t CSize, uint32_t RSize, Edge *Fc, uint32_t *Fr, Edge * Bc, uint32_t * Br, OBF_vertex * m,
//	bool COL_enforced, uint32_t COL_limit, bool trimm, bool quiet)
//{
//	//-----------GPU initialization---------------------------->
//	uint32_t  * d_Fr, *d_Br;
//	Edge * d_Fc, *d_Bc;
//	OBF_vertex * d_m;
//	float temp;
//
//	uint32_t * d_pivot;
//	uint32_t * d_Temps;//[0]-d_FWD, [1]-d_OWCTY, [2]-d_BWD
//	uint32_t Temps[Temp_count];
//	uint32_t terminate = 1;
//	int interruptions = 0;
//	pair <uint32_t, float> result;
//
//
//#ifdef _DEBUG
//int FWD_ints = 0;
//int OWCTY_ints = 0;
//int BWD_ints = 0;
//StopWatchInterface* KernelTime = NULL;
//StopWatchInterface* IntTime = NULL;
//sdkCreateTimer(&KernelTime);
//sdkCreateTimer(&IntTime);
//#endif
//
//	if( !_DeviceSet ) {
//		_DeviceSet = true;
//		checkCudaErrors(cudaSetDevice(1));
//	}
//
//	if ( COL_enforced )
//		COL_limit = RSize / 100;
//
//	cudaError_t e1, e2, e3, e4, e5, e6, e7;
//	checkCudaErrors(e1 = cudaMalloc((void**)&d_Fc, CSize * sizeof(Edge)));
//	checkCudaErrors(e2 = cudaMalloc((void**)&d_Fr, RSize * sizeof(uint32_t)));
//	checkCudaErrors(e3 = cudaMalloc((void**)&d_Bc, CSize * sizeof(Edge)));
//	checkCudaErrors(e4 = cudaMalloc((void**)&d_Br, RSize * sizeof(uint32_t)));
//	checkCudaErrors(e5 = cudaMalloc((void**)&d_m, (RSize - 1) * sizeof(OBF_vertex)));
//	checkCudaErrors(e6 = cudaMalloc((void**)&d_Temps, Temp_count * sizeof(uint32_t)));
//	checkCudaErrors(e7 = cudaMalloc((void**)&d_pivot, sizeof(uint32_t)));
//
//	if (e1 == cudaErrorMemoryAllocation || e2 == cudaErrorMemoryAllocation ||
//		e3 == cudaErrorMemoryAllocation || e4 == cudaErrorMemoryAllocation ||
//		e5 == cudaErrorMemoryAllocation || e6 == cudaErrorMemoryAllocation ||
//		e7 == cudaErrorMemoryAllocation ) {
//		throw "Error: Not enough memory on GPU\n";
//	}
//
////col
//	StopWatchInterface* COLTime = 0;
//	sdkCreateTimer(&COLTime);
//	uint32_t * d_temp_COL;
//	uint32_t * d_temp_COL2;
//	uint32_t * d_COL_OldRange;
//	bool COL_used = false;
//	COL_vertex * d_cm;
//
//	checkCudaErrors(e1 = cudaMalloc((void**)&d_temp_COL, (RSize - 2) * sizeof(uint32_t)));
//	checkCudaErrors(e2 = cudaMalloc((void**)&d_temp_COL2, 64 * sizeof(uint32_t)));
//	checkCudaErrors(e3 = cudaMalloc((void**)&d_COL_OldRange, sizeof(uint32_t)));
//	if (e1 == cudaErrorMemoryAllocation || e2 == cudaErrorMemoryAllocation || e3 == cudaErrorMemoryAllocation ) {
//		throw "Error: Not enough memory on GPU\n";
//	}
//
//	checkCudaErrors(cudaMemcpy(d_Fc, Fc, CSize * sizeof(Edge), cudaMemcpyHostToDevice));
//	checkCudaErrors(cudaMemcpy(d_Fr, Fr, RSize * sizeof(uint32_t), cudaMemcpyHostToDevice));
//	checkCudaErrors(cudaMemcpy(d_Bc, Bc, CSize * sizeof(Edge), cudaMemcpyHostToDevice));
//	checkCudaErrors(cudaMemcpy(d_Br, Br, RSize * sizeof(uint32_t), cudaMemcpyHostToDevice));
//	checkCudaErrors(cudaMemcpy(d_m, m, (RSize - 1) * sizeof(OBF_vertex), cudaMemcpyHostToDevice));
//
//	//unsigned int SCCTime = 0;
//	StopWatchInterface* SCCTime = NULL;
//	sdkCreateTimer(&SCCTime);
//	sdkStartTimer(&SCCTime);
//
//	dim3 grid( ((RSize + 510) / 512), 1, 1 );  // (RSize-2) valid vertecis the first one is undefined
//	dim3 threads(512, 1, 1);
//	dim3 grid1(64, 1, 1);
//	dim3 threads1(blockSize, 1, 1);
//	dim3 grid2(1, 1, 1);
//	dim3 threads2(32, 1, 1);
////<----------GPU initialization-----------------------------
//	if ( !quiet ) {
// 		printf("Computing OBF decomposition \n");
// 		printf("Size: %u, blockdim.x: %u, threaddim.x: %u\n", RSize - 1, grid.x, threads.x);
//		printf( "Vertices: %u\n", RSize - 2 );
//		printf( "Edges: %u\n", CSize );
//	}
////-----------Main algorithm-------------------------------->
//		OBF_vertex initial;
//		initial.setInSCC();
//		checkCudaErrors(cudaMemcpy(&(d_m[0]), &initial, sizeof(OBF_vertex), cudaMemcpyHostToDevice));
////-----------Trimming-------------------------------------->
//	if ( trimm ) {
//		uint32_t tits = 1;
//		do {
//			if ( !(tits++ & 0x7FF) ) {
//				temp = sdkGetTimerValue( &SCCTime );
//				if ( temp > Time_limit )
//					break;
//			}
//			checkCudaErrors(cudaMemset(d_pivot, 0, sizeof(uint32_t)));
//		//	_Trimming_c(m,Bc,Br,RSize - 1,&terminate);
//			_Trimming<<<grid, threads>>>( d_m, d_Bc, d_Br, RSize - 1, d_pivot );
//			checkCudaErrors(cudaMemcpy(&terminate, d_pivot, sizeof(uint32_t), cudaMemcpyDeviceToHost));
//		} while ( terminate );
//		checkCudaErrors(cudaMemset(d_pivot, 0, sizeof(uint32_t)));
//		CheckTerminateSetKernel << <grid, threads >> >(d_m, d_pivot, RSize - 1);
//		//CheckTerminateSetKernel_c (m, &terminate, RSize - 1);
//		checkCudaErrors( cudaMemcpy( &terminate, d_pivot, sizeof(uint32_t), cudaMemcpyDeviceToHost ));
//	}
////<----------Trimming---------------------------------------
//
//
//
////-----------Initial FWD pivot setup----------------------->
//	if ( terminate ) {
//		initial.setInFWD();
//		initial.setRange( terminate );
//		initial.setFWDVisited();
//		initial.setDone2();
//		checkCudaErrors( cudaMemcpy( &(d_m[ terminate ]), &initial, sizeof(OBF_vertex), cudaMemcpyHostToDevice ));
//		terminate = 0;
//	}
//	else
//		terminate = 1;
////<----------Initial FWD pivot setup------------------------
//
//#ifdef _DEBUG
//long int its = 0;
//#endif
//
//checkCudaErrors(cudaMemcpy(Fc, d_Fc, CSize * sizeof(Edge), cudaMemcpyDeviceToHost));
//checkCudaErrors(cudaMemcpy(Fr, d_Fr, RSize * sizeof(uint32_t), cudaMemcpyDeviceToHost));
//checkCudaErrors(cudaMemcpy(Bc, d_Bc, CSize * sizeof(Edge), cudaMemcpyDeviceToHost));
//checkCudaErrors(cudaMemcpy(Br, d_Br, RSize * sizeof(uint32_t), cudaMemcpyDeviceToHost));
//checkCudaErrors(cudaMemcpy(m, d_m, (RSize - 1) * sizeof(OBF_vertex), cudaMemcpyDeviceToHost));
//
//	int i;
//	while ( !terminate ) {
//
//#ifdef _DEBUG
//its++;
////printf("its: %d\n", its);
//cudaThreadSynchronize();
//(sdkStartTimer(&KernelTime));
//#endif
//
//
//		checkCudaErrors( cudaMemset( d_Temps, 0, Temp_count * sizeof(uint32_t) ));
//		//OBFKernel_c( Fc, Fr, Bc, Br, m, Temps, RSize - 1 );
//		OBFKernel<<<grid, threads>>>( d_Fc, d_Fr, d_Bc, d_Br, d_m, d_Temps, RSize - 1 );
//		checkCudaErrors(cudaMemcpy(Temps, d_Temps, Temp_count * sizeof(uint32_t), cudaMemcpyDeviceToHost));
//
//#ifdef _DEBUG
//cudaThreadSynchronize();
//( sdkStopTimer( &KernelTime ));
//( sdkStartTimer( &IntTime ));
//#endif
//
////-----------Interrupt handling---------------------------->
//		for ( i = 0; i < Temp_count; i++ )
//			if ( Temps[ i ] )
//				break;
//
//		if ( i == Temp_count ) {
//			continue;
//		}
//
//interruptions++;
////printf("inters: %d\n", i);
//
////time_limit
//		if ( !(interruptions & 0x7FF) ) {
//			temp = sdkGetTimerValue( &SCCTime );
//			if ( temp > Time_limit )
//				break;
//		}
//
//		checkCudaErrors( cudaMemset( d_pivot, 0, sizeof(uint32_t) ));
//		switch ( i ) {
//			case 0:
//
//#ifdef _DEBUG
//FWD_ints++;
//#endif
//
//				Temps[ 0 ]--;
//
//				if ( COL_enforced )
//					FWDSynchKernel_FWDpivot1<<<grid, threads>>>( d_m, Temps[ 0 ], d_Bc, d_Br, d_pivot, RSize - 1, d_COL_OldRange );
//				else
//					FWDSynchKernel_FWDpivot1<<<grid, threads>>>( d_m, Temps[ 0 ], d_Bc, d_Br, d_pivot, RSize - 1 );
//				FWD_pivot2<<<grid, threads>>>( d_m, Temps[ 0 ], d_pivot, RSize - 1 );
//				break;
//			case 1:
//
//#ifdef _DEBUG
//OWCTY_ints++;
//#endif
//				OWCTYSynchKernel_BWDpivot1<<<grid, threads>>>( d_m, d_Bc, d_Br, Temps[ 1 ], d_pivot, RSize - 1 );
//				BWD_pivot2<<<grid, threads>>>( d_m, Temps[ 1 ], d_pivot, RSize - 1 );
//				if ( COL_enforced )
//					checkCudaErrors( cudaMemcpy( d_COL_OldRange, &(d_Temps[ 1 ]), sizeof(uint32_t), cudaMemcpyDeviceToDevice ));
//				break;
//			case 2:
//
//#ifdef _DEBUG
//BWD_ints++;
//#endif
//				BWDTerminateKernel_OWCTY1<<<grid, threads>>>( d_m, d_Fc, d_Fr, Temps[ 2 ], RSize - 1 );
//				BWDSynchKernel_FWDOWCTYpivot1<<<grid, threads>>>( d_m, Temps[ 2 ], &(d_Temps[ 2 ]), d_pivot, &(d_Temps[ 0 ]), RSize - 1 );
//				FWDOWCTY_pivot2<<<grid, threads>>>( d_m, &(d_Temps[ 2 ]), Temps[ 2 ], d_pivot, &(d_Temps[ 0 ]), RSize - 1 );
//				if ( COL_enforced )
//					checkCudaErrors( cudaMemcpy( d_COL_OldRange, &(d_Temps[ 2 ]), sizeof(uint32_t), cudaMemcpyDeviceToDevice ));
//				break;
//		}
//		if ( COL_enforced ) {
//
//#ifdef _DEBUG
//cudaThreadSynchronize();
//( sdkStartTimer( &COLTime ));
//( sdkStopTimer( &IntTime ));
//#endif
//
//			checkCudaErrors( cudaMemset( d_temp_COL, 0, (RSize - 2) * sizeof(uint32_t) ));
//			OBFCOLcompute<<<grid, threads>>>( d_m, d_temp_COL, d_COL_OldRange, RSize - 1 );
//			cuReduce<<<grid1, threads1>>>( d_temp_COL, d_temp_COL2, RSize - 2 );
//			cuReduce1<<<grid2, threads2>>>( d_temp_COL2, d_pivot );
//			checkCudaErrors( cudaMemcpy( &terminate, d_pivot, sizeof(uint32_t), cudaMemcpyDeviceToHost ));
//			if ( terminate < COL_limit && terminate != 0 ) {
//				COL_used = true;
//				SetCOL<<<grid, threads>>>( d_m, d_COL_OldRange, RSize - 1 );
//			}
//
//#ifdef _DEBUG
//cudaThreadSynchronize();
//checkCudaErrors( sdkStopTimer( &COLTime ));
//checkCudaErrors( sdkStartTimer( &IntTime ));
//#endif
//
//		}
////<----------Interrupt handling-----------------------------
//
////-----------Termination detection------------------------->
//		checkCudaErrors( cudaMemset( d_pivot, 1, sizeof(uint32_t) ));
//		CheckTerminateKernel<<<grid, threads>>>( d_m, d_pivot, RSize - 1 );
//		checkCudaErrors( cudaMemcpy( &terminate, d_pivot, sizeof(uint32_t), cudaMemcpyDeviceToHost ));
////<----------Termination detection--------------------------
//
//#ifdef _DEBUG
//cudaThreadSynchronize();
//( sdkStopTimer( &IntTime ));
//#endif
//
//	}
//	return result;
//}
