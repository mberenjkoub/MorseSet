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

#include "SCC2.h"
#include "scc_kernel.h"
#define BLOCKSIZE 512
#define sample_seeds 50
//#include "graph_generator.h"
//#include "parallel_fwd.h"
//#include "hash_table.h"

/******************/
/* PLAN STRUCTURE */
/******************/
// --- Async
template<class T>
struct plan {
	T *d_data;
	Boundary* b;
	Dimension*d;
	Point* s;
	Point* v;
	fEdge* eg;
	fFace* fc;
	uint32_t* Fe_Edge;
	uint32_t* Fr_Edge;
	uint32_t* Fe_Face;
	uint32_t* Fr_Face;
	//Stream for asynchronous command execution
	cudaStream_t stream;
};

bool _DeviceSet1;

void
runTest(uint32_t* d_In, uint32_t* d_Out, uint32_t* h_out, int num_elements);
#ifdef _DEBUG
void bbin_printf(uint32_t elem, int N = 32, int end = 0)
{
	for ( int i = N - 1; i >= end; i-- )
	printf("%d", (bool)(elem & ((uint32_t)1 << i)));
}
#endif

int iDivUp1(int a, int b) // Round a / b to nearest higher integer value

		{
	return (a % b != 0) ? (a / b + 1) : (a / b);
}

/*********************/
/* SVD PLAN CREATION */
/*********************/
template<class T>
void createPlan(plan<T>& plan, unsigned int NperGPU, unsigned int NEdgeperGPU,
		unsigned int N, unsigned int gpuID) {

	// --- Device allocation
	(cudaSetDevice(gpuID));
	checkCudaErrors(cudaStreamCreate(&plan.stream));
	checkCudaErrors(
			cudaMalloc(&(plan.d_data), NperGPU * sample_seeds* sizeof(T)));
	checkCudaErrors(cudaMalloc(&(plan.v), N * sizeof(Point)));
	checkCudaErrors(
			cudaMalloc(&(plan.Fe_Edge), NEdgeperGPU * sample_seeds * sizeof(uint32_t)));
	checkCudaErrors(
			cudaMalloc(&(plan.Fr_Edge), NEdgeperGPU * sample_seeds * sizeof(uint32_t)));

	checkCudaErrors(
			cudaMalloc(&(plan.Fe_Face), NEdgeperGPU * sample_seeds * sizeof(uint32_t)));
	checkCudaErrors(
			cudaMalloc(&(plan.Fr_Face), NEdgeperGPU * sample_seeds * sizeof(uint32_t)));

	checkCudaErrors(
			cudaMalloc(&(plan.eg), NEdgeperGPU * sample_seeds * sizeof(fEdge)));
	checkCudaErrors(
			cudaMalloc(&(plan.fc), NEdgeperGPU * sample_seeds * sizeof(fFace)));
	checkCudaErrors(cudaMalloc(&(plan.b), 1 * sizeof(Boundary)));
	checkCudaErrors(cudaMalloc(&(plan.d), 1 * sizeof(Dimension)));
	checkCudaErrors(cudaMalloc(&(plan.s), 1 * sizeof(Point)));

}

//=======================================================================================

__device__ __host__ void error(int errorcode) {

	int temp = errorcode;
	printf("error!!\n");
}

__device__ __host__ void trilinearInterpolation(float p1[3], int idex,
		Boundary*b, Dimension* d, Point* step, float* m_x1, float* m_y1,
		float* m_z1, int primaryXDIMENSION, float& vx, float&vy, float& vz) {
	int xDim;
	int yDim;
	int zDim;

//	float vvx, vvy, vvz;
	if (d->x > primaryXDIMENSION) {
		int coef = primaryXDIMENSION / d->x;
		xDim = d->x / coef;
		yDim = d->y / coef;
		zDim = d->z / coef;

	} else if (d->x < primaryXDIMENSION) {
		int coef = primaryXDIMENSION / d->x;
		xDim = d->x * coef;
		yDim = d->y * coef;
		zDim = d->z * coef;

	} else {
		xDim = d->x;
		yDim = d->y;
		zDim = d->z;
	}

	float step_x = (b->high.x - b->low.x) / (xDim - 1);
	float step_y = (b->high.y - b->low.y) / (yDim - 1);
	float step_z = (b->high.z - b->low.z) / (zDim - 1);

	int kk = (p1[2] /*+ exp(-4.0)*/- b->low.z) / step_z;
	int jj = (p1[1] /*+ exp(-4.0)*/- b->low.y) / step_y;
	int ii = (p1[0] /*+ exp(-4.0)*/- b->low.x) / step_x;

	float p0[3];
	p0[0] = b->low.x + ii * step_x;
	p0[1] = b->low.y + jj * step_y;
	p0[2] = b->low.z + kk * step_z;

	double xd = (p1[0] - p0[0]) / step_x; //(p1[0]-step_x*int(p1[0]/step_x))/step_x;
	double yd = (p1[1] - p0[1]) / step_y; //(p1[1]-step_y*int(p1[1]/step_y))/step_y;
	double zd = (p1[2] - p0[2]) / step_z; //(p1[2]-step_z*int(p1[2]/step_z))/step_z;
	float v1[3];
	float v2[3];
	float v3[3];
	float v4[3];
	float v5[3];
	float v6[3];
	float v7[3];
	float v8[3];

	if (kk >= zDim - 1 || jj >= yDim - 1 || ii >= xDim - 1) {
		vx = m_x1[((kk * yDim + jj) * xDim + ii)];
		vy = m_y1[((kk * yDim + jj) * xDim + ii)];
		vz = m_z1[((kk * yDim + jj) * xDim + ii)];
		return;

	}
	if (xd == 0 && yd == 0 && zd == 0) {
		vx = m_x1[((kk * yDim + jj) * xDim + ii)];
		vy = m_y1[((kk * yDim + jj) * xDim + ii)];
		vz = m_z1[((kk * yDim + jj) * xDim + ii)];
		return;
	}

	{

		/*	int xDim = XDIMENSION;
		 int yDim = YDIMENSION;
		 int zDim = ZDIMENSION;*/
		v1[0] = m_x1[((kk * yDim + jj) * xDim + ii)];
		v1[1] = m_y1[((kk * yDim + jj) * xDim + ii)];
		v1[2] = m_z1[((kk * yDim + jj) * xDim + ii)];

		v2[0] = m_x1[((kk * yDim + jj) * xDim + ii + 1)];
		v2[1] = m_y1[((kk * yDim + jj) * xDim + ii + 1)];
		v2[2] = m_z1[((kk * yDim + jj) * xDim + ii + 1)];

		v3[0] = m_x1[((kk * yDim + (jj + 1)) * xDim + ii)];
		v3[1] = m_y1[((kk * yDim + (jj + 1)) * xDim + ii)];
		v3[2] = m_z1[((kk * yDim + (jj + 1)) * xDim + ii)];
		//

		v4[0] = m_x1[((kk * yDim + (jj + 1)) * xDim + ii + 1)];
		v4[1] = m_y1[((kk * yDim + (jj + 1)) * xDim + ii + 1)];
		v4[2] = m_z1[((kk * yDim + (jj + 1)) * xDim + ii + 1)];

		//int idx2 = (k + 1)*(xDim*yDim) + (j + 1)*yDim + i + 1;

		v5[0] = m_x1[(((kk + 1) * yDim + jj) * xDim + ii)];
		v5[1] = m_y1[(((kk + 1) * yDim + jj) * xDim + ii)];
		v5[2] = m_z1[(((kk + 1) * yDim + jj) * xDim + ii)];

		v6[0] = m_x1[(((kk + 1) * yDim + jj) * xDim + ii + 1)];
		v6[1] = m_y1[(((kk + 1) * yDim + jj) * xDim + ii + 1)];
		v6[2] = m_z1[(((kk + 1) * yDim + jj) * xDim + ii + 1)];

		v7[0] = m_x1[(((kk + 1) * yDim + (jj + 1)) * xDim + ii)];
		v7[1] = m_y1[(((kk + 1) * yDim + (jj + 1)) * xDim + ii)];
		v7[2] = m_z1[(((kk + 1) * yDim + (jj + 1)) * xDim + ii)];

		v8[0] = m_x1[(((kk + 1) * yDim + (jj + 1)) * xDim + ii + 1)];
		v8[1] = m_y1[(((kk + 1) * yDim + (jj + 1)) * xDim + ii + 1)];
		v8[2] = m_z1[(((kk + 1) * yDim + (jj + 1)) * xDim + ii + 1)];
	}

	double c00 = v1[0] * (1 - xd) + v2[0] * xd;
	double c10 = v3[0] * (1 - xd) + v4[0] * xd;
	double c01 = v5[0] * (1 - xd) + v6[0] * xd;
	double c11 = v7[0] * (1 - xd) + v8[0] * xd;

	double c0 = c00 * (1 - yd) + c10 * yd;
	double c1 = c01 * (1 - yd) + c11 * yd;

	vx = c0 * (1 - zd) + c1 * zd;

	c00 = v1[1] * (1 - xd) + v2[1] * xd;
	c10 = v3[1] * (1 - xd) + v4[1] * xd;
	c01 = v5[1] * (1 - xd) + v6[1] * xd;
	c11 = v7[1] * (1 - xd) + v8[1] * xd;

	c0 = c00 * (1 - yd) + c10 * yd;
	c1 = c01 * (1 - yd) + c11 * yd;

	vy = c0 * (1 - zd) + c1 * zd;

	c00 = v1[2] * (1 - xd) + v2[2] * xd;
	c10 = v3[2] * (1 - xd) + v4[2] * xd;
	c01 = v5[2] * (1 - xd) + v6[2] * xd;
	c11 = v7[2] * (1 - xd) + v8[2] * xd;

	c0 = c00 * (1 - yd) + c10 * yd;
	c1 = c01 * (1 - yd) + c11 * yd;

	vz = c0 * (1 - zd) + c1 * zd;

	int cxz = 0;

//	getLorenzField1(p1, vvx, vvy, vvz);
//	if (vvx != vx || vvy != vy || vvz != vz) {
//		vx = vvx;
//		vy = vvy;
//		vz = vvz;
//	}
	return;

}

__device__ __host__ void generalstreamlineTracing_single(float p[3],
		bool bForward, float e[3], float* m_x1, float* m_y1, float* m_z1,
		Boundary*b, Dimension* d, Point* step, int currentDimX, int tau) {

	int start_pixel_id = 0;

	float i_kk = ((p[2] - b->low.z) / step->z);
	float i_jj = ((p[1] - b->low.y) / step->y);
	float i_ii = ((p[0] - b->low.x) / step->x);

	float vx, vy, vz;
	float next_i, next_j, next_k;
	float ii, jj, kk;
	start_pixel_id = i_kk * (d->x * d->y) + i_jj * d->x + i_ii;

//for (int j = 0;j < ndim*ndim*ndim;j++)

	float p2[3];

	next_i = p[0];					//samples_x[j];
	next_j = p[1];					//samples_y[j];
	next_k = p[2];					//samples_z[j];
	for (int k = 0; k < tau; k++) {

		p2[0] = next_i;
		p2[1] = next_j;
		p2[2] = next_k;
		/*	if (m_bTornadoFieldSelected)
		 trilinearInterpolation2(next_i,next_j,next_k, start_pixel_id, vx, vy, vz);

		 else*/
		trilinearInterpolation(p2, start_pixel_id, b, d, step, m_x1, m_y1, m_z1,
				currentDimX, vx, vy, vz);

		//get_ABC_flow(next_i,next_j,next_k,vx,vy,vz);
		/*	if (WhichType == 0)

		 else if (WhichType == 1)
		 get_ABC_flow(next_i, next_j, next_k, vx, vy, vz);*/
		//	get_Lorenz_Field(next_i, next_j, next_k, vx, vy, vz);
		//getLorenzField1(p2, vx, vy, vz);
		float dist = sqrt(vx * vx + vy * vy + vz * vz);
		if (dist < 1.0e-6 || next_i < b->low.x || next_i > b->high.x
				|| next_j < b->low.y || next_j > b->high.y || next_k < b->low.z
				|| next_k > b->high.z) {
			//cout<<"ddd"<<endl;
			break;
		}

		vx = (vx / dist) * (step->x / 4.0);
		vy = (vy / dist) * (step->y / 4.0);
		vz = (vz / dist) * (step->z / 4.0);

		ii = next_i;
		jj = next_j;
		kk = next_k;

		if (bForward) {
			next_i = ii + vx; //RK4
			next_j = jj + vy;
			next_k = kk + vz;
		} else {
			next_i = ii - vx; //RK4
			next_j = jj - vy;
			next_k = kk - vz;
		}

		i_kk = ((next_k - b->low.z) / step->z);
		i_jj = ((next_j - b->low.y) / step->y);
		i_ii = ((next_i - b->low.x) / step->x);

		if (i_ii >= 0 && i_jj >= 0 && i_kk >= 0) {
			// 			next_i = lowBoundary+((floor(i_ii+0.5))*step) /*+ (0.5*step)*/;
			// 			next_j = lowBoundary+((floor(i_jj+0.5))*step) /*+ (0.5*step)*/;
			// 			next_k = lowBoundary_z+((floor(i_kk+0.5))*step) /*+ (0.5*step)*/;
		} else
			break;

	}

	e[0] = next_i;
	e[1] = next_j;
	e[2] = next_k;

	return;

}

__global__ void Tracing(ASF_vertex * m, float* m_x1, float* m_y1, float* m_z1,
		Dimension* d, Boundary* b, Point* step, bool bForward,
		uint32_t whichData, int currentXDim, uint32_t start_rows,
		uint32_t num_rows, uint32_t level, int currtau)

		{
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;

	if (row > num_rows)
		return;

	ASF_vertex vertex = m[row];

	if (row == 8701)
		vertex = m[row];
	if (!vertex.checkInBoundary(b) || !vertex.checkInBoundary_StartPoint(b)) {
		vertex.unsetInBoundary();
		m[row] = vertex;
		return;
	}

	float p1[3];
	if (vertex.type == 1) {

		p1[0] = vertex.e.x;
		p1[1] = vertex.e.y;
		p1[2] = vertex.e.z;

		float ep[3];
		generalstreamlineTracing_single(p1, bForward, ep, m_x1, m_y1, m_z1, b,
				d, step, currentXDim, 1);

		vertex.e.x = ep[0];
		vertex.e.y = ep[1];
		vertex.e.z = ep[2];
		//
		//			if (Trace[row][currtau + 10].x != vertex.e.x
		//					|| Trace[row][currtau + 10].y != vertex.e.y
		//					|| Trace[row][currtau + 10].z != vertex.e.z)
		//				vertex.e.x = ep[0];
		//Trace[row][currtau] = vertex.e;
		//vertex.es[currtau] = vertex.e;

	} else if (vertex.type == 2) {
		if (!m[vertex.left].checkInBoundary(b)
				|| !m[vertex.right].checkInBoundary(b)
				|| !m[vertex.left].checkInBoundary_StartPoint(b)
				|| !m[vertex.right].checkInBoundary_StartPoint(b)) {
			vertex.unsetInBoundary();
			m[row] = vertex;
			return;
		}

		p1[0] = vertex.p.x;
		p1[1] = vertex.p.y;
		p1[2] = vertex.p.z;
		float ep[3];
		generalstreamlineTracing_single(p1, bForward, ep, m_x1, m_y1, m_z1, b,
				d, step, currentXDim, currtau + 1);

		Point e1 = m[vertex.left].e;
		Point e2 = m[vertex.right].e;
		vertex.e.x = (e1.x + e2.x) / 2;
		vertex.e.y = (e1.y + e2.y) / 2;
		vertex.e.z = (e1.z + e2.z) / 2;
		//Trace[row][currtau] = vertex.e;

		//vertex.es[currtau] = vertex.e;
		//Trace[row][currtau] = vertex.e;
	} else if (vertex.type == 3) {

		if (!m[vertex.left].checkInBoundary(b)
				|| !m[vertex.right].checkInBoundary(b)
				|| !m[vertex.up].checkInBoundary(b)
				|| !m[vertex.down].checkInBoundary(b)
				|| !m[vertex.left].checkInBoundary_StartPoint(b)
				|| !m[vertex.right].checkInBoundary_StartPoint(b)
				|| !m[vertex.up].checkInBoundary_StartPoint(b)
				|| !m[vertex.down].checkInBoundary_StartPoint(b)) {

			vertex.unsetInBoundary();
			m[row] = vertex;
			return;
		}
		Point e1 = m[vertex.left].e;
		Point e2 = m[vertex.right].e;
		Point e3 = m[vertex.up].e;
		Point e4 = m[vertex.down].e;
		vertex.e.x = (e1.x + e2.x + e3.x + e4.x) / 4;
		vertex.e.y = (e1.y + e2.y + e3.y + e4.y) / 4;
		vertex.e.z = (e1.z + e2.z + e3.z + e4.z) / 4;
		//	Trace[row][currtau] = vertex.e;
		//vertex.es[currtau] = vertex.e;
		//Trace[row][currtau] = vertex.e;
	}

	{
		if (vertex.checkInBoundary(b)) {		//}(vertex.e, b, currtau)) {
			int io, jo, ko;
			vertex.getIndex(b, step, d, io, jo, ko);
			uint32_t range = io + jo * d->x + ko * (d->x * d->y);//_tv.getRange(&_b, &step, &_d);

			vertex.setRange(range);
			//				uint32_t range = vertex.getRange(b, step, d);
			//
			//				vertex.setRange(range);

		} else {
			vertex.unsetInBoundary();
		}
	}

	m[row] = vertex;
	return;
}

__device__ __host__ bool checkEdge(ASF_vertex vertex1, ASF_vertex vertex2,
		Boundary*b, Dimension* d, Point* step, bool bForward, int tau, int v1i,
		int v2i) {
//	int range1 = .getRange_tau(&Trace[v1i][tau], b, step, d, tau);
//	int range2 = vertex2.getRange_tau(&Trace[v2i][tau], b, step, d, tau);

	int range1 = vertex1.getRange(b, step, d);
	int range2 = vertex2.getRange(b, step, d);
	/*if (!bForward)
	 {
	 range1 = vertex1.getRangeBackward();
	 range2 = vertex2.getRangeBackward();
	 }*/
	int iz1 = range1 / (d->x * d->y);
	int iy1 = (range1 - iz1 * (d->x * d->y)) / d->x;
	int ix1 = (range1 - iz1 * (d->x * d->y)) % d->x;

	int iz2 = range2 / (d->x * d->y);
	int iy2 = (range2 - iz2 * (d->x * d->y)) / d->x;
	int ix2 = (range2 - iz2 * (d->x * d->y)) % d->x;

//	vertex1.getIndex(b, step, d, ix1, iy1, iz1);
//	vertex2.getIndex(b, step, d, ix2, iy2, iz2);
	float dist = 0.;
	dist = sqrt(
			(float) ((iz2 - iz1) * (iz2 - iz1) + (iy2 - iy1) * (iy2 - iy1)
					+ (ix2 - ix1) * (ix2 - ix1)));

	double distance = vertex1.e.dist(vertex2.e);
//	if (distance >= 2 * step->x)
//		return false;
//	return true;
	/*if (iz1 == d->z || iz1 == d->z - 1 || iy1 == d->y || iy1 == d->y - 1 || ix1 == d->x || ix1 == d->x - 1)
	 printf("");*/

	/*if (iz2 == d->z || iz2 == d->z - 1 || iy2 == d->y || iy2 == d->y - 1 || ix2 == d->x || ix2 == d->x - 1)
	 printf("");*/
	float realdist = sqrt(
			(vertex1.e.x - vertex2.e.x) * (vertex1.e.x - vertex2.e.x)
					+ (vertex1.e.y - vertex2.e.y) * (vertex1.e.y - vertex2.e.y)
					+ (vertex1.e.z - vertex2.e.z)
							* (vertex1.e.z - vertex2.e.z));
	if (realdist > 4)
		printf("");
	if (realdist <= step->x
			|| (abs(iz1 - iz2) <= 1 && abs(iy1 - iy2) <= 1 && ix1 == ix2)
			|| (abs(iy1 - iy2) <= 1 && abs(ix1 - ix2) <= 1 && iz1 == iz2)
			|| (abs(iz1 - iz2) <= 1 && abs(ix1 - ix2) <= 1 && iy1 == iy2)) //
		//if(dist <= 2.0 || abs(range1 - range2) == (2 * d->x + 1) || abs(range1 - range2) == (2 * d->x*d->y + 1) || abs(range1 - range2) == (2 * d->x*d->y + d->x) || abs(range1 - range2) == (2 * d->x*d->y - d->x) || abs(range1 - range2) == (1 * d->x*d->y + d->x) || abs(range1 - range2) == (1 * d->x*d->y - d->x))//&& realdist <= 2.0)//|| (abs(iz1 - iz2) <= 2 && abs(iy2 - iy1) <= 2 && abs(ix2 - ix1) <= 2))

		return true;

	return false;

}

__device__ __host__ bool DivideEdges(ASF_vertex* v1, ASF_vertex* v2,
		uint32_t range1, Boundary* _b, Dimension* _d, Point* step, float*m_x1,
		float* m_y1, float*m_z1, int currentXDim, uint32_t faceNum,
		bool bForward, int divider, int tau, ASF_vertex* otv, int curvertexid,
		int v1i, int v2i) {

	Point p_left_c;
	int n = divider;
	//ASF_vertex* _tv = new ASF_vertex[n];			// = new ASF_vertex[n-1];
	ASF_vertex vcenter;
	vcenter.type = 1;
	double s = 1.0 / n;
	for (int i = 1; i < n; i++) {

		double t = i * s;
		vcenter.p.x = (t) * v1->p.x + (1 - t) * v2->p.x;
		vcenter.p.y = (t) * v1->p.y + (1 - t) * v2->p.y;
		vcenter.p.z = (t) * v1->p.z + (1 - t) * v2->p.z;

		vcenter.e.x = (t) * v1->e.x + (1 - t) * v2->e.x;
		vcenter.e.y = (t) * v1->e.y + (1 - t) * v2->e.y;
		vcenter.e.z = (t) * v1->e.z + (1 - t) * v2->e.z;

		//if (vcenter.type == 2)
		{
//			for (int j = 0; j < tau + 1; j++) {
//				Trace[curvertexid][j].x = (t * Trace[v1i][j].x)
//						+ (1 - t) * Trace[v2i][j].x;
//				Trace[curvertexid][j].y = (t * Trace[v1i][j].y)
//						+ (1 - t) * Trace[v2i][j].y;
//				Trace[curvertexid][j].z = (t * Trace[v1i][j].z)
//						+ (1 - t) * Trace[v2i][j].z;
//
//			}
		}
		//	vcenter.e = Trace[curvertexid][tau];

		//Trace[curvertexid][tau] = vcenter.e;
		//		AdvectParticle_estimated(&vcenter, m_x1, m_y1, m_z1, _b, _d, step,
		//				currentXDim, tau + 1, bForward);

		if (!vcenter.checkInBoundary(_b)) {
			*otv = vcenter;
			//return false;
			error(0);
			return false;
		}
		uint32_t oldrange = v1->getOldRange();// _tv.getOldRange(&_b, &step, &_d);
		vcenter.setOldRange(oldrange);
		vcenter.setInBoundary();

		int io, jo, ko;
		io = (vcenter.e.x - _b->low.x) / (step->x);
		io = (vcenter.e.y - _b->low.y) / (step->y);
		io = (vcenter.e.z - _b->low.z) / (step->z);
		vcenter.getIndex(_b, step, _d, io, jo, ko);
		uint32_t range = io + jo * _d->x + ko * (_d->x * _d->y);//_tv.getRange(&_b, &step, &_d);

		//vcenter.setRange(range);
		vcenter.range = range;
		//printf("vcenter: %d range:%d",vcenter.range,range);
		//if(v1->getRange() == 4430 )
		if (checkEdge(vcenter, *v1, _b, _d, step, bForward, tau, v1i, v2i)
				&& checkEdge(vcenter, *v2, _b, _d, step, bForward, tau, v1i,
						v2i)) {
			*otv = vcenter;

			return true;
		}

	}

	if (n == 5) {
		*otv = vcenter;
	}
	*otv = vcenter;
	return false;

}

__device__ __host__ void Split_One_Edge(ASF_vertex*m, fFace* fc, fEdge* eg,
		uint32_t* Fe_Edge, uint32_t* Fr_Edge, Dimension* d, Boundary* b,
		Point* step, float*m_x1, float* m_y1, float*m_z1, int currentXDim,
		int tau, bool bForward, uint32_t num_vertex, uint32_t num_edges,
		uint32_t row) {

	int curVertexId = num_vertex + Fr_Edge[row];
	int curEdgeId = num_edges + Fr_Edge[row] * 2;

	fEdge edge1 = eg[row];
	fEdge edge2 = eg[row];

	int v1i = edge1.v1;
	int v2i = edge1.v2;

	if (edge1.bsplit)
		error(2);

	ASF_vertex vertex1 = m[edge1.v1];
	ASF_vertex vertex2 = m[edge1.v2];

//		if(curVertexId == 5846 )
//				error(1);

	if (!vertex1.checkInBoundary_StartPoint(b)
			|| !vertex2.checkInBoundary_StartPoint(b)) {
		eg[row].unsetInBoundary();
		return;
	}

	if ((vertex1.p.x - vertex2.p.x) > step->x
			|| (vertex1.p.y - vertex2.p.y) > step->y
			|| (vertex1.p.z - vertex2.p.z) > step->z)
		error(0);
	if (curVertexId == 61996)
		Fe_Edge[row] = 1;

	ASF_vertex allseeds[5];
	ASF_vertex vedge;

	bool bFound = false;
	int n = 2;
	//for (n = 2; n <= 5; n++)
	{
		//int n = 10;
		//ASF_vertex* vedge_array = new ASF_vertex[n];
		DivideEdges(&vertex1, &vertex2, vertex1.getOldRange(), b, d, step, m_x1,
				m_y1, m_z1, currentXDim, 0, bForward, n, tau, &vedge,
				curVertexId, v1i, v2i);
		{
			bFound = true;
			//allseeds[n - 2] = vedge;
			//			if (tau == 50)
			//				error(1);
//			ASF_vertex* tempp = new ASF_vertex[3];
//			tempp[0] = vertex1;
//			tempp[1] = vertex2;
//			tempp[2] = vedge;
//			Point temppoint[3];
//			temppoint[0] = vertex1.p;
//			temppoint[1] = vertex2.p;
//			temppoint[2] = vedge.p;
//
//			int ii = (vertex1.e.x - b->low.x) / step->x;
//			int jj = (vertex1.e.y - b->low.y) / step->y;
//			int kk = (vertex1.e.z - b->low.z) / step->z;
//
//			int ii1 = (vertex2.e.x - b->low.x) / step->x;
//			int jj1 = (vertex2.e.y - b->low.y) / step->y;
//			int kk1 = (vertex2.e.z - b->low.z) / step->z;
//			int intarray[3];
////			intarray[0] = vertex1.getRange_tau(&Trace[v1i][tau], b, step, d,
////					tau);
////			intarray[1] = vertex2.getRange_tau(&Trace[v2i][tau], b, step, d,
////					tau);
////			intarray[2] = vedge.getRange_tau(&Trace[curVertexId][tau], b, step,
////					d, tau);
//
//			intarray[0] = vertex1.getRange(b, step, d);
//			intarray[1] = vertex2.getRange(b, step, d);
//			intarray[2] = vedge.getRange(b, step, d);
//			//
//			m[curVertexId] = vedge;
//			if (false && tau == 90) {
//
//				ASF_vertex tempm[3];
//				tempm[0] = vertex1;
//				tempm[1] = vertex2;
//				tempm[2] = vedge;
//				int indexarray[3];
//				indexarray[0] = v1i;
//				indexarray[1] = v2i;
//				indexarray[2] = curVertexId;
////				fs->Save_Streamlines_estimated(m, Trace, indexarray, bForward,
////						"st1", m_x1, m_y1, m_z1, b, d, step, currentXDim,
////						tau + 1, 3);
//
//				fs->Save_Streamlines_EndAdvection(m, indexarray, bForward,
//						"st1", m_x1, m_y1, m_z1, b, d, step, currentXDim,
//						tau + 1, 3);
//				fs->save_Quad_FaceWithSeedPoits(m, fc, tempp, eg, Trace, *d,
//						"face85", num_vertex, 1, 3, tau, 1);
//				if (tau == 99)
//					fs->save_Voxel(m, fc, intarray, eg, *d, "face85",
//							num_vertex, 3, tau, 1, 10);
//			}

			//break;
		}

		//allseeds[n - 2] = vedge;
	}

	//	string dataname = "edge";
	//	fs->start_save_Quad_One_Face(m, fc[edge1.E2F[0]], eg, *d, dataname,
	//			num_vertex, 1, tau, 1);
	//
	//	generate_streamlines_sparsely(allseeds, bForward, 1, m_x1, m_y1, m_z1, b, d,
	//			step, currentXDim, tau + 2, n - 1);
	//
	//	fs->save_Quad_FaceWithSeedPoits(m, fc, allseeds, eg, *d, dataname,
	//			num_vertex, 1, 1, tau, 1);
	//
	//	fs->start_save_Quad_FaceWithSeedPoits(m, fc, allseeds, eg, *d, dataname,
	//				num_vertex, 1, 1, tau, 1);

	//	if (!bFound) {
	//		DivideEdges(&vertex1, &vertex2, vertex1.getOldRange(), b, d, step, m_x1,
	//				m_y1, m_z1, currentXDim, 0, bForward, 2, tau, &vedge);
	//
	//	}
	vedge.type = 2;
	vedge.oldrange = vertex1.oldrange;
	vedge.left = edge1.v1;
	vedge.right = edge1.v2;

	edge1.v2 = curVertexId;
	edge2.v1 = curVertexId;
	//edge1.bsplit = true;
	edge1.next = curEdgeId;
	edge1.level = edge1.level + 1;
	edge2.level = edge1.level;
	edge2.Prev = row;
	vedge.level = edge2.level;

	//vedge.fxy = 0;
	//vedge.fyz = 0;
	//vedge.fxz = 0;

	float rgb[3];
	rgb[0] = 0.;
	rgb[1] = 0.;
	rgb[2] = 0.;

	eg[row].bsplit = true;
	eg[row].subedge[0] = curEdgeId;
	eg[row].subedge[1] = curEdgeId + 1;

	edge1.bsplit = false;
	edge2.bsplit = false;
	edge1.parent = row;
	edge2.parent = row;

	float dist1 = vertex1.e.dist(vertex2.e);
	float dist1_1 = vertex1.e.dist(vedge.e);
	float dist1_2 = vertex2.e.dist(vedge.e);

	if (dist1_2 >= dist1 || dist1_2 >= dist1) {
		error(1);
		ASF_vertex oa[3];
		oa[0] = vertex1;
		oa[1] = vertex2;
		oa[2] = vedge;

//		generate_streamlines_sparsely(oa, bForward, 1, m_x1, m_y1, m_z1, b, d,
//				step, currentXDim, 100 + 2, 3);
//
//		fs->save_Quad_FaceWithSeedPoits(m, fc, oa, eg, Trace, *d, "face85",
//				num_vertex, 1, 3, tau, 1);
		//error(1);
	}

	eg[curEdgeId] = edge1;

	eg[curEdgeId + 1] = edge2;

	m[curVertexId] = vedge;

}

__global__ void CheckNeighborhood(ASF_vertex*m, fEdge* eg, uint32_t*Fr,
		Dimension* d, Boundary* b, Point* step, bool bForward,
		uint32_t original_num_rows, uint32_t num_rows, uint32_t level,
		int tau) {
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;

//uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
//	for (uint32_t row = 0; row < num_rows; row++)

	if (row >= num_rows)
		return;
	Fr[row] = 1;
	return;
	int co = 0;
	if (row == 11705)
		co = 0;

	fEdge edge = eg[row];
	ASF_vertex vertex1 = m[edge.v1];
	ASF_vertex vertex2 = m[edge.v2];
	//printf("---- %f%f%f \n", vertex1.p.x, vertex1.p.y, vertex1.p.z);

//	if (vertex1.e.dist(vertex2.e) > step->x) {
//		co = 0;
//	}
	if (!edge.isInBoundary())
		return;
	if (edge.bsplit)
		return;
	//	edge.bsplit = false;

	if (!vertex1.isInBoundary() || !vertex2.isInBoundary()) {
		//eg[row].unsetInBoundary();
		edge.unsetInBoundary();
		eg[row] = edge;
		return;
	}

	if (!vertex1.checkInBoundary_StartPoint(b)
			|| !vertex2.checkInBoundary_StartPoint(b)) {
		//	printf("%f%f%f \n", vertex1.p.x, vertex1.p.y, vertex1.p.z);
		//eg[row].unsetInBoundary();
		edge.unsetInBoundary();
		eg[row] = edge;

		return;
	}

	if (!vertex1.checkInBoundary(b) || !vertex2.checkInBoundary(b)) {
		edge.unsetInBoundary();
		eg[row] = edge;

		return;

	}

	if (abs((int) vertex1.getOldRange() - (int) vertex2.getOldRange()) > 1
			&& abs((int) vertex1.getOldRange() - (int) vertex2.getOldRange())
					!= d->x
			&& abs((int) vertex1.getOldRange() - (int) vertex2.getOldRange())
					!= (d->x * d->y))
		error(0);

	/*	if (vertex1.getOldRange() % (d->x - 1) == 0 && (vertex2.getOldRange() % d->x) == 0)
	 continue;
	 if (vertex1.getOldRange() % (d->x*d->y - 1 )==0 && (vertex2.getOldRange() % d->x*d->y) == 0)
	 continue;*/

	float dist = sqrt(
			(vertex1.e.x - vertex2.e.x) * (vertex1.e.x - vertex2.e.x)
					+ (vertex1.e.y - vertex2.e.y) * (vertex1.e.y - vertex2.e.y)
					+ (vertex1.e.z - vertex2.e.z)
							* (vertex1.e.z - vertex2.e.z));

	if (fabs(vertex1.p.x - vertex2.p.x) > step->x + exp(-6.0)
			|| fabs(vertex1.p.y - vertex2.p.y) > step->y + exp(-6.0)
			|| fabs(vertex1.p.z - vertex2.p.z) > step->z + exp(-6.0))
		error(1);
	if ((row == 493 || row == 518 || row == 685 || row == 494))
		Fr[row] = 0;
	//if (dist >  step->x)// || !checkEdge(vertex1, vertex2, d, bForward))

	if (!checkEdge(vertex1, vertex2, b, d, step, bForward, tau, edge.v1,
			edge.v2)) {
		//vertex1.setInNextLevel_xy();
		//	edge.bsplit = true;

		Fr[row] = 1;
	}

	//			continue;

	eg[row] = edge;

}

__global__ void CheckFace(ASF_vertex*m, fEdge* eg, fFace* fc, uint32_t*Fe_Edge,
		uint32_t* Fe_Face, Boundary* b, Dimension* d, Point* step,
		bool bForward, uint32_t num_face, int tau) {

	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;

	if (row >= num_face)
		return;
	fFace face = fc[row];
	if (row == 442)
		face = fc[row];

	if (face.bsplit) {
		face = fc[row];
		return;
	}

	if (!face.isInBoundary()) {
		return;
	}

	if (!m[face.vertexes[0]].isInBoundary()
			|| !m[face.vertexes[1]].isInBoundary()
			|| !m[face.vertexes[2]].isInBoundary()
			|| !m[face.vertexes[3]].isInBoundary()) {
		fc[row].unsetInBoundary();
		return;

	}

	if (!m[face.vertexes[0]].checkInBoundary_StartPoint(b)
			|| !m[face.vertexes[1]].checkInBoundary_StartPoint(b)
			|| !m[face.vertexes[2]].checkInBoundary_StartPoint(b)
			|| !m[face.vertexes[3]].checkInBoundary_StartPoint(b)) {
		fc[row].unsetInBoundary();
		return;

	}

	int v1i = face.vertexes[0];
	int v2i = face.vertexes[1];
	int v3i = face.vertexes[2];
	int v4i = face.vertexes[3];

	if (!m[face.vertexes[0]].checkInBoundary(b)
			|| !m[face.vertexes[1]].checkInBoundary(b)
			|| !m[face.vertexes[2]].checkInBoundary(b)
			|| !m[face.vertexes[3]].checkInBoundary(b)) {
		fc[row].unsetInBoundary();
		return;

	}

	ASF_vertex v1 = m[face.vertexes[0]];
	ASF_vertex v2 = m[face.vertexes[1]];
	ASF_vertex v3 = m[face.vertexes[2]];
	ASF_vertex v4 = m[face.vertexes[3]];

	if (!checkEdge(v1, v2, b, d, step, bForward, tau, v1i, v2i)
			|| !checkEdge(v2, v3, b, d, step, bForward, tau, v2i, v3i)
			|| !checkEdge(v3, v4, b, d, step, bForward, tau, v3i, v4i)
			|| !checkEdge(v1, v4, b, d, step, bForward, tau, v1i, v4i))

		Fe_Face[row] = 1;

	fc[row] = face;

}

__global__ void EdgeReduction(fEdge*eg, ASF_vertex*m, fFace* Ff,
		uint32_t* Fe_Face, uint32_t* Fe_Edge, uint32_t num_edges) {

	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;

	if (row >= num_edges)
		return;
}

__global__ void SplitEdge(ASF_vertex*m, fFace* fc, fEdge* eg, uint32_t* Fe_Edge,
		uint32_t* Fr_Edge, Dimension* d, Boundary* b, Point* step, float*m_x1,
		float* m_y1, float*m_z1, int currentXDim, bool bForward,
		uint32_t num_vertex, uint32_t num_edges, uint32_t level, int tau) {

	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
//	for (uint32_t row = 0; row < num_rows; row++)

	if (row > num_edges)
		return;

	if (Fe_Edge[row] == 1) {

//			if (row == 3272)
//				error(0);
//				Split_One_Edge(m, fc, eg, Fe_Edge, Fr_Edge, d, b, step, m_x1, m_y1,
//						m_z1, currentXDim, tau, bForward, num_vertex, num_edges,
//						row);

		Split_One_Edge(m, fc, eg, Fe_Edge, Fr_Edge, d, b, step, m_x1, m_y1,
				m_z1, currentXDim, tau, bForward, num_vertex, num_edges, row);
	}

}

pair<uint32_t, float> OBF_Decomposition2(ASF_vertex* _a, fFace *fc, fEdge*eg,
		Boundary* _b, Dimension* _d, Point* step, int _tau, float* m_x1,
		float* m_y1, float* m_z1, ASF_vertex** oVertex, Edge ** oFc,
		uint32_t ** oFr, fFace** oFace, fEdge** oEdge, uint32_t * oRSize,
		uint32_t * oFaceSize, uint32_t * oEdgeSize, uint32_t whichData,
		int CSize, uint32_t _curEdgeSize, uint32_t _curFaceSize,
		int curXDimension, int sampleSeeds, int MorseLevel, bool bForward) {

//-----------CPU initialization---------------------------->

	uint32_t * To, *From;
	uint32_t * To2;

//	Trace = new Point*[CSize * sampleSeeds];

//	for (int j = 0; j < CSize * sampleSeeds; j++) {
//		//if(_a[j].isInBoundary())
//		Trace[j] = new Point[_tau + 1];
//	}

	printf("Vertices: %u\n", CSize - 2);

	uint32_t currTau = _tau;

	uint32_t currSize = CSize;
	uint32_t icycle = currTau;

	uint32_t curEdgeSize = _curEdgeSize;
	uint32_t curFaceSize = _curFaceSize;

	uint32_t* Fe_Edge = new uint32_t[curEdgeSize * sampleSeeds];
	uint32_t* Fr_Edge = new uint32_t[curEdgeSize * sampleSeeds];

	uint32_t* Fe_Edge2 = new uint32_t[curEdgeSize * sampleSeeds];
	uint32_t* Fr_Edge2 = new uint32_t[curEdgeSize * sampleSeeds];

	uint32_t* Fe_Face = new uint32_t[curFaceSize * sampleSeeds];
	uint32_t* Fr_Face = new uint32_t[curFaceSize * sampleSeeds];

	ASF_vertex* currA = new ASF_vertex[currSize * sampleSeeds];

	memcpy(currA, _a, currSize * sizeof(ASF_vertex));

	uint32_t oldsize = CSize;
	int level = 0;

	memset(Fe_Edge, 0, curEdgeSize * sampleSeeds * sizeof(uint32_t));
	memset(Fr_Edge, 0, curEdgeSize * sampleSeeds * sizeof(uint32_t));

	memset(Fe_Edge2, 0, curEdgeSize * sampleSeeds * sizeof(uint32_t));
	memset(Fr_Edge2, 0, curEdgeSize * sampleSeeds * sizeof(uint32_t));

	memset(Fe_Face, 0, curFaceSize * sampleSeeds * sizeof(uint32_t));
	memset(Fr_Face, 0, curFaceSize * sampleSeeds * sizeof(uint32_t));

	int interval = 1;
	int oldFaceSize = curFaceSize;
	int co = 0;

	int c = currSize;
	int nf = curFaceSize;

//================================================================

//-----------GPU initialization---------------------------->

	if (!_DeviceSet1) {
		_DeviceSet1 = true;
		checkCudaErrors(cudaSetDevice(2));
	}
	Dimension* d_d;
	Boundary* d_b;
	Point* d_step;
	float* d_vx;
	float* d_vy;
	float* d_vz;
	ASF_vertex* d_currA;
	uint32_t* d_Fe_Edge;
	uint32_t* d_Fe_Face;
	uint32_t* d_Fr_Edge;
	uint32_t* d_Fr_Face;
	fEdge* d_eg;
	fFace* d_fc;

	cudaError_t e1, e2, e3, e4, e5, e6, e7;
	checkCudaErrors(
			e1 = cudaMalloc((void** )&d_fc,
					curFaceSize * sampleSeeds * sizeof(fFace)));
	checkCudaErrors(
			e2 = cudaMalloc((void** )&d_eg,
					curEdgeSize * sampleSeeds * sizeof(fEdge)));

	checkCudaErrors(
			e5 = cudaMalloc((void** )&d_currA,
					currSize * sampleSeeds * sizeof(ASF_vertex)));
	checkCudaErrors(e3 = cudaMalloc((void** )&d_d, sizeof(Dimension)));
	checkCudaErrors(e4 = cudaMalloc((void** )&d_b, sizeof(Boundary)));
	checkCudaErrors(e4 = cudaMalloc((void** )&d_step, sizeof(Point)));

	checkCudaErrors(
			e6 = cudaMalloc((void** )&d_Fe_Edge,
					curEdgeSize * sampleSeeds * sizeof(uint32_t)));
	checkCudaErrors(
			e6 = cudaMalloc((void** )&d_Fr_Edge,
					curEdgeSize * sampleSeeds * sizeof(uint32_t)));

	checkCudaErrors(
			e7 = cudaMalloc((void** )&d_Fe_Face,
					curFaceSize * sampleSeeds * sizeof(uint32_t)));
	checkCudaErrors(
			e7 = cudaMalloc((void** )&d_Fr_Face,
					curFaceSize * sampleSeeds * sizeof(uint32_t)));

	checkCudaErrors(
			e6 = cudaMalloc((void** )&d_vx,
					currSize  * sizeof(float)));
	checkCudaErrors(
			e6 = cudaMalloc((void** )&d_vy,
					curEdgeSize  * sizeof(float)));

	checkCudaErrors(
			e7 = cudaMalloc((void** )&d_vz,
					curFaceSize * sizeof(float)));


	if (e1 == cudaErrorMemoryAllocation || e2 == cudaErrorMemoryAllocation
			|| e3 == cudaErrorMemoryAllocation
			|| e4 == cudaErrorMemoryAllocation
			|| e5 == cudaErrorMemoryAllocation
			|| e6 == cudaErrorMemoryAllocation
			|| e7 == cudaErrorMemoryAllocation) {
		throw "Error: Not enough memory on GPU\n";
	}
//===============================================================
//	(cudaMallocHost(&inputMatrices, N * sample_seeds * sizeof(ASF_vertex)));
//	fEdge *inputEdges;
//	(cudaMallocHost(&inputEdges, NEdgeperGPU * sample_seeds * sizeof(fEdge)));
//	uint32_t *inputFe_Edge;
//	(cudaMallocHost(&inputFe_Edge,
//			NEdgeperGPU * sample_seeds * sizeof(uint32_t)));
//	uint32_t *inputFr_Edge;
//	(cudaMallocHost(&inputFr_Edge,
//			NEdgeperGPU * sample_seeds * sizeof(uint32_t)));
//
//	uint32_t *inputFe_Face;
//	(cudaMallocHost(&inputFe_Face,
//			NEdgeperGPU * sample_seeds * sizeof(uint32_t)));
//	uint32_t *inputFr_Face;
//	(cudaMallocHost(&inputFr_Face,
//			NEdgeperGPU * sample_seeds * sizeof(uint32_t)));

//================================================================

	dim3 grid(((CSize + 510) / 512), 1, 1); // (RSize-2) valid vertecis the first one is undefined
	dim3 threads(512, 1, 1);
	dim3 grid1(((CSize + 255) / 256), 1, 1);
	dim3 threads1(blockSize, 1, 1);
	dim3 grid2(1, 1, 1);
	dim3 threads2(32, 1, 1);

//===========================================MGPU===============================================
	int GPU_N;
	checkCudaErrors(cudaGetDeviceCount(&GPU_N));
	//cudaError_t e1, e2, e3, e4, e5, e6, e7;

	GPU_N = 1;
	const int numGPUs = 1;
	const int NperGPU = CSize / numGPUs;
	const int NEdgeperGPU = (NperGPU * 3) / numGPUs;
	const int NFaceperGPU = (NperGPU * 3) / numGPUs;
	int tempEdgeNum = 0;
	const int N = NperGPU * numGPUs;
	uint nextpower = 1;

//	plan<ASF_vertex> plan[numGPUs];
//	for (int k = 0; k < numGPUs; k++)
//		createPlan(plan[k], NperGPU, NEdgeperGPU, N * sample_seeds, k);
// --- "Breadth-first" approach - async
//	ASF_vertex *inputMatrices;
//	(cudaMallocHost(&inputMatrices, N * sample_seeds * sizeof(ASF_vertex)));
//	fEdge *inputEdges;
//	(cudaMallocHost(&inputEdges, NEdgeperGPU * sample_seeds * sizeof(fEdge)));
//	uint32_t *inputFe_Edge;
//	(cudaMallocHost(&inputFe_Edge,
//			NEdgeperGPU * sample_seeds * sizeof(uint32_t)));
//	uint32_t *inputFr_Edge;
//	(cudaMallocHost(&inputFr_Edge,
//			NEdgeperGPU * sample_seeds * sizeof(uint32_t)));
//
//	uint32_t *inputFe_Face;
//	(cudaMallocHost(&inputFe_Face,
//			NEdgeperGPU * sample_seeds * sizeof(uint32_t)));
//	uint32_t *inputFr_Face;
//	(cudaMallocHost(&inputFr_Face,
//			NEdgeperGPU * sample_seeds * sizeof(uint32_t)));

//uint32_t* d_Fi;
//	uint32_t* d_Fr;

	//ASF_vertex* h_Out = new ASF_vertex[N * sample_seeds];

	(cudaMemcpyAsync(d_currA, _a, currSize * sizeof(ASF_vertex),
			cudaMemcpyHostToDevice));
	(cudaMemcpyAsync(d_fc, fc, curFaceSize * sizeof(fFace),
			cudaMemcpyHostToDevice));
	(cudaMemcpyAsync(d_vx, m_x1, CSize * sizeof(float), cudaMemcpyHostToDevice));
	(cudaMemcpyAsync(d_vy, m_y1, CSize * sizeof(float), cudaMemcpyHostToDevice));
	(cudaMemcpyAsync(d_vz, m_z1, CSize * sizeof(float), cudaMemcpyHostToDevice));
	(cudaMemcpyAsync(d_eg, eg, curEdgeSize * sizeof(fEdge),
			cudaMemcpyHostToDevice));
	(cudaMemcpyAsync(d_step, step, sizeof(Point), cudaMemcpyHostToDevice));
	(cudaMemcpyAsync(d_b, _b, sizeof(Boundary), cudaMemcpyHostToDevice));
	(cudaMemcpyAsync(d_d, _d, sizeof(Dimension), cudaMemcpyHostToDevice));

	cudaDeviceSynchronize();

	for (int ii = 0; ii < icycle; ii++) {

		Tracing<< <iDivUp1(currSize, BLOCKSIZE), BLOCKSIZE >> >(d_currA, d_vx, d_vy, d_vz, d_d, d_b, d_step, bForward, whichData,
				curXDimension, 0, currSize, level, ii);
		/*CheckNeighborhood_c3 << <iDivUp1(NEdgeperGPU, BLOCKSIZE), BLOCKSIZE >> >(plan[k].d_data, plan[k].eg, plan[k].Fe_Edge, plan[k].d, plan[k].b, plan[k].s, bForward, currSize);
		 if (k < numGPUs)
		 runTest(plan[k].Fe_Edge, plan[k].Fr_Edge, inputFr_Edge, NEdgeperGPU);*/

//	cudaDeviceSynchronize();
		CheckNeighborhood<< <iDivUp1(curEdgeSize, BLOCKSIZE), BLOCKSIZE >> >(d_currA, d_eg, d_Fe_Edge, d_d, d_b, d_step, bForward, currSize,
				curEdgeSize, level, ii);

		runTest(d_Fe_Edge, d_Fr_Edge, Fr_Edge,
				curEdgeSize);

		memset(Fe_Face, 0, curFaceSize * sizeof(uint32_t));
		memset(Fr_Face, 0, curFaceSize * sizeof(uint32_t));
		CheckFace<< <iDivUp1(curFaceSize, BLOCKSIZE), BLOCKSIZE >> >(currA, eg, fc, Fe_Edge, Fe_Face, _b, _d, step,
				bForward, curFaceSize, ii);

		runTest(d_Fe_Face, d_Fr_Face, Fr_Face,
				curFaceSize);
		printf(" Faces = %d \n", Fr_Face[curFaceSize]);

		//======================================================================

		//while (Fr_Face[curFaceSize] > 0) {
		while (Fr_Edge[curEdgeSize] > 0) {
			//ct++;
			//printf("Edges = %d \n", Fr_Edge[curEdgeSize]);
			printf("Edges = %d \n", Fr_Face[curFaceSize]);
			int oldEdgeSize = curEdgeSize;

			SplitEdge << <iDivUp1(curEdgeSize, BLOCKSIZE), BLOCKSIZE >> >(d_currA, d_fc, d_eg, d_Fe_Edge, d_Fr_Edge, d_d, d_b, d_step, d_vx,
					d_vy, d_vz, curXDimension, bForward, currSize, curEdgeSize,
					level, ii);
			cudaDeviceSynchronize();

			currSize = currSize + Fr_Edge[curEdgeSize];
			curEdgeSize = curEdgeSize + Fr_Edge[curEdgeSize] * 2;
//
//			memset(Fr_Edge2, 0, (oldEdgeSize) * sizeof(uint32_t));
//			memset(Fe_Edge2, 0, (oldEdgeSize) * sizeof(uint32_t));
//
//			CheckRemainingEdge(currA, fc, eg, Fe_Edge, Fe_Face, Fe_Edge2, _d,
//					_b, step, bForward, curFaceSize);
//
//			Fr_Edge2[0] = 0;
//			for (int i = 0; i <= oldEdgeSize; i++) {
//				Fr_Edge2[i + 1] = Fr_Edge2[i] + Fe_Edge2[i];
//			}
//
//			printf(" Edges = %d \n", Fr_Edge2[oldEdgeSize]);
//
//			EdgeReduction << <iDivUp1(curEdgeSize, BLOCKSIZE), BLOCKSIZE >> >(plan[k].eg, plan[k].d_data, plan[k].fc, plan[k].Fe_Face, plan[k].Fe_Edge, curEdgeSize);
//			cudaDeviceSynchronize();

			memset(Fr_Edge, 0, curEdgeSize*sizeof(uint32_t));
			cudaMemset(d_Fr_Edge, 0, curEdgeSize*sizeof(uint32_t));

			CheckNeighborhood<< <iDivUp1(curEdgeSize, BLOCKSIZE), BLOCKSIZE >> >(d_currA, d_eg, d_Fe_Edge, d_d, d_b, d_step, bForward, currSize,
					curEdgeSize, level, ii);

			runTest(d_Fe_Edge, d_Fr_Edge, Fr_Edge,
					curEdgeSize);

			//printf("%d \n", inputFr_Edge[curEdgeSize - 1]);

//			SplitEdge << <iDivUp1(curEdgeSize, BLOCKSIZE), BLOCKSIZE >> >(plan[k].eg, plan[k].d_data, plan[k].fc, plan[k].Fe_Edge, plan[k].Fr_Edge, plan[k].d, plan[k].b, plan[k].s, bForward, currSize, curEdgeSize);
//			cudaDeviceSynchronize();
//
//			currSize = currSize + inputFr_Edge[curEdgeSize - 1] + 1;
//
//			curEdgeSize = curEdgeSize + inputFr_Edge[curEdgeSize - 1] + 1;
//
//			CheckNeighborhood << <iDivUp1(curEdgeSize, BLOCKSIZE), BLOCKSIZE >> >(plan[k].eg, plan[k].d_data, plan[k].Fe_Edge, plan[k].d, plan[k].b, plan[k].s, i, bForward, curEdgeSize);
//			cudaMemset(plan[k].Fr_Edge, 0, curEdgeSize*sizeof(uint32_t));
//			//cudaMemset(inputFr_Edge, 0, curEdgeSize*sizeof(uint32_t));
//			memset(inputFr_Edge, 0, curEdgeSize*sizeof(uint32_t));
//
//			runTest(plan[k].Fe_Edge, plan[k].Fr_Edge, inputFr_Edge, curEdgeSize);

		}

		cudaDeviceSynchronize();

		(cudaMemcpyAsync(currA, d_currA, currSize * sizeof(ASF_vertex),
						cudaMemcpyDeviceToHost));
	}

	printf("current size = %d \n", currSize);

	printf("current size = %d \n", currSize);

	Edge* Fc = new Edge[curEdgeSize];

	memset(Fe_Edge, 0, currSize * sizeof(uint32_t));

	for (int i = 0; i < currSize; i++) {
		ASF_vertex vertex = currA[i];
		if (vertex.checkInBoundary(_b)) {

			Fe_Edge[vertex.getOldRange()]++;

		}
	}

//CheckRangeSetKernel << <grid_2, threads >> >(d_m, d_To, currSize - 1);
	To = new uint32_t[CSize];
	memset(Fr_Edge, 0, CSize * sizeof(uint32_t));
	memset(To, 0, CSize * sizeof(uint32_t));

	uint32_t * Fi = new uint32_t[currSize];
	Fr_Edge[0] = 0;
	for (int i = 1; i <= CSize; i++)
		Fr_Edge[i] = Fr_Edge[i - 1] + Fe_Edge[i];

	memset(Fi, 0, currSize * sizeof(uint32_t));
//checkCudaErrors(cudaMemcpy(d_Fr, Fr, currSize * sizeof(uint32_t), cudaMemcpyHostToDevice));
	for (int row = 0; row < currSize; row++) {
		ASF_vertex vertex = currA[row];
		if (vertex.checkInBoundary(_b)
				&& vertex.getOldRange() > 0/*&& vertex.getOldRange() == row*/) {
			uint32_t i = Fr_Edge[vertex.getOldRange() - 1]
					+ Fi[vertex.getOldRange()];
			//	if (bForward)
			Fc[i].setValue(vertex.getRange());
			/*else
			 Fc[i].setValue(vertex.getRangeBackward());
			 */
			//Fc[i].setValidBit();
			Fi[vertex.getOldRange()]++;
		}

	}

	float rgb[3];
	rgb[0] = 0.8;
	rgb[1] = 0.0;
	rgb[2] = 0.0;

//display_voxel(1551, rgb, _d->x, _d->y, _d->z);

//////////////////////////////////////////////////////

	rgb[0] = 0.0;
	rgb[1] = 0.0;
	rgb[2] = 0.8;

	float p[3], ep[3];
//glNewList(index_StreamLine_Lorenz, GL_COMPILE);
	{
		for (int row = 0; row < currSize; row++) {
			uint32_t _index = 0;
			if (currA[row].getOldRange() == 1551) {

				p[0] = currA[row].p.x;
				p[1] = currA[row].p.y;
				p[2] = currA[row].p.z;
				_index = currA[row].getRange();

				//	display_voxel(_index,rgb , _d->x, _d->y, _d->z);
				//	drawPoint(p, rgb);
				//	generalstreamlineTracing_single(p, bForward, ep, false);

			}
		}
	}
//	glEndList();
	uint32_t r = currA[14731].getRange();

	/*if (!bForward)
	 r = currA[14731].getRangeBackward();*/

//currA = 0;
//delete[] _aa;
	*oVertex = currA;
//	*oFace = fc;
//	*oEdge = eg;
	*oFc = Fc;
	*oFr = Fr_Edge;
	*oFaceSize = curFaceSize;
	*oEdgeSize = curEdgeSize;
	*oRSize = currSize;			// Fr_Edge[CSize - 1];

}

