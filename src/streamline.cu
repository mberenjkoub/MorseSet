

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

#include "streamline.h"
//#include "streamline_kernel.h"
#include "graph_generator.h"
#include "scan_common.h"
//#include "parallel_fwd.h"
//#include "hash_table.h"
# define M_PI           3.14159265358979323846  /* pi */
#define OUTBOUNDARY 0xFFFF
void
runTest(uint32_t* d_In, uint32_t* d_Out, uint32_t* h_out, int num_elements);
//void
//runTest(uint32_t* d_in, uint32_t* d_out, uint32_t* h_out, int num_elements);
__host__ __device__ uint32_t getDimRange(Dimension d, Dimension i){ return (i.x + i.y*d.x + i.z*d.x*d.y); }


__host__ __device__ double gauss_fun(double x, double sigma)
{

	return (exp(-pow(x, 2) / (2 * pow(sigma, 2)))) / (sigma*sqrt(2 * M_PI));

}


int iDivUp1(int a, int b) // Round a / b to nearest higher integer value

{
	return (a % b != 0) ? (a / b + 1) : (a / b);
}

__device__ __host__ void  Sort_Endpoints_C(ASF_vertex * m, Boundary b, Dimension d, Point step, uint32_t* output, uint32_t num_rows)
{


	ASF_vertex vertex = m[0];

	int index = 0;
	int index2 = 0;
	int x = 0;
	int y = 0;
	int z = 0;

	float xx = 0.;
	float yy = 0.;
	float zz = 0.;

	uint32_t* out1 = new uint32_t[num_rows];
	uint32_t* out2 = new uint32_t[num_rows];
	memset(out1, 0, num_rows*sizeof(uint32_t));
	memset(out2, 0, num_rows*sizeof(uint32_t));
	for (uint32_t row = 0; row < num_rows; row++)
	{
		vertex = m[row];

		if (!vertex.isInBoundary())
		{
			continue;

		}

		xx = (vertex.e.x - b.low.x) / step.x;
		yy = (vertex.e.y - b.low.y) / step.y;
		zz = (vertex.e.z - b.low.z) / step.z;

		x = (int)xx;
		y = (int)yy;
		z = (int)zz;

		index = z*(d.x*d.y) + y*d.x + x;
		index2 = zz*(d.x*d.y) + yy*d.x + xx;
		/*if (index2 != vertex.getRange())

		output[vertex.getRange()]++;;*/
		out1[index]++;;
		out2[index2]++;;
		output[vertex.getRange()]++;;



	}


}



__global__ void  Sort_Endpoints(ASF_vertex * m, Point* output, uint32_t num_rows)
{
	int row = blockIdx.x * blockDim.x + threadIdx.x;

	if (row >= num_rows)
		return;


}


__global__ void  Assign_Ball(ASF_vertex * m, float radius, uint32_t num_rows)
{
	int row = blockIdx.x * blockDim.x + threadIdx.x;

	if (row >= num_rows)
		return;
	dim3 grid(((num_rows + 510) / 512), 1, 1);  // (RSize-2) valid vertecis the first one is undefined 
	dim3 threads(512, 1, 1);
	float d_output[32678];

	ASF_vertex vertex = m[row];




}


__host__ __device__ void get_Lorenz_Field(float x, float y, float z, float &vxp, float &vyp, float 	&vzp)
{
	// 	float alpha = -1.28805;
	// 	float gamma = -0.502655;
	// 	float beta = 0.0314159;
	// 	vxp = gamma* cos(y) + alpha * sin(z);
	// 	vyp = alpha * cos(x) + beta* sin(x);
	// 	vzp = beta*cos(x) + gamma* sin(y);

	float sigma = 10.0;
	float ro = 28.0;
	float beta = 8.0 / 3.0;
	vxp = sigma * (y - x);
	vyp = (x*(ro - z)) - y;
	vzp = x*y - beta*z;
}



__host__ __device__ void trilinearInterpolation(ASF_vertex*m, Point p1, Point* v, uint32_t i, Dimension* _d, Boundary* _b, Point* step, float&vx, float& vy, float& vz, bool bForward)
{

	uint32_t xDim = _d->x;
	uint32_t yDim = _d->y;
	uint32_t zDim = _d->z;



	//if (!bForward)
	//	p1 = m[i].eb;

	Point highBoundary = _b->high;
	Point lowBoundary = _b->low;



	//	Point temp = (p1 - lowBoundary);
	//	Dimension index = temp.divide(*step);

	int indexx = (p1.x - lowBoundary.x) / step->x;
	int indexy = (p1.y - lowBoundary.y) / step->y;
	int indexz = (p1.z - lowBoundary.z) / step->z;


	uint32_t vi = indexx + indexy* _d->x + indexz*_d->x*_d->y; //getDimRange(*_d, index);
	//	float vx, vy, vz;
	get_Lorenz_Field(p1.x, p1.y, p1.z, vx, vy, vz);

	Point p0 = m[vi].p;


	double xd = (p1.x - p0.x) / step->x;//(p1[0]-step_x*int(p1[0]/step_x))/step_x;
	double yd = (p1.y - p0.y) / step->y;//(p1[1]-step_y*int(p1[1]/step_y))/step_y;
	double zd = (p1.z - p0.z) / step->z;//(p1[2]-step_z*int(p1[2]/step_z))/step_z;

	int kk = indexz;// (p1.z - lowBoundary.z) / step->z;;
	int jj = indexy;// (p1.y - lowBoundary.y) / step->y;;
	int ii = indexx;// (p1.x - lowBoundary.x) / step->x;;

	if (kk >= zDim - 1 || jj >= yDim - 1 || ii >= xDim - 1)
	{
		vx = v[((kk  * yDim + jj) * xDim + ii)].x;
		vy = v[((kk  * yDim + jj) * xDim + ii)].y;
		vz = v[((kk  * yDim + jj) * xDim + ii)].z;
		//return  v[((kk  * yDim + jj) * xDim + ii)];
	}

	if (p0 == p1)
	{
		vx = v[vi].x;
		vy = v[vi].y;
		vz = v[vi].z;

	}


	float v1[3];
	float v2[3];
	float v3[3];
	float v4[3];
	float v5[3];
	float v6[3];
	float v7[3];
	float v8[3];

	//vx = m_x[((kk * yDim + jj) * xDim + ii)];
	//vy = m_y[((kk * yDim + jj) * xDim + ii)];
	//vz = m_z[((kk * yDim + jj) * xDim + ii)];




	v1[0] = v[((kk * yDim + jj) * xDim + ii)].x;
	v1[1] = v[((kk * yDim + jj) * xDim + ii)].y;
	v1[2] = v[((kk * yDim + jj) * xDim + ii)].z;


	v2[0] = v[((kk * yDim + jj) * xDim + ii + 1)].x;
	v2[1] = v[((kk * yDim + jj) * xDim + ii + 1)].y;
	v2[2] = v[((kk * yDim + jj) * xDim + ii + 1)].z;

	v3[0] = v[((kk * yDim + (jj + 1)) * xDim + ii)].x;
	v3[1] = v[((kk * yDim + (jj + 1)) * xDim + ii)].y;
	v3[2] = v[((kk * yDim + (jj + 1)) * xDim + ii)].z;
	//

	v4[0] = v[((kk * yDim + (jj + 1)) * xDim + ii + 1)].x;
	v4[1] = v[((kk * yDim + (jj + 1)) * xDim + ii + 1)].y;
	v4[2] = v[((kk * yDim + (jj + 1)) * xDim + ii + 1)].z;

	//int idx2 = (k + 1)*(xDim*yDim) + (j + 1)*yDim + i + 1;

	v5[0] = v[(((kk + 1) * yDim + jj) * xDim + ii)].x;
	v5[1] = v[(((kk + 1) * yDim + jj) * xDim + ii)].y;
	v5[2] = v[(((kk + 1) * yDim + jj) * xDim + ii)].z;


	v6[0] = v[(((kk + 1) * yDim + jj) * xDim + ii + 1)].x;
	v6[1] = v[(((kk + 1) * yDim + jj) * xDim + ii + 1)].y;
	v6[2] = v[(((kk + 1) * yDim + jj) * xDim + ii + 1)].z;


	v7[0] = v[(((kk + 1) * yDim + (jj + 1)) * xDim + ii)].x;
	v7[1] = v[(((kk + 1) * yDim + (jj + 1)) * xDim + ii)].y;
	v7[2] = v[(((kk + 1) * yDim + (jj + 1)) * xDim + ii)].z;


	v8[0] = v[(((kk + 1) * yDim + (jj + 1)) * xDim + ii + 1)].x;
	v8[1] = v[(((kk + 1) * yDim + (jj + 1)) * xDim + ii + 1)].y;
	v8[2] = v[(((kk + 1) * yDim + (jj + 1)) * xDim + ii + 1)].z;


	double c00 = v1[0] * (1 - xd) + v2[0] * xd;
	double c10 = v3[0] * (1 - xd) + v4[0] * xd;
	double c01 = v5[0] * (1 - xd) + v6[0] * xd;
	double c11 = v7[0] * (1 - xd) + v8[0] * xd;

	double c0 = c00*(1 - yd) + c10*yd;
	double c1 = c01*(1 - yd) + c11*yd;

	//Point tempv = p1;
	vx = c0*(1 - zd) + c1*zd;


	c00 = v1[1] * (1 - xd) + v2[1] * xd;
	c10 = v3[1] * (1 - xd) + v4[1] * xd;
	c01 = v5[1] * (1 - xd) + v6[1] * xd;
	c11 = v7[1] * (1 - xd) + v8[1] * xd;

	c0 = c00*(1 - yd) + c10*yd;
	c1 = c01*(1 - yd) + c11*yd;

	vy = c0*(1 - zd) + c1*zd;

	c00 = v1[2] * (1 - xd) + v2[2] * xd;
	c10 = v3[2] * (1 - xd) + v4[2] * xd;
	c01 = v5[2] * (1 - xd) + v6[2] * xd;
	c11 = v7[2] * (1 - xd) + v8[2] * xd;

	c0 = c00*(1 - yd) + c10*yd;
	c1 = c01*(1 - yd) + c11*yd;

	vz = c0*(1 - zd) + c1*zd;
	float vxx, vyy, vzz;
	//	get_Lorenz_Field(p1.x, p1.y, p1.z, vxx, vyy, vzz);
	if (abs(vx - vxx) > exp(-6.0) || abs(vy - vyy) > exp(-6.0) || abs(vz - vzz) > exp(-6.0))
		printf("");
	//tempv.x = vx;
	//tempv.y = vy;
	//tempv.z = vz;
	//	return tempv;

	return;



}



void _Tracing_c(ASF_vertex * m, Point*v, Dimension* d, Boundary* b, Point* step, bool bForward, uint32_t whichData, uint32_t originalnum_rows, uint32_t num_rows, uint32_t level)
{
	//uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
	for (uint32_t row = 0; row < num_rows; row++)
	{


		if (level % 2 == 0)
		{
			if (row >= num_rows)
				continue;
		}
		else
		{
			if (row >= num_rows || row < originalnum_rows)
				continue;
		}
		/*if (row == 32818)
		printf("");*/
		ASF_vertex vertex = m[row];

		if (!vertex.isInBoundary())
		{
			continue;
		}

		Point e = vertex.e;
		float vx, vy, vz;
		trilinearInterpolation(m, e, v, vertex.getOldRange(), d, b, step, vx, vy, vz, bForward);

		/*if (row == 4591)
		printf("v = %f,%f,%f \n", vx, vy, vz);*/
		float dist = sqrt(vx*vx + vy*vy + vz*vz);

		vx = (vx / (dist*4.0))*step->x;
		vy = (vy / (dist*4.0))*step->y;
		vz = (vz / (dist*4.0))*step->z;
		/*_v = _v / 4.0;
		_v *= *step;*/
		//glColor3f(rgb[0], rgb[1], rgb[2]);    
		//		Point e = vertex.e;
		if (bForward)
		{
			vertex.e.x += vx;
			vertex.e.y += vy;
			vertex.e.z += vz;
		}
		//vertex.e += _v;
		else
		{
			vertex.e.x -= vx;
			vertex.e.y -= vy;
			vertex.e.z -= vz;
		}
		//vertex.eb -= _v;

		float ep[3];
		//generalstreamlineTracing_single(p1, bForward, ep, false);
	//	if (bForward)
		{
			if (ep[0] != vertex.e.x && ep[1] != vertex.e.y && ep[2] != vertex.e.z)
				printf("");
		}
		/*else
		{
			if (ep[0] != vertex.eb.x && ep[1] != vertex.eb.y && ep[2] != vertex.eb.z)
				printf("");

		}*/




		bool xy = false;
		bool yz = false;
		bool xz = false;
		//if (bForward)
		{
			if (vertex.checkInBoundary(b))
			{
				uint32_t range = vertex.getRange(b, step, d);
				if (vertex.isInNextLevel_yz())
					yz = true;
				if (vertex.isInNextLevel_xy())xy = true;
				if (vertex.isInNextLevel_xz())xz = true;
				vertex.setRange(range);
				if (xy)	vertex.setInNextLevel_xy();
				if (xz)vertex.setInNextLevel_xz();
				if (yz)vertex.setInNextLevel_yz();

			}
			else
			{
				vertex.unsetInBoundary();
			}
		}
		/*else
		{
			if (vertex.checkInBoundaryBackward(b))
			{
				if (vertex.isInNextLevel_yz())
					yz = true;
				if (vertex.isInNextLevel_xy())xy = true;
				if (vertex.isInNextLevel_xz())xz = true;
				uint32_t range = vertex.getRangeBackward(b, step, d);
				vertex.setRangeBackward(range);
				if (xy)	vertex.setInNextLevel_xy();
				if (xz)vertex.setInNextLevel_xz();
				if (yz)vertex.setInNextLevel_yz();

			}
			else
			{
				vertex.unsetInBoundary();
			}
		}*/




		m[row] = vertex;
	}
	return;
}


__host__ __device__ ASF_vertex*  DivideEdges(ASF_vertex* v1, ASF_vertex* v2, uint32_t range, Dimension*d, Point* step, Boundary*b, uint32_t faceNum, bool bForward)
{

	Point p_left_c = v1->p;
	ASF_vertex _tv = *v1;
	//return &_tv;
	p_left_c.x = (v1->p.x + v2->p.x) / 2;
	p_left_c.y = (v1->p.y + v2->p.y) / 2;
	p_left_c.z = (v1->p.z + v2->p.z) / 2;
	_tv.p = p_left_c;


	p_left_c.x = (v1->e.x + v2->e.x) / 2;
	p_left_c.y = (v1->e.y + v2->e.y) / 2;
	p_left_c.z = (v1->e.z + v2->e.z) / 2;

	//_tv.p = p_left_c;
	/*if (!bForward)
	{
		p_left_c.x = (v1->eb.x + v2->eb.x) / 2;
		p_left_c.y = (v1->eb.y + v2->eb.y) / 2;
		p_left_c.z = (v1->eb.z + v2->eb.z) / 2;
		_tv.eb = p_left_c;
		uint32_t range = _tv.getRangeBackward(b, step, d);

		_tv.setRangeBackward(range);
	}
	else*/
	{
		_tv.e = p_left_c;

		uint32_t range = _tv.getRange(b, step, d);

		_tv.setRange(range);
	}
	uint32_t oldrange = v1->getOldRange();   // _tv.getOldRange(&_b, &step, &_d);
	_tv.setOldRange(oldrange);
	_tv.setInBoundary();


	return &_tv;

}


__global__ void CheckRangeSetKernel(const ASF_vertex * m, uint32_t * Fr, const uint32_t num_rows)
{
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
	if (row >= num_rows)
		return;
	ASF_vertex vertex = m[row];

	if (vertex.isInBoundary())
		Fr[vertex.getOldRange()] = 1;
	else
		Fr[vertex.getOldRange()] = 0;
	return;

}

__global__ void CheckConnectivityKernel(const ASF_vertex * m, uint32_t * Fr, Dimension*d, const uint32_t num_rows)
{
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
	if (row >= num_rows)
		return;
	int iz = row / (d->x*d->y);
	int iy = (row - iz*(d->x*d->y)) / d->x;
	int ix = (row - iz*(d->x*d->y)) % d->x;

	Fr[row] = 0;

	if (ix == d->x - 1 || iy == d->y - 1 || iz == d->z - 1)
		return;
	ASF_vertex vertex = m[row];
	if (vertex.isInNextLevel_xy())
		Fr[row] += 1;
	if (vertex.isInNextLevel_xz())
		Fr[row] += 1;
	if (vertex.isInNextLevel_yz())
		Fr[row] += 1;

}



__host__ __device__ void CheckConnectivityKernel_c(const ASF_vertex * m, uint32_t * Fr, Dimension*d, const uint32_t num_rows)
{
	//	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
	for (uint32_t row = 0; row < num_rows; row++)
	{
		if (row >= num_rows)
			return;
		int iz = row / (d->x*d->y);
		int iy = (row - iz*(d->x*d->y)) / d->x;
		int ix = (row - iz*(d->x*d->y)) % d->x;

		Fr[row] = 0;

		if (ix == d->x - 1 || iy == d->y - 1 || iz == d->z - 1)
			continue;
		ASF_vertex vertex = m[row];
		if (vertex.isInNextLevel_xy())
			Fr[row] += 1;
		if (vertex.isInNextLevel_xz())
			Fr[row] += 1;
		if (vertex.isInNextLevel_yz())
			Fr[row] += 1;

	}

	return;

}

//__global__ void _CheckFacing(ASF_vertex*m, ASF_vertex*d_a, Dimension* d, Boundary* b, Point* step, bool bForward, uint32_t num_rows)
//{
//	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
//
//	if (row >= num_rows)
//		return;
//
//	ASF_vertex vertex1 = m[row];
//
//	if (!vertex1.isInSCC())
//	{
//		return;
//
//	}
//	ASF_vertex vertex2 = m[vertex1.Fr[0]];
//	ASF_vertex vertex3 = m[vertex1.Fr[1]];
//	ASF_vertex vertex4 = m[vertex1.Fr[2]];
//	ASF_vertex pa = d_a[row];
//	//Point* p = DivideFaces(&vertex1, &vertex2, &vertex3, &vertex4);
//	/* a[row * 5].p = p[0];
//	a[row * 5].p = p[0];
//	a[row * 5].p = p[0];
//	a[row * 5].p = p[0];*/
//	d_a[row * 2] = vertex1;
////	pa.p = p[0];
//	pa.setInBoundary();
//	pa.setOldRange(vertex1.getOldRange());
//	d_a[row * 2 + 1] = pa;
//
//
//
//
//	return;
//}


//void _CheckFacing_c(ASF_vertex*m, ASF_vertex*d_a, Dimension* d, Boundary* b, Point* step, bool bForward, uint32_t num_rows)
//{
//	// uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
//	for (uint32_t row = 0; row < num_rows; row++)
//	{
//
//		if (row >= num_rows)
//			return;
//
//		ASF_vertex vertex1 = m[row];
//		bool b = m[14831].isInSCC();
//
//		if (!vertex1.isInSCC())
//		{
//			continue;
//
//		}
//		ASF_vertex vertex2 = m[vertex1.Fr[0]];
//		ASF_vertex vertex3 = m[vertex1.Fr[1]];
//		ASF_vertex vertex4 = m[vertex1.Fr[2]];
//		ASF_vertex pa = d_a[row];
//		//Point* p = DivideFaces(&vertex1, &vertex2, &vertex3, &vertex4);
//		/* a[row * 5].p = p[0];
//		a[row * 5].p = p[0];
//		a[row * 5].p = p[0];
//		a[row * 5].p = p[0];*/
//		d_a[row * 2] = vertex1;
////		pa.p = p[0];
//		pa.setInBoundary();
//		pa.setOldRange(vertex1.getOldRange());
//		d_a[row * 2 + 1] = pa;
//
//
//	}
//
//	return;
//}


__global__ void Initialize(ASF_vertex * m, Dimension* d, Boundary* b, Point* step, bool bForward, uint32_t num_rows)
{
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;

	if (row >= num_rows)
		return;

	ASF_vertex vertex = m[row];
	vertex.setInBoundary();

	/*if (row < 10)
	printf("%d \n", row);*/

	//vertex.Fr_xy.x = row + 1;
	//vertex.Fr_xy.y = row + 1 + d->x;
	//vertex.Fr_xy.z = row + d->x;

	//vertex.Fr_xz.x = row + 1;
	//vertex.Fr_xz.y = row + d->y * d->x + 1;
	//vertex.Fr_xz.z = row + d->x*d->y;

	//vertex.Fr_yz.x = row + d->x;
	//vertex.Fr_yz.y = row + (1 + d->y) * d->x;
	//vertex.Fr_yz.z = row + d->x*d->y;


	//vertex.setInNextLevel_xy();
	//vertex.setInNextLevel_xz();
	//vertex.setInNextLevel_yz();
	m[row] = vertex;


}

__global__ void Initialize2(ASF_vertex * m, Dimension* d, Boundary* b, Point* step, bool bForward, uint32_t num_rows)
{
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;

	if (row >= num_rows)
		return;

	ASF_vertex vertex = m[row];
	if (vertex.isInBoundary())
	{
		vertex.setInNextLevel_xy();
		vertex.setInNextLevel_xz();
		vertex.setInNextLevel_yz();
	}

	m[row] = vertex;


}



__global__ void _TracingMGPU(ASF_vertex * m, Point*v, Boundary* b, Dimension* d, Point* step, bool bForward, uint32_t num_rows)//Dimension* d, Boundary* b, Point* step, bool bForward, uint32_t whichdata, uint32_t originalnum_rows, uint32_t num_rows, uint32_t level)
{
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;



	if (row >= num_rows)
		printf("");
	ASF_vertex vertex = m[row];

	if (!vertex.isInBoundary())
	{
		return;
	}

	float p1[3];
	//if (bForward)
	{
		p1[0] = vertex.e.x;
		p1[1] = vertex.e.y;
		p1[2] = vertex.e.z;
	}
	/*else
	{
		p1[0] = vertex.eb.x;
		p1[1] = vertex.eb.y;
		p1[2] = vertex.eb.z;
	}*/
	if (row == 191)
		printf("");
	Point e = vertex.e;
	float vx, vy, vz;
	trilinearInterpolation(m, e, v, vertex.getOldRange(), d, b, step, vx, vy, vz, bForward);

	/*if (row == 4591)
	printf("v = %f,%f,%f \n", vx, vy, vz);*/
	float dist = sqrt(vx*vx + vy*vy + vz*vz);

	vx = (vx / (dist*4.0))*step->x;
	vy = (vy / (dist*4.0))*step->y;
	vz = (vz / (dist*4.0))*step->z;
	/*_v = _v / 4.0;
	_v *= *step;*/
	//glColor3f(rgb[0], rgb[1], rgb[2]);    
	//	Point e = vertex.e;
	//if (bForward)
	{
		vertex.e.x += vx;
		vertex.e.y += vy;
		vertex.e.z += vz;
	}
	//vertex.e += _v;
	/*else
	{
		vertex.eb.x -= vx;
		vertex.eb.y -= vy;
		vertex.eb.z -= vz;
	}*/
	//vertex.eb -= _v;

	float ep[3];
	//generalstreamlineTracing_single(p1, bForward, ep, false);
	


	bool xy = false;
	bool yz = false;
	bool xz = false;
//	if (bForward)
	{
		if (vertex.checkInBoundary(b))
		{
			uint32_t range = vertex.getRange(b, step, d);
			if (vertex.isInNextLevel_yz())
				yz = true;
			if (vertex.isInNextLevel_xy())xy = true;
			if (vertex.isInNextLevel_xz())xz = true;
			vertex.setRange(range);
			if (xy)	vertex.setInNextLevel_xy();
			if (xz)vertex.setInNextLevel_xz();
			if (yz)vertex.setInNextLevel_yz();

		}
		else
		{
			vertex.unsetInBoundary();
		}
	}
	/*else
	{
		if (vertex.checkInBoundaryBackward(b))
		{
			if (vertex.isInNextLevel_yz())
				yz = true;
			if (vertex.isInNextLevel_xy())xy = true;
			if (vertex.isInNextLevel_xz())xz = true;
			uint32_t range = vertex.getRangeBackward(b, step, d);
			vertex.setRangeBackward(range);
			if (xy)	vertex.setInNextLevel_xy();
			if (xz)vertex.setInNextLevel_xz();
			if (yz)vertex.setInNextLevel_yz();

		}
		else
		{
			vertex.unsetInBoundary();
		}
	}*/




	m[row] = vertex;

	return;
}





__global__ void ADPSynchKernel_X(ASF_vertex * m, const uint32_t OldRange, const Edge * Bc, const uint32_t * Br, uint32_t * pivot,
	const uint32_t num_rows, uint32_t * COL_pivot = NULL)
{
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
	if (row >= num_rows)
		return;
	ASF_vertex vertex = m[row];


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
				ASF_vertex p_vertex = m[index];

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


__global__ void cuReduceEdges1(uint32_t *g_idata, uint32_t * g_odata)
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


__host__ __device__ void cuReduceEndpoints_C(ASF_vertex* m, uint32_t *Fr, Edge *Fc, bool bForward, const uint32_t n)
{
	/*uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;

	if (row > n)
	return;*/
	ASF_vertex vertex = m[0];

	for (uint32_t row = 0; row < n; row++)
	{
		vertex = m[row];

		//Edge edge ;
		if (vertex.isInBoundary())
		{
		//	if (bForward)
				Fc[Fr[row]].setValue(vertex.getRange());
			/*else
				Fc[Fr[row]].setValue(vertex.getRangeBackward());*/


			Fc[Fr[row]].setValidBit();

			//Fc[Fr[row]] = edge;
		}

	}

}

__global__ void cuReduceEndpoints(ASF_vertex* m, uint32_t *Fr, Edge *Fc, bool bForward, const uint32_t n)
{
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;

	if (row > n)
		return;
	ASF_vertex vertex = m[row];
	//Edge edge ;
	if (vertex.isInBoundary())
	{
		//if (bForward)
			Fc[Fr[row]].setValue(vertex.getRange());
		/*else
			Fc[Fr[row]].setValue(vertex.getRangeBackward());*/


		Fc[Fr[row]].setValidBit();

		//Fc[Fr[row]] = edge;
	}



}


__global__ void cuReduceEdges_shared(ASF_vertex* m, uint32_t *Fr, Edge *Fc, bool bForward, const uint32_t n)
{
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;

	if (row > n)
		return;
	__shared__ unsigned char _sa[32786];
	ASF_vertex vertex = m[row];
	//Edge edge ;
	if (vertex.isInBoundary())
	{
		//if (bForward)
			Fc[Fr[row]].setValue(vertex.getRange());
		/*else
			Fc[Fr[row]].setValue(vertex.getRangeBackward());*/


		Fc[Fr[row]].setValidBit();

		//Fc[Fr[row]] = edge;
	}




}



__global__ void cuReduceEdges(ASF_vertex* m, uint32_t *Fr, Edge *Fc, bool bForward, const uint32_t n)
{
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;

	if (row > n)
		return;
	ASF_vertex vertex = m[row];
	//Edge edge ;
	if (vertex.isInBoundary())
	{
	//	if (bForward)
			Fc[Fr[row]].setValue(vertex.getRange());
	/*	else
			Fc[Fr[row]].setValue(vertex.getRangeBackward());*/


		Fc[Fr[row]].setValidBit();

		//Fc[Fr[row]] = edge;
	}




}



__host__ __device__ bool checkEdge(ASF_vertex vertex1, ASF_vertex vertex2, Dimension* d, bool bForward)
{
	int range1 = vertex1.getRange();
	int range2 = vertex2.getRange();


	/*if (!bForward)
	{
		range1 = vertex1.getRangeBackward();
		range2 = vertex2.getRangeBackward();
	}
*/


	int iz1 = range1 / (d->x*d->y);
	int iy1 = (range1 - iz1*(d->x*d->y)) / d->x;
	int ix1 = (range1 - iz1*(d->x*d->y)) % d->x;

	int iz2 = range2 / (d->x*d->y);
	int iy2 = (range2 - iz2*(d->x*d->y)) / d->x;
	int ix2 = (range2 - iz2*(d->x*d->y)) % d->x;

	float dist = 0.;
	dist = sqrt((float)((iz2 - iz1)*(iz2 - iz1) + (iy2 - iy1)*(iy2 - iy1) + (ix2 - ix1)*(ix2 - ix1)));

	/*if (iz1 == d->z || iz1 == d->z - 1 || iy1 == d->y || iy1 == d->y - 1 || ix1 == d->x || ix1 == d->x - 1)
	printf("");*/

	/*if (iz2 == d->z || iz2 == d->z - 1 || iy2 == d->y || iy2 == d->y - 1 || ix2 == d->x || ix2 == d->x - 1)
	printf("");*/


	//printf("%f --", dist);

	//return false;

	if (dist <= 2) //abs(iz1 - iz2) <= 2 && abs(iy2 - iy1) <= 2 && abs(ix2 - ix1) <= 2)

		return true;
	//if (!bForward)
	//	printf("%d ... %d,%d,..., %d,%d,%d,...,%d,%d,%d \n",vertex1.getOldRange() ,range1,range2, iz1, iy1, ix1,iz2,iy2,ix2);

	return false;




}





//pair <uint32_t, float> Flow_Combinatorialization(ASF_vertex* _a,Point* _v, Boundary* _b, Dimension* _d, uint32_t  _tau, Edge ** oFc, uint32_t ** oFr, uint32_t * oRSize,uint32_t whichData, bool bForward)
//{
//
//	//-----------GPU initialization---------------------------->
//	uint32_t  * d_Fr;
//	uint32_t  * d_To;
//
//	Edge * d_Fc;
//
//	ASF_vertex * d_m;
//	Boundary *d_b;
//	Dimension* d_d;
//	Point* d_step;
//	Point* d_v;
//	float temp;
//
//	
//	uint32_t terminate = 1;
//	int interruptions = 0;
//	uint32_t CSize = _d->x*_d->y*_d->z;
//
//	uint32_t currentSize = _d->x*_d->y*_d->z;
//#ifdef _DEBUG
//	int FWD_ints = 0;
//	int OWCTY_ints = 0;
//	int BWD_ints = 0;
//	StopWatchInterface* KernelTime = 0;
//	StopWatchInterface* IntTime = 0;
//	(sdkCreateTimer(&KernelTime));
//	(sdkCreateTimer(&IntTime));
//#endif
//
//	//if (!_DeviceSet) {
//	//	_DeviceSet = true;
//	//	checkCudaErrors(cudaSetDevice(1));
//	//}
//
//
//	cudaError_t e1, e2, e3, e4, e5, e6, e7;
//	checkCudaErrors(e1 = cudaMalloc((void**)&d_b, sizeof(Boundary)));
//	checkCudaErrors(e2 = cudaMalloc((void**)&d_d, sizeof(Dimension)));
//	checkCudaErrors(e3 = cudaMalloc((void**)&d_v, currentSize* sizeof(Point)));
//	checkCudaErrors(e4 = cudaMalloc((void**)&d_step, sizeof(Point)));
//
//
//	//checkCudaErrors(e4 = cudaMalloc((void**)&d_Br, RSize * sizeof(uint32_t)));
//	checkCudaErrors(e4 = cudaMalloc((void**)&d_Fr, CSize * sizeof(uint32_t)));
//	checkCudaErrors(e4 = cudaMalloc((void**)&d_To, CSize * sizeof(uint32_t)));
//
//
//	checkCudaErrors(e1 = cudaMalloc((void**)&d_b, sizeof(Boundary)));
//
//
//
//
//	if (e1 == cudaErrorMemoryAllocation || e2 == cudaErrorMemoryAllocation ||
//		e3 == cudaErrorMemoryAllocation || e4 == cudaErrorMemoryAllocation ||
//		e5 == cudaErrorMemoryAllocation || e6 == cudaErrorMemoryAllocation ||
//		e7 == cudaErrorMemoryAllocation) {
//		throw "Error: Not enough memory on GPU\n";
//	}
//
//	uint32_t* To = new uint32_t[CSize];
//	uint32_t* From = new uint32_t[CSize];
//	uint32_t* Fr;// = new uint32_t[CSize];
//	//To2 = new uint32_t[CSize];
//	//col
//	//unsigned int COLTime = 0;
//	StopWatchInterface* COLTime = 0;
//	(sdkCreateTimer(&COLTime));
//	uint32_t * d_temp_COL;
//	uint32_t * d_temp_COL2;
//	uint32_t * d_COL_OldRange;
//	bool COL_used = false;
//	COL_vertex * d_cm;
//
//
//
//	Point step = (_b->high - _b->low);
//	step /= (*_d);;
//
//	float radius = step.x;
//	checkCudaErrors(cudaMemcpy(d_b, _b, sizeof(Boundary), cudaMemcpyHostToDevice));
//	checkCudaErrors(cudaMemcpy(d_d, _d, sizeof(Dimension), cudaMemcpyHostToDevice));
//	checkCudaErrors(cudaMemcpy(d_step, &step, sizeof(Dimension), cudaMemcpyHostToDevice));
//
//	//checkCudaErrors(cudaMemcpy(d_m, _a, (CSize)* sizeof(ASF_vertex), cudaMemcpyHostToDevice));
//	checkCudaErrors(cudaMemcpy(d_v, _v, (currentSize)* sizeof(Point), cudaMemcpyHostToDevice));
//
//	//unsigned int SCCTime = 0;
//	StopWatchInterface* SCCTime = 0;
//	(sdkCreateTimer(&SCCTime));
//	(sdkStartTimer(&SCCTime));
//
//	dim3 grid(((CSize + 510) / 512), 1, 1);  // (RSize-2) valid vertecis the first one is undefined 
//	dim3 threads(512, 1, 1);
//	dim3 grid1(((CSize + 255) / 256), 1, 1);
//	dim3 threads1(blockSize, 1, 1);
//	dim3 grid2(1, 1, 1);
//	dim3 threads2(32, 1, 1);
//
//	{
//
//		printf("Vertices: %u\n", CSize - 2);
//		//printf("Edges: %u\n", CSize);
//	}
//	//-----------Main algorithm-------------------------------->
//
//	//-----------Trimming-------------------------------------->
//	//if (trimm) 
//
//	
//	//Point _v;
//	uint32_t currTau = _tau;
//	uint32_t tits = 1;
//	uint32_t currSize = CSize;
//	uint32_t icycle = currTau / tits;
//	ASF_vertex * currA;
//	currA = _a;
//	//checkCudaErrors(e5 = cudaMalloc((void**)&d_m, (CSize)* sizeof(ASF_vertex)));
//	
//	
//	checkCudaErrors(e5 = cudaMalloc((void**)&d_To, (currSize)* 30*sizeof(uint32_t)));
//	checkCudaErrors(e6 = cudaMalloc((void**)&d_Fr, (currSize)* 30*sizeof(uint32_t)));
//	checkCudaErrors(e5 = cudaMalloc((void**)&d_m, (currSize)*30* sizeof(ASF_vertex)));
//	checkCudaErrors(cudaMemcpy(d_m, _a, (currSize)* sizeof(ASF_vertex), cudaMemcpyHostToDevice));
//	Initialize << <grid, threads >> >(d_m, d_d, d_b, d_step, bForward, currSize);
//	//checkCudaErrors(cudaMemcpy(_a, d_m, (CSize)* sizeof(ASF_vertex), cudaMemcpyDeviceToHost));
//
//
//	uint32_t* Fr2 = new uint32_t[CSize * 30];
//	Fr = new uint32_t[CSize * 30];
//	uint32_t* Fi = new uint32_t[CSize*30];
//	ASF_vertex* _aa;
//	//ASF_vertex* tempa = new ASF_vertex[41898];;
//	_aa = new ASF_vertex[CSize*30];
//	currA = _aa;
//	memcpy(currA, _a, currSize*sizeof(ASF_vertex));
//
//	uint32_t oldsize = CSize;
//	int level = 0;
//	//To2 = To;
//	//Fr2 = Fi ;
//
//	do{
//		 tits = 1;
//
//		 dim3 grid(((currSize + 510) / 512), 1, 1);  // (RSize-2) valid vertecis the first one is undefined 
//
//		do {
//
//			//checkCudaErrors(cudaMemset(d_pivot, 0, sizeof(uint32_t)));
//
//			//if (!bForward)
//			//	_Tracing_c (currA,_v, _d, _b, &step, bForward, whichData, CSize, currSize,level );
//			Tracing << <grid, threads >> >(d_m, d_v, d_d, d_b, d_step, bForward, whichData,oldsize, currSize,level );
//			
//			
//
//
//			tits--;
//		} while (tits > 0);// && terminate);
//		tits = 1;
//
//		//if (Fr[CSize - 1] > 0)
//		{
//		//	checkCudaErrors(cudaMemcpy(currA, d_m, (currSize)* sizeof(ASF_vertex), cudaMemcpyDeviceToHost));
//			/*if (icycle == 70)
//				_CheckNeighborhood_c3(currA, Fr, Fi, _d, _b, &step, bForward, oldsize, currSize, level);*/
//			//_CheckNeighborhood_c3(currA, Fr, Fi, _d, _b, &step, bForward, oldsize, currSize, level);
////			_CheckNeighborhood << <grid, threads >> >(d_m, d_To,  d_d, d_b, d_step, bForward, oldsize, currSize, level);
//
//		}
//
//		//_CheckNeighborhood_c(currA, Fr, _d, _b, &step, bForward, CSize, currSize);
//	
//		
//		//CheckConnectivityKernel_c(currA, To, _d, CSize);
//
//		//checkCudaErrors(cudaMemcpy(d_m, currA, (currSize)* sizeof(ASF_vertex), cudaMemcpyHostToDevice));
//		//checkCudaErrors(cudaMemcpy(d_To, Fi, (currSize)* sizeof(uint32_t), cudaMemcpyHostToDevice));
//
//	/*	if (currA[32772].Fr_xy.z > currSize)
//			printf("");*/
//		//CheckConnectivityKernel << <grid, threads >> >(d_m, d_To,d_d, currSize );
//
//		//checkCudaErrors(cudaMemcpy(currA, d_m, (currSize)* sizeof(ASF_vertex), cudaMemcpyDeviceToHost));
//		
//		/*CheckConnectivityKernel_c (_a, To, CSize - 1);*/
//		checkCudaErrors(cudaMemcpy(Fi, d_To, (currSize) * sizeof(uint32_t), cudaMemcpyDeviceToHost));
////		checkCudaErrors(cudaMemcpy(d_To, To, (CSize) * sizeof(uint32_t), cudaMemcpyHostToDevice));
//		//To2 = Fi;
////		reduction(d_To, d_Fr, Fi, Fr2, currSize);
//		checkCudaErrors(cudaMemcpy(d_Fr, Fr2, (currSize)* sizeof(uint32_t), cudaMemcpyHostToDevice));
//		//if (icycle == 5)
//		//	printf("");
//		if (Fr2[currSize - 1] > 0)
//		{
//				
//			printf("%d \n", Fr2[currSize - 1]);
//			checkCudaErrors(cudaMemcpy(currA, d_m, (currSize)* sizeof(ASF_vertex), cudaMemcpyDeviceToHost));
////			InsertSeedPoint << <grid, threads >> > (d_m, d_Fr, d_d, d_step, (currTau / tits) - icycle, currSize, bForward);
//			
//			oldsize = currSize;
//			currSize = currSize + Fr2[currSize - 1] * 5;
//
//			//checkCudaErrors(cudaMemcpy(d_m, currA, (currSize)* sizeof(ASF_vertex), cudaMemcpyHostToDevice));
//			//memset(Fi, 0, currSize*sizeof(uint32_t));
//
//			icycle++;
//			level = 1;
//
//
//		}
//		else
//			level = 0;
//
//	/*	else if (currSize != oldsize)
//			currA = _aa;*/
//		//	Initialize2 << <grid, threads >> >(d_m, d_d, d_b, d_step, bForward, currSize);*/
//		printf(" \n %d \n", icycle);
//		icycle--;
//
//
//	} while (icycle > 0 );// && terminate);
//
//	dim3 grid_2(((currSize + 510) / 512), 1, 1);  // (RSize-2) valid vertecis the first one is undefined 
//
//		//checkCudaErrors(cudaMemset(d_pivot, 0, sizeof(uint32_t)));
//		checkCudaErrors(e3 = cudaMalloc((void**)&d_Fc, currSize * sizeof(Edge)));
//		Edge* Fc = new Edge[currSize];
//		checkCudaErrors(cudaMemcpy(d_m, currA, (currSize)* sizeof(ASF_vertex), cudaMemcpyHostToDevice));
//
//		checkCudaErrors(cudaFree(d_To));
//		checkCudaErrors(cudaFree(d_Fr));
//		checkCudaErrors(e5 = cudaMalloc((void**)&d_To, (currSize)* sizeof(uint32_t)));
//		checkCudaErrors(e6 = cudaMalloc((void**)&d_Fr, (currSize)* sizeof(uint32_t)));
//		checkCudaErrors(cudaMemcpy(currA, d_m, (currSize)* sizeof(ASF_vertex), cudaMemcpyDeviceToHost));
//
//		memset(Fi, 0, currSize*sizeof(uint32_t));
//		for (int i = 0; i < currSize; i++)
//		{
//			ASF_vertex vertex = currA[i];
//			if (vertex.isInBoundary())
//				Fi[vertex.getOldRange()]++;
//		}
//		//CheckRangeSetKernel << <grid_2, threads >> >(d_m, d_To, currSize - 1);
//		To = new uint32_t[currSize];
//		//checkCudaErrors(cudaMemcpy(&terminate, d_pivot, sizeof(uint32_t), cudaMemcpyDeviceToHost));
//	
//
//		//checkCudaErrors(cudaMemcpy(To, d_To, currSize * sizeof(uint32_t), cudaMemcpyDeviceToHost));
//
//		//checkCudaErrors(cudaMemcpy(d_To, Fi, (currSize)* sizeof(uint32_t), cudaMemcpyHostToDevice));
//
////		reduction(d_To, d_Fr, Fi, Fr, CSize);
//		memset(Fi, 0, CSize*sizeof(uint32_t));
//		checkCudaErrors(cudaMemcpy(d_Fr, Fr, currSize * sizeof(uint32_t), cudaMemcpyHostToDevice));
//		for (int row = 0; row < currSize;row++)
//		{
//			ASF_vertex vertex = currA[row];
//			if (vertex.isInBoundary())
//			{
//				uint32_t i = Fr[vertex.getOldRange()] + Fi[vertex.getOldRange()] - 1;
//				if (bForward)
//					Fc[i].setValue(vertex.getRange());
//				else
//					Fc[i].setValue(vertex.getRangeBackward());
//
//				Fc[i].setValidBit();
//				Fi[vertex.getOldRange()]++;
//			}
//					
//		}
//	//	cuReduceEdges << <grid_2, threads >> >(d_m, d_Fr, d_Fc, bForward, currSize);
//		//cuReduceEdges_shared << <grid_2, threads >> >(d_m, d_Fr, d_Fc, bForward, currSize);
//
//	//checkCudaErrors(cudaMemcpy(Fc, d_Fc, CSize * sizeof(Edge), cudaMemcpyDeviceToHost));
//
//
//	
//
//	//===============================================================================
//	//	checkCudaErrors(cudaMemcpy(_a, d_m, (currSize)* sizeof(ASF_vertex), cudaMemcpyDeviceToHost));
//
//
//	//	checkCudaErrors(cudaMemcpy(To, d_To, CSize * sizeof(uint32_t), cudaMemcpyDeviceToHost));
//	//	checkCudaErrors(cudaMemcpy(From, d_From, CSize * sizeof(uint32_t), cudaMemcpyDeviceToHost));
//	//	checkCudaErrors(cudaMemcpy(Fr, d_Fr, currSize * sizeof(uint32_t), cudaMemcpyDeviceToHost));
//	//	checkCudaErrors(cudaMemcpy(Fc, d_Fc, currSize * sizeof(Edge), cudaMemcpyDeviceToHost));
//
//		uint32_t r = currA[14731].getRange();
//
//		if (!bForward)
//			r = currA[14731].getRangeBackward();
//	
////	checkCudaErrors(cudaMemcpy(To, d_To, (1) * sizeof(uint32_t), cudaMemcpyDeviceToHost));
////	checkCudaErrors(cudaMemcpy(_a, d_m, (CSize)* sizeof(ASF_vertex), cudaMemcpyDeviceToHost));
//
//
//
//
//	*oFc = Fc;
//	*oFr = Fr;
//	*oRSize = Fr[CSize - 1];
//
//	//<----------Trimming---------------------------------------
//
//
//	pair <uint32_t, float> result;
//	return result;
//}
////{


//__global__ InsertSeedPoint(ASF_vertex* m, fEdge*eg, uint32_t * Fr, Dimension*d, Point* step,  uint32_t num_vertexes, uint32_t numEdges, bool bForward)
//{
//
//
//	ASF_vertex vertex1 = m[0];
//	//	om[0] = vertex1;
//
//	//for (uint32_t row = 1; row < num_rows; row++)
//	{
//
//		/* if (row < 2*originalnumrows)
//		om[row] = vertex1;*/
//		if ((row >= 1 && (Fr[row] - Fr[row - 1]) == 0))
//			continue;
//
//		//if ((Fr[row] - Fr[row - 1]) == 0)
//		{
//
//			int curVertexIdx = originalnumrows + Fr[row];
//			int curEdgeIdx = num_rows + Fr[row] * 3;
//			int curtriangleIdx = trianglenum + Fr[row] * 2;
//
//			fEdge edge1 = eg[row];
//			fEdge edge2 = eg[row];
//			fEdge edge3 = eg[row];
//			fEdge edge4 = eg[row];
//			ASF_vertex vertex1 = m[edge1.v1];
//			ASF_vertex vertex2 = m[edge1.v2];
//			ASF_vertex vedge = DivideEdges(&vertex1, &vertex2, vertex1.getOldRange(), 0, bForward);
//			m[curVertexIdx] = vedge;
//
//			int triangleIdx1 = edge1.E2T[0];
//			int triangleIdx2 = edge1.E2T[1];
//
//			if (vertex1.getOldRange() == 4529 || vertex2.getOldRange() == 4529)
//				printf("");
//
//			fTriangle t1 = tr[triangleIdx1];
//			if (row == 179098)
//				printf("");
//			if (curtriangleIdx + 1 == 131306 || curtriangleIdx == 131306)// || triangleIdx2 == 14393 || curtriangleIdx == 14393 || curtriangleIdx == 14392)
//				printf("");
//			fTriangle to = tr[triangleIdx1];
//			fTriangle t2 = tr[triangleIdx2];
//			fTriangle to2 = tr[triangleIdx2];
//			fTriangle t3 = tr[triangleIdx1];
//			fTriangle t4 = tr[triangleIdx2];
//
//
//			if (to.edge[0] != row && to.edge[1] != row && to.edge[2] != row)
//				printf("");
//
//
//			if (t1.edge[2] == row)
//			{
//				t1.edge[0] = to.edge[2];
//				t1.edge[1] = curEdgeIdx + 1;
//
//				if (eg[to.edge[1]].v1 == edge1.v1 || eg[to.edge[1]].v2 == edge1.v1)
//				{
//					t1.edge[2] = to.edge[1];
//					t3.edge[1] = to.edge[0];
//
//				}
//				else if (eg[to.edge[0]].v1 == edge1.v1 || eg[to.edge[0]].v2 == edge1.v1)
//				{
//					t1.edge[2] = to.edge[0];
//					t3.edge[1] = to.edge[1];
//
//				}
//
//				t3.edge[0] = curEdgeIdx;
//				t3.edge[2] = curEdgeIdx + 1;
//
//
//			}
//
//			else if (t1.edge[1] == row)
//			{
//				t1.edge[0] = to.edge[1];
//				t1.edge[1] = curEdgeIdx + 1;
//				//t1.edge[2] = to.edge[0];
//
//				if (eg[to.edge[0]].v1 == edge1.v1 || eg[to.edge[0]].v2 == edge1.v1)
//				{
//					t1.edge[2] = to.edge[0];
//					t3.edge[1] = to.edge[2];
//
//				}
//				else if (eg[to.edge[2]].v1 == edge1.v1 || eg[to.edge[2]].v2 == edge1.v1)
//				{
//					t1.edge[2] = to.edge[2];
//					t3.edge[1] = to.edge[0];
//
//				}
//
//				t3.edge[0] = curEdgeIdx;
//				//t3.edge[1] = to.edge[2];
//				t3.edge[2] = curEdgeIdx + 1;
//
//			}
//
//			else if (t1.edge[0] == row)
//			{
//				t1.edge[0] = to.edge[0];
//				t1.edge[1] = curEdgeIdx + 1;
//				//	t1.edge[2] = to.edge[2];
//
//
//				if (eg[to.edge[1]].v1 == edge1.v1 || eg[to.edge[1]].v2 == edge1.v1)
//				{
//					t1.edge[2] = to.edge[1];
//					t3.edge[1] = to.edge[2];
//
//				}
//				else if (eg[to.edge[2]].v1 == edge1.v1 || eg[to.edge[2]].v2 == edge1.v1)
//				{
//					t1.edge[2] = to.edge[2];
//					t3.edge[1] = to.edge[1];
//
//				}
//
//				t3.edge[0] = curEdgeIdx;
//				//	t3.edge[1] = to.edge[1];
//				t3.edge[2] = curEdgeIdx + 1;
//
//
//
//			}
//
//			//=====================================
//
//			edge1.v2 = curVertexIdx;
//			edge2.v1 = curVertexIdx;
//
//			edge3.v1 = curVertexIdx;
//			edge4.v1 = curVertexIdx;
//
//			if (row == 178834)
//				printf("");
//
//			if (t2.edge[2] == row)
//			{
//				t2.edge[0] = to2.edge[2];
//				t2.edge[1] = curEdgeIdx + 2;
//				if (eg[to2.edge[0]].v1 == edge1.v1 || eg[to2.edge[0]].v2 == edge1.v1)
//				{
//					t2.edge[2] = to2.edge[0];
//					t4.edge[1] = to2.edge[1];
//				}
//				else if (eg[to2.edge[1]].v1 == edge1.v1 || eg[to2.edge[1]].v2 == edge1.v1)
//				{
//					t2.edge[0] = to2.edge[1];
//					t4.edge[1] = to2.edge[0];
//
//				}
//
//
//				t4.edge[0] = curEdgeIdx;
//				t4.edge[2] = curEdgeIdx + 2;
//
//
//			}
//
//			else if (t2.edge[1] == row)
//			{
//				t2.edge[0] = to2.edge[1];
//				t2.edge[1] = curEdgeIdx + 2;
//				//t2.edge[2] = to2.edge[0];
//
//				if (eg[to2.edge[0]].v1 == edge1.v1 || eg[to2.edge[0]].v2 == edge1.v1)
//				{
//					t2.edge[2] = to2.edge[0];
//					t4.edge[1] = to2.edge[2];
//				}
//				else if (eg[to2.edge[2]].v1 == edge1.v1 || eg[to2.edge[2]].v2 == edge1.v1)
//				{
//					t2.edge[2] = to2.edge[2];
//					t4.edge[1] = to2.edge[0];
//
//				}
//
//				t4.edge[0] = curEdgeIdx;
//				//t4.edge[1] = to2.edge[2];
//				t4.edge[2] = curEdgeIdx + 2;
//
//			}
//
//			else if (t2.edge[0] == row)
//			{
//				t2.edge[0] = to2.edge[0];
//				t2.edge[1] = curEdgeIdx + 2;
//
//				if (eg[to2.edge[1]].v1 == edge1.v1 || eg[to2.edge[1]].v2 == edge1.v1)
//				{
//					t2.edge[2] = to2.edge[1];
//					t4.edge[1] = to2.edge[2];
//				}
//				else if (eg[to2.edge[2]].v1 == edge1.v1 || eg[to2.edge[2]].v2 == edge1.v1)
//				{
//					t2.edge[2] = to2.edge[2];
//					t4.edge[1] = to2.edge[1];
//
//				}
//				//t2.edge[2] = to2.edge[0];
//
//				t4.edge[0] = curEdgeIdx;
//				//t4.edge[1] = to2.edge[2];
//				t4.edge[2] = curEdgeIdx + 2;
//
//			}
//
//
//
//
//			if (to.edge[0] == row)
//			{
//				if ((eg[to.edge[2]].v1 == edge1.v1) /*|| (eg[to.edge[2]].v2 == edge2.v2)*/)
//					edge3.v2 = eg[to.edge[2]].v2;
//				else if ((eg[to.edge[2]].v2 == edge1.v1))
//					edge3.v2 = eg[to.edge[2]].v1;
//
//
//				else if ((eg[to.edge[1]].v1 == edge1.v1) /*|| (eg[to.edge[2]].v2 == edge2.v2)*/)
//					edge3.v2 = eg[to.edge[1]].v2;
//				else if ((eg[to.edge[1]].v2 == edge1.v1))
//					edge3.v2 = eg[to.edge[1]].v1;
//
//				else if (t1.edge[0] == row)
//					printf("");
//
//			}
//
//
//			else if (to.edge[1] == row)
//			{
//				if (eg[to.edge[2]].v2 == edge1.v1 /*|| eg[to.edge[2]].v2 == edge1.v2*/)
//					edge3.v2 = eg[to.edge[2]].v1;
//				else if ((eg[to.edge[2]].v1 == edge1.v1 /*|| eg[to.edge[2]].v1 == edge1.v2*/))
//					edge3.v2 = eg[to.edge[2]].v2;
//
//				else if (eg[to.edge[0]].v2 == edge1.v1 /*|| eg[to.edge[2]].v2 == edge1.v2*/)
//					edge3.v2 = eg[to.edge[0]].v1;
//				else if ((eg[to.edge[0]].v1 == edge1.v1 /*|| eg[to.edge[2]].v1 == edge1.v2*/))
//					edge3.v2 = eg[to.edge[0]].v2;
//
//
//				else if (t1.edge[1] == row)
//					printf("");
//			}
//
//
//
//			else if (to.edge[2] == row)
//			{
//				if (eg[to.edge[0]].v2 == edge1.v1/* || eg[to.edge[0]].v2 == edge1.v2*/)
//					edge3.v2 = eg[to.edge[0]].v1;
//				else if (eg[to.edge[0]].v1 == edge1.v1 /*|| eg[to.edge[0]].v1 == edge1.v2)*/)
//					edge3.v2 = eg[to.edge[0]].v2;
//
//				else if (eg[to.edge[1]].v2 == edge1.v1/* || eg[to.edge[0]].v2 == edge1.v2*/)
//					edge3.v2 = eg[to.edge[1]].v1;
//				else if (eg[to.edge[1]].v1 == edge1.v1 /*|| eg[to.edge[0]].v1 == edge1.v2)*/)
//					edge3.v2 = eg[to.edge[1]].v2;
//
//				else if (t1.edge[2] == row)
//					printf("");
//
//			}
//
//			//==========================================================================================
//
//
//
//			if (to2.edge[0] == row)
//			{
//				if (eg[to2.edge[1]].v2 == edge1.v1 /*|| eg[to2.edge[1]].v2 == edge1.v2*/)
//					edge4.v2 = eg[to2.edge[1]].v1;
//				else if (eg[to2.edge[1]].v1 == edge1.v1 /*|| eg[to2.edge[1]].v1 == edge1.v2)*/)
//					edge4.v2 = eg[to2.edge[1]].v2;
//
//				else if (eg[to2.edge[2]].v2 == edge1.v1 /*|| eg[to2.edge[1]].v2 == edge1.v2*/)
//					edge4.v2 = eg[to2.edge[2]].v1;
//				else if (eg[to2.edge[2]].v1 == edge1.v1 /*|| eg[to2.edge[1]].v1 == edge1.v2)*/)
//					edge4.v2 = eg[to2.edge[2]].v2;
//
//
//				else if (t2.edge[0] == row)
//					printf("");
//			}
//
//
//
//			else if (to2.edge[1] == row)
//			{
//				if (eg[to2.edge[2]].v2 == edge1.v1 /*|| eg[to2.edge[2]].v1 == edge1.v1*/)
//					edge4.v2 = eg[to2.edge[2]].v1;
//				else if (eg[to2.edge[2]].v1 == edge1.v1 /*|| eg[to2.edge[2]].v1 == edge1.v2)*/)
//					edge4.v2 = eg[to2.edge[2]].v2;
//
//				else if (eg[to2.edge[0]].v2 == edge1.v1 /*|| eg[to2.edge[2]].v1 == edge1.v1*/)
//					edge4.v2 = eg[to2.edge[0]].v1;
//				else if (eg[to2.edge[0]].v1 == edge1.v1 /*|| eg[to2.edge[2]].v1 == edge1.v2)*/)
//					edge4.v2 = eg[to2.edge[0]].v2;
//
//				else if (t2.edge[1] == row)
//					printf("");
//			}
//
//
//			else if (to2.edge[2] == row)
//			{
//				if (eg[to2.edge[1]].v2 == edge1.v1 /*|| eg[to2.edge[1]].v2 == edge1.v1*/)
//					edge4.v2 = eg[to2.edge[1]].v1;
//				else if (eg[to2.edge[1]].v1 == edge1.v1/* || eg[to2.edge[0]].v2 == edge1.v1*/)
//					edge4.v2 = eg[to2.edge[1]].v2;
//
//				else if (eg[to2.edge[0]].v2 == edge1.v1 /*|| eg[to2.edge[1]].v2 == edge1.v1*/)
//					edge4.v2 = eg[to2.edge[0]].v1;
//				else if (eg[to2.edge[0]].v1 == edge1.v1/* || eg[to2.edge[0]].v2 == edge1.v1*/)
//					edge4.v2 = eg[to2.edge[0]].v2;
//
//
//				else if (t2.edge[2] == row)
//					printf("");
//			}
//
//			edge1.v2 = curVertexIdx;
//			edge2.v1 = curVertexIdx;
//
//			edge3.v1 = curVertexIdx;
//			edge4.v1 = curVertexIdx;
//
//			edge2.E2T[0] = curtriangleIdx;
//			edge2.E2T[1] = curtriangleIdx + 1;
//
//			edge3.E2T[0] = triangleIdx1;;
//			edge3.E2T[1] = curtriangleIdx;
//
//			edge4.E2T[0] = triangleIdx2;;
//			edge4.E2T[1] = curtriangleIdx + 1;
//			eg[row] = edge1;
//			eg[curEdgeIdx] = edge2;
//			eg[curEdgeIdx + 1] = edge3;
//			eg[curEdgeIdx + 2] = edge4;
//
//			t3.T2F = to.T2F;
//			t4.T2F = to2.T2F;
//
//
//			if (eg[t4.edge[1]].E2T[0] == triangleIdx2)
//				eg[t4.edge[1]].E2T[0] = curtriangleIdx + 1;
//			else if (eg[t4.edge[1]].E2T[1] == triangleIdx2)
//				eg[t4.edge[1]].E2T[1] = curtriangleIdx + 1;
//
//			if (eg[t3.edge[1]].E2T[0] == triangleIdx1)
//				eg[t3.edge[1]].E2T[0] = curtriangleIdx;
//			else if (eg[t3.edge[1]].E2T[1] == triangleIdx1)
//				eg[t3.edge[1]].E2T[1] = curtriangleIdx;
//			if (curtriangleIdx == 131272 || curtriangleIdx + 1 == 131272)
//				printf("");
//
//			if (to.T2F != to2.T2F)
//				printf("");
//			tr[triangleIdx1] = t1;
//			tr[triangleIdx2] = t2;
//
//			tr[curtriangleIdx] = t3;
//			tr[curtriangleIdx + 1] = t4;
//
//			float rgb[3];
//			rgb[0] = 1.0;
//			rgb[0] = 0.0;
//			rgb[0] = 0.0;
//
//			if (t1.edge[0] != row)
//				printf("");
//
//			if (t1.edge[0] == t1.edge[1] || t1.edge[0] == t1.edge[2] || t1.edge[2] == t1.edge[1])
//				printf("");
//
//			if (t2.edge[0] == t2.edge[1] || t2.edge[0] == t2.edge[2] || t2.edge[2] == t2.edge[1])
//				printf("");
//
//			if (t3.edge[0] == t3.edge[1] || t3.edge[0] == t3.edge[2] || t3.edge[2] == t3.edge[1])
//				printf("");
//
//
//			if (t4.edge[0] == t4.edge[1] || t4.edge[0] == t4.edge[2] || t4.edge[2] == t4.edge[1])
//				printf("");
//
//
//
//			if (eg[t1.edge[0]].v2 != eg[t1.edge[1]].v1 && eg[t1.edge[0]].v2 != eg[t1.edge[1]].v2)
//			{
//				printf("");
//				if (eg[t1.edge[0]].v1 != eg[t1.edge[1]].v1 && eg[t1.edge[0]].v1 != eg[t1.edge[1]].v2)
//					printf("");
//			}
//
//
//			if (eg[t1.edge[1]].v2 != eg[t1.edge[2]].v1 && eg[t1.edge[1]].v2 != eg[t1.edge[2]].v2)
//			{
//
//				printf("");
//				if (eg[t1.edge[1]].v1 != eg[t1.edge[2]].v1 && eg[t1.edge[1]].v1 != eg[t1.edge[2]].v2)
//					printf("");
//			}
//
//
//			if (eg[t1.edge[2]].v2 != eg[t1.edge[0]].v1 && eg[t1.edge[2]].v2 != eg[t1.edge[0]].v2)
//			{
//				printf("");
//				if (eg[t1.edge[2]].v1 != eg[t1.edge[0]].v1 && eg[t1.edge[2]].v1 != eg[t1.edge[0]].v2)
//					printf("");
//			}
//
//
//			//====================================================================================
//
//			if (eg[t2.edge[0]].v2 != eg[t2.edge[1]].v1 && eg[t2.edge[0]].v2 != eg[t2.edge[1]].v2)
//			{
//
//				printf("");
//				if (eg[t2.edge[0]].v1 != eg[t2.edge[1]].v1 && eg[t2.edge[0]].v1 != eg[t2.edge[1]].v2)
//					printf("");
//			}
//
//			if (eg[t2.edge[1]].v2 != eg[t2.edge[2]].v1 && eg[t2.edge[1]].v2 != eg[t2.edge[2]].v2)
//			{
//				printf("");
//				if (eg[t2.edge[1]].v1 != eg[t2.edge[2]].v1 && eg[t2.edge[1]].v1 != eg[t2.edge[2]].v2)
//					printf("");
//			}
//
//
//			if (eg[t2.edge[2]].v2 != eg[t2.edge[0]].v1 && eg[t2.edge[2]].v2 != eg[t2.edge[0]].v2)
//			{
//				printf("");
//				if (eg[t2.edge[2]].v1 != eg[t2.edge[0]].v1 && eg[t2.edge[2]].v1 != eg[t2.edge[0]].v2)
//					printf("");
//			}
//
//
//			//======================================================================================
//
//			if (eg[t3.edge[0]].v2 != eg[t3.edge[1]].v1 && eg[t3.edge[0]].v2 != eg[t3.edge[1]].v2)
//			{
//
//				printf("");
//				if (eg[t3.edge[0]].v1 != eg[t3.edge[1]].v1 && eg[t3.edge[0]].v1 != eg[t3.edge[1]].v2)
//					printf("");
//			}
//			if (eg[t3.edge[1]].v2 != eg[t3.edge[2]].v1 && eg[t3.edge[1]].v2 != eg[t3.edge[2]].v2)
//			{
//				printf("");
//				if (eg[t3.edge[1]].v1 != eg[t3.edge[2]].v1 && eg[t3.edge[1]].v1 != eg[t3.edge[2]].v2)
//					printf("");
//			}
//
//
//			if (eg[t3.edge[2]].v2 != eg[t3.edge[0]].v1 && eg[t3.edge[2]].v2 != eg[t3.edge[0]].v2)
//			{
//				printf("");
//				if (eg[t3.edge[2]].v1 != eg[t3.edge[0]].v1 && eg[t3.edge[2]].v1 != eg[t3.edge[0]].v2)
//					printf("");
//			}
//
//
//			//===========================================================================================
//
//			if (eg[t4.edge[0]].v2 != eg[t4.edge[1]].v1 && eg[t4.edge[0]].v2 != eg[t4.edge[1]].v2)
//			{
//
//				printf("");
//				if (eg[t4.edge[0]].v1 != eg[t4.edge[1]].v1 && eg[t4.edge[0]].v1 != eg[t4.edge[1]].v2)
//					printf("");
//			}
//
//			if (eg[t4.edge[1]].v2 != eg[t4.edge[2]].v1 && eg[t4.edge[1]].v2 != eg[t4.edge[2]].v2)
//			{
//				printf("");
//				if (eg[t4.edge[1]].v1 != eg[t4.edge[2]].v1 && eg[t4.edge[1]].v1 != eg[t4.edge[2]].v2)
//					printf("");
//			}
//
//
//			if (eg[t4.edge[2]].v2 != eg[t4.edge[0]].v1 && eg[t4.edge[2]].v2 != eg[t4.edge[0]].v2)
//			{
//				printf("");
//				if (eg[t4.edge[2]].v1 != eg[t4.edge[0]].v1 && eg[t4.edge[2]].v1 != eg[t4.edge[0]].v2)
//					printf("");
//			}
//
//
//			//============================================================================================
//			continue;
//			if (eg[t1.edge[0]].v2 != eg[t1.edge[1]].v1 || eg[t1.edge[1]].v2 != eg[t1.edge[2]].v1 || eg[t1.edge[2]].v2 != eg[t1.edge[0]].v1)
//				printf("%d,%d,   %d,%d,    %d,%d \n", eg[t1.edge[0]].v2, eg[t1.edge[0]].v1, eg[t1.edge[1]].v1, eg[t1.edge[1]].v2, eg[t1.edge[2]].v2, eg[t1.edge[2]].v1);
//
//
//			if (eg[t2.edge[0]].v2 != eg[t2.edge[1]].v1 || eg[t2.edge[1]].v2 != eg[t2.edge[2]].v1 || eg[t2.edge[2]].v2 != eg[t2.edge[0]].v1)
//				printf("%d,%d,   %d,%d,    %d,%d \n", eg[t2.edge[0]].v2, eg[t2.edge[0]].v1, eg[t2.edge[1]].v1, eg[t2.edge[1]].v2, eg[t2.edge[2]].v2, eg[t2.edge[2]].v1);
//
//
//
//			if (eg[t3.edge[0]].v2 != eg[t3.edge[1]].v2 || eg[t3.edge[1]].v1 != eg[t3.edge[2]].v2 || eg[t3.edge[2]].v1 != eg[t3.edge[0]].v1)
//				printf("%d,%d,   %d,%d,    %d,%d \n", eg[t3.edge[0]].v2, eg[t3.edge[0]].v1, eg[t3.edge[1]].v1, eg[t3.edge[1]].v2, eg[t3.edge[2]].v2, eg[t3.edge[2]].v1);
//
//
//			if (eg[t4.edge[0]].v2 != eg[t4.edge[1]].v1 || eg[t4.edge[1]].v2 != eg[t4.edge[2]].v1 || eg[t4.edge[2]].v2 != eg[t4.edge[0]].v1)
//				printf("%d,%d,   %d,%d,    %d,%d \n", eg[t4.edge[0]].v2, eg[t4.edge[0]].v1, eg[t4.edge[1]].v1, eg[t4.edge[1]].v2, eg[t4.edge[2]].v2, eg[t4.edge[2]].v1);
//
//			continue;
//			if (t1.edge[0] == row && eg[t1.edge[1]].v1 == edge1.v2)
//				edge3.v2 = eg[t1.edge[1]].v2;
//			else if (t1.edge[0] == row && eg[t1.edge[1]].v2 == edge1.v2)
//				edge3.v2 = eg[t1.edge[1]].v1;
//
//
//			if (t1.edge[1] == row && eg[t1.edge[2]].v1 == edge1.v2)
//				edge3.v2 = eg[t1.edge[2]].v2;
//			else if (t1.edge[1] == row && eg[t1.edge[2]].v2 == edge1.v2)
//				edge3.v2 = eg[t1.edge[2]].v1;
//
//			if (t1.edge[2] == row && eg[t1.edge[0]].v1 == edge1.v2)
//				edge3.v2 = eg[t1.edge[0]].v2;
//			else if (t1.edge[2] == row && eg[t1.edge[0]].v2 == edge2.v2)
//				edge3.v2 = eg[t1.edge[1]].v1;
//
//
//			edge1.v2 = curVertexIdx;
//			edge2.v1 = curVertexIdx;
//
//			edge3.v1 = curVertexIdx;
//
//			t3.edge[0] = curVertexIdx;
//			if (t1.edge[0] == row)
//			{
//				if (eg[t1.edge[1]].v2 == edge3.v2)
//					t1.edge[1] = curEdgeIdx + 1;
//
//				if (eg[t3.edge[1]].v1 == edge3.v2)
//					t3.edge[2] = curEdgeIdx + 1;
//
//
//			}
//
//			if (t1.edge[1] == row)
//			{
//				if (eg[t1.edge[2]].v1 == curVertexIdx)
//					t1.edge[1] = curEdgeIdx + 1;
//
//				if (eg[t3.edge[1]].v2 == edge3.v2)
//					t1.edge[2] = curEdgeIdx + 1;
//
//
//			}
//
//			if (t1.edge[2] == row)
//			{
//				if (eg[t1.edge[2]].v1 == curVertexIdx)
//					t1.edge[1] = curEdgeIdx + 1;
//
//				if (eg[t3.edge[1]].v2 == edge3.v2)
//					t1.edge[2] = curEdgeIdx + 1;
//
//
//			}
//
//
//
//			if (t1.edge[0] == row)
//			{
//				t3.edge[0] = curEdgeIdx;
//				t3.edge[2] = curEdgeIdx + 1;
//
//				if (eg[t1.edge[2]].v1 == edge1.v1)
//					t3.edge[1] = t1.edge[1];
//				else
//					t3.edge[1] = t1.edge[2];
//
//
//				edge3.v2 = eg[t1.edge[1]].v2;
//				t1.edge[1] = curEdgeIdx + 1;
//			}
//
//			else if (t1.edge[1] == row)
//			{
//				t3.edge[0] = curEdgeIdx;
//				if (eg[t1.edge[2]].v1 == edge1.v1)
//					t3.edge[1] = t1.edge[0];
//				else
//					t3.edge[1] = t1.edge[2];
//
//				t3.edge[2] = curEdgeIdx + 1;
//
//				edge3.v2 = eg[t1.edge[2]].v2;
//				t1.edge[2] = curEdgeIdx + 1;
//
//			}
//
//			else if (t1.edge[2] == row)
//			{
//				t3.edge[0] = curEdgeIdx;
//				//t3.edge[1] = t1.edge[0];
//				t3.edge[2] = curEdgeIdx + 1;
//				if (eg[t1.edge[1]].v1 == edge1.v2)
//					t3.edge[1] = t1.edge[0];
//				else
//					t3.edge[1] = t1.edge[1];
//
//				edge3.v2 = eg[t1.edge[0]].v2;
//				t1.edge[1] = curEdgeIdx + 1;
//
//			}
//
//			edge4.v1 = curVertexIdx;
//
//			if (t2.edge[0] == row)
//			{
//				t4.edge[0] = curEdgeIdx;
//				t4.edge[2] = curEdgeIdx + 2;
//
//				if (eg[t2.edge[2]].v1 == edge1.v1)
//					t4.edge[1] = t2.edge[1];
//				else
//					t4.edge[1] = t2.edge[0];
//				edge4.v2 = eg[t2.edge[1]].v2;
//				t2.edge[1] = curEdgeIdx + 2;
//			}
//			else if (t2.edge[1] == row)
//			{
//				t4.edge[0] = curEdgeIdx;
//				t4.edge[1] = t2.edge[2];
//				t4.edge[2] = curEdgeIdx + 2;
//
//				if (eg[t4.edge[2]].v1 == edge2.v1)
//					t4.edge[1] = t2.edge[0];
//				else
//					t4.edge[1] = t2.edge[2];
//
//				edge4.v2 = eg[t2.edge[2]].v2;
//				t2.edge[2] = curEdgeIdx + 2;
//
//			}
//			else if (t2.edge[2] == row)
//			{
//				t4.edge[0] = curEdgeIdx;
//				t4.edge[1] = t2.edge[0];
//				t4.edge[2] = curEdgeIdx + 2;
//
//				if (eg[t2.edge[1]].v2 == edge2.v2 || eg[t2.edge[1]].v2 == edge2.v1)
//					t4.edge[1] = t2.edge[1];
//				else
//					t4.edge[1] = t2.edge[0];
//
//				edge4.v2 = eg[t2.edge[0]].v2;
//				t2.edge[1] = curEdgeIdx + 2;
//
//			}
//
//
//			edge2.E2T[0] = curtriangleIdx;
//			edge2.E2T[1] = curtriangleIdx + 1;
//
//			edge3.E2T[0] = triangleIdx1;;
//			edge3.E2T[1] = curtriangleIdx;
//
//			edge4.E2T[0] = triangleIdx2;;
//			edge4.E2T[1] = curtriangleIdx + 1;
//			eg[row] = edge1;
//			eg[curEdgeIdx] = edge2;
//			eg[curEdgeIdx + 1] = edge3;
//			eg[curEdgeIdx + 2] = edge4;
//
//			tr[triangleIdx1] = t1;
//			tr[triangleIdx2] = t2;
//
//			tr[curtriangleIdx] = t3;
//			tr[curtriangleIdx + 1] = t4;
//
//		}
//		continue;
//		co = 0;
//		//fTriangle _t = tr[row];
//
//		int curidx = num_rows;
//		if (vertex1.getOldRange() == 6671)
//			printf("");
//		if (vertex1.isInNextLevel_xy())
//		{
//			ASF_vertex vertex2 = m[vertex1.fx[1]];
//			ASF_vertex vedge = DivideEdges(&vertex1, &vertex2, vertex1.getOldRange(), 0, bForward);
//			curidx = num_rows + (Fr[row - 1] + co);
//			vedge.fy[1] = vertex1.fy[1];
//
//			if (vedge.getOldRange() != vertex1.getOldRange())
//				printf("");
//			vedge.fx[1] = vertex1.fx[1];
//			vedge.fcxy[1] = vertex1.fcxy[1];
//			vertex1.fx[1] = curidx;
//
//			m[curidx] = vedge;
//			vertex1.unsetInFace_xy();
//			co++;
//		}
//
//		if (curidx == 65549)
//			printf("");
//
//		if (vertex1.isInNextLevel_yz())
//		{
//			ASF_vertex vertex2 = m[vertex1.fy[1]];
//			ASF_vertex vedge = DivideEdges(&vertex1, &vertex2, vertex1.getOldRange(), 0, bForward);
//			if (vedge.getOldRange() != vertex1.getOldRange())
//				printf("");
//			curidx = num_rows + (Fr[row - 1] + co);
//			vedge.fx[1] = vertex1.fx[1];
//
//			vedge.fy[1] = vertex1.fy[1];
//			vedge.fcxy[1] = vertex1.fcxy[1];
//
//			vertex1.fy[1] = curidx;
//			m[curidx] = vedge;
//			vertex1.unsetInFace_yz();
//			co++;
//
//		}
//
//
//		if (vertex1.isInNextLevel_xz())
//		{
//			ASF_vertex vertex2 = m[vertex1.fcxy[1]];
//			ASF_vertex vedge = DivideEdges(&vertex1, &vertex2, vertex1.getOldRange(), 0, bForward);
//			curidx = num_rows + (Fr[row - 1] + co);
//			vedge.fcxy[1] = vertex1.fcxy[1];
//
//			if (vedge.getOldRange() != vertex1.getOldRange())
//				printf("");
//			vedge.fy[1] = vertex1.fy[1];
//			vedge.fx[1] = vertex1.fx[1];
//			vertex1.fcxy[1] = curidx;
//			m[curidx] = vedge;
//			vertex1.unsetInFace_xz();
//			co++;
//
//		}
//
//
//		////m[row] = vertex1;
//		//if ((Fr[row] - Fr[row - 1]) < 0)
//		//	printf("");
//		//else
//		m[row] = vertex1;
//	}
//
//}



//=====================================================================================================
__global__ void CheckNeighborhood_c3(fEdge* eg, ASF_vertex*m, uint32_t* Fr, Dimension* d, Boundary* b, Point* step, int i, bool bForward, uint32_t num_rows)
{
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;


	//uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
	//	for (uint32_t row = 0; row < num_rows; row++)

	{

		if (row >= num_rows)
			return;

		Fr[row] = 0;


		int co = 0;
		fEdge edge = eg[row];
		// if (row == 0)


		if (!edge.isInBoundary())
			return;
//		edge.bsplit = false;
		//  printf("%d--", row);
		ASF_vertex vertex1 = m[edge.v1];
		ASF_vertex vertex2 = m[edge.v2];
		if (!vertex1.isInBoundary() || !vertex2.isInBoundary())
		{
			edge.unsetInBoundary();
			eg[row] = edge;

			return;

		}
		/*if (abs(vertex1.p.x - vertex2.p.x) > step->x + exp(-6.0) || abs(vertex1.p.y - vertex2.p.y) > step->y + exp(-6.0) || abs(vertex1.p.z - vertex2.p.z) > step->z + exp(-6.0))
		{

		continue;
		}*/

		if (abs((int)vertex1.getOldRange() - (int)vertex2.getOldRange()) > 1 && abs((int)vertex1.getOldRange() - (int)vertex2.getOldRange()) != d->x && abs((int)vertex1.getOldRange() - (int)vertex2.getOldRange()) != (d->x*d->y))
			printf("");


		float dist = sqrt((vertex1.e.x - vertex2.e.x)*(vertex1.e.x - vertex2.e.x) + (vertex1.e.y - vertex2.e.y)*(vertex1.e.y - vertex2.e.y) + (vertex1.e.z - vertex2.e.z)*(vertex1.e.z - vertex2.e.z));

		//if (abs(vertex1.p.x - vertex2.p.x)  > step->x + exp(-6.0) || abs(vertex1.p.y - vertex2.p.y)  > step->y + exp(-6.0) || abs(vertex1.p.z - vertex2.p.z)  > step->z + exp(-6.0))
		//	printf(" %d-- ", row);
		//// printf("");

		//if (dist > 2 * step->x && (row == 98304 || row == 13778 || row == 98328 || row == 98327))
		//	printf(" %d-- ", row);
		// printf("");
		// printf(" %d-- ", row);

		if (!checkEdge(vertex1, vertex2, d, bForward))
		{
			/*if (row == 98312)
			printf("\n v1 = %d,%d --- range = %d,%d \n",edge.v1,edge.v2, vertex1.getRange(), vertex2.getRange());*/
			Fr[row] = 1;
			/*if (i == 1)
			printf(" %d-- ", row);*/
		}

		eg[row] = edge;



	}





	//if (row > num_rows)
	//	return;
	//fEdge edge = eg[row];
	//Fr[row] = 0;
	//if (!eg[row].isInBoundary())
	//	return;
	//	ASF_vertex vertex1 = m[edge.v1];
	//	ASF_vertex vertex2 = m[edge.v2];
	//	if (!vertex1.isInBoundary() || !vertex2.isInBoundary())
	//	{
	//		edge.unsetInBoundary();
	//		eg[row] = edge;

	//		return;

	//	}
	//	//printf("%d,%d -- ", d->x, d->y);
	//	if (!checkEdge(vertex1, vertex2, d, bForward))
	//	{

	//		Fr[row] = 1;
	//		//printf(" %d-- ", row);
	//	}
	//	//Fr[row] = 10;

	//	//printf("%d -- ", row);

}


//__global__ void CheckNeighborhood_c3(ASF_vertex*m, fEdge* eg, uint32_t* Fr, Dimension* d, Boundary* b, Point* step, bool bForward, uint32_t num_rows)
//{
//
//
//	 uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
////	for (uint32_t row = 0; row < num_rows; row++)
//	
//	{
//
//		 if (row >= num_rows)
//			 return;
//
//		 Fr[row] = 0;
//
//
//		 int co = 0;
//		 fEdge edge = eg[row];
//		// if (row == 0)
//		 if (row>32768)
//		 printf("%d--", row);
//
//		 if (!edge.isInBoundary())
//			 return;
//		 edge.bsplit = false;
//		//  printf("%d--", row);
//		 ASF_vertex vertex1 = m[edge.v1];
//		 ASF_vertex vertex2 = m[edge.v2];
//		 if (!vertex1.isInBoundary() || !vertex2.isInBoundary())
//		 {
//			 edge.unsetInBoundary();
//			 eg[row] = edge;
//
//			 return;
//
//		 }
//		 /*if (abs(vertex1.p.x - vertex2.p.x) > step->x + exp(-6.0) || abs(vertex1.p.y - vertex2.p.y) > step->y + exp(-6.0) || abs(vertex1.p.z - vertex2.p.z) > step->z + exp(-6.0))
//		 {
//
//		 continue;
//		 }*/
//
//		 if (abs((int)vertex1.getOldRange() - (int)vertex2.getOldRange()) > 1 && abs((int)vertex1.getOldRange() - (int)vertex2.getOldRange()) != d->x && abs((int)vertex1.getOldRange() - (int)vertex2.getOldRange()) != (d->x*d->y))
//			 printf("");
//		
//		 
//		 float dist = sqrt((vertex1.e.x - vertex2.e.x)*(vertex1.e.x - vertex2.e.x) + (vertex1.e.y - vertex2.e.y)*(vertex1.e.y - vertex2.e.y) + (vertex1.e.z - vertex2.e.z)*(vertex1.e.z - vertex2.e.z));
//
//		 if (abs(vertex1.p.x - vertex2.p.x)  > step->x + exp(-6.0) || abs(vertex1.p.y - vertex2.p.y)  > step->y + exp(-6.0) || abs(vertex1.p.z - vertex2.p.z)  > step->z + exp(-6.0))
//			 printf(" %d-- ", row);
//		 // printf("");
//
//		 if (dist > 2 * step->x && (row == 98304 || row == 13778 || row == 98328 || row == 98327))
//			 printf(" %d-- ", row);
//		 // printf("");
//		// printf(" %d-- ", row);
//
//		 if (!checkEdge(vertex1, vertex2, d, bForward))
//		 {
//
//			 Fr[row] = 1;
//			 printf(" %d-- ", row);
//		 }
//
//		 eg[row] = edge;
//
//
//
//	}
//	
//}


__global__ void EdgeReduction(fEdge*eg, ASF_vertex*m, fFace* Ff, uint32_t* Fe_Face, uint32_t* Fe_Edge, uint32_t num_edges)
{


	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;

	{

		if (row >= num_edges)
			return;

		if (Fe_Edge[row] == 0)
		{

			fEdge edge = eg[row];
			if (!edge.isInBoundary())
				return;
			//face.bSplit = false;
			for (int i = 0; i < 4; i++)
			if (Fe_Face[edge.E2F[i]] == 1 && edge.level <= Ff[edge.E2F[i]].level)// && edge.length >= Ff[edge.E2F[i]].dx2)
			{
				if (row > 98304)
					printf("");
				Fe_Edge[row] = 1;
				break;
			}
		}

	}

}


__global__ void CheckFace_c3(fFace* fc, ASF_vertex*m, fEdge* eg, uint32_t*Fe_Edge, uint32_t* Fe_Face, uint32_t num_face)
{


	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
	//	for (uint32_t row = 0; row < num_rows; row++)

	{
		if (row >= num_face)
			return;
		if (row >= num_face)
			return;
		if (row == 98317)
			printf("");


		fFace face = fc[row];
//		face.bSplit = false;





		//}

		if (!face.isInBoundary())
			return;
		if (!eg[face.edge[0]].isInBoundary() || !eg[face.edge[1]].isInBoundary() || !eg[face.edge[2]].isInBoundary() || !eg[face.edge[3]].isInBoundary())
			return;
		if (Fe_Edge[face.edge[0]] == 1 || Fe_Edge[face.edge[1]] == 1 || Fe_Edge[face.edge[2]] == 1 || Fe_Edge[face.edge[3]] == 1)
		{
			if (!m[eg[face.edge[0]].v1].isInBoundary() || !m[eg[face.edge[0]].v2].isInBoundary() || !m[eg[face.edge[1]].v2].isInBoundary() || !m[eg[face.edge[3]].v2].isInBoundary())

				return;



			Fe_Face[row] = 1;
			fc[row].level = fc[row].level + 1;
			//face.bSplit = true;
		}
		fc[row] = face;

	}




	return;

}




__global__ void SplitEdge(fEdge* eg, ASF_vertex*m, fFace* fc, uint32_t* Fe_Edge, uint32_t* Fr_Edge, Dimension* d, Boundary* b, Point* step, bool bForward, uint32_t num_vertex, uint32_t num_edges)
{


	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
	//	for (uint32_t row = 0; row < num_rows; row++)

	{

		if (row > num_edges)
			return;

		if (Fe_Edge[row] == 1)
		{
			int curVertexId = num_vertex + Fr_Edge[row];
			int curEdgeId = num_edges + Fr_Edge[row];

			int co = 0;



			fEdge edge1 = eg[row];
			fEdge edge2 = eg[row];
			ASF_vertex vertex1 = m[edge1.v1];
			ASF_vertex vertex2 = m[edge1.v2];


			if ((vertex1.p.x - vertex2.p.x) > step->x || (vertex1.p.y - vertex2.p.y) > step->y || (vertex1.p.z - vertex2.p.z) > step->z)
				printf("");
			ASF_vertex* vedge = DivideEdges(&vertex1, &vertex2, vertex1.getOldRange(), d, step, b, 0, bForward);

			//return;;

	//		edge1.length = edge1.length / 2;
			edge1.v2 = curVertexId;
			edge2.v1 = curVertexId;
//			edge1.bsplit = true;
			edge1.next = curEdgeId;
			edge1.level = edge1.level * 2;
			edge2.level = edge2.level * 2;
			edge2.Prev = row;
			vedge->level = edge2.level;

			float rgb[3];
			rgb[0] = 0.;
			rgb[1] = 0.;
			rgb[2] = 0.;

			//	display_voxel(vertex1.getOldRange(), rgb, d->x, d->y, d->z);

			eg[row] = edge1;
			eg[curEdgeId] = edge2;;
			m[curVertexId] = *vedge;
			//printf("%d,%d,%d range  = %d,%d \n", edge1.v1, edge1.v2, edge2.v2, vertex1.getRange(), vertex2.getRange());



			/*	if (CheckEdgeCondition(m, edge1) == false)
			printf("");*/
			/*if (m[edge1.v1].getOldRange() == 1551)
			drawLine(edge1.v1, edge1.v2, rgb);*/
			rgb[0] = 1.;
			rgb[1] = 0.;
			rgb[2] = 0.;
			/*	if (m[edge2.v1].getOldRange() == 1551)
			drawLine(edge2.v1, edge2.v2, rgb);*/

		}

	}
}

#define BLOCKSIZE 512

template<class T>
__global__ void kernelFunction1(T * __restrict__ d_data, const unsigned int NperGPU) {

	const int tid = threadIdx.x + blockIdx.x * blockDim.x;
	ASF_vertex vertex = d_data[tid];
	if (tid < NperGPU)
	{
		uint32_t rng = vertex.range;
		for (int k = 0; k < 1000; k++) vertex.range = 2 * rng;// *3.0 *d_data[tid];
		d_data[tid] = vertex;
	}


}

/*******************/
/* KERNEL FUNCTION */
/*******************/
template<class T>
__global__ void Initialize_(T * __restrict__ m, fEdge* eg, fFace* fr, Boundary* b, Dimension* d, Point* step, bool bForward, const unsigned int num_rows)
{

	const int row = threadIdx.x + blockIdx.x * blockDim.x;

	if (row >= num_rows)
		return;
	//printf(" %d --", row);
	ASF_vertex vertex = m[row];
	int xdim = d->x;
	int ydim = d->y;
	int zdim = d->z;

	//if (row + d->x < num_rows /*&& row + d->x*d->y < num_rows*/)
	{

		/*if ((row % d->x - 1)==0 || (row % ((d->x*d->y) - 1))==0)
		continue;*/

		vertex.setInBoundary();

		int curEdgeidx = row * 3;


		/*vertex.fxy = curEdgeidx;
		vertex.fyz = curEdgeidx + 1;
		vertex.fxz = curEdgeidx + 2;*/
		m[row] = vertex;
		int iz = row / (d->x*d->y);
		int iy = (row - iz*(d->x*d->y)) / (d->x);
		int ix = (row - iz*(d->x*d->y)) % (d->x);
		/*if (iz == d->z - 1 || iy == d->y - 1 || ix == d->x - 1)
		{
		eg[curEdgeidx].unsetInBoundary();
		eg[curEdgeidx + 1].unsetInBoundary();
		eg[curEdgeidx + 2].unsetInBoundary();
		fr[curEdgeidx].unsetInBoundary();
		fr[curEdgeidx + 1].unsetInBoundary();
		fr[curEdgeidx + 2].unsetInBoundary();
		continue;
		}*/

		if (curEdgeidx == 2979)
			printf("");

		if (row + 1 < num_rows && ix < d->x - 1)
		{

			ASF_vertex vertex2 = m[row + 1];
			if (abs(vertex.p.x - vertex2.p.x) > step->x + exp(-6.0) || abs(vertex.p.y - vertex2.p.y) > step->y + exp(-6.0) || abs(vertex.p.z - vertex2.p.z) > step->z + exp(-6.0))
				eg[curEdgeidx].unsetInBoundary();
			else
			{
				eg[curEdgeidx].v1 = row;
				eg[curEdgeidx].v2 = row + 1;
//				eg[curEdgeidx].bsplit = false;
				eg[curEdgeidx].E2V = row;
				eg[curEdgeidx].E2F[0] = curEdgeidx;
				eg[curEdgeidx].E2F[1] = curEdgeidx + 1;
		//		eg[curEdgeidx].length = step->x;;
				eg[curEdgeidx].level = 1;
				eg[curEdgeidx].Prev = 0;
				eg[curEdgeidx].setInBoundary();
			}



			if (curEdgeidx >= d->x * 3)
			{
				eg[curEdgeidx].E2F[2] = (row - d->x) * 3;

			}
			else
				eg[curEdgeidx].E2F[2] = OUTBOUNDARY;

			if (row >= d->x*d->y)
				eg[curEdgeidx].E2F[3] = (row - d->x*d->y) * 3 + 1;
			else
				eg[curEdgeidx].E2F[3] = OUTBOUNDARY;

		}


		//eg[curEdgeidx ].E2T[0] = (row ) * 4;

		/*if (row > d->x)
		eg[curEdgeidx ].E2T[1] = (row - d->x) * 4 + 2;
		else
		eg[curEdgeidx ].E2T[1] = 0;*/
		if (row + d->x < num_rows && iy < d->y - 1)
		{

			ASF_vertex vertex2 = m[row + d->x];
			if (abs(vertex.p.x - vertex2.p.x) > step->x + exp(-6.0) || abs(vertex.p.y - vertex2.p.y) > step->y + exp(-6.0) || abs(vertex.p.z - vertex2.p.z) > step->z + exp(-6.0))
				eg[curEdgeidx + 1].unsetInBoundary();
			else
			{
				eg[curEdgeidx + 1].v1 = row;
				eg[curEdgeidx + 1].v2 = row + d->x;
//				eg[curEdgeidx + 1].bsplit = false;
				eg[curEdgeidx + 1].E2V = row;

				eg[curEdgeidx + 1].E2F[0] = curEdgeidx;
				eg[curEdgeidx + 1].E2F[1] = curEdgeidx + 2;
//				eg[curEdgeidx + 1].length = step->y;;
				eg[curEdgeidx + 1].level = 1;
				eg[curEdgeidx + 1].Prev = 0;
				eg[curEdgeidx + 1].setInBoundary();





				if (curEdgeidx >= 3)
				{
					eg[curEdgeidx + 1].E2F[2] = (row - 1) * 3;

				}

				else
					eg[curEdgeidx + 1].E2F[2] = OUTBOUNDARY;

				if (curEdgeidx >= d->x*d->y * 3)
					eg[curEdgeidx + 1].E2F[3] = (row - d->x*d->y) * 3 + 2;
				else
					eg[curEdgeidx + 1].E2F[3] = OUTBOUNDARY;
			}
		}



		if (row + d->x*d->y < num_rows && iz < d->z - 1)
		{

			ASF_vertex vertex2 = m[row + d->x*d->y];
			if (abs(vertex.p.x - vertex2.p.x) > step->x + exp(-6.0) || abs(vertex.p.y - vertex2.p.y) > step->y + exp(-6.0) || abs(vertex.p.z - vertex2.p.z) > step->z + exp(-6.0))
				eg[curEdgeidx + 2].unsetInBoundary();
			else
			{

				eg[curEdgeidx + 2].v1 = row;
				eg[curEdgeidx + 2].v2 = row + (d->x*d->y);
//				eg[curEdgeidx + 2].bsplit = false;
				eg[curEdgeidx + 2].E2V = curEdgeidx;

				eg[curEdgeidx + 2].E2F[0] = curEdgeidx + 1;
				eg[curEdgeidx + 2].E2F[1] = curEdgeidx + 2;

//				eg[curEdgeidx + 2].length = step->z;;
				eg[curEdgeidx + 2].level = 1;
				eg[curEdgeidx + 2].Prev = 0;
				eg[curEdgeidx + 2].setInBoundary();


				if (curEdgeidx >= 3 && row%d->x != 0)
					eg[curEdgeidx + 2].E2F[2] = (row - 1) * 3 + 1;
				else
					eg[curEdgeidx + 2].E2F[2] = OUTBOUNDARY;


				if (row >= d->x)
				{
					eg[curEdgeidx + 2].E2F[3] = (row - d->x) * 3 + 2;

				}
				else
					eg[curEdgeidx + 2].E2F[3] = OUTBOUNDARY;
			}

		}
		//eg[curEdgeidx + 1].E2T[0] = (row)* 4 +3;

		/*if (row > 1)
		eg[curEdgeidx + 1].E2T[1] = (row - 1) * 4 +1;
		else
		eg[curEdgeidx + 1].E2T[1] = 0;*/



		if (eg[curEdgeidx].isInBoundary() && eg[curEdgeidx + 1].isInBoundary())
		{


			if ((row + 1) % d->x != 0 && row + (d->x) < num_rows)
			{
				fr[curEdgeidx].edge[0] = curEdgeidx;
				fr[curEdgeidx].edge[1] = (row + 1) * 3 + 1;
				fr[curEdgeidx].edge[2] = (row + d->x) * 3;
				fr[curEdgeidx].edge[3] = curEdgeidx + 1;
				fr[curEdgeidx].F2V = row;
				fr[curEdgeidx].cornerId = row + d->x + 1;
//				fr[curEdgeidx].bSplit = false;
//				fr[curEdgeidx].dx1 = step->x;
//				fr[curEdgeidx].dx2 = step->y;
				fr[curEdgeidx].setInBoundary();
				fr[curEdgeidx].level = 1;;

				if (row + (d->x*d->y) < num_rows)
				{
					fr[curEdgeidx + 1].edge[0] = curEdgeidx;
					fr[curEdgeidx + 1].edge[1] = (row + 1) * 3 + 2;
					fr[curEdgeidx + 1].edge[2] = (row + (d->x*d->y)) * 3;
					fr[curEdgeidx + 1].edge[3] = curEdgeidx + 2;
					fr[curEdgeidx + 1].cornerId = row + (d->x*d->y) + 1;
					fr[curEdgeidx + 1].F2V = row;
//					fr[curEdgeidx + 1].bSplit = false;
//					fr[curEdgeidx + 1].dx1 = step->x;
//					fr[curEdgeidx + 1].dx2 = step->z;
					fr[curEdgeidx + 1].setInBoundary();
					fr[curEdgeidx + 1].level = 1;;


				}
				else
				{
					fr[curEdgeidx + 1].unsetInBoundary();

				}
			}

		}
		else
		{
			fr[curEdgeidx].unsetInBoundary();
			fr[curEdgeidx + 1].unsetInBoundary();

		}


		if (eg[curEdgeidx + 2].isInBoundary() && eg[curEdgeidx + 1].isInBoundary() && (row == 0 || ((row) % d->x*d->y != 0 && row + (d->x*d->y) < num_rows)))
		{
			fr[curEdgeidx + 2].edge[0] = curEdgeidx + 1;
			fr[curEdgeidx + 2].edge[1] = (row + (d->x)) * 3 + 2;
			fr[curEdgeidx + 2].edge[2] = (row + (d->x*d->y)) * 3 + 1;
			fr[curEdgeidx + 2].edge[3] = curEdgeidx + 2;
			fr[curEdgeidx + 2].F2V = row;
			fr[curEdgeidx + 2].cornerId = row + (d->x*(d->y + 1));

//			fr[curEdgeidx + 2].bSplit = false;
//			fr[curEdgeidx + 2].dx1 = step->y;
//			fr[curEdgeidx + 2].dx2 = step->z;
			fr[curEdgeidx + 2].setInBoundary();
			fr[curEdgeidx + 2].level = 1;;


		}
		else
		{
			fr[curEdgeidx + 2].unsetInBoundary();
		}

		m[row] = vertex;
		/*if (curEdgeidx > num_rows && curEdgeidx < num_rows + 10)
		printf("%d,%d --", eg[curEdgeidx].v1, eg[curEdgeidx].v2);*/


	}





	//if (tid < NperGPU) for (int k = 0; k < 1000; k++) d_data[tid].range =2;// d_data[tid] * 3.0 *d_data[tid];

}
/******************/
/* PLAN STRUCTURE */
/******************/
// --- Async
template<class T>
struct plan {
	T               *d_data;
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

#define sample_seeds 50
/*********************/
/* SVD PLAN CREATION */
/*********************/
template<class T>
void createPlan(plan<T>& plan, unsigned int NperGPU, unsigned int NEdgeperGPU, unsigned int N, unsigned int gpuID) {

	// --- Device allocation
	(cudaSetDevice(gpuID));
	checkCudaErrors(cudaStreamCreate(&plan.stream));
	checkCudaErrors(cudaMalloc(&(plan.d_data), NperGPU * sample_seeds* sizeof(T)));
	checkCudaErrors(cudaMalloc(&(plan.v), N * sizeof(Point)));
	checkCudaErrors(cudaMalloc(&(plan.Fe_Edge), NEdgeperGPU * sample_seeds * sizeof(uint32_t)));
	checkCudaErrors(cudaMalloc(&(plan.Fr_Edge), NEdgeperGPU * sample_seeds * sizeof(uint32_t)));

	checkCudaErrors(cudaMalloc(&(plan.Fe_Face), NEdgeperGPU * sample_seeds * sizeof(uint32_t)));
	checkCudaErrors(cudaMalloc(&(plan.Fr_Face), NEdgeperGPU * sample_seeds * sizeof(uint32_t)));

	checkCudaErrors(cudaMalloc(&(plan.eg), NEdgeperGPU * sample_seeds * sizeof(fEdge)));
	checkCudaErrors(cudaMalloc(&(plan.fc), NEdgeperGPU * sample_seeds * sizeof(fFace)));
	checkCudaErrors(cudaMalloc(&(plan.b), 1 * sizeof(Boundary)));
	checkCudaErrors(cudaMalloc(&(plan.d), 1 * sizeof(Dimension)));
	checkCudaErrors(cudaMalloc(&(plan.s), 1 * sizeof(Point)));

}

template<class T>
__global__ void kernelFunction(T * __restrict__ d_data, int i, const unsigned int NperGPU) {

	const int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if (tid < NperGPU) /*for (int k = 0; k < 10; k++)*/
	{
		ASF_vertex d = d_data[tid];


	}

}

//template<class T>
//__global__ void Tracing(T * __restrict__ d_data, Point*v, Boundary* b, Dimension* _d, Point* step, bool bForward,int i, uint32_t NperGPU) {
//
//	const int tid = threadIdx.x + blockIdx.x * blockDim.x;
//
//	/*if (tid == 1)
//		printf(" num rows = %d  ", i);
//*/
//	if (tid < NperGPU) /*for (int k = 0; k < 10; k++)*/
//	{
//	
//
//		ASF_vertex vertex = d_data[tid];
//
//		if (!vertex.isInBoundary())
//		{
//			return;
//
//		}
//
//		/*if (i == 0)
//			printf("%f,%f \n",v[32000].x, v[32000].y);
//		
//		d_data[tid].e.x =  v[vertex.getOldRange()].x;
//		d_data[tid].e.y = vertex.e.y + 0.1;
//		d_data[tid].e.z = vertex.e.z + 0.1;
//
//		return;*/
//		Point e = vertex.e;
//		Point _v = trilinearInterpolation(e, v, vertex.getOldRange(), _d, b, step, bForward);
//		/*vertex.e.x = b->high.x;
//		vertex.e.y = b->high.y;
//		vertex.e.z = b->high.z;
//		if (row < 10)
//		printf("%f,%f \n", _v.x, _v.y);*/
//
//		float dist = _v.getDist();
//		_v.x = (_v.x / (dist*4.0));
//		_v.y = (_v.y / (dist*4.0));
//		_v.z = (_v.z / (dist*4.0));
//
//		_v.x = _v.x*step->x;
//		_v.y = _v.y*step->y;
//		_v.z = _v.z*step->z;
//
//		if (bForward)
//		{
//			vertex.e.x = vertex.e.x + _v.x;
//			vertex.e.y = vertex.e.y + _v.y;
//			vertex.e.z = vertex.e.z + _v.z;
//
//		}
//		else
//		{
//			vertex.eb.x = vertex.eb.x - _v.x;
//			vertex.eb.y = vertex.eb.y - _v.y;
//			vertex.eb.z = vertex.eb.z - _v.z;
//		}
//		if (tid == 4591 )
//			printf(" \n i = %d, %d, %f,%f,%f \n",i,tid, vertex.e.x, vertex.e.y, vertex.e.z);
//	
//		
//
//		bool xy = false;
//		bool yz = false;
//		bool xz = false;
//
//
//		if (bForward)
//		{
//			if (vertex.checkInBoundary(b))
//			{
//				uint32_t range = vertex.getRange(b, step, _d);
//			/*	if (tid == 4591 || tid == 4592)
//					printf("%d raaaaaaaaaaaaangeeeeeeeeeeeee\n", range);*/
//				/*if (vertex.isInNextLevel_yz())
//					yz = true;
//				if (vertex.isInNextLevel_xy())
//					xy = true;
//				if (vertex.isInNextLevel_xz())
//					xz = true;*/
//
//				if (bForward)
//					vertex.setRange(range);
//
//
//				/*if (xy)
//					vertex.setInNextLevel_xy();
//				if (xz)
//					vertex.setInNextLevel_xz();
//				if (yz)
//					vertex.setInNextLevel_yz();*/
//
//			}
//			else
//			{
//				vertex.unsetInBoundary();
//				vertex.unsetInFace_xy();
//				vertex.unsetInFace_yz();
//				vertex.unsetInFace_xz();
//			}
//		}
//		else
//		{
//			if (vertex.checkInBoundaryBackward(b))
//			{
//				uint32_t range = vertex.getRangeBackward(b, step, _d);
//				vertex.setRangeBackward(range);
//
//
//				/*if (xy)
//					vertex.setInNextLevel_xy();
//				if (xz)
//					vertex.setInNextLevel_xz();
//				if (yz)
//					vertex.setInNextLevel_yz();*/
//
//			}
//			else
//			{
//				vertex.unsetInBoundary();
//			/*	vertex.unsetInFace_xy();
//				vertex.unsetInFace_yz();
//				vertex.unsetInFace_xz();*/
//			}
//		}
//		
//		d_data[tid] = vertex;
//	}
//
//}


int isPowerOfTwo(unsigned int x, unsigned int& ipower)
{
	unsigned int powerOfTwo = 1;
	ipower = 0;
	while (powerOfTwo < x && powerOfTwo < 2147483648)
	{
		powerOfTwo *= 2;
		ipower++;

	}
	return (x == powerOfTwo);
}


void scanExclusiveHost(
	uint *dst,
	uint *src,
	uint batchSize,
	uint arrayLength
	)
{
	//for (uint i = 0; i < batchSize; i++, src += arrayLength, dst += arrayLength)
	{
		dst[0] = 0;

		for (uint j = 1; j < arrayLength; j++)
			dst[j] = src[j - 1] + dst[j - 1];
	}
}



__global__ void Tracing(ASF_vertex * m, Point*v, Boundary* b, Dimension* d, Point* step, bool bForward, uint32_t i, uint32_t num_rows)

{
	uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;

	if (row > num_rows)
		return;

	ASF_vertex vertex = m[row];

	if (!vertex.isInBoundary())
	{
		return;

	}

	Point e = vertex.e;
	//if (!bForward)
	//	e = vertex.eb;
	float vx, vy, vz;
	trilinearInterpolation(m, e, v, vertex.getOldRange(), d, b, step, vx, vy, vz, bForward);

	float dist = sqrt(vx*vx + vy*vy + vz*vz);

	vx = (vx / (dist*4.0))*step->x;
	vy = (vy / (dist*4.0))*step->y;
	vz = (vz / (dist*4.0))*step->z;
	/*_v = _v / 4.0;
	_v *= *step;*/
	//glColor3f(rgb[0], rgb[1], rgb[2]);    
	//	Point e = vertex.e;
//	if (bForward)
	{
		vertex.e.x += vx;
		vertex.e.y += vy;
		vertex.e.z += vz;
	}
	//vertex.e += _v;
	//else
	//{
	//	vertex.eb.x -= vx;
	//	vertex.eb.y -= vy;
	//	vertex.eb.z -= vz;
	//}
	//vertex.eb -= _v;
	/*if (row == 4591)
	printf("%f,%f,%f  \n", vertex.e.x, vertex.e.y, vertex.e.z);*/
	float ep[3];
	//generalstreamlineTracing_single(p1, bForward, ep, false);





	bool xy = false;
	bool yz = false;
	bool xz = false;
//	if (bForward)
	{
		if (vertex.checkInBoundary(b))
		{
			uint32_t range = vertex.getRange(b, step, d);
			if (vertex.isInNextLevel_yz())
				yz = true;
			if (vertex.isInNextLevel_xy())xy = true;
			if (vertex.isInNextLevel_xz())xz = true;
			vertex.setRange(range);
			if (xy)	vertex.setInNextLevel_xy();
			if (xz)vertex.setInNextLevel_xz();
			if (yz)vertex.setInNextLevel_yz();

		}
		else
		{
			vertex.unsetInBoundary();
		}
	}
	/*else
	{
		if (vertex.checkInBoundaryBackward(b))
		{
			if (vertex.isInNextLevel_yz())
				yz = true;
			if (vertex.isInNextLevel_xy())xy = true;
			if (vertex.isInNextLevel_xz())xz = true;
			uint32_t range = vertex.getRangeBackward(b, step, d);
			vertex.setRangeBackward(range);
			if (xy)	vertex.setInNextLevel_xy();
			if (xz)vertex.setInNextLevel_xz();
			if (yz)vertex.setInNextLevel_yz();

		}
		else
		{
			vertex.unsetInBoundary();
		}
	}
*/



	m[row] = vertex;


	return;
}



