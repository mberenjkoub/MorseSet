/*
 * ADPAlgorithm2.cpp
 *
 *  Created on: Dec 24, 2016
 *      Author: marzieh
 */

#include "ADPAlgorithm2.h"

namespace std {

ADPAlgorithm2::ADPAlgorithm2() {
	// TODO Auto-generated constructor stub
	fs = new File_Saver();
	IntegrationStep = new Point();
	numJ = 0;
	minvx = 0.;
	minvy = 0.;
	maxvx = 0.;
	maxvy = 0.;

}

ADPAlgorithm2::~ADPAlgorithm2() {
	// TODO Auto-generated destructor stub
}

void ADPAlgorithm2::SetIntegrationStep(Point _is) {
	*IntegrationStep = _is;

}

void ADPAlgorithm2::Initialize(ASF_vertex * m, fEdge** _eg, fFace** _fr,
		Dimension* d, Boundary* b, Point* step, int whichDataset, int currSize,
		int curEdgeSize, int curFaceSize, int sampleSeeds, uint32_t num_rows) {

	int lcuredgeid = 0;
	int lcurfaceid = 0;

	//	if (whichDataset ==)
	//		sampleSeeds = 50;
	//	else

	Trace = new Point*[num_rows * 20 * sampleSeeds];

	for (int j = 0; j < num_rows * sampleSeeds; j++) {
		//if(_a[j].isInBoundary())
		Trace[j] = new Point[300 + 1];
	}

	fEdge* eg = new fEdge[curEdgeSize * sampleSeeds];
	memset(eg, 0, curEdgeSize * sampleSeeds * sizeof(fEdge));

	fFace* fr = new fFace[curFaceSize * sampleSeeds];
	memset(fr, 0, curFaceSize * sampleSeeds * sizeof(fFace));

	//uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
	for (uint32_t row = 0; row < num_rows; row++) {

		if (row >= num_rows)
			return;

		ASF_vertex vertex = m[row];

		int xdim = d->x;
		int ydim = d->y;
		int zdim = d->z;
		lcuredgeid = row * 3;
		lcurfaceid = row * 3;
		m[row] = vertex;

		//		if (vertex.p.x == b->low.x || vertex.p.x == b->high.x
		//				|| vertex.p.y == b->low.y || vertex.p.y == b->high.y
		//				|| vertex.p.z == b->low.z || vertex.p.z == b->high.z)
		//			vertex.unsetInBoundary();
		//		else

		{

			vertex.setInBoundary();

			vertex.type = 1;
			int iz = row / (d->x * d->y);
			int iy = (row - iz * (d->x * d->y)) / (d->x);
			int ix = (row - iz * (d->x * d->y)) % (d->x);

			if (row + 1 < num_rows && ix < d->x - 1) {

				ASF_vertex vertex2 = m[row + 1];
				if (fabs(vertex.p.x - vertex2.p.x) > step->x + exp(-6.0)
						|| fabs(vertex.p.y - vertex2.p.y) > step->y + exp(-6.0)
						|| fabs(vertex.p.z - vertex2.p.z) > step->z + exp(-6.0))
					eg[lcuredgeid].unsetInBoundary();
				else {
					eg[lcuredgeid].v1 = row;
					eg[lcuredgeid].v2 = row + 1;
					//eg[curEdgeidx].bsplit = false;
					eg[lcuredgeid].E2V = row;
					eg[lcuredgeid].E2F[0] = lcuredgeid;
					eg[lcuredgeid].E2F[1] = lcuredgeid + 1;
					//eg[curEdgeidx].length = step->x;;
					eg[lcuredgeid].level = 1;
					eg[lcuredgeid].Prev = 0;
					eg[lcuredgeid].setInBoundary();
				}

				if (lcuredgeid >= d->x * 3) {
					eg[lcuredgeid].E2F[2] = (row - d->x) * 3;

				} else
					eg[lcuredgeid].E2F[2] = OUTBOUNDARY;

				if (row >= d->x * d->y)
					eg[lcuredgeid].E2F[3] = (row - d->x * d->y) * 3 + 1;
				else
					eg[lcuredgeid].E2F[3] = OUTBOUNDARY;

			}

			//eg[curEdgeidx ].E2T[0] = (row ) * 4;

			/*if (row > d->x)
			 eg[curEdgeidx ].E2T[1] = (row - d->x) * 4 + 2;
			 else
			 eg[curEdgeidx ].E2T[1] = 0;*/
			if (row + d->x < num_rows && iy < d->y - 1) {

				ASF_vertex vertex2 = m[row + d->x];
				if (fabs(vertex.p.x - vertex2.p.x) > step->x + exp(-6.0)
						|| fabs(vertex.p.y - vertex2.p.y) > step->y + exp(-6.0)
						|| fabs(vertex.p.z - vertex2.p.z) > step->z + exp(-6.0))
					eg[lcuredgeid + 1].unsetInBoundary();
				else {
					eg[lcuredgeid + 1].v1 = row;
					eg[lcuredgeid + 1].v2 = row + d->x;
					//eg[curEdgeidx + 1].bsplit = false;
					eg[lcuredgeid + 1].E2V = row;

					eg[lcuredgeid + 1].E2F[0] = lcuredgeid;
					eg[lcuredgeid + 1].E2F[1] = lcuredgeid + 2;
					//eg[curEdgeidx + 1].length = step->y;;
					eg[lcuredgeid + 1].level = 1;
					eg[lcuredgeid + 1].Prev = 0;
					eg[lcuredgeid + 1].setInBoundary();

					if (lcuredgeid >= 3) {
						eg[lcuredgeid + 1].E2F[2] = (row - 1) * 3;

					}

					else
						eg[lcuredgeid + 1].E2F[2] = OUTBOUNDARY;

					if (lcuredgeid >= d->x * d->y * 3)
						eg[lcuredgeid + 1].E2F[3] = (row - d->x * d->y) * 3 + 2;
					else
						eg[lcuredgeid + 1].E2F[3] = OUTBOUNDARY;
				}
			}

			if (row + d->x * d->y < num_rows && iz < d->z - 1) {

				ASF_vertex vertex2 = m[row + d->x * d->y];
				if (fabs(vertex.p.x - vertex2.p.x) > step->x + exp(-6.0)
						|| fabs(vertex.p.y - vertex2.p.y) > step->y + exp(-6.0)
						|| fabs(vertex.p.z - vertex2.p.z) > step->z + exp(-6.0))
					eg[lcuredgeid + 2].unsetInBoundary();
				else {

					eg[lcuredgeid + 2].v1 = row;
					eg[lcuredgeid + 2].v2 = row + (d->x * d->y);
					//eg[curEdgeidx + 2].bsplit = false;
					eg[lcuredgeid + 2].E2V = lcuredgeid;

					eg[lcuredgeid + 2].E2F[0] = lcuredgeid + 1;
					eg[lcuredgeid + 2].E2F[1] = lcuredgeid + 2;

					//eg[curEdgeidx + 2].length = step->z;;
					eg[lcuredgeid + 2].level = 1;
					eg[lcuredgeid + 2].Prev = 0;
					eg[lcuredgeid + 2].setInBoundary();

					if (lcuredgeid >= 3 && row % d->x != 0)
						eg[lcuredgeid + 2].E2F[2] = (row - 1) * 3 + 1;
					else
						eg[lcuredgeid + 2].E2F[2] = OUTBOUNDARY;

					if (row >= d->x) {
						eg[lcuredgeid + 2].E2F[3] = (row - d->x) * 3 + 2;

					} else
						eg[lcuredgeid + 2].E2F[3] = OUTBOUNDARY;
				}

			}
			//eg[curEdgeidx + 1].E2T[0] = (row)* 4 +3;

			/*if (row > 1)
			 eg[curEdgeidx + 1].E2T[1] = (row - 1) * 4 +1;
			 else
			 eg[curEdgeidx + 1].E2T[1] = 0;*/

			if (eg[lcuredgeid].isInBoundary()
					&& eg[lcuredgeid + 1].isInBoundary()) {

				if ((row + 1) % d->x != 0 && row + (d->x) < num_rows) {

//					fr[lcuredgeid].eg[0].v1 = row;
//					fr[lcuredgeid].eg[0].v2 = row + 1;
//					fr[lcuredgeid].eg[0].next = 0;
//					fr[lcuredgeid].eg[0].level = 1;
//					fr[lcuredgeid].eg[0].E2V = row;
//					fr[lcuredgeid].eg[0].bsplit = false;
//
//					fr[lcuredgeid].eg[1].v1 = row + 1;
//					fr[lcuredgeid].eg[1].v2 = row + 1 + d->x;
//					fr[lcuredgeid].eg[1].next = 0;
//					fr[lcuredgeid].eg[1].bsplit = false;
//
//					fr[lcuredgeid].eg[2].v1 = row + d->x;
//					fr[lcuredgeid].eg[2].v2 = row + d->x + 1;
//					fr[lcuredgeid].eg[2].next = 0;
//					fr[lcuredgeid].eg[2].bsplit = false;
//
//					fr[lcuredgeid].eg[3].v1 = row;
//					fr[lcuredgeid].eg[3].v2 = row + d->x;
//					fr[lcuredgeid].eg[3].next = 0;
//					fr[lcuredgeid].eg[3].bsplit = false;

					fr[lcuredgeid].edge[0] = lcuredgeid;
					fr[lcuredgeid].edge[1] = (row + 1) * 3 + 1;
					fr[lcuredgeid].edge[2] = (row + d->x) * 3;
					fr[lcuredgeid].edge[3] = lcuredgeid + 1;
					fr[lcuredgeid].F2V = row;
					//fr[lcuredgeid].cornerId = row + d->x + 1;
					//					fr[curEdgeidx].bSplit = false;
					//					fr[curEdgeidx].dx1 = step->x;
					//					fr[curEdgeidx].dx2 = step->y;
					fr[lcuredgeid].setInBoundary();
					fr[lcuredgeid].level = 1;
					;
					fr[lcuredgeid].dir = 1;
					;
					fr[lcuredgeid].vertexes[0] = row;
					;
					fr[lcuredgeid].vertexes[1] = row + 1;
					;
					fr[lcuredgeid].vertexes[2] = row + d->x + 1;
					;
					fr[lcuredgeid].vertexes[3] = row + d->x;
					;
					fr[lcuredgeid].center = 0;

					//					fr[lcuredgeid].state = 0;

					if (row + (d->x * d->y) < num_rows) {
						fr[lcuredgeid + 1].edge[0] = lcuredgeid;
						fr[lcuredgeid + 1].edge[1] = (row + 1) * 3 + 2;
						fr[lcuredgeid + 1].edge[2] = (row + (d->x * d->y)) * 3;
						fr[lcuredgeid + 1].edge[3] = lcuredgeid + 2;

//						fr[lcuredgeid + 1].eg[0].v1 = row;
//						fr[lcuredgeid + 1].eg[0].v2 = row + 1;
//						fr[lcuredgeid + 1].eg[0].next = 0;
//						fr[lcuredgeid + 1].eg[0].bsplit = false;
//
//						fr[lcuredgeid + 1].eg[1].v1 = row + 1;
//						fr[lcuredgeid + 1].eg[1].v2 = row + d->x * d->y + 1;
//						fr[lcuredgeid + 1].eg[1].next = 0;
//						fr[lcuredgeid + 1].eg[1].bsplit = false;
//
//						fr[lcuredgeid + 1].eg[2].v1 = row + (d->x * d->y);
//						fr[lcuredgeid + 1].eg[2].v2 = row + (d->x * d->y) + 1;
//						fr[lcuredgeid + 1].eg[2].next = 0;
//						fr[lcuredgeid + 1].eg[2].bsplit = false;
//
//						fr[lcuredgeid + 1].eg[3].v1 = row;
//						fr[lcuredgeid + 1].eg[3].v2 = row + (d->x * d->y);
//						fr[lcuredgeid + 1].eg[3].next = 0;
//						fr[lcuredgeid + 1].eg[3].bsplit = false;

						//fr[lcuredgeid + 1].cornerId = row + (d->x * d->y) + 1;
						fr[lcuredgeid + 1].F2V = row;
						//						fr[lcuredgeid + 1].bSplit = false;
						//						fr[lcuredgeid + 1].dx1 = step->x;
						//						fr[lcuredgeid + 1].dx2 = step->z;
						fr[lcuredgeid + 1].setInBoundary();
						fr[lcuredgeid + 1].level = 1;
						;
						fr[lcuredgeid + 1].dir = 2;
						;
						fr[lcuredgeid + 1].vertexes[0] = row;
						;
						fr[lcuredgeid + 1].vertexes[1] = row + 1;
						;
						fr[lcuredgeid + 1].vertexes[2] = row + d->x * d->y + 1;
						;
						fr[lcuredgeid + 1].vertexes[3] = row + d->x * d->y;
						;
						fr[lcuredgeid + 1].center = 0;
						//						fr[lcuredgeid + 1].state = 1;

					} else {
						fr[lcuredgeid + 1].unsetInBoundary();

					}
				}

			} else {
				fr[lcuredgeid].unsetInBoundary();
				fr[lcuredgeid + 1].unsetInBoundary();

			}

			if (eg[lcuredgeid + 2].isInBoundary()
					&& eg[lcuredgeid + 1].isInBoundary()
					&& (row == 0
							|| ((row) % d->x * d->y != 0
									&& row + (d->x * d->y) < num_rows))) {
				fr[lcuredgeid + 2].edge[0] = lcuredgeid + 1;
				fr[lcuredgeid + 2].edge[1] = (row + (d->x)) * 3 + 2;
				fr[lcuredgeid + 2].edge[2] = (row + (d->x * d->y)) * 3 + 1;
				fr[lcuredgeid + 2].edge[3] = lcuredgeid + 2;

//				fr[lcuredgeid + 2].eg[0].v1 = row;
//				fr[lcuredgeid + 2].eg[0].v2 = row + d->x;
//				fr[lcuredgeid + 2].eg[0].next = 0;
//				fr[lcuredgeid + 2].eg[0].bsplit = false;
//
//				fr[lcuredgeid + 2].eg[1].v1 = row + d->x;
//				fr[lcuredgeid + 2].eg[1].v2 = row + d->x * (d->y + 1);
//				fr[lcuredgeid + 2].eg[1].next = 0;
//				fr[lcuredgeid + 2].eg[1].bsplit = false;
//
//				fr[lcuredgeid + 2].eg[2].v1 = row + (d->x * d->y);
//				fr[lcuredgeid + 2].eg[2].v2 = row + (d->x * (d->y + 1));
//				fr[lcuredgeid + 2].eg[2].next = 0;
//				fr[lcuredgeid + 2].eg[1].bsplit = false;
//
//				fr[lcuredgeid + 2].eg[3].v1 = row;
//				fr[lcuredgeid + 2].eg[3].v2 = row + (d->x * d->y);
//				fr[lcuredgeid + 2].eg[3].next = 0;
//				fr[lcuredgeid + 2].eg[3].bsplit = false;

				fr[lcuredgeid + 2].F2V = row;
				//fr[lcuredgeid + 2].cornerId = row + (d->x * (d->y + 1));

				//				fr[lcuredgeid + 2].bSplit = false;
				//				fr[lcuredgeid + 2].dx1 = step->y;
				//				fr[lcuredgeid + 2].dx2 = step->z;
				fr[lcuredgeid + 2].setInBoundary();
				fr[lcuredgeid + 2].level = 1;
				;
				fr[lcuredgeid + 2].dir = 3;
				;
				fr[lcuredgeid + 2].vertexes[0] = row;
				;
				fr[lcuredgeid + 2].vertexes[1] = row + d->x;
				;
				fr[lcuredgeid + 2].vertexes[2] = row + d->x * (d->y + 1);
				;
				fr[lcuredgeid + 2].vertexes[3] = row + d->x * (d->y);
				;
				fr[lcuredgeid + 2].center = 0;
				//				fr[lcuredgeid + 2].state = 2;

			} else {
				fr[lcuredgeid + 2].unsetInBoundary();
			}

		}

		m[row] = vertex;

	}
	*_eg = eg;
	*_fr = fr;

	return;

}
void ADPAlgorithm2::error(int errorcode) {

	int temp = errorcode;
	printf("error!!\n");
}

void ADPAlgorithm2::generate_streamlines_sparsely_Ocean(ASF_vertex* m,
		bool bForward, int currData, float* m_x1, float* m_y1, float* m_z1,
		Boundary*b, Dimension* d, Point* step, int currentDimX, int tau,
		int n) {
	Point* op = new Point[n];
	int currn = 0;
	for (int k = 0; k < d->z; k += 4) {
		for (int j = 0; j < d->y; j += 4) {
			//#pragma omp parallel for firstprivate(k,j)
			//	#pragma omp parallel for firstprivate(k,j, ndim) //shared(temp_decomp1,lowBoundary,lowBoundary_z, step,step_z,rgb)
			for (int i = 0; i < d->x; i += 4) {
				op[currn].x = b->low.x + i * step->x;
				op[currn].y = b->low.y + j * step->y;
				op[currn].z = b->low.z + k * step->z;

				currn++;
			}
		}
	}
	//op.push_back(new Point(m->p));
	stringstream ss;
	ss << currData;
	string str = ss.str();
	if (!bForward)
		str = str + "Backward";
	fs->Save_Streamlines2(op, bForward, str, m_x1, m_y1, m_z1, b, d, step,
			currentDimX, tau, currn);

}

void ADPAlgorithm2::generate_streamlines_sparsely(ASF_vertex* m, bool bForward,
		int currData, float* m_x1, float* m_y1, float* m_z1, Boundary*b,
		Dimension* d, Point* step, int currentDimX, int tau, int n) {
	Point* op = new Point[n];
	int currn = 0;
	for (int i = 0; i < n; i++) {
		op[currn] = m[i].p;
		currn++;
	}
	//op.push_back(new Point(m->p));
	stringstream ss;
	ss << currData;
	string str = ss.str();
	if (!bForward)
		str = str + "Backward";
	fs->Save_Streamlines2(op, bForward, str, m_x1, m_y1, m_z1, b, d, step,
			currentDimX, tau, currn);

}

void ADPAlgorithm2::Compute_Seed_Point(ASF_vertex *m, fFace face, int sampleNum,
		Point** outPoints) {
	Point* op = new Point[(sampleNum - 1) * (sampleNum - 1)];
	float dx, dy, dz;
	ASF_vertex v1 = m[face.vertexes[0]];
	ASF_vertex v2 = m[face.vertexes[1]];
	ASF_vertex v3 = m[face.vertexes[2]];
	ASF_vertex v4 = m[face.vertexes[3]];
	float l1 = v1.p.dist(v2.p);
	float l2 = v2.p.dist(v3.p);
	float l3 = v3.p.dist(v4.p);
	float l4 = v1.p.dist(v4.p);
	dx = std::min(l1, l3);
	dy = std::min(l2, l4);

	dx = dx / sampleNum;
	dy = dy / sampleNum;

	if (face.dir == 1) {

		op[0].x = v1.p.x + dx;
		op[0].z = v1.p.z;
		op[0].y = v1.p.y + dy;
		int cnt = 0;
		for (int i = 0; i < sampleNum - 1; i++)
			for (int j = 0; j < sampleNum - 1; j++) {
				op[cnt].x = op[0].x + i * dx;
				op[cnt].y = op[0].y + j * dy;
				op[cnt].z = op[0].z;
				cnt++;
			}
	}

	else if (face.dir == 2) {

		op[0].x = v1.p.x + dx;
		op[0].y = v1.p.y;
		op[0].z = v1.p.z + dy;
		int cnt = 0;
		for (int i = 0; i < sampleNum - 1; i++)
			for (int j = 0; j < sampleNum - 1; j++) {
				op[cnt].x = op[0].x + i * dx;
				op[cnt].z = op[0].z + j * dy;
				op[cnt].y = op[0].y;
				cnt++;
			}
	}

	else if (face.dir == 3) {

		op[0].x = v1.p.x;
		op[0].y = v1.p.y + dx;
		op[0].z = v1.p.z + dy;
		int cnt = 1;
		for (int i = 0; i < sampleNum - 1; i++)
			for (int j = 0; j < sampleNum - 1; j++) {
				op[cnt].y = op[0].y + i * dx;
				op[cnt].z = op[0].z + j * dy;
				op[cnt].x = op[0].x;
				cnt++;
			}
	}

	*outPoints = op;
}

bool ADPAlgorithm2::checkPointInFaceBoundary(ASF_vertex*m, fFace face,
		ASF_vertex vcenter) {
	//	if(!((vcenter.e > m[face.vertexes[0]].e && vcenter.e < m[face.vertexes[1]].e) || (vcenter.e > m[face.vertexes[1]].e && vcenter.e <  m[face.vertexes[0]].e)))
	//			return false

	if (!(vcenter.e.x >= m[face.vertexes[0]].e.x
			&& vcenter.e.x <= m[face.vertexes[1]].e.x)
			|| !(vcenter.e.x >= m[face.vertexes[0]].e.x
					&& vcenter.e.x <= m[face.vertexes[1]].e.x))
		error(0);
	if (!(vcenter.e.y >= m[face.vertexes[0]].e.y
			&& vcenter.e.y <= m[face.vertexes[1]].e.y)
			|| !(vcenter.e.y <= m[face.vertexes[0]].e.y
					&& vcenter.e.y >= m[face.vertexes[1]].e.y))
		error(0);
	if (!(vcenter.e.z >= m[face.vertexes[0]].e.z
			&& vcenter.e.z <= m[face.vertexes[1]].e.z)
			|| !(vcenter.e.z <= m[face.vertexes[0]].e.z
					&& vcenter.e.z >= m[face.vertexes[1]].e.z))
		error(0);

	//==========================================================
	if (!(vcenter.e.x >= m[face.vertexes[1]].e.x
			&& vcenter.e.x <= m[face.vertexes[2]].e.x)
			|| !(vcenter.e.x >= m[face.vertexes[1]].e.x
					&& vcenter.e.x <= m[face.vertexes[2]].e.x))
		error(0);
	if (!(vcenter.e.y >= m[face.vertexes[1]].e.y
			&& vcenter.e.y <= m[face.vertexes[2]].e.y)
			|| !(vcenter.e.y <= m[face.vertexes[1]].e.y
					&& vcenter.e.y >= m[face.vertexes[2]].e.y))
		error(0);
	if (!(vcenter.e.z >= m[face.vertexes[1]].e.z
			&& vcenter.e.z <= m[face.vertexes[2]].e.z)
			|| !(vcenter.e.z <= m[face.vertexes[1]].e.z
					&& vcenter.e.z >= m[face.vertexes[2]].e.z))
		error(0);

	//=============================================================

	if (!(vcenter.e.x >= m[face.vertexes[2]].e.x
			&& vcenter.e.x <= m[face.vertexes[3]].e.x)
			|| !(vcenter.e.x >= m[face.vertexes[2]].e.x
					&& vcenter.e.x <= m[face.vertexes[3]].e.x))
		error(0);
	if (!(vcenter.e.y >= m[face.vertexes[2]].e.y
			&& vcenter.e.y <= m[face.vertexes[3]].e.y)
			|| !(vcenter.e.y <= m[face.vertexes[2]].e.y
					&& vcenter.e.y >= m[face.vertexes[3]].e.y))
		error(0);
	if (!(vcenter.e.z >= m[face.vertexes[2]].e.z
			&& vcenter.e.z <= m[face.vertexes[3]].e.z)
			|| !(vcenter.e.z <= m[face.vertexes[2]].e.z
					&& vcenter.e.z >= m[face.vertexes[3]].e.z))
		error(0);

	//==============================================================

	if (!(vcenter.e.x >= m[face.vertexes[0]].e.x
			&& vcenter.e.x <= m[face.vertexes[3]].e.x)
			|| !(vcenter.e.x >= m[face.vertexes[0]].e.x
					&& vcenter.e.x <= m[face.vertexes[3]].e.x))
		error(0);
	if (!(vcenter.e.y >= m[face.vertexes[0]].e.y
			&& vcenter.e.y <= m[face.vertexes[3]].e.y)
			|| !(vcenter.e.y <= m[face.vertexes[0]].e.y
					&& vcenter.e.y >= m[face.vertexes[3]].e.y))
		error(0);
	if (!(vcenter.e.z >= m[face.vertexes[0]].e.z
			&& vcenter.e.z <= m[face.vertexes[3]].e.z)
			|| !(vcenter.e.z <= m[face.vertexes[0]].e.z
					&& vcenter.e.z >= m[face.vertexes[3]].e.z))
		error(0);

	return true;
}

void ADPAlgorithm2::FindBoundingBox(ASF_vertex*m, ASF_vertex* oa,
		int pointnumber, fFace face, ASF_vertex* bbox, int bbpointnumber) {
	//	ASF_vertex minPoint = oa[0];
	//	ASF_vertex maxPoint = oa[1];
	Point mine = oa[0].e;
	Point maxe = oa[1].e;
	for (int i = 2; i < pointnumber; i++) {
		if (oa[i].e.x <= mine.x) {
			mine.x = oa[i].e.x;
		}
		if (oa[i].e.y <= mine.y) {
			mine.y = oa[i].e.y;
		}

		if (oa[i].e.z <= mine.z) {
			mine.z = oa[i].e.z;
		}

		if (oa[i].e.x >= maxe.x) {
			maxe.x = oa[i].e.x;
		}
		if (oa[i].e.y >= maxe.y) {
			maxe.y = oa[i].e.y;
		}

		if (oa[i].e.z >= maxe.z) {
			maxe.z = oa[i].e.z;
		}

//		if (oa[i].e.x >= maxPoint.e.x && oa[i].e.y >= maxPoint.e.y && oa[i].e.z >= maxPoint.e.z ) {
//			maxPoint = oa[i];
//		}

	}

	Point center;
	center.x = (mine.x + maxe.x) / 2.;
	center.y = (mine.y + maxe.y) / 2.;
	center.z = (mine.z + maxe.z) / 2.;

	float mindist = oa[0].e.dist(center);
	int minIndex = 0;
	for (int i = 1; i < pointnumber - 4; i++) {
		float tempDist = oa[i].e.dist(center);
		if (tempDist < mindist) {
			mindist = tempDist;
			minIndex = i;
		}
	}

	bbox[0].e = mine;
	bbox[1].e = maxe;
	bbox[2] = oa[minIndex];

}

bool ADPAlgorithm2::EstimateBestSeed(fFace* face, ASF_vertex*m, fFace* fc,
		fEdge* eg, int currSize, ASF_vertex* v1, ASF_vertex* v2, ASF_vertex* v3,
		ASF_vertex* v4, ASF_vertex* vm1, ASF_vertex* vm2, ASF_vertex* vm3,
		ASF_vertex* vm4, Boundary*_b, Dimension*_d, Point*_step, float* m_x1,
		float* m_y1, float* m_z1, int currentXDim, int tau, bool bForward,
		ASF_vertex*ov, int sampleNum, int curVertexId, bool bEstimated)

		{

	float curstep[3];

	int n = (sampleNum) * (sampleNum);

	float dx, dy, dz;
	int* hitnumber = new int[n];
	memset(hitnumber, 0, n * sizeof(int));
	ASF_vertex* oa = new ASF_vertex[n + 4];
	ASF_vertex* oa_n = new ASF_vertex[n + 4];
	memset(oa, 0, n * sizeof(ASF_vertex));
	memset(oa_n, 0, n * sizeof(ASF_vertex));
	float l1 = v1->p.dist(v2->p);
	float l2 = v2->p.dist(v3->p);
	float l3 = v3->p.dist(v4->p);
	float l4 = v1->p.dist(v4->p);
	dx = std::min(l1, l3);
	dy = std::min(l2, l4);

	dx = dx / (sampleNum + 1);
	dy = dy / (sampleNum + 1);

	int cnt = 0;
	double s = 1.0 / (sampleNum + 1);
	ASF_vertex temp1;
	ASF_vertex temp2;
	ASF_vertex vedge = *v1;

	int ncnt = 0;
	//	if (sampleNum == 3)
	//		error(1);

	int v1i = face->vertexes[0];
	int v2i = face->vertexes[1];
	int v3i = face->vertexes[2];
	int v4i = face->vertexes[3];

	vedge.p.x = (v1->p.x + v2->p.x + v3->p.x + v4->p.x) / 4;
	vedge.p.y = (v1->p.y + v2->p.y + v3->p.y + v4->p.y) / 4;
	vedge.p.z = (v1->p.z + v2->p.z + v3->p.z + v4->p.z) / 4;

//	vedge.e.x = (v1->e.x + v2->e.x + v3->e.x + v4->e.x) / 4;
//	vedge.e.y = (v1->e.y + v2->e.y + v3->e.y + v4->e.y) / 4;
//	vedge.e.z = (v1->e.z + v2->e.z + v3->e.z + v4->e.z) / 4;

//	vedge.p.x = (vm1->p.x + vm2->p.x + vm3->p.x + vm4->p.x) / 4;
//	vedge.p.y = (vm1->p.y + vm2->p.y + vm3->p.y + vm4->p.y) / 4;
//	vedge.p.z = (vm1->p.z + vm2->p.z + vm3->p.z + vm4->p.z) / 4;
//
//	vedge.e.x = (vm1->e.x + vm2->e.x + vm3->e.x + vm4->e.x) / 4;
//	vedge.e.y = (vm1->e.y + vm2->e.y + vm3->e.y + vm4->e.y) / 4;
//	vedge.e.z = (vm1->e.z + vm2->e.z + vm3->e.z + vm4->e.z) / 4;

	if (bEstimated) {
		vedge.e.x = (v1->e.x + v2->e.x + v3->e.x + v4->e.x) / 4;
		vedge.e.y = (v1->e.y + v2->e.y + v3->e.y + v4->e.y) / 4;
		vedge.e.z = (v1->e.z + v2->e.z + v3->e.z + v4->e.z) / 4;

//		for (int j = 0; j < tau + 1; j++) {
//			Trace[curVertexId][j].x = (Trace[v1i][j].x + Trace[v2i][j].x
//					+ Trace[v3i][j].x + Trace[v4i][j].x) / 4;
//			Trace[curVertexId][j].y = (Trace[v1i][j].y + Trace[v2i][j].y
//					+ Trace[v3i][j].y + Trace[v4i][j].y) / 4;
//			Trace[curVertexId][j].z = (Trace[v1i][j].z + Trace[v2i][j].z
//					+ Trace[v3i][j].z + Trace[v4i][j].z) / 4;
//		}

	} else {
		float p1[3];
		p1[0] = vedge.p.x;
		p1[1] = vedge.p.y;
		p1[2] = vedge.p.z;

		float ep[3];
		for (int k = 0; k < tau + 1; k++) {
			generalstreamlineTracing_single(p1, bForward, ep, m_x1, m_y1, m_z1,
					_b, _d, _step, currentXDim, 1);
			p1[0] = Trace[curVertexId][k].x = ep[0];
			p1[1] = Trace[curVertexId][k].y = ep[1];
			p1[2] = Trace[curVertexId][k].z = ep[2];

		}

		vedge.e.x = ep[0];
		vedge.e.y = ep[1];
		vedge.e.z = ep[2];
	}

//	vedge.e.x = (Trace[v1i][tau].x + Trace[v2i][tau].x + Trace[v3i][tau].x
//			+ Trace[v4i][tau].x) / 4;
//	vedge.e.y = (Trace[v1i][tau].y + Trace[v2i][tau].y + Trace[v3i][tau].y
//			+ Trace[v4i][tau].y) / 4;
//	vedge.e.z = (Trace[v1i][tau].z + Trace[v2i][tau].z + Trace[v3i][tau].z
//			+ Trace[v4i][tau].z) / 4;

	//			AdvectParticle(&vedge, m_x1, m_y1, m_z1, _b, _d, _step, currentXDim,
	//					tau + 1, bForward);
	//	memset(vedge.es,0,100*sizeof(Point));

	if (vedge.checkInBoundary(_b)) {
		//checkPointInFaceBoundary(m, face, vedge);
		int io, jo, ko;
		vedge.getIndex(_b, _step, _d, io, jo, ko);
		uint32_t range = io + jo * _d->x + ko * (_d->x * _d->y);//_tv.getRange(&_b, &step, &_d);

		vedge.setRange(range);
		//if (ESTIMATED)
		{

		}

	} else
		error(1);

	if (tau >= 30) {

		m[curVertexId] = vedge;
		fs->save_Quad_One_Face(m, *face, eg, Trace, *_d, "face85", currSize, 1,
				tau, 1, tau * 2);

		fs->start_save_Quad_One_Face(m, *face, eg, *_d, "face85", currSize, 1,
				tau, 1);

		fs->start_save_Quad_FaceWithSeedPoits(m, fc, oa, eg, *_d, "face85",
				currSize, 1, cnt, tau, 1);

		oa[cnt] = m[face->vertexes[0]];
		oa[cnt + 1] = m[face->vertexes[1]];
		oa[cnt + 2] = m[face->vertexes[2]];
		oa[cnt + 3] = m[face->vertexes[3]];
		oa[cnt + 4] = vedge;

		int intarray[5];
		intarray[0] = v1->getRange_tau(&Trace[v1i][tau], _b, _step, _d, tau);
		intarray[1] = v2->getRange_tau(&Trace[v2i][tau], _b, _step, _d, tau);
		intarray[2] = v3->getRange_tau(&Trace[v3i][tau], _b, _step, _d, tau);
		intarray[3] = v4->getRange_tau(&Trace[v4i][tau], _b, _step, _d, tau);
		intarray[4] = vedge.getRange_tau(&Trace[curVertexId][tau], _b, _step,
				_d, tau);

		Point op[5];
		op[0] = v1->p;
		op[1] = v2->p;
		op[2] = v3->p;
		op[3] = v4->p;
		op[4] = vedge.p;

		int indexarray[5];
		indexarray[0] = v1i;
		indexarray[1] = v2i;
		indexarray[2] = v3i;
		indexarray[3] = v4i;
		indexarray[4] = curVertexId;
		//
		fs->Save_Streamlines_estimated(m, Trace, indexarray, bForward, "st1",
				m_x1, m_y1, m_z1, _b, _d, _step, currentXDim, tau + 2, 5);

		fs->Save_Streamlines2(op, bForward, "st1", m_x1, m_y1, m_z1, _b, _d,
				_step, currentXDim, tau + 2, 5);
		fs->save_Quad_FaceWithSeedPoits(m, fc, oa, eg, Trace, *_d, "face85",
				currSize, 1, 5, tau, 1);
		//if (tau == 99)
		fs->save_Voxel(m, fc, intarray, eg, *_d, "face85", currSize, 5, tau, 1,
				10);

		ASF_vertex test[2];
		test[0] = oa[0];
		fs->save_Quad_FaceWithSeedPoits(m, fc, oa, eg, Trace, *_d, "face85",
				currSize, 1, cnt + 4, tau, 1);

		//		generate_streamlines_sparsely(oa, bForward, 1, m_x1, m_y1, m_z1, _b, _d,
		//				_step, currentXDim, tau + 2, cnt + 4);

	}
	*ov = vedge;

	return false;

}

//bool ADP_Algorithm::FindBestSeed(fFace* face, ASF_vertex*m, fFace* fc,
//		fEdge* eg, int currSize, ASF_vertex* v1, ASF_vertex* v2, ASF_vertex* v3,
//		ASF_vertex* v4, ASF_vertex* vm1, ASF_vertex* vm2, ASF_vertex* vm3,
//		ASF_vertex* vm4, Boundary*_b, Dimension*_d, Point*_step, float* m_x1,
//		float* m_y1, float* m_z1, int currentXDim, int tau, bool bForward,
//		ASF_vertex*ov, int sampleNum)
//
//		{
//
//	float curstep[3];
//
//	int n = (sampleNum) * (sampleNum);
//
//	float dx, dy, dz;
//	int* hitnumber = new int[n];
//	memset(hitnumber, 0, n * sizeof(int));
//	ASF_vertex* oa = new ASF_vertex[n + 4];
//	ASF_vertex* oa_n = new ASF_vertex[n + 4];
//	memset(oa, 0, n * sizeof(ASF_vertex));
//	memset(oa_n, 0, n * sizeof(ASF_vertex));
//	float l1 = v1->p.dist(v2->p);
//	float l2 = v2->p.dist(v3->p);
//	float l3 = v3->p.dist(v4->p);
//	float l4 = v1->p.dist(v4->p);
//
//	int v1i = face->vertexes[0];
//	int v1i = face->vertexes[0];
//	int v1i = face->vertexes[0];
//	int v1i = face->vertexes[0];
//
//	dx = std::min(l1, l3);
//	dy = std::min(l2, l4);
//
//	dx = dx / (sampleNum + 1);
//	dy = dy / (sampleNum + 1);
//
//	int cnt = 0;
//	double s = 1.0 / (sampleNum + 1);
//	ASF_vertex temp1;
//	ASF_vertex temp2;
//	ASF_vertex vedge = *v1;
//
//	int ncnt = 0;
//	//	if (sampleNum == 3)
//	//		error(1);
//	for (int i = 1; i < sampleNum + 1; i++) {
//		for (int j = 1; j < sampleNum + 1; j++) {
//			double t = i * s;
//			double t1 = j * s;
//			temp1.p.x = (t) * v1->p.x + (1 - t) * v2->p.x;
//			temp1.p.y = (t) * v1->p.y + (1 - t) * v2->p.y;
//			temp1.p.z = (t) * v1->p.z + (1 - t) * v2->p.z;
//
////			temp1.e.x = (t) * v1->e.x + (1 - t) * v2->e.x;
////			temp1.e.y = (t) * v1->e.y + (1 - t) * v2->e.y;
////			temp1.e.z = (t) * v1->e.z + (1 - t) * v2->e.z;
//
//			temp1.e.x = (t) * Trace[v1] + (1 - t) * v2->e.x;
//			temp1.e.y = (t) * v1->e.y + (1 - t) * v2->e.y;
//			temp1.e.z = (t) * v1->e.z + (1 - t) * v2->e.z;
//
//			temp2.p.x = (t) * v4->p.x + (1 - t) * v3->p.x;
//			temp2.p.y = (t) * v4->p.y + (1 - t) * v3->p.y;
//			temp2.p.z = (t) * v4->p.z + (1 - t) * v3->p.z;
//
//			temp2.e.x = (t) * v3->e.x + (1 - t) * v4->e.x;
//			temp2.e.y = (t) * v3->e.y + (1 - t) * v4->e.y;
//			temp2.e.z = (t) * v3->e.z + (1 - t) * v4->e.z;
//
//			vedge.p.x = (t1) * temp1.p.x + (1 - t1) * temp2.p.x;
//			vedge.p.y = (t1) * temp1.p.y + (1 - t1) * temp2.p.y;
//			vedge.p.z = (t1) * temp1.p.z + (1 - t1) * temp2.p.z;
//
////			vedge.e.x = (t1) * temp1.e.x + (1 - t1) * temp2.e.x;
////			vedge.e.y = (t1) * temp1.e.y + (1 - t1) * temp2.e.y;
////			vedge.e.z = (t1) * temp1.e.z + (1 - t1) * temp2.e.z;
//
////			AdvectParticle(&vedge, m_x1, m_y1, m_z1, _b, _d, _step, currentXDim,
////					tau + 1, bForward);
//
//			if (vedge.checkInBoundary(_b)) {
//				//checkPointInFaceBoundary(m, face, vedge);
//				int io, jo, ko;
//				vedge.getIndex(_b, _step, _d, io, jo, ko);
//				uint32_t range = io + jo * _d->x + ko * (_d->x * _d->y);//_tv.getRange(&_b, &step, &_d);
//
//				vedge.setRange(range);
//				for (int j = 0; j < tau + 1; j++) {
//
//					temp1.e.x = (t) * v1->e.x + (1 - t) * v2->e.x;
//					temp1.e.y = (t) * v1->e.y + (1 - t) * v2->e.y;
//					temp1.e.z = (t) * v1->e.z + (1 - t) * v2->e.z;
//
//					temp2.e.x = (t) * v3->e.x + (1 - t) * v4->e.x;
//					temp2.e.y = (t) * v3->e.y + (1 - t) * v4->e.y;
//					temp2.e.z = (t) * v3->e.z + (1 - t) * v4->e.z;
//
//					vedge.es[j].x = (t1) * temp1.e.x + (1 - t1) * temp2.e.x;
//					vedge.es[j].y = (t1) * temp1.e.y + (1 - t1) * temp2.e.y;
//					vedge.es[j].z = (t1) * temp1.e.z + (1 - t1) * temp2.e.z;
//				}
//
//			}
////				else {
////				//	error(0);
////				int io, jo, ko;
////				vedge.getIndex(_b, _step, _d, io, jo, ko);
////				uint32_t range = io + jo * _d->x + ko * (_d->x * _d->y);//_tv.getRange(&_b, &step, &_d);
////
////				vedge.setRange(range);
////				oa_n[ncnt++] = vedge;
////				continue;
////
////			}
//			oa[cnt] = vedge;
//			int hitnumber1 = CheckSeedInFace(vedge, *vm1, *vm2, *vm3, *vm4, _b,
//					_d, _step, bForward, tau);
//
//			int hitnumber2 = CheckSeedInFace(vedge, *v1, *v2, *v3, *v4, _b, _d,
//					_step, bForward, tau);
//
//			if (hitnumber1 == 0 || hitnumber2 == 0) {
//				*ov = oa[cnt];
//				return true;
//			}
//			cnt++;
//		}
//	}
//
//	string dataName = "1";
//	if (sampleNum == 3 && cnt == 0) {
//
//		//if (face.parent > 0)
//		fs->save_Quad_One_Face_parent(m, fc[face->parent], eg, *_d, dataName,
//				currSize, 1, tau, 2, tau * 2);
//
//		fs->save_Quad_One_Face(m, *face, eg, *_d, dataName, currSize, 1, tau, 1,
//				tau * 2);
//
//		fs->start_save_Quad_One_Face(m, *face, eg, *_d, dataName, currSize, 1,
//				tau, 1);
//		fs->start_save_Quad_FaceWithSeedPoits(m, fc, oa, eg, *_d, dataName,
//				currSize, 1, cnt, tau, 1);
//
//		oa_n[ncnt] = m[face->vertexes[0]];
//		oa_n[ncnt + 1] = m[face->vertexes[1]];
//		oa_n[ncnt + 2] = m[face->vertexes[2]];
//		oa_n[ncnt + 3] = m[face->vertexes[3]];
//
//		ASF_vertex test[2];
//		test[0] = oa_n[0];
//		fs->save_Quad_FaceWithSeedPoits(m, fc, oa_n, eg, *_d, dataName,
//				currSize, 1, cnt + 4, tau, 1);
//
//		generate_streamlines_sparsely(oa_n, bForward, 1, m_x1, m_y1, m_z1, _b,
//				_d, _step, currentXDim, tau + 2, ncnt + 4);
//
//		int* va = new int[n + 4];
//		for (int i = 0; i < n + 4; i++)
//			va[i] = oa_n[i].getRange();
//
//		fs->save_Voxel(m, fc, va, eg, *_d, "face85", currSize, n + 4, tau, 1,
//				tau * 5);
//
//		//		fs->save_Voxel(m, fc, va, eg, *_d, "face85", currSize, cnt + 4, tau, 1,
//		//				tau * 5);
//
//	}
//
//	if (cnt == 0 && sampleNum == 3) {
//		face->unsetInBoundary();
//		return false;
//	}
//	if (sampleNum == 3) {
//		face->unsetInBoundary();
//
//		return false;
//		//if (face.parent > 0)
//		fs->save_Quad_One_Face_parent(m, *face, eg, *_d, dataName, currSize, 1,
//				tau, 2, tau * 2);
//
//		fs->save_Quad_One_Face_parent(m, fc[face->parent], eg, *_d, dataName,
//				currSize, 3, tau, 3, tau * 3);
//
//		fs->save_Quad_One_Face(m, *face, eg, *_d, dataName, currSize, 1, tau, 1,
//				tau * 2);
//
//		fs->start_save_Quad_One_Face(m, *face, eg, *_d, dataName, currSize, 1,
//				tau, 1);
//
//		fs->start_save_Quad_FaceWithSeedPoits(m, fc, oa, eg, *_d, dataName,
//				currSize, 1, cnt, tau, 1);
//
//		oa[cnt] = m[face->vertexes[0]];
//		oa[cnt + 1] = m[face->vertexes[1]];
//		oa[cnt + 2] = m[face->vertexes[2]];
//		oa[cnt + 3] = m[face->vertexes[3]];
//
//		ASF_vertex test[2];
//		test[0] = oa[0];
//		fs->save_Quad_FaceWithSeedPoits(m, fc, oa, eg, *_d, dataName, currSize,
//				1, cnt + 4, tau, 1);
//
//		generate_streamlines_sparsely(oa, bForward, 1, m_x1, m_y1, m_z1, _b, _d,
//				_step, currentXDim, tau + 2, cnt + 4);
//
//		int* va = new int[cnt + 4];
//		for (int i = 0; i < cnt + 4; i++)
//			va[i] = oa[i].getRange();
//
//		fs->save_Voxel(m, fc, va, eg, *_d, "face85", currSize, n + 4, tau, 1,
//				tau * 5);
//
//		ASF_vertex* bbox = new ASF_vertex[3];
//		FindBoundingBox(m, oa, cnt + 4, *face, bbox, 3);
//
//		fs->save_Quad_FaceWithSeedPoits(m, fc, bbox, eg, *_d, dataName,
//				currSize, 1, 3, tau, 1);
//
//		face->Jacobian = ComputeJacobian(m, *face, m_x1, m_y1, m_z1, _step);
//
//		Jacobianarray[numJ++] = face->Jacobian;
//		//error(0);
////		int min = 0;
////		int minvalue = hitnumber[0];
////		for (int j = 0; j < cnt; j++)
////			if (hitnumber[j] < minvalue) {
////				minvalue = hitnumber[j];
////				min = j;
////			}
//		*ov = bbox[2];
//		return false;
//	}
//
//	return false;
//
//}

void ADPAlgorithm2::AdvectParticle_estimated(ASF_vertex* vedge, float* m_x1,
		float* m_y1, float* m_z1, Boundary* b, Dimension* d, Point* step,
		int currentXDim, int tau, bool bForward) {
	float p[3], e[3];
	p[0] = vedge->p.x;
	p[1] = vedge->p.y;
	p[2] = vedge->p.z;

	if (vedge->checkInBoundary(b)) {
		int io, jo, ko;
		vedge->getIndex(b, step, d, io, jo, ko);
		uint32_t range = io + jo * d->x + ko * (d->x * d->y); //_tv.getRange(&_b, &step, &_d);

		vedge->setRange(range);
//		if (vedge->type == 2) {
//			for(int j = 0; j < tau; j++)
//			{
//				vedge->es[j] =
//			}
//		}

	} else {
		int io, jo, ko;
		vedge->getIndex(b, step, d, io, jo, ko);
		uint32_t range = io + jo * d->x + ko * (d->x * d->y); //_tv.getRange(&_b, &step, &_d);

		vedge->setRange(range);
		vedge->unsetInBoundary();

	}
}

void ADPAlgorithm2::AdvectParticle(ASF_vertex* vedge, float* m_x1, float* m_y1,
		float* m_z1, Boundary* b, Dimension* d, Point* step, int currentXDim,
		int tau, bool bForward) {
	float p[3], e[3];
	p[0] = vedge->p.x;
	p[1] = vedge->p.y;
	p[2] = vedge->p.z;

	generalstreamlineTracing_single(p, bForward, e, m_x1, m_y1, m_z1, b, d,
			step, currentXDim, tau);

	vedge->e.x = e[0];
	vedge->e.y = e[1];
	vedge->e.z = e[2];

	if (vedge->checkInBoundary(b)) {
		int io, jo, ko;
		vedge->getIndex(b, step, d, io, jo, ko);
		uint32_t range = io + jo * d->x + ko * (d->x * d->y); //_tv.getRange(&_b, &step, &_d);

		vedge->setRange(range);

	} else {
		int io, jo, ko;
		vedge->getIndex(b, step, d, io, jo, ko);
		uint32_t range = io + jo * d->x + ko * (d->x * d->y); //_tv.getRange(&_b, &step, &_d);

		vedge->setRange(range);
		vedge->unsetInBoundary();

	}
}

int ADPAlgorithm2::CheckSeedInFace(fFace* face, ASF_vertex vcenter,
		ASF_vertex v1, ASF_vertex v2, ASF_vertex v3, ASF_vertex v4, Boundary*_b,
		Dimension* _d, Point* _step, bool bForward, int tau) {
	int n = 0;

	int v1i = face->vertexes[0];
	int v2i = face->vertexes[1];
	int v3i = face->vertexes[2];
	int v4i = face->vertexes[3];

	if (!checkEdge(vcenter, v1, _b, _d, _step, bForward, tau, v1i, v2i))
		n++;
	if (!checkEdge(vcenter, v2, _b, _d, _step, bForward, tau, v1i, v2i))
		n++;
	if (!checkEdge(vcenter, v3, _b, _d, _step, bForward, tau, v1i, v2i))
		n++;
	if (!checkEdge(vcenter, v4, _b, _d, _step, bForward, tau, v1i, v2i))
		n++;
	return n;
}

void ADPAlgorithm2::TestInsertFaceCenterWithMultipleSeeds(ASF_vertex *m,
		fFace* fc, fEdge* eg, uint32_t* Fr_Face, uint32_t* Fe_Face,
		uint32_t* Fe, Dimension* d, Boundary* b, Point* step, float* m_x1,
		float* m_y1, float* m_z1, int tau, int currentXDim, bool bForward,
		uint32_t num_vertex, uint32_t num_edges, uint32_t num_faces) {

	for (uint32_t row = 0; row < num_faces; row++) {

		if (row >= num_faces)
			return;
		fFace face = fc[row];
		if (Fe_Face[row] == 1 && !face.bsplit && face.isInBoundary()) {

			if (row == 2127)
				face = fc[row];
			//face.level = face.level + 1;

			int curVertexId = num_vertex + Fr_Face[row];
			//int curFaceId = num_faces + (Fr_Face[row] * 4);

			fEdge edge1 = eg[eg[face.edge[0]].subedge[0]];
			fEdge edge2 = eg[eg[face.edge[1]].subedge[0]];
			fEdge edge3 = eg[eg[face.edge[2]].subedge[0]];
			fEdge edge4 = eg[eg[face.edge[3]].subedge[0]];

			ASF_vertex vertex1 = m[face.vertexes[0]];
			ASF_vertex vertex2 = m[face.vertexes[1]];
			ASF_vertex vertex3 = m[face.vertexes[2]];
			ASF_vertex vertex4 = m[face.vertexes[3]];

			//			if (!vertex1.isInBoundary() || !vertex2.isInBoundary()
			//					|| !vertex3.isInBoundary() || !vertex4.isInBoundary())
			//				error(1);

			ASF_vertex v1_m = m[edge1.v2];
			ASF_vertex v2_m = m[edge2.v2];
			ASF_vertex v3_m = m[edge3.v2];
			ASF_vertex v4_m = m[edge4.v2];
			if (!v1_m.isInBoundary() || !v2_m.isInBoundary()
					|| !v3_m.isInBoundary() || !v4_m.isInBoundary()) {
				face.unsetInBoundary();
				fc[row] = face;
				continue;
				//	error(1);
			}

			ASF_vertex seeds;
			ASF_vertex vcenter = vertex1;
			//			vcenter.p.x = (v1_m.p.x + v2_m.p.x + v3_m.p.x + v4_m.p.x) / 4.0;
			//			vcenter.p.y = (v1_m.p.y + v2_m.p.y + v3_m.p.y + v4_m.p.y) / 4.0;
			//			vcenter.p.z = (v1_m.p.z + v2_m.p.z + v3_m.p.z + v4_m.p.z) / 4.0;
			//
			//			vcenter.e.x = (v1_m.e.x + v2_m.e.x + v3_m.e.x + v4_m.e.x) / 4.0;
			//			vcenter.e.y = (v1_m.e.y + v2_m.e.y + v3_m.e.y + v4_m.e.y) / 4.0;
			//			vcenter.e.z = (v1_m.e.z + v2_m.e.z + v3_m.e.z + v4_m.e.z) / 4.0;

			vcenter.p.x =
					(vertex1.p.x + vertex2.p.x + vertex3.p.x + vertex4.p.x)
							/ 4.0;
			vcenter.p.y =
					(vertex1.p.y + vertex2.p.y + vertex3.p.y + vertex4.p.y)
							/ 4.0;
			vcenter.p.z =
					(vertex1.p.z + vertex2.p.z + vertex3.p.z + vertex4.p.z)
							/ 4.0;

			vcenter.e.x =
					(vertex1.e.x + vertex2.e.x + vertex3.e.x + vertex4.e.x)
							/ 4.0;
			vcenter.e.y =
					(vertex1.e.y + vertex2.e.y + vertex3.e.y + vertex4.e.y)
							/ 4.0;
			vcenter.e.z =
					(vertex1.e.z + vertex2.e.z + vertex3.e.z + vertex4.e.z)
							/ 4.0;

			if (vcenter.checkInBoundary(b)) {
				int io, jo, ko;
				vcenter.setInBoundary();
				vcenter.getIndex(b, step, d, io, jo, ko);
				uint32_t range = io + jo * d->x + ko * (d->x * d->y); //_tv.getRange(&_b, &step, &_d);

				vcenter.setRange(range);

			} else
				error(1);

			//			AdvectParticle(&vcenter, m_x1, m_y1, m_z1, b, d, step, currentXDim,
			//					tau + 1, bForward);
			bool bfound = false;
			int* hitnumber = new int[16];
			int c = 0;

			//			for (int i = 1; i < 5; i++) {
			//				if (FindBestSeed(face, &vertex1, &vertex2, &vertex3, &vertex4,
			//						&v1_m, &v2_m, &v3_m, &v4_m, b, d, step, m_x1, m_y1,
			//						m_z1, currentXDim, tau, bForward, &seeds, i) == true) {
			//					bfound = true;
			//					vcenter = seeds;
			//					break;
			//				}
			//			}

			//				for (int j = 0; j < i * i; j++) {
			//					vcenter = seeds[j];
			//					AdvectParticle(&vcenter, m_x1, m_y1, m_z1, b, d, step,
			//							currentXDim, tau + 1, bForward);
			//					hitnumber[c] = CheckSeedInFace(vcenter, vertex1, vertex2,
			//							vertex3, vertex4, b, d, step, bForward);
			//					if (hitnumber[c] == 0) {
			//						bfound = true;
			//						break;
			//					}
			//				}
			//				if (bfound)
			//					break;
			//
			//			}
			//			if (!bfound)
			//				c = 0;
			face.center = curVertexId;
			fc[row] = face;
			m[curVertexId] = vcenter;

		}
	}

}

void ADPAlgorithm2::InsertFaceCenterWithMultipleSeeds(ASF_vertex *m, fFace* fc,
		fEdge* eg, uint32_t* Fr_Face, uint32_t* Fe_Face, uint32_t* Fe,
		Dimension* d, Boundary* b, Point* step, float* m_x1, float* m_y1,
		float* m_z1, int tau, int currentXDim, bool bForward,
		uint32_t num_vertex, uint32_t num_edges, uint32_t num_faces) {

	for (uint32_t row = 0; row < num_faces; row++) {

		if (row >= num_faces)
			return;

		if (Fe_Face[row] == 1) {
			fFace face = fc[row];

			if (row == 633)
				face = fc[row];
			//face.level = face.level + 1;

			int curVertexId = num_vertex + Fr_Face[row];
			//int curFaceId = num_faces + (Fr_Face[row] * 4);

			fEdge edge1 = eg[eg[face.edge[0]].subedge[0]];
			fEdge edge2 = eg[eg[face.edge[1]].subedge[0]];
			fEdge edge3 = eg[eg[face.edge[2]].subedge[0]];
			fEdge edge4 = eg[eg[face.edge[3]].subedge[0]];

			if (!face_correctnessCheck(face, eg))
				error(1);
			ASF_vertex vertex1 = m[face.vertexes[0]];
			ASF_vertex vertex2 = m[face.vertexes[1]];
			ASF_vertex vertex3 = m[face.vertexes[2]];
			ASF_vertex vertex4 = m[face.vertexes[3]];

			ASF_vertex v1_m = m[edge1.v2];
			ASF_vertex v2_m = m[edge2.v2];
			ASF_vertex v3_m = m[edge3.v2];
			ASF_vertex v4_m = m[edge4.v2];

//			if (curVertexId == 6120)
//				error(0);

			ASF_vertex seeds;
			ASF_vertex vcenter = vertex1;
			bool bfound = false;
			int* hitnumber = new int[16];
			int c = 0;

			//for (int i = 1; i < 4; i++)
			{
				int i = 1;
				EstimateBestSeed(&face, m, fc, eg, num_vertex, &vertex1,
						&vertex2, &vertex3, &vertex4, &v1_m, &v2_m, &v3_m,
						&v4_m, b, d, step, m_x1, m_y1, m_z1, currentXDim, tau,
						bForward, &seeds, i, curVertexId, true);

				bfound = true;
				vcenter = seeds;
				//		break;

			}

			seeds.type = 1;

//			seeds.left = face.vertexes[0];
//			seeds.right = face.vertexes[2];
//			seeds.up = face.vertexes[1];
//			seeds.down = face.vertexes[3];

			seeds.left = edge1.v2;
			seeds.right = edge3.v2;
			seeds.up = edge2.v2;
			seeds.down = edge4.v2;

			vcenter = seeds;

			if (!face.isInBoundary()) {
				error(1);
				face.center = curVertexId;
				fc[row] = face;
				vcenter.unsetInBoundary();
				m[curVertexId] = vcenter;

				return;
			}

			if (!bfound) {
				//	Fe_moresplit[row] = 1;
				vcenter = seeds;
				//	error(1);
			}
			if (!vcenter.checkInBoundary(b)) {
				face.unsetInBoundary();
				error(1);
			}
			face.center = curVertexId;
			fc[row] = face;
			m[curVertexId] = vcenter;

		}
	}

}

bool ADPAlgorithm2::face_correctnessCheck(fFace face, fEdge* eg) {
	fEdge eg1 = eg[face.edge[0]];
	fEdge eg2 = eg[face.edge[1]];
	fEdge eg3 = eg[face.edge[2]];
	fEdge eg4 = eg[face.edge[3]];

	bool berror = true;
	int facelevel = face.level;
	if (facelevel != eg1.level || face.level != eg2.level
			|| face.level != eg3.level || face.level != eg4.level) {
		error(1);
		berror = false;
		//return false;
	}

	if (!(face.vertexes[0] == eg1.v1 && face.vertexes[1] == eg1.v2)) {

		berror = false;
		error(1);
	}

	if (!(face.vertexes[1] == eg2.v1 && face.vertexes[2] == eg2.v2)) {
		berror = false;
		error(1);
	}
	if (!(face.vertexes[3] == eg3.v1 && face.vertexes[2] == eg3.v2)) {
		berror = false;
		error(1);
	}

	if (!(face.vertexes[0] == eg4.v1 && face.vertexes[3] == eg4.v2)) {
		berror = false;
		error(1);
	}
	return berror;

}

void ADPAlgorithm2::UpdatesubEdge(fEdge* eg, fEdge child1, int row,
		int curFaceId1, int curFaceId2) {

	if (child1.subedge[0] > 0 && child1.subedge[1] > 0) {
		fEdge temp = eg[child1.subedge[0]];
		bool bfound = false;
		for (int i = 0; i < 4; i++) {
			if (temp.E2F[i] == row) {
				eg[child1.subedge[0]].E2F[i] = curFaceId1;
				bfound = true;
				break;
			}

		}
		if (!bfound)
			error(0);
		else {
			//eg[child1.subedge[0]] = temp;
			UpdatesubEdge(eg, temp, row, curFaceId1, curFaceId1);
		}
		bfound = false;
		temp = eg[child1.subedge[1]];
		for (int i = 0; i < 4; i++) {
			if (temp.E2F[i] == row) {
				eg[child1.subedge[1]].E2F[i] = curFaceId2;
				bfound = true;
				break;
			}

		}
		if (!bfound)
			error(0);
		else { //eg[child1.subedge[1]] = temp;
			UpdatesubEdge(eg, temp, row, curFaceId2, curFaceId2);
		}
	}

}
void ADPAlgorithm2::FindSubEdge(fEdge* eg, fFace face, fEdge* child1,
		fEdge* child2, int whichedge) {
	fEdge edge1 = eg[face.edge[whichedge]];
	fEdge edge1_1 = eg[face.edge[whichedge]];
	fEdge edge1_2 = eg[face.edge[whichedge]];

	if (edge1.bsplit) {
		edge1_1 = eg[edge1.subedge[0]];
		edge1_2 = eg[edge1.subedge[1]];
	} else {
		error(0);
		//printf("error!!");
	}
	*child1 = edge1_1;
	*child2 = edge1_2;
}

void ADPAlgorithm2::SplitFace2(ASF_vertex *m, fFace* fc, fEdge* eg,
		uint32_t* Fr_Face, uint32_t* Fe_Face, uint32_t* Fe, Dimension* d,
		Boundary* b, Point* step, bool bForward, uint32_t num_vertex,
		uint32_t num_edges, uint32_t num_faces, int curtau,
		int splitnumInCurtau, int whichData, float* m_x1, float* m_y1,
		float* m_z1, int tau) {

	fFace* disconnectedFace_o = new fFace[Fr_Face[num_faces]];
	fFace* originalFace = new fFace[Fr_Face[num_faces]];
	fFace* wholeface = new fFace[num_faces];
	fFace* subFace = new fFace[Fr_Face[num_faces] * 4];

	int cfnum = 0;
	int cfwholenum = 0;
	for (uint32_t row = 0; row < num_faces; row++)

	{

		if (row >= num_faces)
			return;

		fFace face = fc[row];

		if (Fe_Face[row] == 1 && face.isInBoundary()) {
			//	face.level = face.level + 1;

//			if (row == 252 || row == 4251)
//				error(1);

			if (face.center == 0) {
				face = fc[row];
				continue;
			}
			int curVertexId = face.center;

			int curEdgeId = num_edges + (Fr_Face[row] * 4);
			int curFaceId = num_faces + (Fr_Face[row] * 4);
			if (fabs(curEdgeId - 15840) < 4) //|| row == 101217)
				face = fc[row];

			fEdge edge1 = eg[face.edge[0]];
			fEdge edge1_1 = eg[face.edge[0]];
			fEdge edge1_2 = eg[face.edge[0]];
			FindSubEdge(eg, face, &edge1_1, &edge1_2, 0);

			fEdge edge2 = eg[face.edge[1]];
			fEdge edge2_1 = eg[face.edge[1]];
			fEdge edge2_2 = eg[face.edge[1]];

			FindSubEdge(eg, face, &edge2_1, &edge2_2, 1);

			fEdge edge3 = eg[face.edge[2]];
			fEdge edge3_1 = edge3;
			fEdge edge3_2 = eg[face.edge[2]];
			FindSubEdge(eg, face, &edge3_1, &edge3_2, 2);

			fEdge edge4 = eg[face.edge[3]];
			fEdge edge4_1 = eg[face.edge[3]];
			fEdge edge4_2 = eg[face.edge[3]];
			FindSubEdge(eg, face, &edge4_1, &edge4_2, 3);

			int co = 0;

			fFace face1 = face;
			fFace face2 = face;
			fFace face3 = face;
			fFace face4 = face;
			fEdge edge1_c = edge1_1;
			edge1_c.v1 = edge1_1.v2;
			edge1_c.v2 = curVertexId;
			edge1_c.subedge[0] = -1;
			edge1_c.subedge[1] = -1;
			edge1_c.E2F[0] = curFaceId;
			edge1_c.E2F[1] = curFaceId + 1;
			edge1_c.E2F[2] = -1;
			edge1_c.E2F[3] = -1;
			edge1_c.bsplit = false;
			edge1_c.next = 0;
			edge1_c.level = face.level + 1;
			edge1_c.parent = -1;
			eg[curEdgeId + co] = edge1_c;

//			bool bfound = false;
//			for (int k = 0; k < 4; k++) {
//				if (edge1_1.E2F[k] == row) {
//					eg[edge1.subedge[0]].E2F[k] = curFaceId;
//					bfound = true;
//
//				}
//			}
//			if (!bfound)
//				error(0);
//
//			if (edge1_1.subedge[0] > 0)
//				UpdatesubEdge(eg, edge1_1, row, curFaceId, curFaceId);
//
//			bfound = false;
//			for (int k = 0; k < 4; k++) {
//				if (edge1_2.E2F[k] == row) {
//					eg[edge1.subedge[1]].E2F[k] = curFaceId + 1;
//					bfound = true;
//					break;
//				}
//			}
//			if (!bfound)
//				error(0);
//
//			if (edge1_2.subedge[0] > 0)
//				UpdatesubEdge(eg, edge1_2, row, curFaceId + 1, curFaceId + 1);

			//edge1_1 = eg[edge1.subedge[0]];
			//edge1_2 = eg[edge1.subedge[1]];
			//eg[edge1.subedge[0]] = edge1_1;

			int c1 = curEdgeId + co;

			co++;
			//======================================================
			fEdge edge2_c = edge2_1;
			edge2_c.v2 = edge2_1.v2;
			edge2_c.v1 = curVertexId;

			edge2_c.E2F[0] = curFaceId + 1;
			edge2_c.E2F[1] = curFaceId + 2;
			edge2_c.E2F[2] = -1;
			edge2_c.E2F[3] = -1;
			edge2_c.next = 0;
			edge2_c.bsplit = false;
			edge2_c.subedge[0] = -1;
			edge2_c.subedge[1] = -1;
			edge2_c.parent = -1;
//			bfound = false;
//			for (int k = 0; k < 4; k++) {
//				if (edge2_1.E2F[k] == row) {
//					eg[edge2.subedge[0]].E2F[k] = curFaceId + 1;
//					bfound = true;
//
//				}
//			}
//			if (!bfound)
//				error(0);
//
//			if (edge2_1.subedge[0] > 0)
//				UpdatesubEdge(eg, edge2_1, row, curFaceId + 1, curFaceId + 1);
//
//			bfound = false;
//			for (int k = 0; k < 4; k++) {
//				if (edge2_2.E2F[k] == row) {
//					eg[edge2.subedge[1]].E2F[k] = curFaceId + 2;
//					bfound = true;
//
//				}
//			}
//			if (!bfound)
//				error(0);
//
//			if (edge2_2.subedge[0] > 0)
//				UpdatesubEdge(eg, edge2_2, row, curFaceId + 2, curFaceId + 2);

			//	UpdatesubEdge( eg, edge2,  row, curFaceId + 1, curFaceId + 2);
			//			edge2_1 = eg[edge2.subedge[0]];
			//			edge2_2 = eg[edge2.subedge[1]];

			int c2 = curEdgeId + co;
			eg[curEdgeId + co] = edge2_c;

			co++;
			//========================================================
			//bool bfound = false;
			fEdge edge3_c = edge3_1;
			edge3_c.v2 = edge3_1.v2;
			edge3_c.v1 = curVertexId;
			edge3_c.subedge[0] = -1;
			edge3_c.subedge[1] = -1;
			edge3_c.E2F[0] = curFaceId + 2;
			edge3_c.E2F[1] = curFaceId + 3;
			edge3_c.E2F[2] = -1;
			edge3_c.E2F[3] = -1;
			edge3_c.next = 0;
			edge3_c.bsplit = false;
			edge3_c.parent = -1;
			eg[curEdgeId + co] = edge3_c;

//			for (int k = 0; k < 4; k++) {
//				if (edge3_1.E2F[k] == row) {
//					eg[edge3.subedge[0]].E2F[k] = curFaceId + 3;
//					bfound = true;
//					break;
//				}
//			}
//			if (!bfound)
//				error(0);
//
//			if (edge3_1.subedge[0] > 0)
//				UpdatesubEdge(eg, edge3_1, row, curFaceId + 3, curFaceId + 3);
//
//			bfound = false;
//
//			for (int k = 0; k < 4; k++) {
//				if (edge3_2.E2F[k] == row) {
//					eg[edge3.subedge[1]].E2F[k] = curFaceId + 2;
//					bfound = true;
//					break;
//				}
//			}
//			if (!bfound)
//				error(0);
//
//			if (edge3_2.subedge[0] > 0)
//				UpdatesubEdge(eg, edge3_2, row, curFaceId + 2, curFaceId + 2);

			//bfound = false;
			//eg[edge3.subedge[0]] = edge3_1;
			//eg[edge3.subedge[1]] = edge3_2;

			int c3 = curEdgeId + co;
			eg[curEdgeId + co] = edge3_c;

			co++;
			//=======================================================

			fEdge edge4_c = edge4_1;
			edge4_c.v1 = edge4_1.v2;
			edge4_c.v2 = curVertexId;
			edge4_c.subedge[0] = -1;
			edge4_c.subedge[1] = -1;
			edge4_c.E2F[0] = curFaceId;
			edge4_c.E2F[1] = curFaceId + 3;
			edge4_c.E2F[2] = -1;
			edge4_c.E2F[3] = -1;
			edge4_c.next = 0;
			edge4_c.bsplit = false;
			edge4_c.parent = -1;
			eg[curEdgeId + co] = edge4_c;
//			for (int k = 0; k < 4; k++) {
//				if (edge4_1.E2F[k] == row) {
//					eg[edge4.subedge[0]].E2F[k] = curFaceId;
//					bfound = true;
//					break;
//				}
//			}
//			if (!bfound)
//				error(0);
//
//			if (edge4_1.subedge[0] > 0)
//				UpdatesubEdge(eg, edge4_1, row, curFaceId, curFaceId);
//
//			bfound = false;
//
//			for (int k = 0; k < 4; k++) {
//				if (edge4_2.E2F[k] == row) {
//					eg[edge4.subedge[1]].E2F[k] = curFaceId + 3;
//					bfound = true;
//					break;
//				}
//			}
//			if (!bfound)
//				error(0);
//
//			if (edge4_2.subedge[0] > 0)
//				UpdatesubEdge(eg, edge4_2, row, curFaceId + 3, curFaceId + 3);
//
//			bfound = false;

			//eg[edge4.subedge[0]] = edge4_1;
			//eg[edge4.subedge[1]] = edge4_2;

			int c4 = curEdgeId + co;
			eg[curEdgeId + co] = edge4_c;

			co++;

			//=====================================================
			face1.edge[0] = edge1.subedge[0];
			face1.edge[1] = curEdgeId;
			face1.edge[2] = curEdgeId + 3;
			face1.edge[3] = edge4.subedge[0];
			face1.vertexes[0] = face.vertexes[0];
			face1.vertexes[1] = edge1_1.v2;
			face1.vertexes[2] = curVertexId;
			face1.vertexes[3] = edge4_1.v2;
			face1.center = 0;
			face1.level = face.level + 1;
			face1.parent = row;
			face1.F2V = face.F2V;
			face1.isInBoundary();
			//==================================
			face2.edge[0] = edge1.subedge[1];
			face2.edge[1] = edge2.subedge[0];
			face2.edge[2] = curEdgeId + 1;
			face2.edge[3] = curEdgeId;
			face2.vertexes[0] = edge1_1.v2;
			face2.vertexes[1] = face.vertexes[1];
			face2.vertexes[2] = edge2_1.v2;
			face2.vertexes[3] = curVertexId;
			face2.center = 0;
			face2.level = face.level + 1;
			face2.parent = row;
			//====================================
			face3.edge[0] = curEdgeId + 1;
			face3.edge[1] = edge2.subedge[1];
			face3.edge[2] = edge3.subedge[1];
			face3.edge[3] = curEdgeId + 2;
			face3.vertexes[0] = curVertexId;
			face3.vertexes[1] = edge2_1.v2;
			;
			face3.vertexes[2] = face.vertexes[2];
			face3.vertexes[3] = edge3_1.v2;
			face3.center = 0;
			face3.level = face.level + 1;
			face3.parent = row;
			//====================================

			face4.edge[0] = curEdgeId + 3;
			face4.edge[1] = curEdgeId + 2;
			;
			face4.edge[2] = edge3.subedge[0];
			face4.edge[3] = edge4.subedge[1];
			face4.vertexes[0] = edge4_1.v2;
			face4.vertexes[1] = curVertexId;
			face4.vertexes[2] = edge3_1.v2;
			face4.vertexes[3] = face.vertexes[3];
			face4.center = 0;
			face4.level = face.level + 1;
			face4.parent = row;

			fc[row].bsplit = true;
			fc[row].subface[0] = curFaceId;
			fc[row].subface[1] = curFaceId + 1;
			fc[row].subface[2] = curFaceId + 2;
			fc[row].subface[3] = curFaceId + 3;

			face_correctnessCheck(face, eg);
			face_correctnessCheck(face1, eg);
			face_correctnessCheck(face2, eg);
			face_correctnessCheck(face3, eg);
			face_correctnessCheck(face4, eg);

			//fc[row] = face;
			face1.subface[0] = face1.subface[1] = face1.subface[2] =
					face1.subface[3] = -1;
			face1.bsplit = false;
			fc[curFaceId] = face1;
			face2.bsplit = false;
			face2.subface[0] = face2.subface[1] = face2.subface[2] =
					face2.subface[3] = -1;
			fc[curFaceId + 1] = face2;
			face3.bsplit = false;
			face3.subface[0] = face3.subface[1] = face3.subface[2] =
					face3.subface[3] = -1;
			fc[curFaceId + 2] = face3;
			face4.bsplit = false;
			face4.subface[0] = face4.subface[1] = face4.subface[2] =
					face4.subface[3] = -1;
			fc[curFaceId + 3] = face4;
			//====================================

			if (false
					&& (!CheckOneFace_vertices(m, eg, &face1, b, d, step,
							bForward, tau)
							|| !CheckOneFace_vertices(m, eg, &face2, b, d, step,
									bForward, tau)
							|| !CheckOneFace_vertices(m, eg, &face3, b, d, step,
									bForward, tau)
							|| !CheckOneFace_vertices(m, eg, &face4, b, d, step,
									bForward, tau))) {
				int tau = curtau;
				int currSize = num_vertex;
				string dataName = "face";
				fs->save_Quad_One_Face_parent(m, face, eg, *d, dataName,
						num_vertex, 1, curtau, 2, curtau * 2);

				//				fs->save_Quad_One_Face_parent(m, fc[fc[face->parent].parent],
				//						eg, *_d, dataName, currSize, 3, tau, 3, tau * 3);

				fs->save_Quad_One_Face(m, face1, eg, Trace, *d, dataName,
						num_vertex, 1, tau, 1, tau * 2);
				fs->save_Quad_One_Face(m, face2, eg, Trace, *d, dataName,
						currSize, 2, tau, 2, tau * 2);
				fs->save_Quad_One_Face(m, face3, eg, Trace, *d, dataName,
						currSize, 3, tau, 3, tau * 2);
				fs->save_Quad_One_Face(m, face4, eg, Trace, *d, dataName,
						currSize, 4, tau, 4, tau * 2);

				fs->start_save_Quad_One_Face(m, face, eg, *d, dataName,
						currSize, 1, tau, 1);

				fs->start_save_Quad_One_Face(m, face1, eg, *d, dataName,
						currSize, 1, tau, 1);
				fs->start_save_Quad_One_Face(m, face2, eg, *d, dataName,
						currSize, 2, tau, 2);
				fs->start_save_Quad_One_Face(m, face3, eg, *d, dataName,
						currSize, 3, tau, 3);
				fs->start_save_Quad_One_Face(m, face4, eg, *d, dataName,
						currSize, 4, tau, 4);

				int cnt = 0;

				ASF_vertex oa[9];
				oa[cnt] = m[face.vertexes[0]];
				oa[cnt + 1] = m[face.vertexes[1]];
				oa[cnt + 2] = m[face.vertexes[2]];
				oa[cnt + 3] = m[face.vertexes[3]];
				oa[cnt + 4] = m[face.center];
				oa[cnt + 5] = m[edge1_1.v2];
				oa[cnt + 6] = m[edge2_1.v2];
				oa[cnt + 7] = m[edge3_1.v2];
				oa[cnt + 8] = m[edge4_1.v2];

				fs->start_save_Quad_FaceWithSeedPoits(m, fc, oa, eg, *d,
						dataName, currSize, 1, cnt, tau, 1);

				//				ASF_vertex test[2];
				//				test[0] = oa[0];
				//				fs->save_Quad_FaceWithSeedPoits(m, fc, oa, eg, *_d, dataName,
				//						currSize, 1, cnt + 4, tau, 1);

				generate_streamlines_sparsely(oa, bForward, 1, m_x1, m_y1, m_z1,
						b, d, step, d->x, tau + 2, 9);

				int* va = new int[9];
				for (int i = 0; i < 9; i++)
					va[i] = oa[i].getRange();

				fs->save_Voxel(m, fc, va, eg, *d, "face85", currSize, 9, tau, 1,
						tau * 5);
				//			int temp = 0;
				//			/*if(eg[105348].E2F[0] == 101217 || eg[105348].E2F[0] == 63417)//*/if (eg[101385].E2F[0]
				//					== 101217)
				//				temp = 1;
				//			subFace[cfnum * 4] = face1;
				//			subFace[cfnum * 4 + 1] = face2;
				//			subFace[cfnum * 4 + 2] = face3;
				//			subFace[cfnum * 4 + 3] = face4;
				//			originalFace[cfnum] = face;
				//			disconnectedFace_o[cfnum++] = face;
				//
				//			if (row == 252) {
				//				int tau = 1;
				//				string dataName = "face85";
				//				fs->save_Quad_One_Face(m, face1, eg, *d, dataName, num_vertex,
				//						1, tau, 2, tau * 2);
				//				fs->save_Quad_One_Face(m, face2, eg, *d, dataName, num_vertex,
				//						1, tau, 2, tau * 2);
				//				fs->save_Quad_One_Face(m, face3, eg, *d, dataName, num_vertex,
				//						1, tau, 2, tau * 2);
				//				fs->save_Quad_One_Face(m, face4, eg, *d, dataName, num_vertex,
				//						1, tau, 2, tau * 2);
			}
			//		} else
			//			wholeface[cfwholenum++] = face;
			int temp = 0;
			//		if (eg[105348].E2F[0] == 101217 || eg[105348].E2F[0] == 63417) //if(eg[101385].E2F[0] == 101217)
			//			temp = 1;

		}
	}

	stringstream ss;
	ss << whichData;
//	string str = ss.str() + "original";
//	fs->save_Quad_Face_original(m, originalFace, eg, *d, str, num_vertex, cfnum,
//			curtau, splitnumInCurtau);
//	str = ss.str() + "original_children";
//	fs->save_Quad_Face_original(m, subFace, eg, *d, str, num_vertex, cfnum * 4,
//			curtau, splitnumInCurtau);
	string str = ss.str();
//	fs->save_Quad_Face(m, disconnectedFace_o, eg, *d, str, num_vertex, cfnum,
//			curtau, splitnumInCurtau);
//	str = ss.str();			// + "childern";
//	fs->save_Quad_Face(m, subFace, eg, *d, str, num_vertex, cfnum * 4, curtau,
//			splitnumInCurtau);
//	str = ss.str() + "rest";
//	fs->save_Quad_Face(m, wholeface, eg, *d, str, num_vertex, cfwholenum,
//			curtau, splitnumInCurtau);

}

void ADPAlgorithm2::SplitFace(ASF_vertex *m, fFace* fc, fEdge* eg,
		uint32_t* Fr_Face, uint32_t* Fe_Face, uint32_t* Fe, Dimension* d,
		Boundary* b, Point* step, bool bForward, uint32_t num_vertex,
		uint32_t num_edges, uint32_t num_faces, int curtau,
		int splitnumInCurtau, int whichData, float* m_x1, float* m_y1,
		float* m_z1, int tau) {

	fFace* disconnectedFace_o = new fFace[Fr_Face[num_faces]];
	fFace* originalFace = new fFace[Fr_Face[num_faces]];
	fFace* wholeface = new fFace[num_faces];
	fFace* subFace = new fFace[Fr_Face[num_faces] * 4];

	int cfnum = 0;
	int cfwholenum = 0;
	for (uint32_t row = 0; row < num_faces; row++)

	{

		if (row >= num_faces)
			return;

		fFace face = fc[row];

		if (Fe_Face[row] == 1 && face.isInBoundary()) {

			if (!eg[face.edge[0]].isInBoundary()
					|| !eg[face.edge[1]].isInBoundary()
					|| !eg[face.edge[2]].isInBoundary()
					|| !eg[face.edge[3]].isInBoundary()) {
				fc[row].unsetInBoundary();
				continue;
			}
			//	face.level = face.level + 1;

//			if (row == 252 || row == 4251)
//				error(1);

			if (face.center == 0) {
				face = fc[row];
				continue;
			}
			int curVertexId = face.center;

			int curEdgeId = num_edges + (Fr_Face[row] * 4);
			int curFaceId = num_faces + (Fr_Face[row] * 4);
			if (fabs(curEdgeId - 15840) < 4) //|| row == 101217)
				face = fc[row];

			fEdge edge1 = eg[face.edge[0]];
			fEdge edge1_1 = eg[face.edge[0]];
			fEdge edge1_2 = eg[face.edge[0]];
			FindSubEdge(eg, face, &edge1_1, &edge1_2, 0);

			fEdge edge2 = eg[face.edge[1]];
			fEdge edge2_1 = eg[face.edge[1]];
			fEdge edge2_2 = eg[face.edge[1]];

			FindSubEdge(eg, face, &edge2_1, &edge2_2, 1);

			fEdge edge3 = eg[face.edge[2]];
			fEdge edge3_1 = edge3;
			fEdge edge3_2 = eg[face.edge[2]];
			FindSubEdge(eg, face, &edge3_1, &edge3_2, 2);

			fEdge edge4 = eg[face.edge[3]];
			fEdge edge4_1 = eg[face.edge[3]];
			fEdge edge4_2 = eg[face.edge[3]];
			FindSubEdge(eg, face, &edge4_1, &edge4_2, 3);

			int co = 0;

			fFace face1 = face;
			fFace face2 = face;
			fFace face3 = face;
			fFace face4 = face;
			fEdge edge1_c = edge1_1;
			edge1_c.v1 = edge1_1.v2;
			edge1_c.v2 = curVertexId;
			edge1_c.subedge[0] = -1;
			edge1_c.subedge[1] = -1;
			edge1_c.E2F[0] = curFaceId;
			edge1_c.E2F[1] = curFaceId + 1;
			edge1_c.E2F[2] = -1;
			edge1_c.E2F[3] = -1;
			edge1_c.bsplit = false;
			edge1_c.next = 0;
			edge1_c.level = face.level + 1;
			edge1_c.parent = -1;
			eg[curEdgeId + co] = edge1_c;

			bool bfound = false;
			for (int k = 0; k < 4; k++) {
				if (edge1_1.E2F[k] == row) {
					eg[edge1.subedge[0]].E2F[k] = curFaceId;
					bfound = true;

				}
			}
			if (!bfound)
				error(0);

			if (edge1_1.subedge[0] > 0)
				UpdatesubEdge(eg, edge1_1, row, curFaceId, curFaceId);

			bfound = false;
			for (int k = 0; k < 4; k++) {
				if (edge1_2.E2F[k] == row) {
					eg[edge1.subedge[1]].E2F[k] = curFaceId + 1;
					bfound = true;
					break;
				}
			}
			if (!bfound)
				error(0);

			if (edge1_2.subedge[0] > 0)
				UpdatesubEdge(eg, edge1_2, row, curFaceId + 1, curFaceId + 1);

			//edge1_1 = eg[edge1.subedge[0]];
			//edge1_2 = eg[edge1.subedge[1]];
			//eg[edge1.subedge[0]] = edge1_1;

			int c1 = curEdgeId + co;

			co++;
			//======================================================
			fEdge edge2_c = edge2_1;
			edge2_c.v2 = edge2_1.v2;
			edge2_c.v1 = curVertexId;

			edge2_c.E2F[0] = curFaceId + 1;
			edge2_c.E2F[1] = curFaceId + 2;
			edge2_c.E2F[2] = -1;
			edge2_c.E2F[3] = -1;
			edge2_c.next = 0;
			edge2_c.bsplit = false;
			edge2_c.subedge[0] = -1;
			edge2_c.subedge[1] = -1;
			edge2_c.parent = -1;
			bfound = false;
			for (int k = 0; k < 4; k++) {
				if (edge2_1.E2F[k] == row) {
					eg[edge2.subedge[0]].E2F[k] = curFaceId + 1;
					bfound = true;

				}
			}
			if (!bfound)
				error(0);

			if (edge2_1.subedge[0] > 0)
				UpdatesubEdge(eg, edge2_1, row, curFaceId + 1, curFaceId + 1);

			bfound = false;
			for (int k = 0; k < 4; k++) {
				if (edge2_2.E2F[k] == row) {
					eg[edge2.subedge[1]].E2F[k] = curFaceId + 2;
					bfound = true;

				}
			}
			if (!bfound)
				error(0);

			if (edge2_2.subedge[0] > 0)
				UpdatesubEdge(eg, edge2_2, row, curFaceId + 2, curFaceId + 2);

			//	UpdatesubEdge( eg, edge2,  row, curFaceId + 1, curFaceId + 2);
			//			edge2_1 = eg[edge2.subedge[0]];
			//			edge2_2 = eg[edge2.subedge[1]];

			int c2 = curEdgeId + co;
			eg[curEdgeId + co] = edge2_c;

			co++;
			//========================================================
			//bool bfound = false;
			fEdge edge3_c = edge3_1;
			edge3_c.v2 = edge3_1.v2;
			edge3_c.v1 = curVertexId;
			edge3_c.subedge[0] = -1;
			edge3_c.subedge[1] = -1;
			edge3_c.E2F[0] = curFaceId + 2;
			edge3_c.E2F[1] = curFaceId + 3;
			edge3_c.E2F[2] = -1;
			edge3_c.E2F[3] = -1;
			edge3_c.next = 0;
			edge3_c.bsplit = false;
			edge3_c.parent = -1;
			eg[curEdgeId + co] = edge3_c;
			for (int k = 0; k < 4; k++) {
				if (edge3_1.E2F[k] == row) {
					eg[edge3.subedge[0]].E2F[k] = curFaceId + 3;
					bfound = true;
					break;
				}
			}
			if (!bfound)
				error(0);

			if (edge3_1.subedge[0] > 0)
				UpdatesubEdge(eg, edge3_1, row, curFaceId + 3, curFaceId + 3);

			bfound = false;

			for (int k = 0; k < 4; k++) {
				if (edge3_2.E2F[k] == row) {
					eg[edge3.subedge[1]].E2F[k] = curFaceId + 2;
					bfound = true;
					break;
				}
			}
			if (!bfound)
				error(0);

			if (edge3_2.subedge[0] > 0)
				UpdatesubEdge(eg, edge3_2, row, curFaceId + 2, curFaceId + 2);
			bfound = false;
			//eg[edge3.subedge[0]] = edge3_1;
			//eg[edge3.subedge[1]] = edge3_2;

			int c3 = curEdgeId + co;
			eg[curEdgeId + co] = edge3_c;

			co++;
			//=======================================================

			fEdge edge4_c = edge4_1;
			edge4_c.v1 = edge4_1.v2;
			edge4_c.v2 = curVertexId;
			edge4_c.subedge[0] = -1;
			edge4_c.subedge[1] = -1;
			edge4_c.E2F[0] = curFaceId;
			edge4_c.E2F[1] = curFaceId + 3;
			edge4_c.E2F[2] = -1;
			edge4_c.E2F[3] = -1;
			edge4_c.next = 0;
			edge4_c.bsplit = false;
			edge4_c.parent = -1;
			eg[curEdgeId + co] = edge4_c;
			for (int k = 0; k < 4; k++) {
				if (edge4_1.E2F[k] == row) {
					eg[edge4.subedge[0]].E2F[k] = curFaceId;
					bfound = true;
					break;
				}
			}
			if (!bfound)
				error(0);

			if (edge4_1.subedge[0] > 0)
				UpdatesubEdge(eg, edge4_1, row, curFaceId, curFaceId);

			bfound = false;

			for (int k = 0; k < 4; k++) {
				if (edge4_2.E2F[k] == row) {
					eg[edge4.subedge[1]].E2F[k] = curFaceId + 3;
					bfound = true;
					break;
				}
			}
			if (!bfound)
				error(0);

			if (edge4_2.subedge[0] > 0)
				UpdatesubEdge(eg, edge4_2, row, curFaceId + 3, curFaceId + 3);

			bfound = false;

			//eg[edge4.subedge[0]] = edge4_1;
			//eg[edge4.subedge[1]] = edge4_2;

			int c4 = curEdgeId + co;
			eg[curEdgeId + co] = edge4_c;

			co++;

			//=====================================================
			face1.edge[0] = edge1.subedge[0];
			face1.edge[1] = curEdgeId;
			face1.edge[2] = curEdgeId + 3;
			face1.edge[3] = edge4.subedge[0];
			face1.vertexes[0] = face.vertexes[0];
			face1.vertexes[1] = edge1_1.v2;
			face1.vertexes[2] = curVertexId;
			face1.vertexes[3] = edge4_1.v2;
			face1.center = 0;
			face1.level = face.level + 1;
			face1.parent = row;
			face1.F2V = face.F2V;
			face1.isInBoundary();
			//==================================
			face2.edge[0] = edge1.subedge[1];
			face2.edge[1] = edge2.subedge[0];
			face2.edge[2] = curEdgeId + 1;
			face2.edge[3] = curEdgeId;
			face2.vertexes[0] = edge1_1.v2;
			face2.vertexes[1] = face.vertexes[1];
			face2.vertexes[2] = edge2_1.v2;
			face2.vertexes[3] = curVertexId;
			face2.center = 0;
			face2.level = face.level + 1;
			face2.parent = row;
			//====================================
			face3.edge[0] = curEdgeId + 1;
			face3.edge[1] = edge2.subedge[1];
			face3.edge[2] = edge3.subedge[1];
			face3.edge[3] = curEdgeId + 2;
			face3.vertexes[0] = curVertexId;
			face3.vertexes[1] = edge2_1.v2;
			;
			face3.vertexes[2] = face.vertexes[2];
			face3.vertexes[3] = edge3_1.v2;
			face3.center = 0;
			face3.level = face.level + 1;
			face3.parent = row;
			//====================================

			face4.edge[0] = curEdgeId + 3;
			face4.edge[1] = curEdgeId + 2;
			;
			face4.edge[2] = edge3.subedge[0];
			face4.edge[3] = edge4.subedge[1];
			face4.vertexes[0] = edge4_1.v2;
			face4.vertexes[1] = curVertexId;
			face4.vertexes[2] = edge3_1.v2;
			face4.vertexes[3] = face.vertexes[3];
			face4.center = 0;
			face4.level = face.level + 1;
			face4.parent = row;

			fc[row].bsplit = true;
			fc[row].subface[0] = curFaceId;
			fc[row].subface[1] = curFaceId + 1;
			fc[row].subface[2] = curFaceId + 2;
			fc[row].subface[3] = curFaceId + 3;

			face_correctnessCheck(face, eg);
			face_correctnessCheck(face1, eg);
			face_correctnessCheck(face2, eg);
			face_correctnessCheck(face3, eg);
			face_correctnessCheck(face4, eg);

			//fc[row] = face;
			face1.subface[0] = face1.subface[1] = face1.subface[2] =
					face1.subface[3] = -1;
			face1.bsplit = false;
			fc[curFaceId] = face1;
			face2.bsplit = false;
			face2.subface[0] = face2.subface[1] = face2.subface[2] =
					face2.subface[3] = -1;
			fc[curFaceId + 1] = face2;
			face3.bsplit = false;
			face3.subface[0] = face3.subface[1] = face3.subface[2] =
					face3.subface[3] = -1;
			fc[curFaceId + 2] = face3;
			face4.bsplit = false;
			face4.subface[0] = face4.subface[1] = face4.subface[2] =
					face4.subface[3] = -1;
			fc[curFaceId + 3] = face4;
			//====================================

			if (false
					&& (!CheckOneFace_vertices(m, eg, &face1, b, d, step,
							bForward, tau)
							|| !CheckOneFace_vertices(m, eg, &face2, b, d, step,
									bForward, tau)
							|| !CheckOneFace_vertices(m, eg, &face3, b, d, step,
									bForward, tau)
							|| !CheckOneFace_vertices(m, eg, &face4, b, d, step,
									bForward, tau))) {
				int tau = curtau;
				int currSize = num_vertex;
				string dataName = "face";
				fs->save_Quad_One_Face_parent(m, face, eg, *d, dataName,
						num_vertex, 1, curtau, 2, curtau * 2);

				//				fs->save_Quad_One_Face_parent(m, fc[fc[face->parent].parent],
				//						eg, *_d, dataName, currSize, 3, tau, 3, tau * 3);

				fs->save_Quad_One_Face(m, face1, eg, Trace, *d, dataName,
						num_vertex, 1, tau, 1, tau * 2);
				fs->save_Quad_One_Face(m, face2, eg, Trace, *d, dataName,
						currSize, 2, tau, 2, tau * 2);
				fs->save_Quad_One_Face(m, face3, eg, Trace, *d, dataName,
						currSize, 3, tau, 3, tau * 2);
				fs->save_Quad_One_Face(m, face4, eg, Trace, *d, dataName,
						currSize, 4, tau, 4, tau * 2);

				fs->start_save_Quad_One_Face(m, face, eg, *d, dataName,
						currSize, 1, tau, 1);

				fs->start_save_Quad_One_Face(m, face1, eg, *d, dataName,
						currSize, 1, tau, 1);
				fs->start_save_Quad_One_Face(m, face2, eg, *d, dataName,
						currSize, 2, tau, 2);
				fs->start_save_Quad_One_Face(m, face3, eg, *d, dataName,
						currSize, 3, tau, 3);
				fs->start_save_Quad_One_Face(m, face4, eg, *d, dataName,
						currSize, 4, tau, 4);

				int cnt = 0;

				ASF_vertex oa[9];
				oa[cnt] = m[face.vertexes[0]];
				oa[cnt + 1] = m[face.vertexes[1]];
				oa[cnt + 2] = m[face.vertexes[2]];
				oa[cnt + 3] = m[face.vertexes[3]];
				oa[cnt + 4] = m[face.center];
				oa[cnt + 5] = m[edge1_1.v2];
				oa[cnt + 6] = m[edge2_1.v2];
				oa[cnt + 7] = m[edge3_1.v2];
				oa[cnt + 8] = m[edge4_1.v2];

				fs->start_save_Quad_FaceWithSeedPoits(m, fc, oa, eg, *d,
						dataName, currSize, 1, cnt, tau, 1);

				//				ASF_vertex test[2];
				//				test[0] = oa[0];
				//				fs->save_Quad_FaceWithSeedPoits(m, fc, oa, eg, *_d, dataName,
				//						currSize, 1, cnt + 4, tau, 1);

				generate_streamlines_sparsely(oa, bForward, 1, m_x1, m_y1, m_z1,
						b, d, step, d->x, tau + 2, 9);

				int* va = new int[9];
				for (int i = 0; i < 9; i++)
					va[i] = oa[i].getRange();

				fs->save_Voxel(m, fc, va, eg, *d, "face85", currSize, 9, tau, 1,
						tau * 5);
				//			int temp = 0;
				//			/*if(eg[105348].E2F[0] == 101217 || eg[105348].E2F[0] == 63417)//*/if (eg[101385].E2F[0]
				//					== 101217)
				//				temp = 1;
				//			subFace[cfnum * 4] = face1;
				//			subFace[cfnum * 4 + 1] = face2;
				//			subFace[cfnum * 4 + 2] = face3;
				//			subFace[cfnum * 4 + 3] = face4;
				//			originalFace[cfnum] = face;
				//			disconnectedFace_o[cfnum++] = face;
				//
				//			if (row == 252) {
				//				int tau = 1;
				//				string dataName = "face85";
				//				fs->save_Quad_One_Face(m, face1, eg, *d, dataName, num_vertex,
				//						1, tau, 2, tau * 2);
				//				fs->save_Quad_One_Face(m, face2, eg, *d, dataName, num_vertex,
				//						1, tau, 2, tau * 2);
				//				fs->save_Quad_One_Face(m, face3, eg, *d, dataName, num_vertex,
				//						1, tau, 2, tau * 2);
				//				fs->save_Quad_One_Face(m, face4, eg, *d, dataName, num_vertex,
				//						1, tau, 2, tau * 2);
			}
			//		} else
			//			wholeface[cfwholenum++] = face;
			int temp = 0;
			//		if (eg[105348].E2F[0] == 101217 || eg[105348].E2F[0] == 63417) //if(eg[101385].E2F[0] == 101217)
			//			temp = 1;

		}
	}

	stringstream ss;
	ss << whichData;
//	string str = ss.str() + "original";
//	fs->save_Quad_Face_original(m, originalFace, eg, *d, str, num_vertex, cfnum,
//			curtau, splitnumInCurtau);
//	str = ss.str() + "original_children";
//	fs->save_Quad_Face_original(m, subFace, eg, *d, str, num_vertex, cfnum * 4,
//			curtau, splitnumInCurtau);
	string str = ss.str();
//	fs->save_Quad_Face(m, disconnectedFace_o, eg, *d, str, num_vertex, cfnum,
//			curtau, splitnumInCurtau);
//	str = ss.str();			// + "childern";
//	fs->save_Quad_Face(m, subFace, eg, *d, str, num_vertex, cfnum * 4, curtau,
//			splitnumInCurtau);
//	str = ss.str() + "rest";
//	fs->save_Quad_Face(m, wholeface, eg, *d, str, num_vertex, cfwholenum,
//			curtau, splitnumInCurtau);

}

bool ADPAlgorithm2::checkEdge(ASF_vertex vertex1, ASF_vertex vertex2,
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

void ADPAlgorithm2::CheckRemainingEdge(ASF_vertex*m, fFace*fc, fEdge* eg,
		uint32_t*Fe_Edge, uint32_t* Fe_Face, uint32_t* Fr, Dimension* d,
		Boundary* b, Point* step, bool bForward, uint32_t num_rows) {

	for (uint32_t row = 0; row < num_rows; row++) {
		if (Fe_Face[row] == 1) {

			fFace face = fc[row];
			int facelevel = face.level;
//			fEdge eg1 = eg[face.edge[0]];
//			fEdge eg2 = eg[face.edge[1]];
//			fEdge eg3 = eg[face.edge[2]];
//			fEdge eg4 = eg[face.edge[3]];

			if (!eg[face.edge[0]].bsplit)
				Fr[face.edge[0]] = 1;
			if (!eg[face.edge[1]].bsplit)
				Fr[face.edge[1]] = 1;
			if (!eg[face.edge[2]].bsplit)
				Fr[face.edge[2]] = 1;
			if (!eg[face.edge[3]].bsplit)
				Fr[face.edge[3]] = 1;
		}
	}
	return;

	bool bcollapse = false;

	for (uint32_t row = 0; row < num_rows; row++) {

		if (row >= num_rows)
			return;

		if (row == 229)
			Fr[row] = 0;
		Fr[row] = 0;
		fEdge edge = eg[row];
		if (Fe_Edge[row] == 0 && edge.isInBoundary() && !edge.bsplit) {

			int co = 0;

			if (edge.bsplit)
				continue;

			if (!edge.isInBoundary())
				continue;
			if (edge.level < fc[edge.E2F[0]].level
					|| edge.level < fc[edge.E2F[1]].level
					|| edge.level < fc[edge.E2F[2]].level
					|| edge.level < fc[edge.E2F[3]].level)
				error(0);
			for (int i = 0; i < 4; i++) {
				if (edge.E2F[i] > 0 && Fe_Face[edge.E2F[i]] == 1) {
					fFace face = fc[edge.E2F[i]];
					if (face.level < edge.level) {
						Fr[row] = 0;

					} else if (face.level >= edge.level) {
						Fr[row] = 1;
						if (face.edge[0] != row && face.edge[1] != row
								&& face.edge[2] != row && face.edge[3] != row)
							error(1);

						break;
					}

				}
			}

			//			if (edge.E2F[1] > 0 && Fe_Face[edge.E2F[1]] == 1)
			//				Fr[row] = 1;
			//
			//			if (edge.E2F[2] > 0 && Fe_Face[edge.E2F[2]] == 1)
			//				Fr[row] = 1;
			//
			//			if (edge.E2F[3] > 0 && Fe_Face[edge.E2F[3]] == 1)
			//				Fr[row] = 1;

		}

	}

	return;
}

void ADPAlgorithm2::CheckNeighborhood(ASF_vertex*m, fEdge* eg, uint32_t*Fr,
		Dimension* d, Boundary* b, Point* step, bool bForward,
		uint32_t original_num_rows, uint32_t num_rows, uint32_t level,
		int tau) {

	for (uint32_t row = 0; row < num_rows; row++) {

		if (row >= num_rows)
			return;

		Fr[row] = 0;
		int co = 0;
		if (row == 11705)
			co = 0;

		fEdge edge = eg[row];
		ASF_vertex vertex1 = m[edge.v1];
		ASF_vertex vertex2 = m[edge.v2];
		//printf("---- %f%f%f \n", vertex1.p.x, vertex1.p.y, vertex1.p.z);

		if (vertex1.e.dist(vertex2.e) > step->x) {
			co = 0;
		}
		if (!edge.isInBoundary())
			continue;
		if (edge.bsplit)
			continue;
		//	edge.bsplit = false;

		if (!vertex1.isInBoundary() || !vertex2.isInBoundary()) {
			//eg[row].unsetInBoundary();
			edge.unsetInBoundary();
			eg[row] = edge;
			continue;
		}

		if (!vertex1.checkInBoundary_StartPoint(b)
				|| !vertex2.checkInBoundary_StartPoint(b)) {
			//	printf("%f%f%f \n", vertex1.p.x, vertex1.p.y, vertex1.p.z);
			//eg[row].unsetInBoundary();
			edge.unsetInBoundary();
			eg[row] = edge;

			continue;
		}

		if (!vertex1.checkInBoundary(b) || !vertex2.checkInBoundary(b)) {
			edge.unsetInBoundary();
			eg[row] = edge;

			continue;

		}

		if (abs((int) vertex1.getOldRange() - (int) vertex2.getOldRange()) > 1
				&& abs(
						(int) vertex1.getOldRange()
								- (int) vertex2.getOldRange()) != d->x
				&& abs(
						(int) vertex1.getOldRange()
								- (int) vertex2.getOldRange()) != (d->x * d->y))
			error(0);

		/*	if (vertex1.getOldRange() % (d->x - 1) == 0 && (vertex2.getOldRange() % d->x) == 0)
		 continue;
		 if (vertex1.getOldRange() % (d->x*d->y - 1 )==0 && (vertex2.getOldRange() % d->x*d->y) == 0)
		 continue;*/

		float dist = sqrt(
				(vertex1.e.x - vertex2.e.x) * (vertex1.e.x - vertex2.e.x)
						+ (vertex1.e.y - vertex2.e.y)
								* (vertex1.e.y - vertex2.e.y)
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
	return;
}

bool ADPAlgorithm2::DivideEdges(ASF_vertex* v1, ASF_vertex* v2, uint32_t range1,
		Boundary* _b, Dimension* _d, Point* step, float*m_x1, float* m_y1,
		float*m_z1, int currentXDim, uint32_t faceNum, bool bForward,
		int divider, int tau, ASF_vertex* otv, int curvertexid, int v1i,
		int v2i, bool bEstimated) {

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

		if (!bEstimated) {
			float p1[3];
			p1[0] = vcenter.p.x;
			p1[1] = vcenter.p.y;
			p1[2] = vcenter.p.z;

			float ep[3];
//			generalstreamlineTracing_single(p1, bForward, ep, m_x1, m_y1, m_z1,
//					_b, _d, step, currentXDim, tau + 1);
//
//			vcenter.e.x = ep[0];
//			vcenter.e.y = ep[1];
//			vcenter.e.z = ep[2];

			for (int k = 0; k < tau + 1; k++) {
				generalstreamlineTracing_single(p1, bForward, ep, m_x1, m_y1,
						m_z1, _b, _d, step, currentXDim, 1);
				p1[0] = Trace[curvertexid][k].x = ep[0];
				p1[1] = Trace[curvertexid][k].y = ep[1];
				p1[2] = Trace[curvertexid][k].z = ep[2];
			}
			vcenter.e.x = ep[0];
			vcenter.e.y = ep[1];
			vcenter.e.z = ep[2];

			generalstreamlineTracing_single(p1, bForward, ep, m_x1, m_y1, m_z1,
					_b, _d, step, currentXDim, tau + 1);

		} else {
			vcenter.e.x = (t) * v1->e.x + (1 - t) * v2->e.x;
			vcenter.e.y = (t) * v1->e.y + (1 - t) * v2->e.y;
			vcenter.e.z = (t) * v1->e.z + (1 - t) * v2->e.z;

			//if (vcenter.type == 2)
			{
//				for (int j = 0; j < tau + 1; j++) {
//					Trace[curvertexid][j].x = (t * Trace[v1i][j].x)
//							+ (1 - t) * Trace[v2i][j].x;
//					Trace[curvertexid][j].y = (t * Trace[v1i][j].y)
//							+ (1 - t) * Trace[v2i][j].y;
//					Trace[curvertexid][j].z = (t * Trace[v1i][j].z)
//							+ (1 - t) * Trace[v2i][j].z;
//////
//				}
			}
		}
		//	vcenter.e = Trace[curvertexid][tau];

		//Trace[curvertexid][tau] = vcenter.e;
		//		AdvectParticle_estimated(&vcenter, m_x1, m_y1, m_z1, _b, _d, step,
		//				currentXDim, tau + 1, bForward);

		if (!vcenter.checkInBoundary(_b)) {
			*otv = vcenter;
			//return false;
			error(0);
			continue;
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

void ADPAlgorithm2::SplitRemainingEdge(ASF_vertex*m, fFace* fc, fEdge* eg,
		uint32_t* Fe_Edge, uint32_t* Fr_Edge, uint32_t* Fe_Face, Dimension* d,
		Boundary* b, Point* step, float*m_x1, float* m_y1, float*m_z1,
		int currentXDim, int tau, bool bForward, uint32_t num_vertex,
		uint32_t oldedgeSize, uint32_t num_edges) {
	int co = 0;
	for (uint32_t row = 0; row < oldedgeSize; row++) {

		if (row == 262680)
			co = 0;
		if (Fe_Edge[row] == 1) {

//			if (row == 3272)
//				error(0);

			Split_One_Edge(m, fc, eg, Fe_Edge, Fr_Edge, d, b, step, m_x1, m_y1,
					m_z1, currentXDim, tau, bForward, num_vertex, num_edges,
					row);
//			Split_One_Edge(m, fc, eg, Fe_Edge, Fr_Edge, d, b, step, m_x1, m_y1,
//					m_z1, currentXDim, tau, bForward, num_vertex, num_edges,
//					row);
		}
	}

}

//
//	int temp = 0;
//	for (uint32_t row = 0; row < num_edges; row++) {
//		if (row == 469)
//			temp = 0;
//		if (Fe_Edge[row] == 1) {
//
//			int curVertexId = num_vertex + Fr_Edge[row];
//			int curEdgeId = num_edges + Fr_Edge[row] * 2;
//
//			fEdge edge1 = eg[row];
//
//			if (edge1.bsplit) {
//				Fe_Edge[row] = 0;
//				continue;
//			}
//
//			if (Fe_Face[edge1.E2F[0]] == 1 || Fe_Face[edge1.E2F[1]] == 1
//					|| Fe_Face[edge1.E2F[2]] == 1
//					|| Fe_Face[edge1.E2F[3]] == 1) {
//
//				int co = 0;
//
//				fEdge edge1 = eg[row];
//				fEdge edge2 = eg[row];
//
//				ASF_vertex vertex1 = m[edge1.v1];
//				ASF_vertex vertex2 = m[edge1.v2];
//
//				if ((vertex1.p.x - vertex2.p.x) > step->x
//						|| (vertex1.p.y - vertex2.p.y) > step->y
//						|| (vertex1.p.z - vertex2.p.z) > step->z)
//					error(1);
//
//				ASF_vertex vedge;
//				bool bFound = false;
//				//for (int n = 2; n <= 5; n++) {
//					int n = 2;
//
//
//
//					//ASF_vertex* vedge_array = new ASF_vertex[n];
////					if (DivideEdges(&vertex1, &vertex2, vertex1.getOldRange(),
////							b, d, step, m_x1, m_y1, m_z1, currentXDim, 0,
////							bForward, n, tau, &vedge) == true) {
////						bFound = true;
////					}
//				//}
//
//				if (!bFound)
//					error(0);
//				vedge.left = edge1.v1;
//				vedge.right = edge1.v2;
//
//				edge1.v2 = curVertexId;
//				edge2.v1 = curVertexId;
//				//edge1.bsplit = true;
//				edge1.next = curEdgeId;
//				edge1.level = edge1.level + 1;
//				edge2.level = edge1.level;
//				edge2.Prev = row;
//				vedge.level = edge2.level;
//
//				float dist1 = vertex1.e.dist(vertex2.e);
//				float dist1_1 = vedge.e.dist(vertex1.e);
//				float dist1_2 = vedge.e.dist(vertex2.e);
//
//				if (dist1_2 >= dist1 || dist1_2 >= dist1)
//					error(1);
//
//				eg[row].bsplit = true;
//				eg[row].subedge[0] = curEdgeId;
//				eg[row].subedge[1] = curEdgeId + 1;
//
//				edge1.bsplit = false;
//				edge2.bsplit = false;
//				edge1.parent = row;
//				edge2.parent = row;
//				eg[curEdgeId] = edge1;
//
//				eg[curEdgeId + 1] = edge2;
//
//				m[curVertexId] = vedge;
//
//			}
//
//		}
//	}
//
//}

void ADPAlgorithm2::Split_One_Edge(ASF_vertex*m, fFace* fc, fEdge* eg,
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
				curVertexId, v1i, v2i, true);
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
	if (!vedge.checkInBoundary(b)) {
		edge1.unsetInBoundary();
		edge2.unsetInBoundary();
		eg[row].unsetInBoundary();
	}
	vedge.type = 1;
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
	if (dist1_1 == 0 || dist1 < exp(-4.0)) {
		edge1.unsetInBoundary();
		edge2.unsetInBoundary();
		eg[row].unsetInBoundary();

		eg[curEdgeId] = edge1;

		eg[curEdgeId + 1] = edge2;

		m[curVertexId] = vedge;
		return;

	}

	if (dist1_1 >= dist1 || dist1_2 >= dist1) {
		error(1);

//		ASF_vertex oa[3];
//
//		DivideEdges(&vertex1, &vertex2, vertex1.getOldRange(), b, d, step, m_x1,
//				m_y1, m_z1, currentXDim, 0, bForward, n, tau, &vedge,
//				curVertexId, v1i, v2i, true);
//		oa[0] = vertex1;
//		oa[1] = vertex2;
//		oa[2] = vedge;
//		generate_streamlines_sparsely(oa, bForward, 1, m_x1, m_y1, m_z1, b, d,
//				step, currentXDim, tau + 2, 3);
//		int indexarray[3];
//		indexarray[0] = v1i;
//		indexarray[1] = v2i;
//		indexarray[2] = curVertexId;
//
//		m[curVertexId] = vedge;
//		fs->Save_Streamlines_estimated(m, Trace, indexarray, bForward, "st1",
//				m_x1, m_y1, m_z1, b, d, step, currentXDim, tau + 2, 3);
//		fs->save_Quad_FaceWithSeedPoits(m, fc, oa, eg, Trace, *d, "face85",
//				num_vertex, 1, 3, tau, 1);
		//error(1);
	}

	eg[curEdgeId] = edge1;

	eg[curEdgeId + 1] = edge2;

	m[curVertexId] = vedge;

}

void ADPAlgorithm2::SplitEndEdge(ASF_vertex*m, fFace* fc, fEdge* eg,
		uint32_t* Fe_Edge, uint32_t* Fr_Edge, Dimension* d, Boundary* b,
		Point* step, float*m_x1, float* m_y1, float*m_z1, int currentXDim,
		bool bForward, uint32_t num_vertex, uint32_t num_edges, uint32_t level,
		int tau) {
	int co = 0;
	for (uint32_t row = 0; row < num_edges; row++) {

		if (tau == 100)
			co = 0;
		if (Fe_Edge[row] == 1) {

//			if (row == 3272)
//				error(0);
//				Split_One_Edge(m, fc, eg, Fe_Edge, Fr_Edge, d, b, step, m_x1, m_y1,
//						m_z1, currentXDim, tau, bForward, num_vertex, num_edges,
//						row);

			Split_One_Edge(m, fc, eg, Fe_Edge, Fr_Edge, d, b, step, m_x1, m_y1,
					m_z1, currentXDim, tau, bForward, num_vertex, num_edges,
					row);
		}
	}

}
bool ADPAlgorithm2::check_child_edge(fEdge edge, fEdge* eg, uint32_t*Fe_Edge) {
	fEdge edge1 = eg[edge.subedge[0]];
	fEdge edge2 = eg[edge.subedge[2]];

	bool b1 = false;
	bool b2 = false;
	if (Fe_Edge[edge.subedge[0]] == 1 || Fe_Edge[edge.subedge[1]] == 1)
		return true;
	else if (edge1.bsplit) {
		b1 = check_child_edge(edge1, eg, Fe_Edge);

	} else if (edge2.bsplit) {
		b2 = check_child_edge(edge2, eg, Fe_Edge);
	}
	if (b1 || b2)
		return true;
	return false;
}

void ADPAlgorithm2::CheckFace(ASF_vertex*m, fEdge* eg, fFace* fc,
		uint32_t*Fe_Edge, uint32_t* Fe_Face, Dimension* d, uint32_t num_face) {

	for (uint32_t row = 0; row < num_face; row++) {

		if (row >= num_face)
			return;

		fFace face = fc[row];
		if (row == 494)
			face = fc[row];

		if (face.bsplit) {
			face = fc[row];
			continue;
		}

		if (!face.isInBoundary()) {
			continue;
		}

		//		if (!m[face.vertexes[0]].checkInBoundary(b)
		//				|| !m[face.vertexes[1]].isInBoundary()
		//				|| !m[face.vertexes[2]].isInBoundary()
		//				|| !m[face.vertexes[3]].isInBoundary()) {
		//			fc[row].unsetInBoundary();
		//			continue;
		//
		//		}

		fEdge edge1 = eg[face.edge[0]];
		fEdge edge2 = eg[face.edge[1]];
		fEdge edge3 = eg[face.edge[2]];
		fEdge edge4 = eg[face.edge[3]];
		//if(eg[face.edge[0]].level > face.level || eg[face.edge[1]].level > face.level || eg[face.edge[2]].level > face.level || eg[face.edge[3]].level > face.level)
		if (Fe_Edge[face.edge[0]] == 1 || Fe_Edge[face.edge[1]] == 1
				|| Fe_Edge[face.edge[2]] == 1 || Fe_Edge[face.edge[3]] == 1) {
			Fe_Face[row] = 1;
		} else {
			if (edge1.bsplit) {
				if (check_child_edge(edge1, eg, Fe_Edge) == true)
					Fe_Face[row] = 1;
			}
			if (edge2.bsplit) {
				if (check_child_edge(edge2, eg, Fe_Edge) == true)
					Fe_Face[row] = 1;
			}
			if (edge3.bsplit) {
				if (check_child_edge(edge3, eg, Fe_Edge) == true)
					Fe_Face[row] = 1;
			}
			if (edge4.bsplit) {
				if (check_child_edge(edge4, eg, Fe_Edge) == true)
					Fe_Face[row] = 1;
			}
		}

		fc[row] = face;

	}

	return;
}

void ADPAlgorithm2::CheckFace_vertices(ASF_vertex*m, fEdge* eg, fFace* fc,
		uint32_t*Fe_Edge, uint32_t* Fe_Face, Boundary* b, Dimension* d,
		Point* step, bool bForward, uint32_t num_face, int tau) {

	for (uint32_t row = 0; row < num_face; row++) {

		if (row >= num_face)
			return;

		fFace face = fc[row];
		if (row == 442)
			face = fc[row];

		if (face.bsplit) {
			face = fc[row];
			continue;
		}

		if (!face.isInBoundary()) {
			continue;
		}

		if (!m[face.vertexes[0]].isInBoundary()
				|| !m[face.vertexes[1]].isInBoundary()
				|| !m[face.vertexes[2]].isInBoundary()
				|| !m[face.vertexes[3]].isInBoundary()) {
			fc[row].unsetInBoundary();
			continue;

		}

		if (!eg[face.edge[0]].isInBoundary() || !eg[face.edge[1]].isInBoundary()
				|| !eg[face.edge[2]].isInBoundary()
				|| !eg[face.edge[3]].isInBoundary()) {
			fc[row].unsetInBoundary();
			continue;

		}

		if (!m[face.vertexes[0]].checkInBoundary_StartPoint(b)
				|| !m[face.vertexes[1]].checkInBoundary_StartPoint(b)
				|| !m[face.vertexes[2]].checkInBoundary_StartPoint(b)
				|| !m[face.vertexes[3]].checkInBoundary_StartPoint(b)) {
			fc[row].unsetInBoundary();
			continue;

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
			continue;

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

	return;
}

bool ADPAlgorithm2::CheckOneFace_vertices(ASF_vertex* m, fEdge* eg, fFace* face,
		Boundary* b, Dimension* d, Point* step, bool bForward, int tau) {
	ASF_vertex v1 = m[face->vertexes[0]];
	ASF_vertex v2 = m[face->vertexes[1]];
	ASF_vertex v3 = m[face->vertexes[2]];
	ASF_vertex v4 = m[face->vertexes[3]];

	int v1i = face->vertexes[0];
	int v2i = face->vertexes[1];
	int v3i = face->vertexes[2];
	int v4i = face->vertexes[3];

	if (!checkEdge(v1, v2, b, d, step, bForward, tau, v1i, v2i)
			|| !checkEdge(v2, v3, b, d, step, bForward, tau, v1i, v2i)
			|| !checkEdge(v3, v4, b, d, step, bForward, tau, v1i, v2i)
			|| !checkEdge(v1, v4, b, d, step, bForward, tau, v1i, v2i))
		return true;

	return false;

}

float ADPAlgorithm2::ComputeJacobian(ASF_vertex*a, fFace face, float* m_x1,
		float* m_y1, float* m_z1, Point* step) {

	if (face.dir == 1)		// xy face
			{
		float vx1 = m_x1[face.vertexes[0]] * 0.5 + m_x1[face.vertexes[3]] * 0.5;
		float vx2 = m_x1[face.vertexes[0]] * 0.5 + m_x1[face.vertexes[1]] * 0.5;
		float vx3 =
				(m_x1[face.vertexes[1]] * 0.5 + m_x1[face.vertexes[2]] * 0.5);
		float vx4 =
				(m_x1[face.vertexes[3]] * 0.5 + m_x1[face.vertexes[2]] * 0.5);

		float vy1 =
				(m_y1[face.vertexes[0]] * 0.5 + m_y1[face.vertexes[3]] * 0.5);
		float vy2 =
				(m_y1[face.vertexes[0]] * 0.5 + m_y1[face.vertexes[1]] * 0.5);
		float vy3 =
				(m_y1[face.vertexes[1]] * 0.5 + m_y1[face.vertexes[2]] * 0.5);
		float vy4 =
				(m_y1[face.vertexes[3]] * 0.5 + m_y1[face.vertexes[2]] * 0.5);

		float J11 = (vx3 - vx1) / step->x;
		float J21 = (vy3 - vy1) / step->x;
		float J12 = (vx4 - vx2) / step->y;
		float J22 = (vy4 - vy2) / step->y;

		return J11 + J22;

	}

	if (face.dir == 2)		// xz face
			{
		float vx1 =
				(m_x1[face.vertexes[0]] * 0.5 + m_x1[face.vertexes[3]] * 0.5);
		float vx2 =
				(m_x1[face.vertexes[0]] * 0.5 + m_x1[face.vertexes[1]] * 0.5);
		float vx3 =
				(m_x1[face.vertexes[1]] * 0.5 + m_x1[face.vertexes[2]] * 0.5);
		float vx4 =
				(m_x1[face.vertexes[3]] * 0.5 + m_x1[face.vertexes[2]] * 0.5);

		float vz1 =
				(m_z1[face.vertexes[0]] * 0.5 + m_z1[face.vertexes[3]] * 0.5);
		float vz2 =
				(m_z1[face.vertexes[0]] * 0.5 + m_z1[face.vertexes[1]] * 0.5);
		float vz3 =
				(m_z1[face.vertexes[1]] * 0.5 + m_z1[face.vertexes[2]] * 0.5);
		float vz4 =
				(m_z1[face.vertexes[3]] * 0.5 + m_z1[face.vertexes[2]] * 0.5);

		float J11 = (vx3 - vx1) / step->x;
		float J21 = (vz3 - vz1) / step->x;
		float J12 = (vx4 - vx2) / step->y;
		float J22 = (vz4 - vz2) / step->y;

		return J11 + J22;

	}

	if (face.dir == 3)		// yz face
			{
		float vy1 =
				(m_y1[face.vertexes[0]] * 0.5 + m_y1[face.vertexes[3]] * 0.5);
		float vy2 =
				(m_y1[face.vertexes[0]] * 0.5 + m_y1[face.vertexes[1]] * 0.5);
		float vy3 =
				(m_y1[face.vertexes[1]] * 0.5 + m_y1[face.vertexes[2]] * 0.5);
		float vy4 =
				(m_y1[face.vertexes[3]] * 0.5 + m_y1[face.vertexes[2]] * 0.5);

		float vz1 =
				(m_z1[face.vertexes[0]] * 0.5 + m_z1[face.vertexes[3]] * 0.5);
		float vz2 =
				(m_z1[face.vertexes[0]] * 0.5 + m_z1[face.vertexes[1]] * 0.5);
		float vz3 =
				(m_z1[face.vertexes[1]] * 0.5 + m_z1[face.vertexes[2]] * 0.5);
		float vz4 =
				(m_z1[face.vertexes[3]] * 0.5 + m_z1[face.vertexes[2]] * 0.5);

		float J11 = (vy3 - vy1) / step->y;
		float J21 = (vz3 - vz1) / step->y;
		float J12 = (vy4 - vy2) / step->z;
		float J22 = (vz4 - vz2) / step->z;

		return J11 + J22;

	}

}

void ADPAlgorithm2::ComputeJacobianForEachFace(ASF_vertex*a, fFace* fc,
		int numface, float* m_x1, float* m_y1, float* m_z1, Boundary*b,
		Point* step, float& minJ, float& maxJ) {
	float minJacobian = 1000.;
	float maxJacobian = -1000.;

	for (int row = 0; row < numface; row++) {

		fFace face = fc[row];
		if (!a[face.vertexes[0]].checkInBoundary(b)
				|| !a[face.vertexes[1]].checkInBoundary(b)
				|| !a[face.vertexes[2]].checkInBoundary(b)
				|| !a[face.vertexes[3]].checkInBoundary(b))
			continue;
		face.Jacobian = ComputeJacobian(a, face, m_x1, m_y1, m_z1, step);
		if (face.Jacobian < minJacobian)
			minJacobian = face.Jacobian;
		if (face.Jacobian > maxJacobian)
			maxJacobian = face.Jacobian;

		fc[row] = face;
	}
	minJ = minJacobian;
	maxJ = maxJacobian;

}

void ADPAlgorithm2::Flow_Combinatorialization_C2(ASF_vertex* _a, fFace *fc,
		fEdge*eg, Boundary* _b, Dimension* _d, Point* step, int _tau,
		float* m_x1, float* m_y1, float* m_z1, ASF_vertex** oVertex,
		Edge ** oFc, uint32_t ** oFr, fFace** oFace, fEdge** oEdge,
		uint32_t * oRSize, uint32_t * oFaceSize, uint32_t * oEdgeSize,
		uint32_t whichData, int CSize, uint32_t _curEdgeSize,
		uint32_t _curFaceSize, int curXDimension, int sampleSeeds,
		int MorseLevel, bool bForward) {

//-----------GPU initialization---------------------------->

	uint32_t * To, *From;
	uint32_t * To2;

	printf("Vertices: %u\n", CSize - 2);

	uint32_t currTau = _tau;

	uint32_t currSize = CSize;
	uint32_t icycle = currTau;
//	if (MorseLevel == 3)
		generate_streamlines_sparsely_Ocean(_a, bForward, whichData, m_x1, m_y1,
				m_z1, _b, _d, step, curXDimension, 100, currSize);

//
//	generate_streamlines_sparsely(_a, bForward, whichData, m_x1, m_y1, m_z1, _b,
//			_d, step, curXDimension, currTau, currSize);

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
//currSize = currSize * 2;

//	generate_streamlines_sparsely(currA, bForward, whichData, m_x1, m_y1, m_z1,
//			_b, _d, step, curXDimension, 200, currSize);

	//Jacobianarray = new float[_d->x * _d->y * _d->z];

//	float minJacobian, maxJacobian;
//
//	ComputeJacobianForEachFace(currA, fc, curFaceSize, m_x1, m_y1, m_z1, _b,
//			step, minJacobian, maxJacobian);

	int c = currSize;
	int nf = curFaceSize;

	memset(Fe_Face, 0, curFaceSize * sizeof(uint32_t));
	memset(Fr_Face, 0, curFaceSize * sizeof(uint32_t));
	CheckFace_vertices(currA, eg, fc, Fe_Edge, Fe_Face, _b, _d, step, bForward,
			curFaceSize, currTau);

	Fr_Face[0] = 0;
	for (int i = 0; i <= curFaceSize; i++) {
		Fr_Face[i + 1] = Fr_Face[i] + Fe_Face[i];

	}

	int* voxelarray = new int[curFaceSize * 4];
	int voxelumber = 0;

	int faceSplittingCounter = 0;
	for (int ii = 0; ii < icycle; ii++) {

//		if (ii == 80 && !bForward)
//			error(0);

		int ct = 0;
		Tracing_c(currA, m_x1, m_y1, m_z1, _d, _b, step, bForward, whichData,
				curXDimension, 0, currSize, level, ii);

//		Tracing_c_estimated(currA, m_x1, m_y1, m_z1, _d, _b, step, bForward,
//				whichData, curXDimension, CSize, currSize, level);

		CheckNeighborhood(currA, eg, Fe_Edge, _d, _b, step, bForward, currSize,
				curEdgeSize, level, ii);

		Fr_Edge[0] = 0;
		for (int i = 0; i <= curEdgeSize; i++) {

			Fr_Edge[i + 1] = Fr_Edge[i] + Fe_Edge[i];
		}

		memset(Fe_Face, 0, curFaceSize * sizeof(uint32_t));
		memset(Fr_Face, 0, curFaceSize * sizeof(uint32_t));
		CheckFace_vertices(currA, eg, fc, Fe_Edge, Fe_Face, _b, _d, step,
				bForward, curFaceSize, ii);

		Fr_Face[0] = 0;
		voxelumber = 0;
		for (int i = 0; i <= curFaceSize; i++) {
			Fr_Face[i + 1] = Fr_Face[i] + Fe_Face[i];

		}

		//	if (Fr_Face[curFaceSize] > 200)
		{

//			fs->save_Quad_Face_Index(currA, fc, Fe_Face, eg, *_d, "difference",
//					currSize, curFaceSize, currTau, ii);
//			fs->save_Voxel(currA, fc, voxelarray, eg, *_d, "outer", currSize,
//					voxelumber, currTau, 10, 20);
//			fs->Save_Streamlines_estimated(currA, Trace, voxelarray, bForward,
//					"streamline", m_x1, m_y1, m_z1, _b, _d, step, curXDimension,
//					ii + 1, voxelumber);
		}

		printf(" Faces = %d \n", Fr_Face[curFaceSize]);

		faceSplittingCounter = 0;
		//======================================================================

		while (Fr_Face[curFaceSize] > 0 && faceSplittingCounter++ < 3) {
			//while (Fr_Edge[curEdgeSize] > 0) {
			//ct++;
			//printf("Edges = %d \n", Fr_Edge[curEdgeSize]);
			printf("Edges = %d \n", Fr_Face[curFaceSize]);
			int oldEdgeSize = curEdgeSize;

			//=========================================================================

			//if (Fr_Edge[curEdgeSize] > 0)
			SplitEndEdge(currA, fc, eg, Fe_Edge, Fr_Edge, _d, _b, step, m_x1,
					m_y1, m_z1, curXDimension, bForward, currSize, curEdgeSize,
					level, ii);

//
//			CheckNeighborhood(currA, eg, Fe_Edge, _d, _b, step, bForward,
//					currSize, curEdgeSize, level, ii);

//			Fr_Edge[0] = 0;
//			for (int i = 0; i <= curEdgeSize; i++) {
//
//				Fr_Edge[i + 1] = Fr_Edge[i] + Fe_Edge[i];
//
//			}

			currSize = currSize + Fr_Edge[curEdgeSize];
			curEdgeSize = curEdgeSize + Fr_Edge[curEdgeSize] * 2;

			memset(Fr_Edge2, 0, (oldEdgeSize) * sizeof(uint32_t));
			memset(Fe_Edge2, 0, (oldEdgeSize) * sizeof(uint32_t));

			CheckRemainingEdge(currA, fc, eg, Fe_Edge, Fe_Face, Fe_Edge2, _d,
					_b, step, bForward, curFaceSize);

			Fr_Edge2[0] = 0;
			for (int i = 0; i <= oldEdgeSize; i++) {
				Fr_Edge2[i + 1] = Fr_Edge2[i] + Fe_Edge2[i];
			}

			printf(" Edges = %d \n", Fr_Edge2[oldEdgeSize]);
			//
			//

			SplitRemainingEdge(currA, fc, eg, Fe_Edge2, Fr_Edge2, Fe_Face, _d,
					_b, step, m_x1, m_y1, m_z1, curXDimension, ii, bForward,
					currSize, oldEdgeSize, curEdgeSize);
			//
			//
			currSize = currSize + Fr_Edge2[oldEdgeSize];
			curEdgeSize = curEdgeSize + Fr_Edge2[oldEdgeSize] * 2;

			InsertFaceCenterWithMultipleSeeds(currA, fc, eg, Fr_Face, Fe_Face,
					Fe_Edge, _d, _b, step, m_x1, m_y1, m_z1, ii, curXDimension,
					bForward, currSize, curEdgeSize, curFaceSize);

			//			InsertFaceCenter(currA, fc, eg, Fr_Face, Fe_Face, Fe_Edge, _d, _b,
			//					step, m_x1, m_y1, m_z1, ii, curXDimension, bForward,
			//					currSize, curEdgeSize, curFaceSize);

			currSize = currSize + Fr_Face[curFaceSize];
			//
			int oldfacesize = curFaceSize;

			//if (MorseLevel == 1)
			SplitFace(currA, fc, eg, Fr_Face, Fe_Face, Fe_Edge, _d, _b, step,
					bForward, currSize, curEdgeSize, curFaceSize, ii, ct,
					whichData, m_x1, m_y1, m_z1, ii);
//			else
//				SplitFace2(currA, fc, eg, Fr_Face, Fe_Face, Fe_Edge, _d, _b,
//						step, bForward, currSize, curEdgeSize, curFaceSize, ii,
//						ct, whichData, m_x1, m_y1, m_z1, ii);

			//printf("\n %d %d \n", eg[105348].E2F[0], eg[101385].E2F[0]);
			curEdgeSize = curEdgeSize + Fr_Face[curFaceSize] * 4;

			curFaceSize = curFaceSize + Fr_Face[curFaceSize] * 4;

			memset(Fe_Face, 0, curFaceSize * sizeof(uint32_t));
			memset(Fr_Face, 0, curFaceSize * sizeof(uint32_t));
			CheckFace_vertices(currA, eg, fc, Fe_Edge, Fe_Face, _b, _d, step,
					bForward, curFaceSize, ii);

			Fr_Face[0] = 0;
			for (int i = 0; i <= curFaceSize; i++) {
				Fr_Face[i + 1] = Fr_Face[i] + Fe_Face[i];

			}

			printf(" Faces = %d \n", Fr_Face[curFaceSize]);

			//============================================================================

			memset(Fr_Edge, 0, (curEdgeSize) * sizeof(uint32_t));
			memset(Fe_Edge, 0, (curEdgeSize) * sizeof(uint32_t));

			CheckNeighborhood(currA, eg, Fe_Edge, _d, _b, step, bForward,
					currSize, curEdgeSize, level, ii);

			Fr_Edge[0] = 0;
			for (int i = 0; i <= curEdgeSize; i++) {

				Fr_Edge[i + 1] = Fr_Edge[i] + Fe_Edge[i];

			}

			//			if (Fr_Face[curFaceSize] > 0 && Fr_Edge[curEdgeSize] == 0)
			//				error(222);

			printf("Edges = %d \n", Fr_Edge[curEdgeSize]);

			int ct = 0;

		}

		printf(" \n %d \n", ii);
	}

	int maxlevel = 1;
	int prevmaxlevel = 1;
	int maxfacelevelindex = 0;
	for (int i = 0; i < curFaceSize; i++) {
		if (fc[i].level > maxlevel) {
			prevmaxlevel = maxfacelevelindex;
			maxlevel = fc[i].level;
			maxfacelevelindex = i;
		}

	}
//
//	int index = 1126;
//	fs->start_save_Quad_One_Face(currA, fc[index], eg, *_d, "face", currSize,
//			curFaceSize, currTau, 1);
//
//
//	ASF_vertex* oa = new ASF_vertex[10000];
//	int* indexarray = new int[10000];
//	int* intarray = new int[10000];
//	int* facearray = new int[10000];
//
//	facearray[0] = index;
//	int cv = 1;
//	//for (int i = 0; i < 4; i++) {
//	facearray[cv++] = fc[index].subface[0];
//	facearray[cv++] = fc[index].subface[1];
//	facearray[cv++] = fc[index].subface[2];
//	facearray[cv++] = fc[index].subface[3];
////	}
//
//	for (int i = 0; i < cv; i++) {
//		if (fc[facearray[i]].bsplit) {
//			facearray[cv++] = fc[facearray[i]].subface[0];
//			facearray[cv++] = fc[facearray[i]].subface[1];
//			facearray[cv++] = fc[facearray[i]].subface[2];
//			facearray[cv++] = fc[facearray[i]].subface[3];
//		}
//	}
//
//	int cr = 0;
//	for (int j = 0; j < cv; j++) {
////		fs->save_Quad_One_Face(currA, fc[facearray[j]], eg, Trace, *_d, "Lorenz",
////						currSize, j, currTau+ j, 1, currTau * 2);
//
//		intarray[cr++] = currA[fc[facearray[j]].vertexes[0]].getRange();
//		intarray[cr++] = currA[fc[facearray[j]].vertexes[1]].getRange();
//		intarray[cr++] = currA[fc[facearray[j]].vertexes[2]].getRange();
//		intarray[cr++] = currA[fc[facearray[j]].vertexes[3]].getRange();
//
//	}
//	fs->save_Voxel(currA, fc, intarray, eg, *_d, "face85", currSize, cr,
//			currTau, 1, 100);
//	cr = 0;
//	for (int j = 0; j < cv; j++) {
//
//		oa[cr] = currA[fc[facearray[j]].vertexes[0]];
//		indexarray[cr++] = fc[facearray[j]].vertexes[0];
//		oa[cr] = currA[fc[facearray[j]].vertexes[1]];
//		indexarray[cr++] = fc[facearray[j]].vertexes[1];
//		oa[cr] = currA[fc[facearray[j]].vertexes[2]];
//		indexarray[cr++] = fc[facearray[j]].vertexes[2];
//		oa[cr] = currA[fc[facearray[j]].vertexes[3]];
//		indexarray[cr++] = fc[facearray[j]].vertexes[3];
//
//	}
//
//	fs->Save_Streamlines_estimated(currA, Trace, indexarray, bForward, "st1",
//			m_x1, m_y1, m_z1, _b, _d, step, curXDimension, currTau, cr);
//	generate_streamlines_sparsely(oa, bForward, 1, m_x1, m_y1, m_z1, _b, _d,
//			step, curXDimension, currTau + 2, cr);
//
//	fs->save_Quad_FaceWithSeedPoits(currA, fc, oa, eg, Trace, *_d, "face85",
//			currSize, 1, cr, currTau, 1);
//	//for (int j = 0; j < 100; j += 10)
//
//	{
//		int j = currTau;
//
//		if (whichData == Lorenz) {
//			oa[0] = currA[fc[1126].vertexes[0]];
//			oa[1] = currA[fc[1126].vertexes[1]];
//			oa[2] = currA[fc[1126].vertexes[2]];
//			oa[3] = currA[fc[1126].vertexes[3]];
//
//			oa[4] = currA[fc[1126].center];
//			indexarray[0] = fc[1126].vertexes[0];
//			indexarray[1] = fc[1126].vertexes[1];
//			indexarray[2] = fc[1126].vertexes[2];
//			indexarray[3] = fc[1126].vertexes[3];
//			indexarray[4] = fc[1126].center;
//			fs->Save_Streamlines_estimated(currA, Trace, indexarray, bForward,
//					"st1", m_x1, m_y1, m_z1, _b, _d, step, curXDimension, j, 4);
//
//			generate_streamlines_sparsely(oa, bForward, 1, m_x1, m_y1, m_z1, _b,
//					_d, step, curXDimension, j + 2, 4);
//			fs->save_Quad_One_Face(currA, fc[1126], eg, Trace, *_d, "Lorenz",
//					currSize, 1, j + 1, 1, currTau * 2);
//
////			intarray[0] = currA[fc[1126].vertexes[0]].getRange_tau(&Trace[fc[1126].vertexes[0]][j], _b, step, _d,
////					j);
////			intarray[1] = currA[fc[1126].vertexes[1]].getRange_tau(&Trace[fc[1126].vertexes[1]][j], _b, step, _d,
////					j);
////			intarray[2] = currA[fc[1126].vertexes[2]].getRange_tau(&Trace[fc[1126].vertexes[2]][j], _b, step, _d,
////					j);
////			intarray[3] = currA[fc[1126].vertexes[3]].getRange_tau(&Trace[fc[1126].vertexes[3]][j], _b, step, _d,
////					j);
//
//			intarray[0] = currA[fc[1126].vertexes[0]].getRange();
//			intarray[1] = currA[fc[1126].vertexes[1]].getRange();
//			intarray[2] = currA[fc[1126].vertexes[2]].getRange();
//			intarray[3] = currA[fc[1126].vertexes[3]].getRange();
//			fs->save_Voxel(currA, fc, intarray, eg, *_d, "face85", currSize, 4,
//					currTau, 1, 10);
//
//			int fi = 13405;
//			intarray[0] = currA[fc[fi].vertexes[0]].getRange();
//			intarray[1] = currA[fc[fi].vertexes[1]].getRange();
//			intarray[2] = currA[fc[fi].vertexes[2]].getRange();
//			intarray[3] = currA[fc[fi].vertexes[3]].getRange();
//			fs->save_Voxel(currA, fc, intarray, eg, *_d, "face85", currSize, 4,
//					currTau, 1, 10);
//
//			fi = 17300;
//			intarray[0] = currA[fc[fi].vertexes[0]].getRange();
//			intarray[1] = currA[fc[fi].vertexes[1]].getRange();
//			intarray[2] = currA[fc[fi].vertexes[2]].getRange();
//			intarray[3] = currA[fc[fi].vertexes[3]].getRange();
//			fs->save_Voxel(currA, fc, intarray, eg, *_d, "face85", currSize, 4,
//					currTau, 1, 10);
//
//			fi = 30565;
//			intarray[0] = currA[fc[fi].vertexes[0]].getRange();
//			intarray[1] = currA[fc[fi].vertexes[1]].getRange();
//			intarray[2] = currA[fc[fi].vertexes[2]].getRange();
//			intarray[3] = currA[fc[fi].vertexes[3]].getRange();
//			fs->save_Voxel(currA, fc, intarray, eg, *_d, "face85", currSize, 4,
//					currTau, 1, 10);
//
//			fi = 59177;
//			intarray[0] = currA[fc[fi].vertexes[0]].getRange();
//			intarray[1] = currA[fc[fi].vertexes[1]].getRange();
//			intarray[2] = currA[fc[fi].vertexes[2]].getRange();
//			intarray[3] = currA[fc[fi].vertexes[3]].getRange();
//			fs->save_Voxel(currA, fc, intarray, eg, *_d, "face85", currSize, 4,
//					currTau, 1, 10);
//
//		} else if (whichData == Benard)
//			fs->save_Quad_One_Face(currA, fc[3309], eg, Trace, *_d, "Benard",
//					currSize, 1, j + 1, 1, currTau * 2);
//		else if (whichData == Tornado)
//			fs->save_Quad_One_Face(currA, fc[1212], eg, Trace, *_d, "Tornado",
//					currSize, 1, j + 1, 1, currTau * 2);
//
//	}
//
//	int* faceVertexMapping = new int[currSize];
//	memset(faceVertexMapping, 0, currSize * sizeof(int));
//	for (int j = 0; j < curFaceSize; j++) {
//		if (!fc[j].bsplit) {
//			faceVertexMapping[fc[j].vertexes[0]] = 1;
//			faceVertexMapping[fc[j].vertexes[1]] = 1;
//			faceVertexMapping[fc[j].vertexes[2]] = 1;
//			faceVertexMapping[fc[j].vertexes[3]] = 1;
//
//		}
//		else
//			fc[j].unsetInBoundary();
//	}
//
//	for (int j = 0; j < currSize; j++) {
//		if (faceVertexMapping[j] == 0)
//			currA[j].unsetInBoundary();
//	}
//	for (int j = 0; j < curEdgeSize; j++) {
//		if (faceVertexMapping[eg[j].v1] == 0 || faceVertexMapping[eg[j].v2] == 0 )
//			eg[j].unsetInBoundary();
//	}

	memset(Fe_Face, 0, curFaceSize * sizeof(uint32_t));
	memset(Fr_Face, 0, curFaceSize * sizeof(uint32_t));
	CheckFace_vertices(currA, eg, fc, Fe_Edge, Fe_Face, _b, _d, step, bForward,
			curFaceSize, currTau);

	Fr_Face[0] = 0;
	for (int i = 0; i <= curFaceSize; i++) {
		Fr_Face[i + 1] = Fr_Face[i] + Fe_Face[i];

	}

	printf("current size = %d \n", currSize);

	Edge* Fc = new Edge[curEdgeSize];

	memset(Fe_Edge, 0, currSize * sizeof(uint32_t));

	for (int i = 0; i < currSize; i++) {
		ASF_vertex vertex = currA[i];
		if (vertex.isInBoundary() && vertex.checkInBoundary(_b)) {

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

//	delete[] Fr_Edge2;
//	delete[] Fr_Face;
//	delete[] Fe_Edge;
//	delete[] Fe_Edge2;
//	delete[] Fe_Face;

//	for (int j = 0; j < CSize * sampleSeeds; j++) {
//		//if(_a[j].isInBoundary())
//		delete[] Trace[j];			// = new Point[_tau + 1];
//	}

//	drawImage(currA, 1300);

//<----------Trimming---------------------------------------

//pair <uint32_t, float> result;
//return result;
}

void ADPAlgorithm2::Flow_Combinatorialization_C(ASF_vertex* _a, fFace *fc,
		fEdge*eg, Boundary* _b, Dimension* _d, Point* step, int _tau,
		float* m_x1, float* m_y1, float* m_z1, ASF_vertex** oVertex,
		Edge ** oFc, uint32_t ** oFr, fFace** oFace, fEdge** oEdge,
		uint32_t * oRSize, uint32_t * oFaceSize, uint32_t * oEdgeSize,
		uint32_t whichData, int CSize, uint32_t _curEdgeSize,
		uint32_t _curFaceSize, int curXDimension, int sampleSeeds,
		int MorseLevel, bool bForward) {

//-----------GPU initialization---------------------------->

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

	for (int ii = 0; ii < icycle; ii++) {

//		if (ii == 80 && !bForward)
//			error(0);

		int ct = 0;
		Tracing_c(currA, m_x1, m_y1, m_z1, _d, _b, step, bForward, whichData,
				curXDimension, 0, currSize, level, ii);

//		Tracing_c_estimated(currA, m_x1, m_y1, m_z1, _d, _b, step, bForward,
//				whichData, curXDimension, CSize, currSize, level);

		CheckNeighborhood(currA, eg, Fe_Edge, _d, _b, step, bForward, currSize,
				curEdgeSize, level, ii);

		Fr_Edge[0] = 0;
		for (int i = 0; i <= curEdgeSize; i++) {

			Fr_Edge[i + 1] = Fr_Edge[i] + Fe_Edge[i];
		}

		//======================================================================

		//while (Fr_Face[curFaceSize] > 0) {
		while (Fr_Edge[curEdgeSize] > 0) {
			//ct++;
			printf("Edges = %d \n", Fr_Edge[curEdgeSize]);
			int oldEdgeSize = curEdgeSize;

			//=========================================================================

			//if (Fr_Edge[curEdgeSize] > 0)
			SplitEndEdge(currA, fc, eg, Fe_Edge, Fr_Edge, _d, _b, step, m_x1,
					m_y1, m_z1, curXDimension, bForward, currSize, curEdgeSize,
					level, ii);

			currSize = currSize + Fr_Edge[curEdgeSize];
			curEdgeSize = curEdgeSize + Fr_Edge[curEdgeSize] * 2;

			CheckNeighborhood(currA, eg, Fe_Edge, _d, _b, step, bForward,
					currSize, curEdgeSize, level, ii);

			Fr_Edge[0] = 0;
			for (int i = 0; i <= curEdgeSize; i++) {

				Fr_Edge[i + 1] = Fr_Edge[i] + Fe_Edge[i];

			}
			int ct = 0;

		}

		printf(" \n %d \n", ii);
	}
	int ct = 0;

	int oldcurrSize = currSize;
//	if (false)
	for (int ii = icycle - 1; ii < icycle; ii++) {
		Tracing_c(currA, m_x1, m_y1, m_z1, _d, _b, step, bForward, whichData,
				curXDimension, oldcurrSize, currSize, level, ii);

		memset(Fe_Face, 0, curFaceSize * sizeof(uint32_t));
		memset(Fr_Face, 0, curFaceSize * sizeof(uint32_t));
		CheckFace_vertices(currA, eg, fc, Fe_Edge, Fe_Face, _b, _d, step,
				bForward, curFaceSize, ii);

		Fr_Face[0] = 0;
		for (int i = 0; i <= curFaceSize; i++) {
			Fr_Face[i + 1] = Fr_Face[i] + Fe_Face[i];
		}

		printf(" Faces = %d \n", Fr_Face[curFaceSize]);

		while (Fr_Face[curFaceSize] > 0) {
			//ct++;
			printf("Faces = %d \n", Fr_Face[curFaceSize]);
			//=======================================================================
			memset(Fe_Edge2, 0, curEdgeSize * sizeof(uint32_t));
			memset(Fr_Edge2, 0, curEdgeSize * sizeof(uint32_t));

			CheckRemainingEdge(currA, fc, eg, Fe_Edge, Fe_Face, Fe_Edge2, _d,
					_b, step, bForward, curFaceSize);

			Fr_Edge2[0] = 0;
			for (int i = 0; i <= curEdgeSize; i++) {
				Fr_Edge2[i + 1] = Fr_Edge2[i] + Fe_Edge2[i];
			}

			printf(" Edges = %d \n", Fr_Edge2[curEdgeSize]);
			//
			//
			SplitRemainingEdge(currA, fc, eg, Fe_Edge2, Fr_Edge2, Fe_Face, _d,
					_b, step, m_x1, m_y1, m_z1, curXDimension, ii, bForward,
					currSize, curEdgeSize, curEdgeSize);
			//
			//
			currSize = currSize + Fr_Edge2[curEdgeSize];
			curEdgeSize = curEdgeSize + Fr_Edge2[curEdgeSize] * 2;

			InsertFaceCenterWithMultipleSeeds(currA, fc, eg, Fr_Face, Fe_Face,
					Fe_Edge, _d, _b, step, m_x1, m_y1, m_z1, ii, curXDimension,
					bForward, currSize, curEdgeSize, curFaceSize);

			//			InsertFaceCenter(currA, fc, eg, Fr_Face, Fe_Face, Fe_Edge, _d, _b,
			//					step, m_x1, m_y1, m_z1, ii, curXDimension, bForward,
			//					currSize, curEdgeSize, curFaceSize);

			currSize = currSize + Fr_Face[curFaceSize];
			//
			int oldfacesize = curFaceSize;

			if (MorseLevel == 1)
				SplitFace(currA, fc, eg, Fr_Face, Fe_Face, Fe_Edge, _d, _b,
						step, bForward, currSize, curEdgeSize, curFaceSize, ii,
						ct, whichData, m_x1, m_y1, m_z1, ii);
			else
				SplitFace2(currA, fc, eg, Fr_Face, Fe_Face, Fe_Edge, _d, _b,
						step, bForward, currSize, curEdgeSize, curFaceSize, ii,
						ct, whichData, m_x1, m_y1, m_z1, ii);

			//printf("\n %d %d \n", eg[105348].E2F[0], eg[101385].E2F[0]);
			curEdgeSize = curEdgeSize + Fr_Face[curFaceSize] * 4;

			curFaceSize = curFaceSize + Fr_Face[curFaceSize] * 4;

			memset(Fe_Face, 0, curFaceSize * sizeof(uint32_t));
			memset(Fr_Face, 0, curFaceSize * sizeof(uint32_t));
			CheckFace_vertices(currA, eg, fc, Fe_Edge, Fe_Face, _b, _d, step,
					bForward, curFaceSize, ii);

			Fr_Face[0] = 0;
			for (int i = 0; i <= curFaceSize; i++) {
				Fr_Face[i + 1] = Fr_Face[i] + Fe_Face[i];

			}

			printf(" Faces = %d \n", Fr_Face[curFaceSize]);

			//============================================================================

			memset(Fr_Edge, 0, (curEdgeSize) * sizeof(uint32_t));
			memset(Fe_Edge, 0, (curEdgeSize) * sizeof(uint32_t));

			CheckNeighborhood(currA, eg, Fe_Edge, _d, _b, step, bForward,
					currSize, curEdgeSize, level, ii);

			Fr_Edge[0] = 0;
			for (int i = 0; i <= curEdgeSize; i++) {

				Fr_Edge[i + 1] = Fr_Edge[i] + Fe_Edge[i];

			}

			//			if (Fr_Face[curFaceSize] > 0 && Fr_Edge[curEdgeSize] == 0)
			//				error(222);

			printf("Edges = %d \n", Fr_Edge[curEdgeSize]);
			ct++;
			//			if (Fr_Face[curFaceSize] > 0 && ct >= 3) {
			//				ct = 0;
			//				break;
			//			}

		}

		/*	else if (currSize != oldsize)
		 currA = _aa;*/
		//	Initialize2 << <grid, threads >> >(d_m, d_d, d_b, d_step, bForward, currSize);*/
		printf(" \n %d \n", ii);

	}	// while (icycle > 1 || Fr_Edge[currSize - 1] > 0);// && terminate);

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

//	delete[] Fr_Edge2;
//	delete[] Fr_Face;
//	delete[] Fe_Edge;
//	delete[] Fe_Edge2;
//	delete[] Fe_Face;

//	for (int j = 0; j < CSize * sampleSeeds; j++) {
//		//if(_a[j].isInBoundary())
//		delete[] Trace[j];			// = new Point[_tau + 1];
//	}

//	drawImage(currA, 1300);

//<----------Trimming---------------------------------------

//pair <uint32_t, float> result;
//return result;
}

void ADPAlgorithm2::Tracing_c_estimated(ASF_vertex * m, float* m_x1,
		float* m_y1, float* m_z1, Dimension* d, Boundary* b, Point* step,
		bool bForward, uint32_t whichData, int currentXDim,
		uint32_t start_num_rows, uint32_t num_rows, uint32_t level) {
//uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
	for (uint32_t row = start_num_rows; row < num_rows; row++) {

		ASF_vertex vertex = m[row];

		if (!vertex.isInBoundary()) {// || !vertex.checkInBoundary_StartPoint(b)) {
			vertex.unsetInBoundary();
			continue;
		}

//		vertex.e.x = ep[0];
//		vertex.e.y = ep[1];
//		vertex.e.z = ep[2];

		bool xy = false;
		bool yz = false;
		bool xz = false;
		//	if (bForward)
		{
			if (vertex.checkInBoundary(b)) {
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
	}

	return;
}

void ADPAlgorithm2::Tracing_c(ASF_vertex * m, float* m_x1, float* m_y1,
		float* m_z1, Dimension* d, Boundary* b, Point* step, bool bForward,
		uint32_t whichData, int currentXDim, uint32_t start_rows,
		uint32_t num_rows, uint32_t level, int currtau) {
	//uint32_t row = blockIdx.x * blockDim.x + threadIdx.x;
	for (uint32_t row = start_rows; row < num_rows; row++) {

		ASF_vertex vertex = m[row];

		if (row == 8701)
			vertex = m[row];
		if (!vertex.checkInBoundary(b)
				|| !vertex.checkInBoundary_StartPoint(b)) {
			vertex.unsetInBoundary();
			m[row] = vertex;
			continue;
		}

		float p1[3];
		if (vertex.type == 1 || vertex.type == 2) {

			p1[0] = vertex.e.x;
			p1[1] = vertex.e.y;
			p1[2] = vertex.e.z;

			float ep[3];
			generalstreamlineTracing_single(p1, bForward, ep, m_x1, m_y1, m_z1,
					b, d, step, currentXDim, 1);

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

		}
//		else if (vertex.type == 2) {
//			if (!m[vertex.left].checkInBoundary(b)
//					|| !m[vertex.right].checkInBoundary(b)
//					|| !m[vertex.left].checkInBoundary_StartPoint(b)
//					|| !m[vertex.right].checkInBoundary_StartPoint(b)) {
//				vertex.unsetInBoundary();
//				m[row] = vertex;
//				continue;
//			}
//
//			p1[0] = vertex.p.x;
//			p1[1] = vertex.p.y;
//			p1[2] = vertex.p.z;
//			float ep[3];
//			generalstreamlineTracing_single(p1, bForward, ep, m_x1, m_y1, m_z1,
//					b, d, step, currentXDim, currtau + 1);
//
//			Point e1 = m[vertex.left].e;
//			Point e2 = m[vertex.right].e;
//			vertex.e.x = (e1.x + e2.x) / 2;
//			vertex.e.y = (e1.y + e2.y) / 2;
//			vertex.e.z = (e1.z + e2.z) / 2;
//			//Trace[row][currtau] = vertex.e;
//			float p1[3];
//
//			//vertex.es[currtau] = vertex.e;
//			Trace[row][currtau] = vertex.e;
//		}
		else if (vertex.type == 3) {

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
				continue;
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
			Trace[row][currtau] = vertex.e;
		}

		bool xy = false;
		bool yz = false;
		bool xz = false;
		//	if (bForward)
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
	}

}

void ADPAlgorithm2::trilinearInterpolation(float p1[3], int idex, Boundary*b,
		Dimension* d, Point* step, float* m_x1, float* m_y1, float* m_z1,
		int primaryXDIMENSION, float& vx, float&vy, float& vz) {
	int xDim;
	int yDim;
	int zDim;

	float vvx, vvy, vvz;
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

	double xd = (p1[0] - p0[0]) / step_x;//(p1[0]-step_x*int(p1[0]/step_x))/step_x;
	double yd = (p1[1] - p0[1]) / step_y;//(p1[1]-step_y*int(p1[1]/step_y))/step_y;
	double zd = (p1[2] - p0[2]) / step_z;//(p1[2]-step_z*int(p1[2]/step_z))/step_z;
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

void ADPAlgorithm2::getLorenzField1(float p[3], float& vx, float& vy,
		float& vz) {

	float sigma = 10.0;
	float ro = 28.0;
	float beta = 8.0 / 3.0;
	vx = sigma * (p[1] - p[0]);
	vy = (p[0] * (ro - p[2])) - p[1];
	vz = p[0] * p[1] - beta * p[2];

}
void ADPAlgorithm2::generalstreamlineTracing_single(float p[3], bool bForward,
		float e[3], float* m_x1, float* m_y1, float* m_z1, Boundary*b,
		Dimension* d, Point* step, int currentDimX, int tau) {

	int start_pixel_id = 0;
	int end_pixel_id = 0;

	int counter_p = 0;

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
}
