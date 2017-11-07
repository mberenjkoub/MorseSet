/*
 * CPUDriver.cpp
 *
 *  Created on: Oct 22, 2016
 *      Author: marzieh
 */

#include "CPUDriver.h"

CPU_Driver::CPU_Driver() {
	// TODO Auto-generated constructor stub
	_d = new Dimension();
	_b = new Boundary();
	_step = new Point();

}

CPU_Driver::~CPU_Driver() {
	// TODO Auto-generated destructor stub
}

void CPU_Driver::Init(int currentStepNumber, int currentDataset) {
	currenttau = currentStepNumber;
	currentData = currentDataset;
	if (currentDataset == Lorenz) {
		InitializeLorenzFieldParameter();

	} else if (currentDataset == Benard) {
		InitializeBenardFieldParameter();

	} else if (currentDataset == Cylinder) {
		InitializeCylinderFieldParameter();
	} else if (currentDataset == Tornado) {
		InitializeTornadoFieldParameter();
	} else if (currentDataset == Ocean) {
		InitializeOceanFieldParameter(); //FieldParameter();
	} else if (currentDataset == Hurricane) {
		InitializeHurricaneFieldParameter(); //FieldParameter();
	}

//    _b->high.x = highBoundary_x;_b->high.y = highBoundary_y;_b->high.z = highBoundary_z;
//    _b->low.x = lowBoundary_x;_b->low.y = lowBoundary_y;_b->low.z = lowBoundary_z;
//    step->x = (highBoundary_x - lowBoundary_x) / (XDIMENSION - 1);
//    step->y = (highBoundary_y - lowBoundary_y) / (YDIMENSION - 1);
//    step->z = (highBoundary_z - lowBoundary_z) / (ZDIMENSION - 1);
//
//    InitArray = InitialArray;
//    _d->x = XDIMENSION;
//    _d->y = YDIMENSION;
//    _d->z = ZDIMENSION;
//
//
//    currentXDimension = currentxdim;
//    Size = XDIMENSION*YDIMENSION*ZDIMENSION;
//    m_x1 = vx;
//    m_y1 = vy;
//    m_z1 = vz;

}

void CPU_Driver::InitializeBenardFieldParameter() {
	_d->x = 128 / 4;
	_d->y = 32 / 4;
	_d->z = 64 / 4;

	_b->low.x = -16.0;
	_b->high.x = 16.0;
	_b->low.y = -4.0;
	_b->high.y = 4.0;
	_b->low.z = -8.0;
	_b->high.z = 8.0;

	_step->x = (_b->high.x - _b->low.x) / (_d->x - 1);
	_step->y = (_b->high.y - _b->low.y) / (_d->y - 1);
	_step->z = (_b->high.z - _b->low.z) / (_d->z - 1);

	currentXDimension = 128;

	scc = new SCC(_d->x, _d->y, _d->z, currentXDimension, _b->high.x, _b->low.x,
			_b->high.y, _b->low.y, _b->high.z, _b->low.z, &m_x1, &m_y1, &m_z1,
			currentData);

}

void CPU_Driver::InitializeOceanFieldParameter() {

	int dcoef = 10;
	_d->x = 1440 / dcoef;
	_d->y = 720 / dcoef;
	_d->z = 50 / dcoef;
	_b->low.x = -36.0;
	_b->high.x = 36.0;
	_b->low.y = -18.0;
	_b->high.y = 18.0;
	_b->low.z = -1.0;
	_b->high.z = 1.5;

	_step->x = (_b->high.x - _b->low.x) / (_d->x - 1);
	_step->y = (_b->high.y - _b->low.y) / (_d->y - 1);
	_step->z = (_b->high.z - _b->low.z) / (_d->z - 1);
	currentXDimension = _d->x * dcoef;

	scc = new SCC(_d->x, _d->y, _d->z, currentXDimension, _b->high.x, _b->low.x,
			_b->high.y, _b->low.y, _b->high.z, _b->low.z, &m_x1, &m_y1, &m_z1,
			currentData);

}

void CPU_Driver::InitializeHurricaneFieldParameter() {

	int dcoef = 10;
	_d->x = 500 / dcoef;
	_d->y = 500 / dcoef;
	_d->z = 100 / dcoef;
	_b->low.x = -30.0;
	_b->high.x = 30.0;
	_b->low.y = -30.0;
	_b->high.y = 30.0;
	_b->low.z = -6.0;
	_b->high.z = 6.0;

	_step->x = (_b->high.x - _b->low.x) / (_d->x - 1);
	_step->y = (_b->high.y - _b->low.y) / (_d->y - 1);
	_step->z = (_b->high.z - _b->low.z) / (_d->z - 1);
	currentXDimension = _d->x * dcoef;

	scc = new SCC(_d->x, _d->y, _d->z, currentXDimension, _b->high.x, _b->low.x,
			_b->high.y, _b->low.y, _b->high.z, _b->low.z, &m_x1, &m_y1, &m_z1,
			currentData);

}

void CPU_Driver::InitializeLorenzFieldParameter() {
	int coef = 2;
	_d->x = 32 / coef;
	_d->y = 32 / coef;
	_d->z = 32 / coef;
	_b->low.x = -30.02;
	_b->high.x = 30.;
	_b->low.y = -30.04;
	_b->high.y = 30.0;
	_b->low.z = -10.01;
	_b->high.z = 50.04;

	_step->x = (_b->high.x - _b->low.x) / (_d->x - 1);
	_step->y = (_b->high.y - _b->low.y) / (_d->y - 1);
	_step->z = (_b->high.z - _b->low.z) / (_d->z - 1);
	currentXDimension = _d->x * coef;

	scc = new SCC(_d->x, _d->y, _d->z, currentXDimension, _b->high.x, _b->low.x,
			_b->high.y, _b->low.y, _b->high.z, _b->low.z, &m_x1, &m_y1, &m_z1,
			currentData);

}

void CPU_Driver::InitializeCylinderFieldParameter() {

	int coef = 4;
	_d->x = 192 / coef;
	_d->y = 64 / coef;
	_d->z = 48 / coef;
	_b->low.x = -12.0;
	_b->high.x = 20.0;
	_b->low.y = -4.0;
	_b->high.y = 4.0;
	_b->low.z = 0.0;
	_b->high.z = 6.0;
	currentXDimension = _d->x * coef;

	_step->x = (_b->high.x - _b->low.x) / (_d->x - 1);
	_step->y = (_b->high.y - _b->low.y) / (_d->y - 1);
	_step->z = (_b->high.z - _b->low.z) / (_d->z - 1);

	scc = new SCC(_d->x, _d->y, _d->z, currentXDimension, _b->high.x, _b->low.x,
			_b->high.y, _b->low.y, _b->high.z, _b->low.z, &m_x1, &m_y1, &m_z1,
			currentData);

}

void CPU_Driver::InitializeTornadoFieldParameter() {

	int coef = 2;
	_d->x = 32 / coef;
	_d->y = 32 / coef;
	_d->z = 32 / coef;
	_b->low.x = -30.0;
	_b->high.x = 30.0;
	_b->low.y = -30.0;
	_b->high.y = 30.0;
	_b->low.z = -10.0;
	_b->high.z = 50.0;
	currentXDimension = _d->x * 2;

	_step->x = (_b->high.x - _b->low.x) / (_d->x - 1);
	_step->y = (_b->high.y - _b->low.y) / (_d->y - 1);
	_step->z = (_b->high.z - _b->low.z) / (_d->z - 1);

	scc = new SCC(_d->x, _d->y, _d->z, currentXDimension, _b->high.x, _b->low.x,
			_b->high.y, _b->low.y, _b->high.z, _b->low.z, &m_x1, &m_y1, &m_z1,
			currentData);

}

void CPU_Driver::save_VTK_File(bool* _MorseSetData) {
	std::string str1;
	;

	if (currentData == Lorenz) {
		//if (XDIMENSION == primaryXDIMENSION)
		str1 = "Lorenz";

	}

	else if (currentData == Cylinder) {
		//if (XDIMENSION == primaryXDIMENSION)
		str1 = "Cylinder";

	}

	else if (currentData == Benard) {
		//if (XDIMENSION == primaryXDIMENSION)
		str1 = "Benard";

	} else if (currentData == Tornado) {
		//if (XDIMENSION == primaryXDIMENSION)
		str1 = "Tornado";

	} else if (currentData == Ocean) {
		//if (XDIMENSION == primaryXDIMENSION)
		str1 = "Ocean";

	}
	std::ofstream fout;
	string frame_str;
	stringstream ss;
	ss << currenttau * MorseLevel;
	ss >> frame_str;
//	fout.open(("frame" + frame_str + ".vtk").c_str(), ios::out);
	string fileName = str1 + frame_str + ".vtk";
	FILE* f_p1 = fopen((str1 + frame_str + ".vtk").c_str(), "wb");
	printf("%s \n", fileName.c_str());

	float step[3];

	////--- the head of vtk file ---, dataset is DATASET STRUCTURED_POINTS
	fprintf(f_p1,
			"# vtk DataFile Version 2.0\nrotationGradient\nASCII\nDATASET STRUCTURED_POINTS\n");
	fprintf(f_p1, "DIMENSIONS %d %d %d\n", _d->x, _d->y, _d->z); //-- the dimension of the strucure GRID
	fprintf(f_p1, "SPACING %f %f %f\nORIGIN %f %f %f\nPOINT_DATA %d\n",
			_step->x, _step->y, _step->z, _b->low.x, _b->low.y, _b->low.z,
			_d->x * _d->y * _d->z);
	fprintf(f_p1, "SCALARS volume_scalars float 1\nLOOKUP_TABLE default\n");

	float value = 0;
	// Define your scalars DATA[x][y][z] and dimensions NX, NY and NZ here!

	//if (XDIMENSION == primaryXDIMENSION)
	{
		//ofstream rawfile("MYRAWFILE_firstlevel.raw", ios::out | ios::binary);

		for (int row = 0; row < _d->x * _d->y * _d->z; row++) {
			if (_MorseSetData[row])
				value = 1.0;
			else
				value = 0.0;
			fprintf(f_p1, "%f\n", value);
			//rawfile.write((char*)&value, sizeof (value));
		}
		/*	}
		 }*/

		//rawfile.close();
	}

	fclose(f_p1);

}

void CPU_Driver::keep_Morse_Data2(fFace* fc, fEdge* eg, ASF_vertex* m,
		fFace** oFc, fEdge** oeg, ASF_vertex** om, bool* MorseSet,
		uint32_t numEdges, uint32_t numFaces, uint32_t num_vertexes,
		int samplesize, int* oFaceSize, int* oVertexSize, int* oEdgeSize) {

//	ASF_vertex* om_ = new ASF_vertex[num_vertexes];
//	fEdge* oeg_ = new fEdge[numEdges * samplesize];
//	fFace* oFc_ = new fFace[numFaces * samplesize];
//
//	memcpy(om_, m, num_vertexes * sizeof(ASF_vertex));
//	memcpy(oeg_, eg, numEdges * sizeof(fEdge));
//	memcpy(oFc_, fc, numFaces * sizeof(fFace));
//
//	*om = om_;
//	*oeg = oeg_;
//	*oFc = oFc_;
//	*oFaceSize = numFaces;
//	*oEdgeSize = numEdges;
//	*oVertexSize = num_vertexes;
//	return;

	int nface = 0;
	int nedge = 0;
	int nvertex = 0;
	bool error = false;

	nreference = new int[num_vertexes];
	nedgereference = new int[numEdges];
	nfacereference = new int[numFaces];

	memset(nreference, 0, num_vertexes * sizeof(int));
	memset(nedgereference, 0, numEdges * sizeof(int));
	memset(nfacereference, 0, numFaces * sizeof(int));

	ASF_vertex* m_ = new ASF_vertex[num_vertexes];
	fEdge* eg_ = new fEdge[numEdges];
	fFace* Fc_ = new fFace[numFaces];

	for (int row = 0; row < num_vertexes; row++) {
		ASF_vertex vertex = m[row];

		if (!MorseSet[m[row].getOldRange()]) {
			vertex.unsetInBoundary();
			m[row] = vertex;
		}
//		if (vertex.isInBoundary()) {
//			//m_.push_back(vertex);
//			if (vertex.type == 3) {
//				if (!m[vertex.left].isInBoundary()
//						|| !m[vertex.right].isInBoundary()
//						|| !m[vertex.down].isInBoundary()
//						|| !m[vertex.up].isInBoundary()) {
//					vertex.unsetInBoundary();
//					m[row] = vertex;
//					continue;
//				}
//			} else if (vertex.type == 2) {
//				if (!m[vertex.left].isInBoundary()
//						|| !m[vertex.right].isInBoundary()) {
//					vertex.unsetInBoundary();
//					m[row] = vertex;
//					continue;
//				}
//			}
//			m_[nvertex] = vertex;
//			nreference[row] = nvertex;
//			nvertex++;
//
//		}
	}

	for (int row = 0; row < numEdges; row++) {
		fEdge edge = eg[row];

		if ((!m[edge.v1].isInBoundary() || !m[edge.v2].isInBoundary())) {

			edge.unsetInBoundary();
			eg[row] = edge;
		}
	}

	for (int row = 0; row < numFaces; row++) {
		fFace face = fc[row];

		if (!m[face.vertexes[0]].isInBoundary()
				|| !m[face.vertexes[1]].isInBoundary()
				|| !m[face.vertexes[2]].isInBoundary()
				|| !m[face.vertexes[3]].isInBoundary()) {

			face.unsetInBoundary();
			fc[row] = face;
		}
	}

	ASF_vertex* om_ = new ASF_vertex[num_vertexes];
	fEdge* oeg_ = new fEdge[numEdges * samplesize];
	fFace* oFc_ = new fFace[numFaces * samplesize];

	memcpy(om_, m, num_vertexes * sizeof(ASF_vertex));
	memcpy(oeg_, eg, numEdges * sizeof(fEdge));
	memcpy(oFc_, fc, numFaces * sizeof(fFace));

	*om = om_;
	*oeg = oeg_;
	*oFc = oFc_;
	*oFaceSize = numFaces;
	*oEdgeSize = numEdges;
	*oVertexSize = num_vertexes;
//	*oFc = fc_;
//	*oeg = &eg_[0];
//	*om = &m_[0];

}

void CPU_Driver::keep_Morse_Data(fFace* fc, fEdge* eg, ASF_vertex* m,
		fFace** oFc, fEdge** oeg, ASF_vertex** om, bool* MorseSet,
		uint32_t numEdges, uint32_t numFaces, uint32_t num_vertexes,
		int samplesize, int* oFaceSize, int* oVertexSize, int* oEdgeSize) {

	int nface = 0;
	int nedge = 0;
	int nvertex = 0;
	bool error = false;

	nreference = new int[num_vertexes];
	nedgereference = new int[numEdges];
	nfacereference = new int[numFaces];

	memset(nreference, 0, num_vertexes * sizeof(int));
	memset(nedgereference, 0, numEdges * sizeof(int));
	memset(nfacereference, 0, numFaces * sizeof(int));

	ASF_vertex* m_ = new ASF_vertex[num_vertexes];
	fEdge* eg_ = new fEdge[numEdges];
	fFace* Fc_ = new fFace[numFaces];

	for (int row = 0; row < num_vertexes; row++) {
		ASF_vertex vertex = m[row];

		if (!MorseSet[m[row].getOldRange()]) {
			vertex.unsetInBoundary();
			m[row] = vertex;
		}
		if (vertex.isInBoundary()) {
			//m_.push_back(vertex);
			if (vertex.type == 3) {
				if (!m[vertex.left].isInBoundary()
						|| !m[vertex.right].isInBoundary()
						|| !m[vertex.down].isInBoundary()
						|| !m[vertex.up].isInBoundary()) {
					vertex.unsetInBoundary();
					m[row] = vertex;
					continue;
				}
			} else if (vertex.type == 2) {
				if (!m[vertex.left].isInBoundary()
						|| !m[vertex.right].isInBoundary()) {
					vertex.unsetInBoundary();
					m[row] = vertex;
					continue;
				}
			}
			m_[nvertex] = vertex;
			nreference[row] = nvertex;
			nvertex++;

		}
	}

	for (int row = 0; row < numEdges; row++) {
		fEdge edge = eg[row];

		if ((m[edge.v1].isInBoundary() && m[edge.v2].isInBoundary())) {

			edge.v1 = nreference[edge.v1];
			edge.v2 = nreference[edge.v2];

			nedgereference[row] = nedge;

			eg_[nedge] = edge;
			//eg_.push_back(edge);
			if (nedge == 2155)
				error = false;
			nedge++;

		} else {
			edge.unsetInBoundary();
			eg[row] = edge;
		}
	}

	for (int row = 0; row < numFaces; row++) {
		fFace face = fc[row];
		int index = face.getFaceVertex();
		if (m[face.vertexes[0]].isInBoundary()
				&& m[face.vertexes[1]].isInBoundary()
				&& m[face.vertexes[2]].isInBoundary()
				&& m[face.vertexes[3]].isInBoundary()) {
			face.vertexes[0] = nreference[face.vertexes[0]];
			face.vertexes[1] = nreference[face.vertexes[1]];
			face.vertexes[2] = nreference[face.vertexes[2]];
			face.vertexes[3] = nreference[face.vertexes[3]];

			face.edge[0] = nedgereference[face.edge[0]];
			face.edge[1] = nedgereference[face.edge[1]];
			face.edge[2] = nedgereference[face.edge[2]];
			face.edge[3] = nedgereference[face.edge[3]];

			nfacereference[row] = nface;

			Fc_[nface] = face;
			//fc_.push_back(face);
			nface++;
		}

		else {
			face.unsetInBoundary();
			fc[row] = face;
		}
	}

	for (int j = 0; j < numFaces; j++) {
		if (fc[j].bsplit) {
			Fc_[nfacereference[j]].subface[0] =
					nfacereference[fc[j].subface[0]];
			Fc_[nfacereference[j]].subface[1] =
					nfacereference[fc[j].subface[1]];
			Fc_[nfacereference[j]].subface[2] =
					nfacereference[fc[j].subface[2]];
			Fc_[nfacereference[j]].subface[3] =
					nfacereference[fc[j].subface[3]];

		}
	}

	for (int j = 0; j < numEdges; j++) {

		if (eg_[nedgereference[j]].v1 == 8146
				&& eg_[nedgereference[j]].v2 == 8351)
			//if (nedgereference[j] == 167)
			error = true;
		if (nedgereference[j] > 0) {
			if (eg[j].bsplit) {
				if (eg[eg[j].subedge[0]].isInBoundary()
						&& eg[eg[j].subedge[1]].isInBoundary()) {
					eg_[nedgereference[j]].subedge[0] =
							nedgereference[eg[j].subedge[0]];
					eg_[nedgereference[j]].subedge[1] =
							nedgereference[eg[j].subedge[1]];
				} else {
					eg_[nedgereference[j]].unsetInBoundary();
					eg[j].unsetInBoundary();
				}
			}

			if (fc[eg[j].E2F[0]].isInBoundary()
					&& fc[eg[j].E2F[1]].isInBoundary() && eg[j].E2F[2] == -1
					&& eg[j].E2F[3] == -1) {
				eg_[nedgereference[j]].E2F[0] = nfacereference[eg[j].E2F[0]];
				eg_[nedgereference[j]].E2F[1] = nfacereference[eg[j].E2F[1]];
//				if (nfacereference[eg[j].E2F[0]] > 0
//						&& fc[eg[j].E2F[0]].isInBoundary())
//					Fc_[nfacereference[eg[j].E2F[0]]].unsetInBoundary();
//				if (nfacereference[eg[j].E2F[1]] > 0
//						&& fc[eg[j].E2F[1]].isInBoundary())
//					Fc_[nfacereference[eg[j].E2F[1]]].unsetInBoundary();

			} else if (fc[eg[j].E2F[0]].isInBoundary()
					&& fc[eg[j].E2F[1]].isInBoundary()
					&& fc[eg[j].E2F[2]].isInBoundary()
					&& fc[eg[j].E2F[3]].isInBoundary()) {
				eg_[nedgereference[j]].E2F[0] = nfacereference[eg[j].E2F[0]];
				eg_[nedgereference[j]].E2F[1] = nfacereference[eg[j].E2F[1]];
				eg_[nedgereference[j]].E2F[2] = nfacereference[eg[j].E2F[2]];
				eg_[nedgereference[j]].E2F[3] = nfacereference[eg[j].E2F[3]];
			}

			else {
				eg_[nedgereference[j]].unsetInBoundary();
				eg[j].unsetInBoundary();
				if (nfacereference[eg[j].E2F[0]] > 0
						&& fc[eg[j].E2F[0]].isInBoundary())
					Fc_[nfacereference[eg[j].E2F[0]]].unsetInBoundary();
				if (nfacereference[eg[j].E2F[1]] > 0
						&& fc[eg[j].E2F[1]].isInBoundary())
					Fc_[nfacereference[eg[j].E2F[1]]].unsetInBoundary();
				if (nfacereference[eg[j].E2F[2]] > 0
						&& fc[eg[j].E2F[2]].isInBoundary())
					Fc_[nfacereference[eg[j].E2F[2]]].unsetInBoundary();
				if (nfacereference[eg[j].E2F[3]] > 0
						&& fc[eg[j].E2F[3]].isInBoundary())
					Fc_[nfacereference[eg[j].E2F[3]]].unsetInBoundary();

			}

		}

	}

	for (int j = 0; j < num_vertexes; j++) {
		if (nreference[j] > 0 && m[j].isInBoundary()) {
			if (m_[nreference[j]].type == 2) {
				m_[nreference[j]].left = nreference[m[j].left];
				m_[nreference[j]].right = nreference[m[j].right];

			} else if (m_[nreference[j]].type == 3) {
				m_[nreference[j]].left = nreference[m[j].left];
				m_[nreference[j]].right = nreference[m[j].right];
				m_[nreference[j]].up = nreference[m[j].up];
				m_[nreference[j]].down = nreference[m[j].down];

			}
		}
	}

	ASF_vertex* om_ = new ASF_vertex[nvertex];
	fEdge* oeg_ = new fEdge[nedge * samplesize];
	fFace* oFc_ = new fFace[nface * samplesize];

	memcpy(om_, m_, nvertex * sizeof(ASF_vertex));
	memcpy(oeg_, eg_, nedge * sizeof(fEdge));
	memcpy(oFc_, Fc_, nface * sizeof(fFace));

	delete[] m_;
	delete[] eg_;
	delete[] Fc_;
//
//	for (int j = 0; j < nvertex; j++) {
//		om_[j] = m_.at(j);
//	}

//	int nface = fc_.size();

////
//
	*om = om_; //om_;
	*oeg = oeg_; //oeg_;
	*oFc = oFc_; //oFc_;
	*oFaceSize = nface;
	*oEdgeSize = nedge;
	*oVertexSize = nvertex;
//	*oFc = fc_;
//	*oeg = &eg_[0];
//	*om = &m_[0];

}

void CPU_Driver::Run() {
	Edge * Fc; // forward columns
	uint32_t * Fr = NULL; // forward ranges
	uint32_t size = _d->x * _d->y * _d->z;
	ADPAlgorithm2* m_ADP_Algorithm = new ADPAlgorithm2();
	//ADP_Algorithm* m_ADP_Algorithm = new ADPAlgorith();
	bool* _MorseSetData;
	_MorseSetData = new bool[size];

	// Backwards arrays
	Edge * Bc; // backward columns
	uint32_t * Br; // backward ranges
	uint32_t RSize = 0;
	ASF_vertex* Initial_a;
	scc->getParameter(&Initial_a);

	ASF_vertex * o_a[2];
	ASF_vertex* oa_F;
	ASF_vertex* oa_B;
	fEdge* oEdge_F;
	fEdge* oEdge_B;
	fFace* oFace_F;
	fFace* oFace_B;
	uint32_t oFaceSize_F;
	uint32_t oFaceSize_B;
	uint32_t oEdgeSize_F;
	uint32_t oEdgeSize_B;
	int level = 1;

	//glNewList(index_StreamLine_Lorenz, GL_COMPILE);

	MorseLevel = 1;
	fFace*fc;
	fFace*fc_b;
	fEdge* eg;
	fEdge* eg_b;
	Edge* oFc, *oBc;
	uint32_t* oFr, *oBr;
	int k = 0;
	int sampleSeeds = 100;

//	if (currentData == Ocean)
//		sampleSeeds = 1;
	uint32_t oRSize_F, oRSize;
	uint32_t oRSize_B, oBSize;

	StopWatchInterface* COLTime = 0;
	(sdkCreateTimer(&COLTime));
	(sdkStartTimer(&COLTime));

	fEdge* edge_f = new fEdge[size * 3 * sampleSeeds];
	fEdge* edge_b = new fEdge[size * 3 * sampleSeeds];

	fFace* face_f = new fFace[size * 3 * sampleSeeds];
	fFace* face_b = new fFace[size * 3 * sampleSeeds];

	m_ADP_Algorithm->Initialize(Initial_a, &eg, &fc, _d, _b, _step, currentData,
			size, size * 3, size * 3, sampleSeeds, size);

	memcpy(edge_f, eg, size * 3 * sampleSeeds);
	memcpy(edge_b, eg, size * 3 * sampleSeeds);

	memcpy(face_f, fc, size * 3 * sampleSeeds);
	memcpy(face_b, fc, size * 3 * sampleSeeds);

	m_ADP_Algorithm->Flow_Combinatorialization_C2(Initial_a, face_f, edge_f, _b,
			_d, _step, currenttau, m_x1, m_y1, m_z1, &oa_F, &oFc, &oFr,
			&oFace_F, &oEdge_F, &oRSize_F, &oFaceSize_F, &oEdgeSize_F,
			currentData, size, size * 3, size * 3, currentXDimension,
			sampleSeeds, MorseLevel, true);

//	memcpy(edge_f, eg, size * 3 * sampleSeeds);
//	memcpy(edge_b, eg, size * 3 * sampleSeeds);
//
//	memcpy(face_f, fc, size * 3 * sampleSeeds);
//	memcpy(face_b, fc, size * 3 * sampleSeeds);
//	m_ADP_Algorithm->Flow_Combinatorialization_C2(Initial_a, face_f, edge_f, _b,
//			_d, _step, currenttau, m_x1, m_y1, m_z1, &oa_F, &oFc, &oFr,
//			&oFace_F, &oEdge_F, &oRSize_F, &oFaceSize_F, &oEdgeSize_F,
//			currentData, size, size * 3, size * 3, currentXDimension,
//			sampleSeeds, MorseLevel, true);

//	m_ADP_Algorithm->Flow_Combinatorialization_C2(oa_F, face_f, edge_f, _b, _d,
//			_step, currenttau, m_x1, m_y1, m_z1, &oa_F, &oFc, &oFr, &oFace_F,
//			&oEdge_F, &oRSize_F, &oFaceSize_F, &oEdgeSize_F, currentData,
//			oRSize_F, oEdgeSize_F, oFaceSize_F, currentXDimension, sampleSeeds,
//			MorseLevel, true);

//	m_ADP_Algorithm->Flow_Combinatorialization_C2(oa_F, face_f, edge_f, _b, _d,
//			_step, currenttau, m_x1, m_y1, m_z1, &oa_F, &oFc, &oFr, &oFace_F,
//			&oEdge_F, &oRSize_F, &oFaceSize_F, &oEdgeSize_F, currentData,
//			oRSize_F, oEdgeSize_F, oFaceSize_F, currentXDimension, sampleSeeds,
//			MorseLevel, true);

	m_ADP_Algorithm->Flow_Combinatorialization_C2(Initial_a, face_b, edge_b, _b,
			_d, _step, currenttau, m_x1, m_y1, m_z1, &oa_B, &oBc, &oBr,
			&oFace_B, &oEdge_B, &oRSize_B, &oFaceSize_B, &oEdgeSize_B,
			currentData, size, size * 3, size * 3, currentXDimension,
			sampleSeeds, MorseLevel, false);

	scc->setParameter(Initial_a, oFc, oFr, oBc, oBr, size, oRSize_F, oRSize_B);
	Load_AllEdges(oFc, oFr, oBc, oBr, size, oRSize_F, oRSize_B, &Fc, &Fr, &Bc,
			&Br, &RSize);
	scc->Run(RSize, size, Fc, Fr, Bc, Br);
	(sdkStopTimer(&COLTime));
	float f = sdkGetTimerValue(&COLTime);
	printf("Running time = %f ms\n", f);

	scc->getMorseSetResult(&_MorseSetData);
	save_VTK_File(_MorseSetData);

	File_Saver fs;
//	fs.save_Quad_Face_Final(oa_F, face_f, edge_f, _MorseSetData, *_d, _b,
//			"face85", oRSize_F, oFaceSize_F, currenttau, currenttau);

	int MorseSetSize = 0;
	int* intarray = new int[size];
	for (int j = 0; j < size; j++)
		if (_MorseSetData[j]) {
			intarray[MorseSetSize] = j;
			MorseSetSize++;
		}


	fs.save_Voxel(Initial_a, face_f, intarray, edge_f, *_d, "Tornado", size,
			MorseSetSize, currenttau, 1, 10);

//	return;
	fFace* nextface_F;
	fEdge* nextedge_F;
	ASF_vertex* nextvertex_F;

	fFace* nextface_B;
	fEdge* nextedge_B;
	ASF_vertex* nextvertex_B;

	int nextFaceSize_F;
	int nextFaceSize_B;

	int nextvertexSize_F;
	int nextvertexSize_B;

	int nextedgeSize_F;
	int nextedgeSize_B;

	while (MorseLevel < 20) {
		//currenttau = currenttau ;

		//sampleSeeds = currenttau;

		File_Saver fs;
		fs.save_Quad_Face(oa_F, face_f, edge_f, *_d, "Lorenz", oRSize_F,
				oFaceSize_F, currenttau * MorseLevel, 1);

		fs.save_Quad_Face(oa_B, face_b, edge_b, *_d, "Backward", oRSize_B,
					oFaceSize_B, currenttau * MorseLevel, 1);

		keep_Morse_Data(face_f, edge_f, oa_F, &nextface_F, &nextedge_F,
				&nextvertex_F, _MorseSetData, oEdgeSize_F, oFaceSize_F,
				oRSize_F, sampleSeeds, &nextFaceSize_F, &nextvertexSize_F,
				&nextedgeSize_F);
		keep_Morse_Data(face_b, edge_b, oa_B, &nextface_B, &nextedge_B,
				&nextvertex_B, _MorseSetData, oEdgeSize_B, oFaceSize_B,
				oRSize_B, sampleSeeds, &nextFaceSize_B, &nextvertexSize_B,
				&nextedgeSize_B);

		MorseLevel++;

		delete[] face_f;
		delete[] face_b;
		delete[] edge_f;
		delete[] edge_b;

		/*

		 delete[] face_f;
		 delete[] face_b;
		 delete[] edge_f;
		 delete[] edge_b;

		 nextFaceSize_F = oFaceSize_F;
		 nextFaceSize_B = oFaceSize_B;

		 nextedgeSize_F = oEdgeSize_F;
		 nextedgeSize_B = oEdgeSize_B;
		 */
//		edge_f = new fEdge[nextedgeSize_F * sampleSeeds];
//		edge_b = new fEdge[nextedgeSize_B * sampleSeeds];
//		face_f = new fFace[nextFaceSize_F * sampleSeeds];
//		face_b = new fFace[nextFaceSize_B * sampleSeeds];
		/*
		 face_f = oFace_F;
		 face_b = oFace_B;

		 edge_f = oEdge_F;
		 edge_b = oEdge_B;
		 //
		 //		memcpy(face_f, nextface_F, nextFaceSize_F * sizeof(fFace));
		 //		memcpy(face_b, nextface_B, nextFaceSize_B * sizeof(fFace));
		 //
		 //		memcpy(edge_f, nextedge_F, nextedgeSize_F * sizeof(fEdge));
		 //		memcpy(edge_b, nextedge_B, nextedgeSize_B * sizeof(fEdge));

		 //		delete[] oa_F;
		 //		delete[] oa_B;
		 //		delete[] oFc;
		 //		delete[] oFr;
		 //		delete[] oBc;
		 //		delete[] oBr;

		 //		delete[] nextedge_F;
		 //		delete[] nextedge_B;
		 //		delete[] nextface_F;
		 //		delete[] nextface_B;


		 */
//
//		uint32_t* Fe_Edge = new uint32_t[nextedgeSize_F];
//		uint32_t* Fe_Face = new uint32_t[nextFaceSize_F];
//		uint32_t* Fr_Face = new uint32_t[nextFaceSize_F];
//
//		m_ADP_Algorithm->CheckFace_vertices(oa_F, edge_f, face_f, Fe_Edge,
//				Fe_Face, _b, _d, _step, true, oFaceSize_F, currenttau);
//
//		Fr_Face[0] = 0;
//		for (int i = 0; i <= nextFaceSize_F; i++) {
//			Fr_Face[i + 1] = Fr_Face[i] + Fe_Face[i];
//		}
//
//		m_ADP_Algorithm->CheckFace_vertices(nextvertex_F, nextedge_F,
//				nextface_F, Fe_Edge, Fe_Face, _b, _d, _step, true,
//				nextFaceSize_F, currenttau);
//
//		Fr_Face[0] = 0;
//		for (int i = 0; i <= nextFaceSize_F; i++) {
//			Fr_Face[i + 1] = Fr_Face[i] + Fe_Face[i];
//		}
//
//		printf(" Faces = %d \n", Fr_Face[nextFaceSize_F]);
		face_f = nextface_F;
		face_b = nextface_B;
		edge_f = nextedge_F;
		edge_b = nextedge_B;
		oa_F = nextvertex_F;
		oa_B = nextvertex_B;
		//sampleSeeds = currenttau*2;

		memset(_MorseSetData, 0, sizeof(bool) * size);

		{
			m_ADP_Algorithm->Flow_Combinatorialization_C2(oa_F, face_f, edge_f,
					_b, _d, _step, currenttau, m_x1, m_y1, m_z1, &oa_F, &oFc,
					&oFr, &oFace_F, &oEdge_F, &oRSize_F, &oFaceSize_F,
					&oEdgeSize_F, currentData, nextvertexSize_F, nextedgeSize_F,
					nextFaceSize_F, currentXDimension, sampleSeeds, MorseLevel,
					true);

			m_ADP_Algorithm->Flow_Combinatorialization_C2(oa_B, face_b, edge_b,
					_b, _d, _step, currenttau, m_x1, m_y1, m_z1, &oa_B, &oBc,
					&oBr, &oFace_B, &oEdge_B, &oRSize_B, &oFaceSize_B,
					&oEdgeSize_B, currentData, nextvertexSize_B, nextedgeSize_B,
					nextFaceSize_B, currentXDimension, sampleSeeds, MorseLevel,
					false);

			SCC* scc1 = new SCC(_d->x, _d->y, _d->z, currentXDimension,
					_b->high.x, _b->low.x, _b->high.y, _b->low.y, _b->high.z,
					_b->low.z, &m_x1, &m_y1, &m_z1, currentData);

			scc1->setParameter(Initial_a, oFc, oFr, oBc, oBr, size, oRSize_F,
					oRSize_B);
			Load_AllEdges(oFc, oFr, oBc, oBr, size, oRSize_F, oRSize_B, &Fc,
					&Fr, &Bc, &Br, &RSize);
			scc1->Run(RSize, size, Fc, Fr, Bc, Br);
			(sdkStopTimer(&COLTime));
			float f = sdkGetTimerValue(&COLTime);
			printf("Running time = %f ms\n", f);

			scc1->getMorseSetResult(&_MorseSetData);
			save_VTK_File(_MorseSetData);
			nextvertexSize_F = oRSize_F;
			nextvertexSize_B = oRSize_B;
			nextedgeSize_F = oEdgeSize_F;
			nextedgeSize_B = oEdgeSize_B;

			nextvertex_F = oa_F;
			nextvertex_B = oa_B;
//			oa_F = nextvertex_F;
//			oa_B = nextvertex_B;

		}
//		File_Saver fs;
//		fs.save_Quad_Face_Final(oa_F, face_f, edge_f, _MorseSetData, *_d, _b,
//				"face85", oRSize_F, oFaceSize_F, currenttau, currenttau);
	}
}

void CPU_Driver::Run_GPU() {
	Edge * Fc; // forward columns
	uint32_t * Fr = NULL; // forward ranges
	uint32_t size = _d->x * _d->y * _d->z;
	ADPAlgorithm2* m_ADP_Algorithm = new ADPAlgorithm2();
	//ADP_Algorithm* m_ADP_Algorithm = new ADPAlgorith();
	bool* _MorseSetData;
	_MorseSetData = new bool[size];

	// Backwards arrays
	Edge * Bc; // backward columns
	uint32_t * Br; // backward ranges
	uint32_t RSize = 0;
	ASF_vertex* Initial_a;
	scc->getParameter(&Initial_a);

	ASF_vertex * o_a[2];
	ASF_vertex* oa_F;
	ASF_vertex* oa_B;
	fEdge* oEdge_F;
	fEdge* oEdge_B;
	fFace* oFace_F;
	fFace* oFace_B;
	uint32_t oFaceSize_F;
	uint32_t oFaceSize_B;
	uint32_t oEdgeSize_F;
	uint32_t oEdgeSize_B;
	int level = 1;

	//glNewList(index_StreamLine_Lorenz, GL_COMPILE);

	MorseLevel = 1;
	fFace*fc;
	fFace*fc_b;
	fEdge* eg;
	fEdge* eg_b;
	Edge* oFc, *oBc;
	uint32_t* oFr, *oBr;
	int k = 0;
	int sampleSeeds = 10;

	if (currentData == Ocean)
		sampleSeeds = 1;
	uint32_t oRSize_F, oRSize;
	uint32_t oRSize_B, oBSize;

	StopWatchInterface* COLTime = 0;
	(sdkCreateTimer(&COLTime));
	(sdkStartTimer(&COLTime));

	fEdge* edge_f = new fEdge[size * 3 * sampleSeeds];
	fEdge* edge_b = new fEdge[size * 3 * sampleSeeds];

	fFace* face_f = new fFace[size * 3 * sampleSeeds];
	fFace* face_b = new fFace[size * 3 * sampleSeeds];

	m_ADP_Algorithm->Initialize(Initial_a, &eg, &fc, _d, _b, _step, currentData,
			size, size * 3, size * 3, sampleSeeds, size);

	memcpy(edge_f, eg, size * 3 * sampleSeeds);
	memcpy(edge_b, eg, size * 3 * sampleSeeds);

	memcpy(face_f, fc, size * 3 * sampleSeeds);
	memcpy(face_b, fc, size * 3 * sampleSeeds);

	SCC2* scc2 = new SCC2();

	scc2->Run(Initial_a, face_f, edge_f, _b, _d, _step, currenttau, m_x1, m_y1,
			m_z1, &oa_F, &oFc, &oFr, &oFace_F, &oEdge_F, &oRSize_F,
			&oFaceSize_F, &oEdgeSize_F, currentData, size, size * 3, size * 3,
			currentXDimension, sampleSeeds, MorseLevel, true);

//	memcpy(edge_f, eg, size * 3 * sampleSeeds);
//	memcpy(edge_b, eg, size * 3 * sampleSeeds);
//
//	memcpy(face_f, fc, size * 3 * sampleSeeds);
//	memcpy(face_b, fc, size * 3 * sampleSeeds);
//	m_ADP_Algorithm->Flow_Combinatorialization_C2(Initial_a, face_f, edge_f, _b,
//			_d, _step, currenttau, m_x1, m_y1, m_z1, &oa_F, &oFc, &oFr,
//			&oFace_F, &oEdge_F, &oRSize_F, &oFaceSize_F, &oEdgeSize_F,
//			currentData, size, size * 3, size * 3, currentXDimension,
//			sampleSeeds, MorseLevel, true);

//	m_ADP_Algorithm->Flow_Combinatorialization_C2(oa_F, face_f, edge_f, _b, _d,
//			_step, currenttau, m_x1, m_y1, m_z1, &oa_F, &oFc, &oFr, &oFace_F,
//			&oEdge_F, &oRSize_F, &oFaceSize_F, &oEdgeSize_F, currentData,
//			oRSize_F, oEdgeSize_F, oFaceSize_F, currentXDimension, sampleSeeds,
//			MorseLevel, true);

//	m_ADP_Algorithm->Flow_Combinatorialization_C2(oa_F, face_f, edge_f, _b, _d,
//			_step, currenttau, m_x1, m_y1, m_z1, &oa_F, &oFc, &oFr, &oFace_F,
//			&oEdge_F, &oRSize_F, &oFaceSize_F, &oEdgeSize_F, currentData,
//			oRSize_F, oEdgeSize_F, oFaceSize_F, currentXDimension, sampleSeeds,
//			MorseLevel, true);

	m_ADP_Algorithm->Flow_Combinatorialization_C2(Initial_a, face_b, edge_b, _b,
			_d, _step, currenttau, m_x1, m_y1, m_z1, &oa_B, &oBc, &oBr,
			&oFace_B, &oEdge_B, &oRSize_B, &oFaceSize_B, &oEdgeSize_B,
			currentData, size, size * 3, size * 3, currentXDimension,
			sampleSeeds, MorseLevel, false);

	scc->setParameter(Initial_a, oFc, oFr, oBc, oBr, size, oRSize_F, oRSize_B);
	Load_AllEdges(oFc, oFr, oBc, oBr, size, oRSize_F, oRSize_B, &Fc, &Fr, &Bc,
			&Br, &RSize);
	scc->Run(RSize, size, Fc, Fr, Bc, Br);
	(sdkStopTimer(&COLTime));
	float f = sdkGetTimerValue(&COLTime);
	printf("Running time = %f ms\n", f);

	scc->getMorseSetResult(&_MorseSetData);
	save_VTK_File(_MorseSetData);

	File_Saver fs;
//	fs.save_Quad_Face_Final(oa_F, face_f, edge_f, _MorseSetData, *_d, _b,
//			"face85", oRSize_F, oFaceSize_F, currenttau, currenttau);

	int MorseSetSize = 0;
	int* intarray = new int[size];
	for (int j = 0; j < size; j++)
		if (_MorseSetData[j]) {
			intarray[MorseSetSize] = j;
			MorseSetSize++;
		}

	fs.save_Voxel(Initial_a, face_f, intarray, edge_f, *_d, "Tornado", size,
			MorseSetSize, currenttau, 1, 10);

//	return;
	fFace* nextface_F;
	fEdge* nextedge_F;
	ASF_vertex* nextvertex_F;

	fFace* nextface_B;
	fEdge* nextedge_B;
	ASF_vertex* nextvertex_B;

	int nextFaceSize_F;
	int nextFaceSize_B;

	int nextvertexSize_F;
	int nextvertexSize_B;

	int nextedgeSize_F;
	int nextedgeSize_B;

	while (MorseLevel < 20) {
		//currenttau = currenttau ;

		//sampleSeeds = currenttau;

		File_Saver fs;
		fs.save_Quad_Face(oa_F, face_f, edge_f, *_d, "Lorenz", oRSize_F,
				oFaceSize_F, currenttau * MorseLevel, 1);
		keep_Morse_Data2(face_f, edge_f, oa_F, &nextface_F, &nextedge_F,
				&nextvertex_F, _MorseSetData, oEdgeSize_F, oFaceSize_F,
				oRSize_F, sampleSeeds, &nextFaceSize_F, &nextvertexSize_F,
				&nextedgeSize_F);
		keep_Morse_Data2(face_b, edge_b, oa_B, &nextface_B, &nextedge_B,
				&nextvertex_B, _MorseSetData, oEdgeSize_B, oFaceSize_B,
				oRSize_B, sampleSeeds, &nextFaceSize_B, &nextvertexSize_B,
				&nextedgeSize_B);

		MorseLevel++;

		delete[] face_f;
		delete[] face_b;
		delete[] edge_f;
		delete[] edge_b;

		/*

		 delete[] face_f;
		 delete[] face_b;
		 delete[] edge_f;
		 delete[] edge_b;

		 nextFaceSize_F = oFaceSize_F;
		 nextFaceSize_B = oFaceSize_B;

		 nextedgeSize_F = oEdgeSize_F;
		 nextedgeSize_B = oEdgeSize_B;
		 */
//		edge_f = new fEdge[nextedgeSize_F * sampleSeeds];
//		edge_b = new fEdge[nextedgeSize_B * sampleSeeds];
//		face_f = new fFace[nextFaceSize_F * sampleSeeds];
//		face_b = new fFace[nextFaceSize_B * sampleSeeds];
		/*
		 face_f = oFace_F;
		 face_b = oFace_B;

		 edge_f = oEdge_F;
		 edge_b = oEdge_B;
		 //
		 //		memcpy(face_f, nextface_F, nextFaceSize_F * sizeof(fFace));
		 //		memcpy(face_b, nextface_B, nextFaceSize_B * sizeof(fFace));
		 //
		 //		memcpy(edge_f, nextedge_F, nextedgeSize_F * sizeof(fEdge));
		 //		memcpy(edge_b, nextedge_B, nextedgeSize_B * sizeof(fEdge));

		 //		delete[] oa_F;
		 //		delete[] oa_B;
		 //		delete[] oFc;
		 //		delete[] oFr;
		 //		delete[] oBc;
		 //		delete[] oBr;

		 //		delete[] nextedge_F;
		 //		delete[] nextedge_B;
		 //		delete[] nextface_F;
		 //		delete[] nextface_B;


		 */
//
//		uint32_t* Fe_Edge = new uint32_t[nextedgeSize_F];
//		uint32_t* Fe_Face = new uint32_t[nextFaceSize_F];
//		uint32_t* Fr_Face = new uint32_t[nextFaceSize_F];
//
//		m_ADP_Algorithm->CheckFace_vertices(oa_F, edge_f, face_f, Fe_Edge,
//				Fe_Face, _b, _d, _step, true, oFaceSize_F, currenttau);
//
//		Fr_Face[0] = 0;
//		for (int i = 0; i <= nextFaceSize_F; i++) {
//			Fr_Face[i + 1] = Fr_Face[i] + Fe_Face[i];
//		}
//
//		m_ADP_Algorithm->CheckFace_vertices(nextvertex_F, nextedge_F,
//				nextface_F, Fe_Edge, Fe_Face, _b, _d, _step, true,
//				nextFaceSize_F, currenttau);
//
//		Fr_Face[0] = 0;
//		for (int i = 0; i <= nextFaceSize_F; i++) {
//			Fr_Face[i + 1] = Fr_Face[i] + Fe_Face[i];
//		}
//
//		printf(" Faces = %d \n", Fr_Face[nextFaceSize_F]);
		face_f = nextface_F;
		face_b = nextface_B;
		edge_f = nextedge_F;
		edge_b = nextedge_B;
		oa_F = nextvertex_F;
		oa_B = nextvertex_B;
		//sampleSeeds = currenttau*2;

		memset(_MorseSetData, 0, sizeof(bool) * size);

		{
			m_ADP_Algorithm->Flow_Combinatorialization_C2(oa_F, face_f, edge_f,
					_b, _d, _step, currenttau, m_x1, m_y1, m_z1, &oa_F, &oFc,
					&oFr, &oFace_F, &oEdge_F, &oRSize_F, &oFaceSize_F,
					&oEdgeSize_F, currentData, nextvertexSize_F, nextedgeSize_F,
					nextFaceSize_F, currentXDimension, sampleSeeds, MorseLevel,
					true);

			m_ADP_Algorithm->Flow_Combinatorialization_C2(oa_B, face_b, edge_b,
					_b, _d, _step, currenttau, m_x1, m_y1, m_z1, &oa_B, &oBc,
					&oBr, &oFace_B, &oEdge_B, &oRSize_B, &oFaceSize_B,
					&oEdgeSize_B, currentData, nextvertexSize_B, nextedgeSize_B,
					nextFaceSize_B, currentXDimension, sampleSeeds, MorseLevel,
					false);

			SCC* scc1 = new SCC(_d->x, _d->y, _d->z, currentXDimension,
					_b->high.x, _b->low.x, _b->high.y, _b->low.y, _b->high.z,
					_b->low.z, &m_x1, &m_y1, &m_z1, currentData);

			scc1->setParameter(Initial_a, oFc, oFr, oBc, oBr, size, oRSize_F,
					oRSize_B);
			Load_AllEdges(oFc, oFr, oBc, oBr, size, oRSize_F, oRSize_B, &Fc,
					&Fr, &Bc, &Br, &RSize);
			scc1->Run(RSize, size, Fc, Fr, Bc, Br);
			(sdkStopTimer(&COLTime));
			float f = sdkGetTimerValue(&COLTime);
			printf("Running time = %f ms\n", f);

			scc1->getMorseSetResult(&_MorseSetData);
			save_VTK_File(_MorseSetData);
			nextvertexSize_F = oRSize_F;
			nextvertexSize_B = oRSize_B;
			nextedgeSize_F = oEdgeSize_F;
			nextedgeSize_B = oEdgeSize_B;

			nextvertex_F = oa_F;
			nextvertex_B = oa_B;
//			oa_F = nextvertex_F;
//			oa_B = nextvertex_B;

		}
//		File_Saver fs;
//		fs.save_Quad_Face_Final(oa_F, face_f, edge_f, _MorseSetData, *_d, _b,
//				"face85", oRSize_F, oFaceSize_F, currenttau, currenttau);
	}
}

//void CPU_Driver::doubleResolution()
//{
//
//}
