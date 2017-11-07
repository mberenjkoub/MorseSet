/*
 * SCC2.cpp
 *
 *  Created on: Jan 4, 2017
 *      Author: marzieh
 */

#include <SCC2.h>



pair<uint32_t, float> OBF_Decomposition2(ASF_vertex* _a, fFace *fc,
fEdge*eg, Boundary* _b, Dimension* _d, Point* step, int _tau,
float* m_x1, float* m_y1, float* m_z1, ASF_vertex** oVertex,
Edge ** oFc, uint32_t ** oFr, fFace** oFace, fEdge** oEdge,
uint32_t * oRSize, uint32_t * oFaceSize, uint32_t * oEdgeSize,
uint32_t whichData, int CSize, uint32_t _curEdgeSize,
uint32_t _curFaceSize, int curXDimension, int sampleSeeds,
int MorseLevel, bool bForward);

SCC2::SCC2(int xdim, int ydim, int zdim, int pxdim, float high_x, float low_x,
		float high_y, float low_y, float high_z, float low_z, float** _vx,
		float** _vy, float** _vz, int whichData) {

	size.x = xdim;
	size.y = ydim;
	size.z = zdim;
	CSize = size.x * size.y * size.z;

	int coef = (pxdim / xdim);
	vx = new float[xdim * ydim * zdim * coef * coef * coef];
	vy = new float[xdim * ydim * zdim * coef * coef * coef];
	vz = new float[xdim * ydim * zdim * coef * coef * coef];

	int origxDim = xdim*coef;
	int origyDim = ydim*coef;
	int origzDim = zdim*coef;

	b.high.x = high_x;
	b.low.x = low_x;
	b.high.y = high_y;
	b.low.y = low_y;
	b.high.z = high_z;
	b.low.z = low_z;
	a = new ASF_vertex[CSize];
	_m = new OBF_vertex[CSize - 1];
//	_a = new ASF_vertex[(RSize - 1) * 2];
	b_MorseSet = new bool[CSize];
	whichDataType = whichData;
	primaryXDim = pxdim;
	if (whichData == Lorenz) {

		generateData(a,origxDim,origyDim,origzDim, &size, &b, vx, vy, vz, &vData);
	} else if (whichData == Ocean) {
		//getOceanRawData(&vx, &vy, &vz);
		getOceanData(&vx,&vy,&vz);
		//printf("\n %f", vz[0][0]);

		generateDataOcean(a, &size, &b, primaryXDim, vx, vy, vz, &vData);
	} else if (whichData == Benard) {
		getBenardData(&vx, &vy, &vz);
		//printf("\n %f", vz[0][0]);

		generateDataBenard(a, &size, &b, primaryXDim, vx, vy, vz, &vData);
	} else if (whichData == Tornado) {
		getTornadoData(origxDim, origyDim, origzDim, 1000, vx, vy, vz);

		generateDataTornado(a, &size, &b, primaryXDim, vx, vy, vz, &vData);
	}

	else if (whichData == Cylinder) {
		getCylinderData(vx, vy, vz, b);
		//printf("\n %f", vz[0][0]);
		generateDataCylinder(a, &size, &b, primaryXDim, vx, vy, vz, &vData);
	}

	else if (whichData == Hurricane) {
			getHuricaneData(&vx, &vy, &vz);
			//printf("\n %f", vz[0][0]);

			generateDataHurricane(a, &size, &b, primaryXDim, vx, vy, vz, &vData);
			//generateDataBenard(a, &size, &b, primaryXDim, vx, vy, vz, &vData);
		}
	*_vx = vx;
	*_vy = vy;
	*_vz = vz;
//	Run(100);

}

SCC2::SCC2()
{

}

void SCC2::getParameter(ASF_vertex** _a) //,  Boundary*_b, Dimension*_d,uint32_t* _seedsize)
		{
	*_a = a;
	//*_v = vData;
	//*_b = b;
	//*_d = size;
	//*_seedsize = seednumber;

}

void SCC2::loadBenardData(float** vx, float** vy, float** vz) {
	getBenardData(vx, vy, vz);
	//printf("\n %f", vz[0][0]);

	generateDataBenard(a, &size, &b, primaryXDim, *vx, *vy, *vz, &vData);

}

void SCC2::loadTornadoData(float** vx, float** vy, float** vz) {
	//getTornadoData( vx, vy, vz);
	////printf("\n %f", vz[0][0]);

	//generateDataTornado(a, &size, &b, primaryXDim, *vx, *vy, *vz, &vData);

}

void SCC2::setParameter(ASF_vertex* _a, Edge* _Fc_f, uint32_t* _Fr_f,
		Edge* _Bc_b, uint32_t* _Br_b, uint32_t _CSize, uint32_t _RSize_F,
		uint32_t _RSize_B) {
	//m_Asf = _a;
	Fc_f = _Fc_f;
	Fr_f = _Fr_f;
	Bc_b = _Bc_b;
	Br_b = _Br_b;
	CSize = _CSize;
	RSize_F = _RSize_F;
	RSize_B = _RSize_B;

}

void SCC2::SetOutput(uint32_t** oFr, uint32_t **oBr, ASF_vertex**oa) {
	*oa = _oa;
	*oBr = Br_b;
	*oFr = Fr_f;

}

void SCC2::Run(ASF_vertex* _a, fFace *fc,
		fEdge*eg, Boundary* _b, Dimension* _d, Point* step, int _tau,
		float* m_x1, float* m_y1, float* m_z1, ASF_vertex** oVertex,
		Edge ** oFc, uint32_t ** oFr, fFace** oFace, fEdge** oEdge,
		uint32_t * oRSize, uint32_t * oFaceSize, uint32_t * oEdgeSize,
		uint32_t whichData, int CSize, uint32_t _curEdgeSize,
		uint32_t _curFaceSize, int curXDimension, int sampleSeeds,
		int MorseLevel, bool bForward) {



			pair<uint32_t, float> result = OBF_Decomposition2(_a,fc,eg,_b,_d,step,
					_tau,m_x1,m_y1,m_z1,oVertex,oFc,oFr,oFace,oEdge, oRSize,oFaceSize,oEdgeSize,whichData,CSize,
					_curEdgeSize,_curFaceSize,curXDimension,sampleSeeds,MorseLevel,bForward);


}


void SCC2::loadOceanData(float** vx, float** vy, float** vz) {

	//getOceanData(vx, vy, vz);
	//printf("\n %f", vz[0][0]);

	//generateDataOcean(a, &size, &b, primaryXDim, *vx, *vy, *vz,&vData);

}

//void SCC::getProbability(float** _fProbability)
//{
//	float* fProb = new float[CSize];
//	for (int i = 0; i < CSize; i++)
//	{
//		fProb[i] = m_Asf[i].probability;
//
//	}
//	*_fProbability = fProb;
//
//}

void SCC2::getMorseSetResult(bool** _bMorseSetData) {

	for (int i = 0; i < CSize; i++) {
		b_MorseSet[i] = false;

		if (m_Asf[i].isInSCC())
			b_MorseSet[i] = true;
	}
	*_bMorseSetData = b_MorseSet;
}
SCC2::~SCC2() {

	if (!dve) {
		delete[] Fr;
		delete[] Fc;
		delete[] Br;
		delete[] Bc;
	}
}


