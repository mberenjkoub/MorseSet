/*
 * SCC.cpp
 *
 *  Created on: Oct 22, 2016
 *      Author: marzieh
 */

#include <SCC.h>

pair<uint32_t, float> OBF_Decomposition(uint32_t, uint32_t, Edge *, uint32_t *,
		Edge *, uint32_t *, OBF_vertex *, ASF_vertex*a, Dimension*d, bool,
		uint32_t, bool trim = false, bool quiet = false);

SCC::SCC(int xdim, int ydim, int zdim, int pxdim, float high_x, float low_x,
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

void SCC::getParameter(ASF_vertex** _a) //,  Boundary*_b, Dimension*_d,uint32_t* _seedsize)
		{
	*_a = a;
	//*_v = vData;
	//*_b = b;
	//*_d = size;
	//*_seedsize = seednumber;

}

void SCC::loadBenardData(float** vx, float** vy, float** vz) {
	getBenardData(vx, vy, vz);
	//printf("\n %f", vz[0][0]);

	generateDataBenard(a, &size, &b, primaryXDim, *vx, *vy, *vz, &vData);

}

void SCC::loadTornadoData(float** vx, float** vy, float** vz) {
	//getTornadoData( vx, vy, vz);
	////printf("\n %f", vz[0][0]);

	//generateDataTornado(a, &size, &b, primaryXDim, *vx, *vy, *vz, &vData);

}

void SCC::setParameter(ASF_vertex* _a, Edge* _Fc_f, uint32_t* _Fr_f,
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

void SCC::SetOutput(uint32_t** oFr, uint32_t **oBr, ASF_vertex**oa) {
	*oa = _oa;
	*oBr = Br_b;
	*oFr = Fr_f;

}

void SCC::Run(uint32_t RSize, uint32_t CSize, Edge* Fc, uint32_t *Fr, Edge* Bc,
		uint32_t* Br) {

	std::vector<int> v_Fr;
	std::vector<Edge> v_Fc;

	std::vector<int> v_Br;
	std::vector<Edge> v_Bc;
	//tau = 100;
	char c = 'o';	// argv[1][0];
	try {
		switch (c) {

		case 'o': {

			/*	StopWatchInterface* COLTime = 0;
			 (sdkCreateTimer(&COLTime));
			 (sdkStartTimer(&COLTime));*/

			pair<uint32_t, float> result = OBF_Decomposition(RSize, CSize, Fc,
					Fr, Bc, Br, _m, a, &size, /*(argc > 4) ? bool(atoi(argv[3])) :*/
					false,
					/*(argc > 5) ? atoi(argv[4]) :*/1000, /*(argc > 3) ? bool(atoi(argv[2])) :*/
					true, true/*(argv[1][1] != '\0')*/);

			/*	(sdkStopTimer(&COLTime));
			 float f = sdkGetTimerValue(&COLTime);
			 printf("Running time = %f ms\n", f);*/
			m_Asf = new ASF_vertex[CSize - 1];
			m_Asf = a;

			//result1 = Adaptive_Sampling(a, _a, &b, &size, &tau, &Bc_b, &Br_b, &RSize_B, false);

			//if (argv[1][1] != '\0')
			printf(" scc = %f", result.second);
			delete[] _m;
		}
			break;

		}
	} catch (const char * e) {
		printf("%s\n", e);
		return;
	}

}

void SCC::Run(uint32_t tau) {

	std::vector<int> v_Fr;
	std::vector<Edge> v_Fc;

	std::vector<int> v_Br;
	std::vector<Edge> v_Bc;
	//tau = 100;
	char c = 'o';	// argv[1][0];
	try {
		switch (c) {

		case 'o': {

			StopWatchInterface* COLTime = 0;
			(sdkCreateTimer(&COLTime));
			(sdkStartTimer(&COLTime));

			//forward
			/*	do
			 {*/
			/*	pair <uint32_t, float> result1 = Flow_Combinatorialization(a, &_oa, vData,vx,vy,vz, &b, &size, tau, &Fc_f, &Fr_f, &RSize_F, primaryXDim, whichDataType, true);
			 seednumber = RSize_F;
			 */
			/*v_Fc.insert(v_Fc.begin(), Fc_f, Fc_f + RSize_F);
			 v_Fr.insert(v_Fr.begin(), Fr_f, Fr_f + RSize_F);*/
			//	ASF_vertex *_oa2;
			//result1 = Flow_Combinatorialization(a, &_oa2, vData,vx,vy,vz, &b, &size, tau, &Bc_b, &Br_b, &RSize_B, primaryXDim, whichDataType, false);
			/*v_Bc.insert(v_Bc.begin(), Bc_b, Bc_b + RSize_B);
			 v_Br.insert(v_Br.begin(), Br_b, Br_b + RSize_B);*/

			//ASF_vertex* n_a = new ASF_vertex[CSize * 8];
			//result1 = Adaptive_Sampling(a,n_a, &b, &size,  &RSize_F, true);
			/*	break;
			 } while (true);*/

			//backward
			//do
			//{
			//	pair <uint32_t, float> result1 = Flow_Combinatorialization(a,vData, &b, &size, tau, &Bc_b, &Br_b, &RSize_B, whichDataType, false);
			//	v_Bc.insert(v_Bc.begin(), Bc_b, Bc_b + RSize_B);
			//	v_Br.insert(v_Br.begin(), Br_b, Br_b + RSize_B);
			//	//result1 = Flow_Combinatorialization(a, &b, &size, &tau, &Bc_b, &Br_b, &RSize_B, whichDataType, false);
			//	break;
			//} while (true);

			//Load_AllEdges(Fc_f, Fr_f, Bc_b, Br_b, &Fc, &Fr, &Bc, &Br, &CSize, &RSize);

			Load_AllEdges(Fc_f, Fr_f, Bc_b, Br_b, CSize, RSize_F, RSize_B, &Fc,
					&Fr, &Bc, &Br, &RSize);

			pair<uint32_t, float> result = OBF_Decomposition(RSize, CSize, Fc,
					Fr, Bc, Br, _m, a, &size, /*(argc > 4) ? bool(atoi(argv[3])) :*/
					false,
					/*(argc > 5) ? atoi(argv[4]) :*/1000, /*(argc > 3) ? bool(atoi(argv[2])) :*/
					true, true/*(argv[1][1] != '\0')*/);

			(sdkStopTimer(&COLTime));
			float f = sdkGetTimerValue(&COLTime);
			printf("Running time = %f ms\n", f);
			m_Asf = new ASF_vertex[CSize - 1];
			m_Asf = a;

			//result1 = Adaptive_Sampling(a, _a, &b, &size, &tau, &Bc_b, &Br_b, &RSize_B, false);

			//if (argv[1][1] != '\0')
			printf("%f", result.second);
			delete[] _m;
		}
			break;

		}
	} catch (const char * e) {
		printf("%s\n", e);
		return;
	}

}

void SCC::loadOceanData(float** vx, float** vy, float** vz) {

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

void SCC::getMorseSetResult(bool** _bMorseSetData) {

	for (int i = 0; i < CSize; i++) {
		b_MorseSet[i] = false;

		if (m_Asf[i].isInSCC())
			b_MorseSet[i] = true;
	}
	*_bMorseSetData = b_MorseSet;
}
SCC::~SCC() {

	if (!dve) {
		delete[] Fr;
		delete[] Fc;
		delete[] Br;
		delete[] Bc;
	}
}


