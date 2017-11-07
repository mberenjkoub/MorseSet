/*
 * SCC.h
 *
 *  Created on: Oct 22, 2016
 *      Author: marzieh
 */

#ifndef SCC_H_
#define SCC_H_


#include <stdint.h>
//#include "graph_generator.h"
#include "scc_kernel.h"
#include "streamline.h"
#include "load.h"

const int Temp_count = 3;
const float Time_limit = 50000.0;

//pair <uint32_t, float> FB_Decomposition(uint32_t, uint32_t, Edge *, uint32_t *, Edge *, uint32_t *, FB_vertex *, bool, bool quiet = false);
//
//pair <uint32_t, float> COL_Decomposition(uint32_t, uint32_t, Edge *, uint32_t *, Edge *, uint32_t *, COL_vertex *, ordering,
//	bool quiet = false);


class SCC
{
public:
	SCC(int xdim, int ydim, int zdim,int pxdim, float high_x, float low_x, float high_y, float low_y, float high_z, float low_z, float** _vx, float** _vy, float** _vz, int whichData);
	void Run(uint32_t tau);
	void getMorseSetResult(bool** _bMorseSetData);
	void loadOceanData(float** vx, float** vy, float** vz);
	void loadBenardData(float** vx, float** vy, float** vz);
	void loadTornadoData(float** vx, float** vy, float** vz);
	void getProbability(float** _fProbability);
	void getParameter(ASF_vertex** _a);//,  Boundary*_b, Dimension*_d,uint32_t* Rsize);
	void setParameter(ASF_vertex* _a, Edge* _Fc_f, uint32_t* _Fr_f, Edge* _Bc_b, uint32_t* _Br_b, uint32_t _CSize, uint32_t _RSize_F, uint32_t _RSize_B);
	void SetOutput(uint32_t** oFr, uint32_t **oBr, ASF_vertex**oa);
	void Run(uint32_t RSize, uint32_t CSize, Edge* Fc, uint32_t *Fr, Edge* Bc, uint32_t* Br);





	~SCC();

private:
	ASF_vertex* m_Asf;
	Edge * Fc_f; // forward columns
	uint32_t * Fr_f; // forward ranges

	Edge * Bc_b; // forward columns
	uint32_t * Br_b; // forward ranges
	uint32_t seednumber = 0;
	Dimension size;
	Boundary b;
	uint32_t tau;
	int primaryXDim;

	Point* vData;
	ASF_vertex* a;
	ASF_vertex* _oa;
	ASF_vertex* _a;
	OBF_vertex * _m;
	// CSR representation
	uint32_t CSize; // column arrays size
	uint32_t RSize; // range arrays size
	uint32_t RSize_F; // range arrays size
	uint32_t RSize_B; // range arrays size

	// Forwards arrays
	Edge * Fc; // forward columns
	uint32_t * Fr = NULL; // forward ranges
	// Backwards arrays
	Edge * Bc; // backward columns
	uint32_t * Br; // backward ranges
	//GraphGenerator * gen = NULL;
	bool dve = true;
	float* vx;
	float* vy;
	float* vz;


	uint32_t iState;
	bool *b_MorseSet;

	uint32_t whichDataType = -1;
};

#endif /* SCC_H_ */
