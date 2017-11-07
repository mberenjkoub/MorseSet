/*
 * SCC2.h
 *
 *  Created on: Jan 4, 2017
 *      Author: marzieh
 */

#ifndef SCC2_H_
#define SCC2_H_


#include <stdint.h>
//#include "graph_generator.h"
#include "scc_kernel.h"
#include "streamline.h"
#include "load.h"



//pair <uint32_t, float> FB_Decomposition(uint32_t, uint32_t, Edge *, uint32_t *, Edge *, uint32_t *, FB_vertex *, bool, bool quiet = false);
//
//pair <uint32_t, float> COL_Decomposition(uint32_t, uint32_t, Edge *, uint32_t *, Edge *, uint32_t *, COL_vertex *, ordering,
//	bool quiet = false);


class SCC2
{
public:
	SCC2();
	SCC2(int xdim, int ydim, int zdim,int pxdim, float high_x, float low_x, float high_y, float low_y, float high_z, float low_z, float** _vx, float** _vy, float** _vz, int whichData);

	void getMorseSetResult(bool** _bMorseSetData);
	void loadOceanData(float** vx, float** vy, float** vz);
	void loadBenardData(float** vx, float** vy, float** vz);
	void loadTornadoData(float** vx, float** vy, float** vz);
	void getProbability(float** _fProbability);
	void getParameter(ASF_vertex** _a);//,  Boundary*_b, Dimension*_d,uint32_t* Rsize);
	void setParameter(ASF_vertex* _a, Edge* _Fc_f, uint32_t* _Fr_f, Edge* _Bc_b, uint32_t* _Br_b, uint32_t _CSize, uint32_t _RSize_F, uint32_t _RSize_B);
	void SetOutput(uint32_t** oFr, uint32_t **oBr, ASF_vertex**oa);
	void Run(ASF_vertex* _a, fFace *fc,
	fEdge*eg, Boundary* _b, Dimension* _d, Point* step, int _tau,
	float* m_x1, float* m_y1, float* m_z1, ASF_vertex** oVertex,
	Edge ** oFc, uint32_t ** oFr, fFace** oFace, fEdge** oEdge,
	uint32_t * oRSize, uint32_t * oFaceSize, uint32_t * oEdgeSize,
	uint32_t whichData, int CSize, uint32_t _curEdgeSize,
	uint32_t _curFaceSize, int curXDimension, int sampleSeeds,
	int MorseLevel, bool bForward);




	~SCC2();

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

#endif /* SCC2_H_ */
