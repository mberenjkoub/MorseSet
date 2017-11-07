/*
 * CPUDriver.h
 *
 *  Created on: Oct 22, 2016
 *      Author: marzieh
 */

#ifndef CPUDRIVER_H_
#define CPUDRIVER_H_
#include "streamline_kernel .h"
#include "SCC.h"

#include "ADPAlgorithm2.h"
#include "SCC2.h"

class CPU_Driver {
	enum Color {
		Lorenz, Benard, Cylinder, Ocean, Tornado
	};
public:
	CPU_Driver();
	virtual ~CPU_Driver();
	void Init(int currentStepNumber, int currentDataset);
	void InitializeLorenzFieldParameter();
	void InitializeTornadoFieldParameter();
	void InitializeBenardFieldParameter();
	void Run();
	void save_VTK_File(bool* _MorseSetData);
	void InitializeHurricaneFieldParameter();

	void InitializeOceanFieldParameter();
	void InitializeCylinderFieldParameter();
	void keep_Morse_Data(fFace* fc, fEdge* eg, ASF_vertex* m, fFace** oFc,
			fEdge** oeg, ASF_vertex** om, bool* MorseSet, uint32_t numEdges,
			uint32_t numFaces, uint32_t num_vertexes, int samplesize,
			int* oFaceSize, int* oVertexSize, int* oEdgeSize);

	void keep_Morse_Data2(fFace* fc, fEdge* eg, ASF_vertex* m, fFace** oFc,
			fEdge** oeg, ASF_vertex** om, bool* MorseSet, uint32_t numEdges,
			uint32_t numFaces, uint32_t num_vertexes, int samplesize,
			int* oFaceSize, int* oVertexSize, int* oEdgeSize);
	void Run_GPU();
	int currenttau;
	int currentData;
	int currentXDimension;
	Boundary *_b;
	Dimension* _d;
	Point* _step;
	float* m_x1;
	float* m_y1;
	float* m_z1;
	SCC* scc;
	int* nreference;
	int* nedgereference;
	int* nfacereference;
	int MorseLevel;

};

#endif /* CPUDRIVER_H_ */
