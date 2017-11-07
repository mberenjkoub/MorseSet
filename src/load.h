/*
 * load.h
 *
 *  Created on: Oct 25, 2016
 *      Author: marzieh
 */

#ifndef LOAD_H_
#define LOAD_H_
#include "scc_kernel.h"
#include "streamline_kernel .h"

#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <vector>
#include <stdint.h>
//#include "graph_generator.h"
#include "scc_kernel.h"
#include "streamline.h"
#include <iostream>
#include <fstream>
#include <iostream>
#include <string>
#include <cassert>
#include <cstdlib>
#include <sstream>
#include <algorithm>



/* This is the name of the data file we will read. */
#define FILE_NAME_U "Ocean/UVEL.1440x720x50.19920102.nc"
#define FILE_NAME_V "Ocean/VVEL.1440x720x50.19920102.nc"
#define FILE_NAME_W "Ocean/WVEL.1440x720x50.19920102.nc"
/* We are reading 2D data, a 6 x 12 grid. */
#define NX 1440
#define NY 720
#define NZ 50
/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

using namespace std;

void loadGraph(const char *, uint32_t *, uint32_t *, Edge **, uint32_t **);

void loadFullGraph(const char *, uint32_t *, uint32_t *, Edge **, uint32_t **,
		Edge **, uint32_t **);

void generateData(ASF_vertex * _a,int oxxdim, int oydim, int ozdim, Dimension *_size, Boundary* _b, float* vx,
		float *vy, float*vz, Point** v);

void Load_AllEdges(Edge * Fc_, uint32_t * Fr_, Edge * Bc_, uint32_t * Br_,
		uint32_t vSize, uint32_t RSize_F, uint32_t RSize_B, Edge ** oFc,
		uint32_t ** oFr, Edge ** oBc, uint32_t ** oBr, uint32_t * oRSize);

void getOceanData(float** uvel, float** vvel, float** wvel);

void getTornadoData(int xs, int ys, int zs, int time, float *m_x1, float* m_y1,
		float* m_z1);

int getCylinderData(float* m_x1, float* m_y1, float* m_z1, Boundary b);

void generateDataCylinder(ASF_vertex* _a, Dimension *_size, Boundary* _b,
		int primxDim, float* vx, float* vy, float* vz, Point** v);

void generateDataOcean(ASF_vertex* _a, Dimension *_size, Boundary* _b,
		int primxDim, float* vx, float* vy, float* vz, Point** v);

void generateDataBenard(ASF_vertex* _a, Dimension *_size, Boundary* _b,
		int primxDim, float* vx, float* vy, float* vz, Point** v);

void generateDataTornado(ASF_vertex* _a, Dimension *_size, Boundary* _b,
		int primxDim, float* vx, float* vy, float* vz, Point** v);

void getBenardData(float** vx, float** vy, float** vz);

void getOceanData(float** uvel, float** vvel, float** wvel);

int getOceanRawData(float** vx, float** vy, float** vz);

void get_Lorenz_Field(Point *p, float& vx, float& vy, float& vz);

void getHuricaneData(float** vx, float** vy, float** vz);

void generateDataHurricane(ASF_vertex* _a, Dimension *_size, Boundary* _b,
		int primxDim, float* vx, float* vy, float* vz, Point** v);

void generateDataOcean(ASF_vertex* _a, Dimension *_size, Boundary* _b,
		int primxDim, float* vx, float* vy, float* vz, Point** v);

class load {
public:
	load();
	virtual ~load();

};

#endif /* LOAD_H_ */
