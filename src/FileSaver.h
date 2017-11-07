/*
 * FileSaver.h
 *
 *  Created on: Nov 9, 2016
 *      Author: marzieh
 */

#ifndef FILESAVER_H_
#define FILESAVER_H_
#include "streamline_kernel .h"
#include <string.h>

class File_Saver {
public:
	File_Saver();
	virtual ~File_Saver();

	void getLorenzField1(float p[3], float& vx, float& vy, float& vz);

	void Save_Streamlines2(Point* op, bool bForward, string dataName,
			float* m_x1, float* m_y1, float* m_z1, Boundary*b, Dimension* d,
			Point* step, int currentDimX, int tau, int n);
	void generalstreamlineTracing_single(float p[3], bool bForward, float e[3],
			float* m_x1, float* m_y1, float* m_z1, Boundary*b, Dimension* d,
			Point* step, int currentDimX, int tau);
	void trilinearInterpolation(float p1[3], int idex, Boundary*b, Dimension* d,
			Point* step, float* m_x1, float* m_y1, float* m_z1,
			int primaryXDIMENSION, float& vx, float&vy, float& vz);
	void save_Quad_Face(ASF_vertex*m, fFace* fc, fEdge* eg, Dimension d,
			string dataName, uint32_t csize, uint32_t face_num, int curtau,
			int curstep);
	void save_Quad_Face_original(ASF_vertex*m, fFace* fc, fEdge* eg,
			Dimension d, string dataName, uint32_t csize, uint32_t face_num,
			int curtau, int curstep);
	void generate_voxelImage(ASF_vertex*m, Dimension* d, const int& frame,
			const int& PARTICLE);

	void save_Quad_FaceWithSeedPoits(ASF_vertex*m, fFace* fc, ASF_vertex*seeds,
			fEdge* eg, Point** Trace, Dimension d, string dataName,
			uint32_t csize, uint32_t face_num, int seednum, int curtau,
			int curstep);

	void save_Quad_One_Face(ASF_vertex*m, fFace fc, fEdge* eg, Point** Trace,
			Dimension d, string dataName, uint32_t csize, uint32_t face_num,
			int curtau, int curstep, int scalarvalue);

	void start_save_Quad_One_Face(ASF_vertex*m, fFace fc, fEdge* eg,
			Dimension d, string dataName, uint32_t csize, uint32_t face_num,
			int curtau, int curstep);

	void start_save_Quad_FaceWithSeedPoits(ASF_vertex*m, fFace* fc,
			ASF_vertex*seeds, fEdge* eg, Dimension d, string dataName,
			uint32_t csize, uint32_t face_num, int seednum, int curtau,
			int curstep);

	void save_Outerapproximation_File(bool* _MorseSetData, Boundary*_b,
			Dimension*_d, Point* _step, int currenttau);

	void save_Voxel(ASF_vertex*m, fFace* fc, int* voxelarray, fEdge* eg,
			Dimension d, string dataName, uint32_t csize, uint32_t voxel_num,
			int curtau, int curstep, int scalarvalue);

	void save_Quad_One_Face_parent(ASF_vertex*m, fFace fc, fEdge* eg,
			Dimension d, string dataName, uint32_t csize, uint32_t face_num,
			int curtau, int curstep, int scalarvalue);

	void SaveWholeEdges(ASF_vertex*m, fFace* fc, ASF_vertex*seeds, fEdge* eg,
			Dimension d, string dataName, uint32_t csize, uint32_t face_num,
			int seednum, int curtau, int curstep);
	void Save_Streamlines_Ocean(Point* op, bool bForward, string dataName,
			float* m_x1, float* m_y1, float* m_z1, Boundary*b, Dimension* d,
			Point* step, int currentDimX, int tau, int n);

	void Save_Streamlines_estimated(ASF_vertex* op, Point** Trace,
			int* intarray, bool bForward, string dataName, float* m_x1,
			float* m_y1, float* m_z1, Boundary*b, Dimension* d, Point* step,
			int currentDimX, int tau, int n);
	void save_Quad_Face_Final(ASF_vertex*m, fFace* fc, fEdge* eg, Point** Trace,
			bool* MorseSet, Dimension d, Boundary*b, string dataName,
			uint32_t csize, uint32_t face_num, int curtau, int curstep);

	void Save_Streamlines_EndAdvection(ASF_vertex* op, int* intarray,
			bool bForward, string dataName, float* m_x1, float* m_y1,
			float* m_z1, Boundary*b, Dimension* d, Point* step, int currentDimX,
			int tau, int n);

	void save_Quad_Face_Index(ASF_vertex*m, fFace* fc,
			uint32_t* index, fEdge* eg, Dimension d, string dataName,
			uint32_t csize, uint32_t face_num, int curtau, int curstep);

};

/* namespace std */

#endif /* FILESAVER_H_ */
