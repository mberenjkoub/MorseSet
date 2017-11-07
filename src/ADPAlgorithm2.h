/*
 * ADPAlgorithm2.h
 *
 *  Created on: Dec 24, 2016
 *      Author: marzieh
 */

#ifndef ADPALGORITHM2_H_
#define ADPALGORITHM2_H_
#include "streamline_kernel .h"
#include "scc_kernel.h"
#include "math.h"
#include "FileSaver.h"
#include <cstdlib>
#define OUTBOUNDARY -1
namespace std {

class ADPAlgorithm2 {
public:
	ADPAlgorithm2();
	virtual ~ADPAlgorithm2();

	void SetIntegrationStep(Point _is);
		void Initialize(ASF_vertex * m, fEdge** _eg, fFace** _fr, Dimension* d,
				Boundary* b, Point* step, int whichDataset, int currSize,
				int curEdgeSize, int curFaceSize, int sampleSeeds,
				uint32_t num_rows);
		void Flow_Combinatorialization_C(ASF_vertex* _a, fFace *fc,
				fEdge*eg, Boundary* _b, Dimension* _d, Point* step, int _tau,
				float* m_x1, float* m_y1, float* m_z1, ASF_vertex** oVertex,
				Edge ** oFc, uint32_t ** oFr, fFace** oFace, fEdge** oEdge,
				uint32_t * oRSize, uint32_t * oFaceSize, uint32_t * oEdgeSize,
				uint32_t whichData, int CSize, uint32_t _curEdgeSize,
				uint32_t _curFaceSize, int curXDimension, int sampleSeeds,
				int MorseLevel, bool bForward);

		bool EstimateBestSeed(fFace* face, ASF_vertex*m, fFace* fc,
				fEdge* eg, int currSize, ASF_vertex* v1, ASF_vertex* v2, ASF_vertex* v3,
				ASF_vertex* v4, ASF_vertex* vm1, ASF_vertex* vm2, ASF_vertex* vm3,
				ASF_vertex* vm4, Boundary*_b, Dimension*_d, Point*_step, float* m_x1,
				float* m_y1, float* m_z1, int currentXDim, int tau, bool bForward,
				ASF_vertex*ov, int sampleNum,int curVertexId, bool bEstimated);

		void Flow_Combinatorialization_C2(ASF_vertex* _a, fFace *fc,
				fEdge*eg, Boundary* _b, Dimension* _d, Point* step, int _tau,
				float* m_x1, float* m_y1, float* m_z1, ASF_vertex** oVertex,
				Edge ** oFc, uint32_t ** oFr, fFace** oFace, fEdge** oEdge,
				uint32_t * oRSize, uint32_t * oFaceSize, uint32_t * oEdgeSize,
				uint32_t whichData, int CSize, uint32_t _curEdgeSize,
				uint32_t _curFaceSize, int curXDimension, int sampleSeeds,
				int MorseLevel, bool bForward);


		void CheckNeighborhood(ASF_vertex*m, fEdge* eg, uint32_t*Fr, Dimension* d,
				Boundary* b, Point* step, bool bForward, uint32_t original_num_rows,
				uint32_t num_rows, uint32_t level,int tau);
		bool checkEdge(ASF_vertex vertex1, ASF_vertex vertex2, Boundary*b,
				Dimension* d, Point* step, bool bForward,int tau,int v1i,int v2i);
		void SplitEndEdge(ASF_vertex*m, fFace* fc, fEdge* eg, uint32_t* Fe_Edge,
				uint32_t* Fr_Edge, Dimension* d, Boundary* b, Point* step,
				float*m_x1, float* m_y1, float*m_z1, int currentXDim, bool bForward,
				uint32_t num_vertex, uint32_t num_edges, uint32_t level, int tau);
		bool DivideEdges(ASF_vertex* v1, ASF_vertex* v2, uint32_t range1,
				Boundary* _b, Dimension* _d, Point* step, float*m_x1, float* m_y1,
				float*m_z1, int currentXDim, uint32_t faceNum, bool bForward, int n,
				int tau, ASF_vertex* otv ,int curvertexid,int v1i,int v2i,bool bEstimated);
		void Tracing_c(ASF_vertex * m, float* m_x1, float* m_y1, float* m_z1,
				Dimension* d, Boundary* b, Point* step, bool bForward,
				uint32_t whichData, int currentXDim,uint32_t start_rows, uint32_t num_rows,
				uint32_t level,int currtau);
		void generalstreamlineTracing_single(float p[3], bool bForward, float e[3],
				float* m_x1, float* m_y1, float* m_z1, Boundary*b, Dimension* d,
				Point* step, int currentDimX, int tau);
		void trilinearInterpolation(float p1[3]/*, float p1[3],*/, int idex,
				Boundary*b, Dimension* d, Point* step, float* m_x1, float* m_y1,
				float* m_z1, int primaryXDIMENSION, float& vx, float&vy, float& vz);
		void CheckFace(ASF_vertex*m, fEdge* eg, fFace* fc, uint32_t*Fe_Edge,
				uint32_t* Fe_Face, Dimension* d, uint32_t num_face);

		void SplitFace(ASF_vertex *m, fFace* fc, fEdge* eg, uint32_t* Fr_Face,
				uint32_t* Fe_Face, uint32_t* Fe, Dimension* d, Boundary* b,
				Point* step, bool bForward, uint32_t num_vertex, uint32_t num_edges,
				uint32_t num_faces, int curtau, int splitnumInCurtau, int whichData,
				float* m_x1, float* m_y1, float* m_z1,int tau);
		void SplitRemainingEdge(ASF_vertex*m, fFace* fc, fEdge* eg,
				uint32_t* Fe_Edge, uint32_t* Fr_Edge, uint32_t* Fe_Face,
				Dimension* d, Boundary* b, Point* step, float*m_x1, float* m_y1,
				float*m_z1, int currentXDim, int tau, bool bForward,
				uint32_t num_vertex,uint32_t oldedgeSize, uint32_t num_edges);
		void CheckRemainingEdge(ASF_vertex*m, fFace*fc, fEdge* eg, uint32_t*Fr,
				uint32_t* Fe_Face, uint32_t* Fe_Edge2, Dimension* d, Boundary* b,
				Point* step, bool bForward, uint32_t num_rows);

		bool face_correctnessCheck(fFace face, fEdge* eg);
		void UpdatesubEdge(fEdge* eg, fEdge child1, int row, int curFaceId1,
				int curFaceId2);
		void FindSubEdge(fEdge* eg, fFace face, fEdge* child1, fEdge* child2,
				int whichedge);
		void error(int errorcode);
		void InsertFaceCenterWithMultipleSeeds(ASF_vertex *m, fFace* fc, fEdge* eg,
				uint32_t* Fr_Face, uint32_t* Fe_Face, uint32_t* Fe, Dimension* d,
				Boundary* b, Point* step, float* m_x1, float* m_y1, float* m_z1,
				int tau, int currentXDim, bool bForward, uint32_t num_vertex,
				uint32_t num_edges, uint32_t num_faces);

		void generate_streamlines_sparsely(ASF_vertex* m, bool bForward,
				int currData, float* m_x1, float* m_y1, float* m_z1, Boundary*b,
				Dimension* d, Point* step, int currentDimX, int tau, int n);

	//	bool FindBestSeed(fFace* face, ASF_vertex*m, fFace* fc, fEdge* eg,
	//			int currSize, ASF_vertex* v1, ASF_vertex* v2, ASF_vertex* v3,
	//			ASF_vertex* v4, ASF_vertex* vm1, ASF_vertex* vm2, ASF_vertex* vm3,
	//			ASF_vertex* vm4, Boundary*_b, Dimension*_d, Point*_step,
	//			float* m_x1, float* m_y1, float* m_z1, int currentXDim, int tau,
	//			bool bForward, ASF_vertex*ov, int sampleNum);

		void AdvectParticle(ASF_vertex* vedge, float* m_x1, float* m_y1,
				float* m_z1, Boundary* b, Dimension* d, Point* step,
				int currentXDim, int tau, bool bForward);

		int CheckSeedInFace(fFace* face, ASF_vertex vcenter, ASF_vertex v1, ASF_vertex v2,
				ASF_vertex v3, ASF_vertex v4, Boundary*_b, Dimension* _d,
				Point* _step, bool bForward,int tau);
		void Compute_Seed_Point(ASF_vertex *m, fFace face, int sampleNum,
				Point** outPoints);

		void TestInsertFaceCenterWithMultipleSeeds(ASF_vertex *m, fFace* fc,
				fEdge* eg, uint32_t* Fr_Face, uint32_t* Fe_Face, uint32_t* Fe,
				Dimension* d, Boundary* b, Point* step, float* m_x1, float* m_y1,
				float* m_z1, int tau, int currentXDim, bool bForward,
				uint32_t num_vertex, uint32_t num_edges, uint32_t num_faces);


		bool checkPointInFaceBoundary(ASF_vertex*m, fFace face, ASF_vertex vcenter);

		bool check_child_edge(fEdge edge, fEdge* eg, uint32_t*Fe_Edge);

		void CheckFace_vertices(ASF_vertex*m, fEdge* eg, fFace* fc,
				uint32_t*Fe_Edge, uint32_t* Fe_Face, Boundary* b, Dimension* d,
				Point* step, bool bForward, uint32_t num_face,int tau);

		void Split_One_Edge(ASF_vertex*m, fFace* fc, fEdge* eg, uint32_t* Fe_Edge,
				uint32_t* Fr_Edge, Dimension* d, Boundary* b, Point* step,
				float*m_x1, float* m_y1, float*m_z1, int currentXDim, int tau,
				bool bForward, uint32_t num_vertex, uint32_t num_edges,
				uint32_t row);


		void getLorenzField1(float p[3], float& vx, float& vy, float& vz) ;






		bool CheckOneFace_vertices(ASF_vertex* currA, fEdge* eg, fFace* face,
				Boundary* b, Dimension* d, Point* step, bool bForward,int tau);

		void FindBoundingBox(ASF_vertex*m, ASF_vertex* oa, int pointnumber,
				fFace face, ASF_vertex* bbox, int bbpointnumber);

		void generate_streamlines_sparsely_Ocean(ASF_vertex* m, bool bForward,
				int currData, float* m_x1, float* m_y1, float* m_z1, Boundary*b,
				Dimension* d, Point* step, int currentDimX, int tau, int n);
		float ComputeJacobian(ASF_vertex*a, fFace face, float* m_x1,
				float* m_y1, float* m_z1, Point* step);
		void ComputeJacobianForEachFace(ASF_vertex*a, fFace* fc, int numface,
				float* m_x1, float* m_y1, float* m_z1,Boundary*b, Point* step,float& minJ, float &maxJ);

		void AdvectParticle_estimated(ASF_vertex* vedge, float* m_x1, float* m_y1,
				float* m_z1, Boundary* b, Dimension* d, Point* step, int currentXDim,
				int tau, bool bForward);

		void Tracing_c_estimated(ASF_vertex * m, float* m_x1,
				float* m_y1, float* m_z1, Dimension* d, Boundary* b, Point* step,
				bool bForward, uint32_t whichData, int currentXDim,
				uint32_t start_num_rows, uint32_t num_rows, uint32_t level) ;

		void SplitFace2(ASF_vertex *m, fFace* fc, fEdge* eg,
				uint32_t* Fr_Face, uint32_t* Fe_Face, uint32_t* Fe, Dimension* d,
				Boundary* b, Point* step, bool bForward, uint32_t num_vertex,
				uint32_t num_edges, uint32_t num_faces, int curtau,
				int splitnumInCurtau, int whichData, float* m_x1, float* m_y1,
				float* m_z1, int tau);

	private:
		File_Saver * fs;
		Point* IntegrationStep;
		float* Jacobianarray;
		int numJ;
		float minvx;
		float minvy;
		float maxvx;
		float maxvy;
		float minvz;
		float maxvz;
		Point** Trace;
		bool ESTIMATED = false;


};

} /* namespace std */

#endif /* ADPALGORITHM2_H_ */
