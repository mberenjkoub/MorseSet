/*
 * FileSaver.cpp
 *
 *  Created on: Nov 9, 2016
 *      Author: marzieh
 */

#include <FileSaver.h>

File_Saver::File_Saver() {
	// TODO Auto-generated constructor stub

}

File_Saver::~File_Saver() {
	// TODO Auto-generated destructor stub
}

void File_Saver::save_Quad_Face_original(ASF_vertex*m, fFace* fc, fEdge* eg,
		Dimension d, string dataName, uint32_t csize, uint32_t face_num,
		int curtau, int curstep) {
	int index[10][30];
	int PARTICLE = csize;
	int frame = 50;
	std::ofstream fout;
	string curtau_s;
	string curstep_s;
	stringstream ss;
	ss << curtau;
	ss >> curtau_s;
	stringstream ss2;
	ss2 << curstep;
	ss2 >> curstep_s;
	string fileName = "QuadMesh" + dataName + curtau_s + "--" + curstep_s
			+ ".vtk";
	printf("%s \n", fileName.c_str());
	fout.open(fileName.c_str(), ios::out);
	fout << "# vtk DataFile Version 2.0" << endl;
	fout << "hex mesh vtk data - converted from .off" << endl;
	fout << "ASCII" << endl;
	fout << "DATASET POLYDATA " << endl;
	fout << "POINTS " << PARTICLE << " float" << endl;
	for (int j = 0; j < PARTICLE; ++j) {
		//fout << data_information[frame-1][j][0] << " " << data_information[frame-1][j][1] << " " << data_information[frame-1][j][2] << endl;
		fout << m[j].p.x << " " << m[j].p.y << " " << m[j].p.z << endl;
	}

	int vcont = 0;
	int ranges[1000];
	//const int& unit_ = 59;
	//ranges[0] = 4560;

//	for (int j = 1; j < face_num; j++)
//	{
//		// VTK CELL_TYPE 11
//		if (fc[j].getFaceVertex() == idx && fc[j].dir == dir /*&& j == 98524*/)
//		{
//			ranges[vcont] = j;
//			vcont++;
//		}
//	}
	const int& cube_ = face_num;

	fout << "POLYGONS  " << cube_ << " " << cube_ * 5 << endl;

	int index_ = 0;	//idx;
	//fout << 4 << " " << fc[idx*3+dir-1].vertexes[0] << " " << fc[idx*3+dir-1].vertexes[1] << " " << fc[idx*3+dir-1].vertexes[2] << " " << fc[idx*3+dir-1].vertexes[3] << endl;

	for (int j = 0; j < cube_; j++) {
		//index_ = fc[j].vertexes[0];
		fout << 4 << " " << fc[j].vertexes[0] << " " << fc[j].vertexes[1] << " "
				<< fc[j].vertexes[2] << " " << fc[j].vertexes[3] << endl;

	}
	fout << "POINT_DATA " << PARTICLE << endl;
	fout << "SCALARS V_Scalars int " << endl;
	fout << "LOOKUP_TABLE V_Table " << endl;

	for (uint32_t row = 0; row < PARTICLE; row++) {
		fout << 10 << endl;
	}
}

void File_Saver::save_Outerapproximation_File(bool* _MorseSetData, Boundary*_b,
		Dimension*_d, Point* _step, int currenttau) {
	std::string str1 = "outer_approximation";

	std::ofstream fout;
	string frame_str;
	stringstream ss;
	ss << currenttau;
	ss >> frame_str;
//	fout.open(("frame" + frame_str + ".vtk").c_str(), ios::out);
	FILE* f_p1 = fopen((str1 + frame_str + ".vtk").c_str(), "wb");

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

void File_Saver::generate_voxelImage(ASF_vertex*m, Dimension* d,
		const int& frame, const int& PARTICLE) {
	std::ofstream fout;
	string frame_str;
	stringstream ss;
	ss << frame;
	ss >> frame_str;
	fout.open(("frame" + frame_str + ".vtk").c_str(), ios::out);
	fout << "# vtk DataFile Version 3.0" << endl;
	fout << "PBF example" << endl;
	fout << "ASCII" << endl;
	fout << "DATASET UNSTRUCTURED_GRID" << endl;
	fout << "POINTS " << PARTICLE << " float" << endl;
	for (int j = 0; j < PARTICLE; ++j) {
		//fout << data_information[frame-1][j][0] << " " << data_information[frame-1][j][1] << " " << data_information[frame-1][j][2] << endl;
		fout << m[j].p.x << " " << m[j].p.y << " " << m[j].p.z << endl;
	}
	const int& unit_ = d->x - 1;
	const int& cube_ = unit_ * unit_ * unit_;
	fout << "CELLS " << cube_ << " " << 9 * cube_ << endl;
	for (int j = 0; j < cube_; j++) {
		//first decide which z it lies on
		const int& z_ = j / (unit_ * unit_);
		const int& y_ = j % (unit_ * unit_) / unit_;
		const int& x_ = j % (unit_ * unit_) % unit_;
		const int& index_ = d->x * d->x * z_ + y_ * d->x + x_;
		// VTK CELL_TYPES 12
		//fout << 8 << " " << index_ << " " << index_+1 << " " << index_+61 << " " << index_+60 << " " << index_+3600 << " " << index_+3600+1 << " " << index_+ 3600+1+60 << " " << index_+3600+60 << endl;

		// VTK CELL_TYPE 11
		fout << 8 << " " << index_ << " " << index_ + 1 << " " << index_ + d->x
				<< " " << index_ + d->x + 1 << " " << index_ + (d->x * d->y)
				<< " " << index_ + (d->x * d->y) + 1 << " "
				<< index_ + (d->x * d->y) + d->x << " "
				<< index_ + (d->x * d->y) + d->x + 1 << endl;
	}
	fout << "CELL_TYPES " << cube_ << endl;
	for (int j = 0; j < cube_; j++) {
		fout << 11 << endl;
	}
	fout << "POINT_DATA " << PARTICLE << endl;
	fout << "SCALARS FTLE float 1" << endl;
	fout << "LOOKUP_TABLE velocity_table" << endl;
	for (int j = 0; j < PARTICLE; j++) {
		fout << 5 << endl;
	}
}

void File_Saver::save_Quad_Face_Index(ASF_vertex*m, fFace* fc,
		uint32_t* Fe_Face, fEdge* eg, Dimension d, string dataName,
		uint32_t csize, uint32_t face_num, int curtau, int curstep) {
	int index[10][30];
	int PARTICLE = csize;
	int frame = 50;
	std::ofstream fout;
	string curtau_s;
	string curstep_s;
	stringstream ss;
	ss << curtau;
	ss >> curtau_s;
	stringstream ss2;
	ss2 << curstep;
	ss2 >> curstep_s;
	string fileName = "/home/marzieh/cuda-workspace/Result/AllQuadMesh"
			+ dataName + curtau_s + "--" + curstep_s + ".vtk";
	printf("%s \n", fileName.c_str());
	fout.open(fileName.c_str(), ios::out);
	fout << "# vtk DataFile Version 2.0" << endl;
	fout << "hex mesh vtk data - converted from .off" << endl;
	fout << "ASCII" << endl;
	fout << "DATASET POLYDATA " << endl;
	fout << "POINTS " << PARTICLE << " float" << endl;
	for (int j = 0; j < PARTICLE; ++j) {
		//fout << data_information[frame-1][j][0] << " " << data_information[frame-1][j][1] << " " << data_information[frame-1][j][2] << endl;
		fout << m[j].e.x << " " << m[j].e.y << " " << m[j].e.z << endl;
	}

	int vcont = 0;
	int ranges[1000];
	//const int& unit_ = 59;
	//ranges[0] = 4560;

	for (int j = 0; j < face_num; j++) {
		// VTK CELL_TYPE 11
		if (Fe_Face[j] == 1 && fc[j].level > 1)	//if (!fc[j].bsplit && fc[j].isInBoundary())
				{

			vcont++;
		}
	}
	const int& cube_ = vcont;

	fout << "POLYGONS  " << cube_ << " " << cube_ * 5 << endl;

	int index_ = 0;	//idx;
	//fout << 4 << " " << fc[idx*3+dir-1].vertexes[0] << " " << fc[idx*3+dir-1].vertexes[1] << " " << fc[idx*3+dir-1].vertexes[2] << " " << fc[idx*3+dir-1].vertexes[3] << endl;

	for (int j = 0; j < face_num; j++) {
		//index_ = fc[j].vertexes[0];
		if (Fe_Face[j] == 1 && fc[j].level > 1)	//!fc[j].bsplit && fc[j].isInBoundary())
			fout << 4 << " " << fc[j].vertexes[0] << " " << fc[j].vertexes[1]
					<< " " << fc[j].vertexes[2] << " " << fc[j].vertexes[3]
					<< endl;

	}
	fout << "POINT_DATA " << PARTICLE << endl;
	fout << "SCALARS V_Scalars int " << endl;
	fout << "LOOKUP_TABLE V_Table " << endl;

	for (uint32_t row = 0; row < PARTICLE; row++) {
		fout << 10 << endl;
	}
}

void File_Saver::save_Quad_Face(ASF_vertex*m, fFace* fc, fEdge* eg, Dimension d,
		string dataName, uint32_t csize, uint32_t face_num, int curtau,
		int curstep) {
	int index[10][30];
	int PARTICLE = csize;
	int frame = 50;
	std::ofstream fout;
	string curtau_s;
	string curstep_s;
	stringstream ss;
	ss << curtau;
	ss >> curtau_s;
	stringstream ss2;
	ss2 << curstep;
	ss2 >> curstep_s;
	string fileName = "/home/marzieh/cuda-workspace/Result/AllQuadMesh"
			+ dataName + curtau_s + "--" + curstep_s + ".vtk";
	printf("%s \n", fileName.c_str());
	fout.open(fileName.c_str(), ios::out);
	fout << "# vtk DataFile Version 2.0" << endl;
	fout << "hex mesh vtk data - converted from .off" << endl;
	fout << "ASCII" << endl;
	fout << "DATASET POLYDATA " << endl;
	fout << "POINTS " << PARTICLE << " float" << endl;
	for (int j = 0; j < PARTICLE; ++j) {
		//fout << data_information[frame-1][j][0] << " " << data_information[frame-1][j][1] << " " << data_information[frame-1][j][2] << endl;
		fout << m[j].e.x << " " << m[j].e.y << " " << m[j].e.z << endl;
	}

	int vcont = 0;
	int ranges[1000];
	//const int& unit_ = 59;
	//ranges[0] = 4560;

	for (int j = 0; j < face_num; j++) {
		// VTK CELL_TYPE 11
		if (!fc[j].bsplit && fc[j].isInBoundary()) {

			vcont++;
		}
	}
	const int& cube_ = vcont;

	fout << "POLYGONS  " << cube_ << " " << cube_ * 5 << endl;

	int index_ = 0;	//idx;
	//fout << 4 << " " << fc[idx*3+dir-1].vertexes[0] << " " << fc[idx*3+dir-1].vertexes[1] << " " << fc[idx*3+dir-1].vertexes[2] << " " << fc[idx*3+dir-1].vertexes[3] << endl;

	for (int j = 0; j < face_num; j++) {
		//index_ = fc[j].vertexes[0];
		if (!fc[j].bsplit && fc[j].isInBoundary())
			fout << 4 << " " << fc[j].vertexes[0] << " " << fc[j].vertexes[1]
					<< " " << fc[j].vertexes[2] << " " << fc[j].vertexes[3]
					<< endl;

	}
	fout << "POINT_DATA " << PARTICLE << endl;
	fout << "SCALARS V_Scalars int " << endl;
	fout << "LOOKUP_TABLE V_Table " << endl;

	for (uint32_t row = 0; row < PARTICLE; row++) {
		fout << 10 << endl;
	}
}

void File_Saver::save_Quad_Face_Final(ASF_vertex*m, fFace* fc, fEdge* eg,
		Point** Trace, bool* MorseSet, Dimension d, Boundary*b, string dataName,
		uint32_t csize, uint32_t face_num, int curtau, int curstep) {
	int index[10][30];
	int PARTICLE = csize;
	int frame = 50;
	std::ofstream fout;
	string curtau_s;
	string curstep_s;
	stringstream ss;
	ss << curtau;
	ss >> curtau_s;
	stringstream ss2;
	ss2 << curstep;
	ss2 >> curstep_s;
	string fileName = "/home/marzieh/cuda-workspace/Test_Code_Face/" + dataName
			+ "/AllQuadMesh" + dataName + curtau_s + "--" + curstep_s + ".vtk";
	printf("%s \n", fileName.c_str());
	fout.open(fileName.c_str(), ios::out);
	fout << "# vtk DataFile Version 2.0" << endl;
	fout << "hex mesh vtk data - converted from .off" << endl;
	fout << "ASCII" << endl;
	fout << "DATASET POLYDATA " << endl;
	fout << "POINTS " << PARTICLE << " float" << endl;
	for (int j = 0; j < PARTICLE; ++j) {
		//fout << data_information[frame-1][j][0] << " " << data_information[frame-1][j][1] << " " << data_information[frame-1][j][2] << endl;
//		fout << m[j].es[curtau - 1].x << " " << m[j].es[curtau - 1].y << " "
//				<< m[j].es[curtau - 1].z << endl;

		fout << Trace[j][curtau].x << " " << Trace[j][curtau].y << " "
				<< Trace[j][curtau].z << endl;
		//fout << m[j].es[10].x << " " << m[j].es[10].y << " " << m[j].es[10].z	<< endl;
	}

	int vcont = 0;
	int ranges[1000];
	//const int& unit_ = 59;
	//ranges[0] = 4560;

//	for (int j = 1; j < face_num; j++)
//	{
//		// VTK CELL_TYPE 11
//		if (fc[j].getFaceVertex() == idx && fc[j].dir == dir /*&& j == 98524*/)
//		{
//			ranges[vcont] = j;
//			vcont++;
//		}
//	}
	int cube_ = 0;	//face_num;
	for (int j = 0; j < face_num; j++) {
		if (m[fc[j].vertexes[0]].checkInBoundary(b)
				|| m[fc[j].vertexes[1]].checkInBoundary(b)
				|| m[fc[j].vertexes[2]].checkInBoundary(b)
				|| m[fc[j].vertexes[3]].checkInBoundary(b)) {
			//index_ = fc[j].vertexes[0];
			if (!fc[j].bsplit
//					&& (MorseSet[m[fc[j].vertexes[0]].getOldRange()]
//							|| MorseSet[m[fc[j].vertexes[1]].getOldRange()]
//							|| MorseSet[m[fc[j].vertexes[2]].getOldRange()]
//							|| MorseSet[m[fc[j].vertexes[3]].getOldRange()])
			)
				cube_++;
		}
	}

	fout << "POLYGONS  " << cube_ << " " << cube_ * 5 << endl;

	int index_ = 0;	//idx;
	//fout << 4 << " " << fc[idx*3+dir-1].vertexes[0] << " " << fc[idx*3+dir-1].vertexes[1] << " " << fc[idx*3+dir-1].vertexes[2] << " " << fc[idx*3+dir-1].vertexes[3] << endl;

	for (int j = 0; j < face_num; j++) {
		//index_ = fc[j].vertexes[0];
		if (m[fc[j].vertexes[0]].checkInBoundary(b)
				|| m[fc[j].vertexes[1]].checkInBoundary(b)
				|| m[fc[j].vertexes[2]].checkInBoundary(b)
				|| m[fc[j].vertexes[3]].checkInBoundary(b)) {
			if (!fc[j].bsplit
//					&& (MorseSet[m[fc[j].vertexes[0]].getOldRange()]
//							|| MorseSet[m[fc[j].vertexes[1]].getOldRange()]
//							|| MorseSet[m[fc[j].vertexes[2]].getOldRange()]
//							|| MorseSet[m[fc[j].vertexes[3]].getOldRange()])
			)
				fout << 4 << " " << fc[j].vertexes[0] << " "
						<< fc[j].vertexes[1] << " " << fc[j].vertexes[2] << " "
						<< fc[j].vertexes[3] << endl;
		}

	}
	fout << "POINT_DATA " << PARTICLE << endl;
	fout << "SCALARS V_Scalars int " << endl;
	fout << "LOOKUP_TABLE V_Table " << endl;

	for (uint32_t row = 0; row < PARTICLE; row++) {
		fout << 10 << endl;
	}
}

void File_Saver::save_Quad_One_Face(ASF_vertex*m, fFace fc, fEdge* eg,
		Point** Trace, Dimension d, string dataName, uint32_t csize,
		uint32_t face_num, int curtau, int curstep, int scalarvalue) {
	int index[10][30];
	int PARTICLE = csize;
	int frame = 50;
	std::ofstream fout;
	string curtau_s;
	string curstep_s;
	stringstream ss;
	ss << curtau;
	ss >> curtau_s;
	stringstream ss2;
	ss2 << curstep;
	ss2 >> curstep_s;
	string fileName = "/home/marzieh/cuda-workspace/Result/QuadMesh" + dataName
			+ curtau_s + "--" + curstep_s + ".vtk";
	printf("%s \n", fileName.c_str());
	fout.open(fileName.c_str(), ios::out);
	fout << "# vtk DataFile Version 2.0" << endl;
	fout << "hex mesh vtk data - converted from .off" << endl;
	fout << "ASCII" << endl;
	fout << "DATASET POLYDATA " << endl;
	fout << "POINTS " << PARTICLE << " float" << endl;
	for (int j = 0; j < PARTICLE; ++j) {
		//fout << data_information[frame-1][j][0] << " " << data_information[frame-1][j][1] << " " << data_information[frame-1][j][2] << endl;
//		fout << m[j].es[curtau].x << " " << m[j].es[curtau].y << " "
//				<< m[j].es[curtau].z << endl;

		fout << Trace[j][curtau].x << " " << Trace[j][curtau].y << " "
				<< Trace[j][curtau].z << endl;
	}

	int vcont = 0;
	int ranges[1000];
	//const int& unit_ = 59;
	//ranges[0] = 4560;

//	for (int j = 1; j < face_num; j++)
//	{
//		// VTK CELL_TYPE 11
//		if (fc[j].getFaceVertex() == idx && fc[j].dir == dir /*&& j == 98524*/)
//		{
//			ranges[vcont] = j;
//			vcont++;
//		}
//	}
	const int& cube_ = face_num;

	fout << "POLYGONS  " << cube_ << " " << cube_ * 5 << endl;

	int index_ = 0;	//idx;
	//fout << 4 << " " << fc[idx*3+dir-1].vertexes[0] << " " << fc[idx*3+dir-1].vertexes[1] << " " << fc[idx*3+dir-1].vertexes[2] << " " << fc[idx*3+dir-1].vertexes[3] << endl;

	for (int j = 0; j < cube_; j++) {
		//index_ = fc[j].vertexes[0];
		fout << 4 << " " << fc.vertexes[0] << " " << fc.vertexes[1] << " "
				<< fc.vertexes[2] << " " << fc.vertexes[3] << endl;

	}
	fout << "POINT_DATA " << PARTICLE << endl;
	fout << "SCALARS V_Scalars int " << endl;
	fout << "LOOKUP_TABLE V_Table " << endl;

	for (uint32_t row = 0; row < PARTICLE; row++) {
		fout << scalarvalue << endl;
	}
}

void File_Saver::save_Quad_One_Face_parent(ASF_vertex*m, fFace fc, fEdge* eg,
		Dimension d, string dataName, uint32_t csize, uint32_t face_num,
		int curtau, int curstep, int scalarvalue) {
	int index[10][30];
	int PARTICLE = csize;
	int frame = 50;
	std::ofstream fout;
	string curtau_s;
	string curstep_s;
	stringstream ss;
	ss << curtau;
	ss >> curtau_s;
	stringstream ss2;
	ss2 << curstep;
	ss2 >> curstep_s;
	string fileName = "/home/marzieh/cuda-workspace/Test_Code_Face/" + dataName
			+ "/QuadMesh" + dataName + curtau_s + "parent--" + curstep_s
			+ ".vtk";
	printf("%s \n", fileName.c_str());
	fout.open(fileName.c_str(), ios::out);
	fout << "# vtk DataFile Version 2.0" << endl;
	fout << "hex mesh vtk data - converted from .off" << endl;
	fout << "ASCII" << endl;
	fout << "DATASET POLYDATA " << endl;
	fout << "POINTS " << PARTICLE << " float" << endl;
	for (int j = 0; j < PARTICLE; ++j) {
		//fout << data_information[frame-1][j][0] << " " << data_information[frame-1][j][1] << " " << data_information[frame-1][j][2] << endl;
		fout << m[j].e.x << " " << m[j].e.y << " " << m[j].e.z << endl;
	}

	int vcont = 0;
	int ranges[1000];
	//const int& unit_ = 59;
	//ranges[0] = 4560;

//	for (int j = 1; j < face_num; j++)
//	{
//		// VTK CELL_TYPE 11
//		if (fc[j].getFaceVertex() == idx && fc[j].dir == dir /*&& j == 98524*/)
//		{
//			ranges[vcont] = j;
//			vcont++;
//		}
//	}
	const int& cube_ = face_num;

	fout << "POLYGONS  " << cube_ << " " << cube_ * 5 << endl;

	int index_ = 0;	//idx;
	//fout << 4 << " " << fc[idx*3+dir-1].vertexes[0] << " " << fc[idx*3+dir-1].vertexes[1] << " " << fc[idx*3+dir-1].vertexes[2] << " " << fc[idx*3+dir-1].vertexes[3] << endl;

	for (int j = 0; j < cube_; j++) {
		//index_ = fc[j].vertexes[0];
		fout << 4 << " " << fc.vertexes[0] << " " << fc.vertexes[1] << " "
				<< fc.vertexes[2] << " " << fc.vertexes[3] << endl;

	}
	fout << "POINT_DATA " << PARTICLE << endl;
	fout << "SCALARS V_Scalars int " << endl;
	fout << "LOOKUP_TABLE V_Table " << endl;

	for (uint32_t row = 0; row < PARTICLE; row++) {
		fout << scalarvalue << endl;
	}
}

void File_Saver::save_Voxel(ASF_vertex*m, fFace* fc, int* voxelarray, fEdge* eg,
		Dimension d, string dataName, uint32_t csize, uint32_t voxel_num,
		int curtau, int curstep, int scalarvalue) {
	//int index[10][30];
	int PARTICLE = csize;
	int frame = 50;
	std::ofstream fout;
	string curtau_s;
	string curstep_s;
	stringstream ss;
	ss << curtau;
	ss >> curtau_s;
	stringstream ss2;
	ss2 << curstep;
	ss2 >> curstep_s;
	string fileName = "/home/marzieh/cuda-workspace/Result/outervoxel"
			+ curtau_s + "--" + curstep_s + ".vtk";
	printf("%s \n", fileName.c_str());
	fout.open(fileName.c_str(), ios::out);
	fout << "# vtk DataFile Version 2.0" << endl;
	fout << "hex mesh vtk data - converted from .off" << endl;
	fout << "ASCII" << endl;
	fout << "DATASET POLYDATA " << endl;
	fout << "POINTS " << PARTICLE << " float" << endl;
	for (int j = 0; j < PARTICLE; ++j) {
		//fout << data_information[frame-1][j][0] << " " << data_information[frame-1][j][1] << " " << data_information[frame-1][j][2] << endl;
		fout << m[j].p.x << " " << m[j].p.y << " " << m[j].p.z << endl;
	}

	int vcont = 0;
	int ranges[1000];
	//const int& unit_ = 59;
	//ranges[0] = 4560;

//	for (int j = 1; j < face_num; j++)
//	{
//		// VTK CELL_TYPE 11
//		if (fc[j].getFaceVertex() == idx && fc[j].dir == dir /*&& j == 98524*/)
//		{
//			ranges[vcont] = j;
//			vcont++;
//		}
//	}
	const int& cube_ = voxel_num * 5;

	fout << "POLYGONS  " << cube_ << " " << cube_ * 5 << endl;

	int index_ = 0;	//idx;
	//fout << 4 << " " << fc[idx*3+dir-1].vertexes[0] << " " << fc[idx*3+dir-1].vertexes[1] << " " << fc[idx*3+dir-1].vertexes[2] << " " << fc[idx*3+dir-1].vertexes[3] << endl;

	for (int j = 0; j < voxel_num; j++) {
		index_ = voxelarray[j];
		fout << 4 << " " << fc[index_ * 3].vertexes[0] << " "
				<< fc[index_ * 3].vertexes[1] << " "
				<< fc[index_ * 3].vertexes[2] << " "
				<< fc[index_ * 3].vertexes[3] << endl;

		fout << 4 << " " << fc[index_ * 3 + 1].vertexes[0] << " "
				<< fc[index_ * 3 + 1].vertexes[1] << " "
				<< fc[index_ * 3 + 1].vertexes[2] << " "
				<< fc[index_ * 3 + 1].vertexes[3] << endl;

		fout << 4 << " " << fc[index_ * 3 + 2].vertexes[0] << " "
				<< fc[index_ * 3 + 2].vertexes[1] << " "
				<< fc[index_ * 3 + 2].vertexes[2] << " "
				<< fc[index_ * 3 + 2].vertexes[3] << endl;

		fout << 4 << " " << fc[(index_ + d.x * d.y) * 3].vertexes[0] << " "
				<< fc[(index_ + d.x * d.y) * 3].vertexes[1] << " "
				<< fc[(index_ + d.x * d.y) * 3].vertexes[2] << " "
				<< fc[(index_ + d.x * d.y) * 3].vertexes[3] << endl;

		fout << 4 << " " << fc[(index_ + d.x) * 3 + 1].vertexes[0] << " "
				<< fc[(index_ + d.x) * 3 + 1].vertexes[1] << " "
				<< fc[(index_ + d.x) * 3 + 1].vertexes[2] << " "
				<< fc[(index_ + d.x) * 3 + 1].vertexes[3] << endl;

	}
	fout << "POINT_DATA " << PARTICLE << endl;
	fout << "SCALARS V_Scalars int " << endl;
	fout << "LOOKUP_TABLE V_Table " << endl;

	for (uint32_t row = 0; row < PARTICLE; row++) {
		fout << scalarvalue << endl;
	}
}

void File_Saver::start_save_Quad_One_Face(ASF_vertex*m, fFace fc, fEdge* eg,
		Dimension d, string dataName, uint32_t csize, uint32_t face_num,
		int curtau, int curstep) {
	int index[10][30];
	int PARTICLE = csize;
	int frame = 50;
	std::ofstream fout;
	string curtau_s;
	string curstep_s;
	stringstream ss;
	ss << curtau;
	ss >> curtau_s;
	stringstream ss2;
	ss2 << curstep;
	ss2 >> curstep_s;
	string fileName = "/home/marzieh/cuda-workspace/Result/start_QuadMesh" + dataName + curtau_s + "--" + curstep_s
			+ ".vtk";
	printf("%s \n", fileName.c_str());
	fout.open(fileName.c_str(), ios::out);
	fout << "# vtk DataFile Version 2.0" << endl;
	fout << "hex mesh vtk data - converted from .off" << endl;
	fout << "ASCII" << endl;
	fout << "DATASET POLYDATA " << endl;
	fout << "POINTS " << PARTICLE << " float" << endl;
	for (int j = 0; j < PARTICLE; ++j) {
		//fout << data_information[frame-1][j][0] << " " << data_information[frame-1][j][1] << " " << data_information[frame-1][j][2] << endl;
		fout << m[j].p.x << " " << m[j].p.y << " " << m[j].p.z << endl;
	}

	int vcont = 0;
	int ranges[1000];
//const int& unit_ = 59;
//ranges[0] = 4560;

//	for (int j = 1; j < face_num; j++)
//	{
//		// VTK CELL_TYPE 11
//		if (fc[j].getFaceVertex() == idx && fc[j].dir == dir /*&& j == 98524*/)
//		{
//			ranges[vcont] = j;
//			vcont++;
//		}
//	}
	const int& cube_ = face_num;

	fout << "POLYGONS  " << cube_ << " " << cube_ * 5 << endl;

	int index_ = 0;	//idx;
//fout << 4 << " " << fc[idx*3+dir-1].vertexes[0] << " " << fc[idx*3+dir-1].vertexes[1] << " " << fc[idx*3+dir-1].vertexes[2] << " " << fc[idx*3+dir-1].vertexes[3] << endl;

	for (int j = 0; j < cube_; j++) {
		//index_ = fc[j].vertexes[0];
		fout << 4 << " " << fc.vertexes[0] << " " << fc.vertexes[1] << " "
				<< fc.vertexes[2] << " " << fc.vertexes[3] << endl;

	}
	fout << "POINT_DATA " << PARTICLE << endl;
	fout << "SCALARS V_Scalars int " << endl;
	fout << "LOOKUP_TABLE V_Table " << endl;

	for (uint32_t row = 0; row < PARTICLE; row++) {
		fout << 10 << endl;
	}
}

void File_Saver::save_Quad_FaceWithSeedPoits(ASF_vertex*m, fFace* fc,
		ASF_vertex*seeds, fEdge* eg, Point** Trace, Dimension d,
		string dataName, uint32_t csize, uint32_t face_num, int seednum,
		int curtau, int curstep) {
	int index[10][30];
	int PARTICLE = seednum;
	int frame = 50;
	std::ofstream fout;
	string curtau_s;
	string curstep_s;
	stringstream ss;
	ss << curtau;
	ss >> curtau_s;
	stringstream ss2;
	ss2 << curstep;
	ss2 >> curstep_s;
	string fileName = "/home/marzieh/cuda-workspace/Result/QuadMeshandSeeds" + dataName + curtau_s + "--" + curstep_s
			+ ".vtk";
	printf("%s \n", fileName.c_str());
	fout.open(fileName.c_str(), ios::out);
	fout << "# vtk DataFile Version 2.0" << endl;
	fout << "hex mesh vtk data - converted from .off" << endl;
	fout << "ASCII" << endl;
	fout << "DATASET POLYDATA " << endl;
	fout << "POINTS " << PARTICLE << " float" << endl;
	for (int j = 0; j < PARTICLE; ++j) {
		//fout << data_information[frame-1][j][0] << " " << data_information[frame-1][j][1] << " " << data_information[frame-1][j][2] << endl;
//		fout << seeds[j].es[curtau].x << " " << seeds[j].es[curtau].y << " "
//				<< seeds[j].es[curtau].z << endl;

//		fout << Trace[j][curtau].x << " " << Trace[j][curtau].y << " "
//						<< Trace[j][curtau].z << endl;

		fout << seeds[j].e.x << " " << seeds[j].e.y << " " << seeds[j].e.z
				<< endl;
	}

	int vcont = 0;
	int ranges[1000];
//const int& unit_ = 59;
//ranges[0] = 4560;

//	for (int j = 1; j < face_num; j++)
//	{
//		// VTK CELL_TYPE 11
//		if (fc[j].getFaceVertex() == idx && fc[j].dir == dir /*&& j == 98524*/)
//		{
//			ranges[vcont] = j;
//			vcont++;
//		}
//	}
	const int& cube_ = face_num;

	fout << "VERTICES  " << seednum << " " << seednum * 2 << endl;

	for (int j = 0; j < seednum; j++) {
		//index_ = fc[j].vertexes[0];
		fout << 1 << " " << j << endl;

	}

	fout << "POINT_DATA " << PARTICLE << endl;
	fout << "SCALARS V_Scalars int " << endl;
	fout << "LOOKUP_TABLE V_Table " << endl;

	for (uint32_t row = 0; row < PARTICLE; row++) {
		fout << 50 << endl;
	}
}

void File_Saver::start_save_Quad_FaceWithSeedPoits(ASF_vertex*m, fFace* fc,
		ASF_vertex*seeds, fEdge* eg, Dimension d, string dataName,
		uint32_t csize, uint32_t face_num, int seednum, int curtau,
		int curstep) {
	int index[10][30];
	int PARTICLE = seednum;
	int frame = 50;
	std::ofstream fout;
	string curtau_s;
	string curstep_s;
	stringstream ss;
	ss << curtau;
	ss >> curtau_s;
	stringstream ss2;
	ss2 << curstep;
	ss2 >> curstep_s;
	string fileName = "/home/marzieh/cuda-workspace/Test_Code_Face/" + dataName
			+ "/start_QuadMeshandSeeds" + dataName + curtau_s + "--" + curstep_s
			+ ".vtk";
	printf("%s \n", fileName.c_str());
	fout.open(fileName.c_str(), ios::out);
	fout << "# vtk DataFile Version 2.0" << endl;
	fout << "hex mesh vtk data - converted from .off" << endl;
	fout << "ASCII" << endl;
	fout << "DATASET POLYDATA " << endl;
	fout << "POINTS " << PARTICLE << " float" << endl;
	for (int j = 0; j < PARTICLE; ++j) {
		//fout << data_information[frame-1][j][0] << " " << data_information[frame-1][j][1] << " " << data_information[frame-1][j][2] << endl;
		fout << seeds[j].p.x << " " << seeds[j].p.y << " " << seeds[j].p.z
				<< endl;
	}

	int vcont = 0;
	int ranges[1000];
//const int& unit_ = 59;
//ranges[0] = 4560;

//	for (int j = 1; j < face_num; j++)
//	{
//		// VTK CELL_TYPE 11
//		if (fc[j].getFaceVertex() == idx && fc[j].dir == dir /*&& j == 98524*/)
//		{
//			ranges[vcont] = j;
//			vcont++;
//		}
//	}
	const int& cube_ = face_num;

	fout << "VERTICES  " << seednum << " " << seednum * 2 << endl;

	for (int j = 0; j < seednum; j++) {
		//index_ = fc[j].vertexes[0];
		fout << 1 << " " << j << endl;

	}

	fout << "POINT_DATA " << PARTICLE << endl;
	fout << "SCALARS V_Scalars int " << endl;
	fout << "LOOKUP_TABLE V_Table " << endl;

	for (uint32_t row = 0; row < PARTICLE; row++) {
		fout << 50 << endl;
	}
}

void File_Saver::trilinearInterpolation(float p1[3], int idex, Boundary*b,
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

void File_Saver::getLorenzField1(float p[3], float& vx, float& vy, float& vz) {

	float sigma = 10.0;
	float ro = 28.0;
	float beta = 8.0 / 3.0;
	vx = sigma * (p[1] - p[0]);
	vy = (p[0] * (ro - p[2])) - p[1];
	vz = p[0] * p[1] - beta * p[2];

}

void File_Saver::SaveWholeEdges(ASF_vertex*m, fFace* fc, ASF_vertex*seeds,
		fEdge* eg, Dimension d, string dataName, uint32_t csize,
		uint32_t face_num, int seednum, int curtau, int curstep) {
	int index[10][30];
	int PARTICLE = seednum;
	int frame = 50;
	std::ofstream fout;
	string curtau_s;
	string curstep_s;
	stringstream ss;
	ss << curtau;
	ss >> curtau_s;
	stringstream ss2;
	ss2 << curstep;
	ss2 >> curstep_s;
	string fileName = "/home/marzieh/cuda-workspace/Test_Code_Face/edge"
			+ dataName + curtau_s + "--" + curstep_s + ".vtk";
	printf("%s \n", fileName.c_str());
	fout.open(fileName.c_str(), ios::out);
	fout << "# vtk DataFile Version 2.0" << endl;
	fout << "hex mesh vtk data - converted from .off" << endl;
	fout << "ASCII" << endl;
	fout << "DATASET POLYDATA " << endl;
	fout << "POINTS " << PARTICLE << " float" << endl;

}

void File_Saver::generalstreamlineTracing_single(float p[3], bool bForward,
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

	next_i = p[0];			//samples_x[j];
	next_j = p[1];			//samples_y[j];
	next_k = p[2];			//samples_z[j];
	for (int k = 0; k < tau; k++) {

		p2[0] = next_i;
		p2[1] = next_j;
		p2[2] = next_k;
		/*	if (m_bTornadoFieldSelected)
		 trilinearInterpolation2(next_i,next_j,next_k, start_pixel_id, vx, vy, vz);

		 else*/
		float vx1, vy1, vz1;
		trilinearInterpolation(p2, start_pixel_id, b, d, step, m_x1, m_y1, m_z1,
				currentDimX, vx, vy, vz);

		//getLorenzField1(p2,vx,vy,vz);

		if (vx != vx1 || vy1 != vy || vz1 != vz)
			counter_p = 0;

		//get_ABC_flow(next_i,next_j,next_k,vx,vy,vz);
		/*	if (WhichType == 0)
		 get_Lorenz_Field(next_i, next_j, next_k, vx, vy, vz);
		 else if (WhichType == 1)
		 get_ABC_flow(next_i, next_j, next_k, vx, vy, vz);*/
		//	get_Lorenz_Field(next_i, next_j, next_k, vx, vy, vz);
		float dist = sqrt(vx * vx + vy * vy + vz * vz);
		if (next_i < b->low.x || next_i > b->high.x || next_j < b->low.y
				|| next_j > b->high.y || next_k < b->low.z
				|| next_k > b->high.z) {
			counter_p = 0;
			//cout<<"ddd"<<endl;
			break;
		}
		if (dist < 1.0e-6) {
			counter_p = 0;
			//cout<<"ddd"<<endl;
			break;
		}

		vx = (vx / dist) * (step->x / 4.0);
		vy = (vy / dist) * (step->y / 4.0);
		vz = (vz / dist) * (step->z / 4.0);

		if (vx != 0 || vy != 0 || vz != 0)
			dist = 0;
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

void File_Saver::Save_Streamlines_Ocean(Point* op, bool bForward,
		string dataName, float* m_x1, float* m_y1, float* m_z1, Boundary*b,
		Dimension* d, Point* step, int currentDimX, int tau, int n) {

//==============================================================================
//==============================================================================

	int** index;
	index = new int*[n];
	for (int j = 0; j < n; j++)
		index[j] = new int[1000];
	int PARTICLE = n * tau;
	string dataType = "Lorenz";
	std::ofstream fout;
	string frame_str;
	stringstream ss;
	ss << tau;
	ss >> frame_str;
	string filename = "streamline" + dataName + "--" + frame_str + ".vtk";
	printf("streamline file = %s \n", filename.c_str());
	fout.open(("streamline" + dataName + "--" + frame_str + ".vtk").c_str(),
			ios::out);
	fout << "# vtk DataFile Version 2.0" << endl;
	fout << "hex mesh vtk data - converted from .off" << endl;
	fout << "ASCII" << endl;
	fout << "DATASET POLYDATA " << endl;
	fout << "POINTS " << PARTICLE << " float" << endl;
//	fFace face = fc[idx*3 + 2];

	for (int k = 0; k < n; k++) {

		float poin[1250][3];

		float p[3], ep[3];
		p[0] = op[k].x;
		p[1] = op[k].y;
		p[2] = op[k].z;

		fout << p[0] << " " << p[1] << " " << p[2] << endl;
		index[k][0] = k * tau;

		for (int j = 1; j < tau; j++) {
			generalstreamlineTracing_single(p, bForward, ep, m_x1, m_y1, m_z1,
					b, d, step, currentDimX, 1);
			poin[j][0] = p[0] = ep[0];
			poin[j][1] = p[1] = ep[1];
			poin[j][2] = p[2] = ep[2];
			index[k][j] = k * tau + j;

			//for (uint32_t row = 0; row < 20; row++)

			//	ASF_vertex v = a[row];
			/*	float p1[3], ep[3];
			 p1[0] = v.e.x;
			 p1[1] = v.e.y;
			 p1[2] = v.e.z;*/
			fout << poin[j][0] << " " << poin[j][1] << " " << poin[j][2]
					<< endl;
			if (poin[j][0] == 0.)
				printf("");
			//	cvertex++;

		}
	}

//for (uint32_t row = 0; row < currfaceSize; row++)
//{
//	fFace face = fc[row];

//	if (face.getFaceVertex() == idx && face.dir == 3)
//	{
//		//for (int i = 0; i < 4; i++)

//		{
//			//	if (a[face.vertexes[i]].getOldRange() == idx)
//			{

//				//	fout << 4 << " " << face.vertexes[0] << " " << face.vertexes[1] << " " << face.vertexes[2] << " " << face.vertexes[3] << endl;
//				cvertex++;
//			}

//		}
//	}
//}

	fout << endl;

	fout << "VERTICES   " << n << " " << n * 2 << endl;
	for (int k = 0; k < n; k++) {
		fout << 1 << " " << index[k][0];
		fout << endl;

	}
	fout << endl;

	fout << "LINES  " << n << " " << n * (tau + 1) << endl;
//	cvertex = 0;

	for (uint32_t row = 0; row < n; row++) {

		//for (int i = 0; i < 4; i++)

		//	if (a[face.vertexes[i]].getOldRange() == idx)
		{

			fout << tau << " " << index[row][0];
			for (int j = 1; j < tau; j++)
				fout << " " << index[row][j];// << " " << face.vertexes[2] << " " << face.vertexes[3] << endl;
			//cvertex++;
		}
		fout << endl;

	}

	fout << "POINT_DATA " << PARTICLE << endl;
	fout << "SCALARS V_Scalars int " << endl;
	fout << "LOOKUP_TABLE V_Table " << endl;

	for (uint32_t row = 0; row < PARTICLE; row++) {
		fout << 2 << endl;
	}

}

void File_Saver::Save_Streamlines2(Point* op, bool bForward, string dataName,
		float* m_x1, float* m_y1, float* m_z1, Boundary*b, Dimension* d,
		Point* step, int currentDimX, int tau, int n) {

//==============================================================================
//==============================================================================

	int** index;
	index = new int*[n];
	for (int j = 0; j < n; j++)
		index[j] = new int[300];
	int PARTICLE = n * tau;
	string dataType = "Lorenz";
	std::ofstream fout;
	string frame_str;
	stringstream ss;
	ss << tau;
	ss >> frame_str;
	string filename = "realstreamline" + dataName + "--" + frame_str + ".vtk";
	printf("streamline file = %s \n", filename.c_str());
	fout.open((filename).c_str(),
			ios::out);
	fout << "# vtk DataFile Version 2.0" << endl;
	fout << "hex mesh vtk data - converted from .off" << endl;
	fout << "ASCII" << endl;
	fout << "DATASET POLYDATA " << endl;
	fout << "POINTS " << PARTICLE << " float" << endl;
//	fFace face = fc[idx*3 + 2];

	for (int k = 0; k < n; k++) {

		float poin[3250][3];

		float p[3], ep[3];
		p[0] = op[k].x;
		p[1] = op[k].y;
		p[2] = op[k].z;

		fout << p[0] << " " << p[1] << " " << p[2] << endl;
		index[k][0] = k * tau;

		for (int j = 1; j < tau; j++) {
			generalstreamlineTracing_single(p, bForward, ep, m_x1, m_y1, m_z1,
					b, d, step, currentDimX, 1);
//			if(ep[0]< b->low.x || ep[0]> b->high.x ||  ep[1]< b->low.y || ep[1]> b->high.y || ep[2]< b->low.z || ep[2]> b->high.z)
//				continue;
			poin[j][0] = p[0] = ep[0];
			poin[j][1] = p[1] = ep[1];
			poin[j][2] = p[2] = ep[2];
			index[k][j] = k * tau + j;

			//for (uint32_t row = 0; row < 20; row++)

			//	ASF_vertex v = a[row];
			/*	float p1[3], ep[3];
			 p1[0] = v.e.x;
			 p1[1] = v.e.y;
			 p1[2] = v.e.z;*/
			fout << poin[j][0] << " " << poin[j][1] << " " << poin[j][2]
					<< endl;
			if (poin[j][0] == 0.)
				printf("");
			//	cvertex++;

		}
	}

//for (uint32_t row = 0; row < currfaceSize; row++)
//{
//	fFace face = fc[row];

//	if (face.getFaceVertex() == idx && face.dir == 3)
//	{
//		//for (int i = 0; i < 4; i++)

//		{
//			//	if (a[face.vertexes[i]].getOldRange() == idx)
//			{

//				//	fout << 4 << " " << face.vertexes[0] << " " << face.vertexes[1] << " " << face.vertexes[2] << " " << face.vertexes[3] << endl;
//				cvertex++;
//			}

//		}
//	}
//}

	fout << endl;

	fout << "VERTICES   " << n << " " << n * 2 << endl;
	for (int k = 0; k < n; k++) {
		fout << 1 << " " << index[k][0];
		fout << endl;

	}
	fout << endl;

	fout << "LINES  " << n << " " << n * (tau + 1) << endl;
//	cvertex = 0;

	for (uint32_t row = 0; row < n; row++) {

		//for (int i = 0; i < 4; i++)

		//	if (a[face.vertexes[i]].getOldRange() == idx)
		{

			fout << tau << " " << index[row][0];
			for (int j = 1; j < tau; j++)
				fout << " " << index[row][j];// << " " << face.vertexes[2] << " " << face.vertexes[3] << endl;
			//cvertex++;
		}
		fout << endl;

	}

	fout << "POINT_DATA " << PARTICLE << endl;
	fout << "SCALARS V_Scalars int " << endl;
	fout << "LOOKUP_TABLE V_Table " << endl;

	for (uint32_t row = 0; row < PARTICLE; row++) {
		fout << 20 << endl;
	}

}

void File_Saver::Save_Streamlines_estimated(ASF_vertex* op, Point** Trace,
		int* intarray, bool bForward, string dataName, float* m_x1, float* m_y1,
		float* m_z1, Boundary*b, Dimension* d, Point* step, int currentDimX,
		int tau, int n) {

//==============================================================================
//==============================================================================

	int** index;
	index = new int*[n];
	for (int j = 0; j < n; j++)
		index[j] = new int[1000];
	int PARTICLE = n * tau;
	string dataType = "Lorenz";
	std::ofstream fout;
	string frame_str;
	stringstream ss;
	ss << tau;
	ss >> frame_str;
	string filename = "streamline" + dataName + "--" + frame_str + ".vtk";
	printf("streamline file = %s \n", filename.c_str());
	fout.open(("streamline" + dataName + "--" + frame_str + ".vtk").c_str(),
			ios::out);
	fout << "# vtk DataFile Version 2.0" << endl;
	fout << "hex mesh vtk data - converted from .off" << endl;
	fout << "ASCII" << endl;
	fout << "DATASET POLYDATA " << endl;
	fout << "POINTS " << PARTICLE << " float" << endl;
	bool error = false;
//	fFace face = fc[idx*3 + 2];

	for (int i = 0; i < n; i++) {

		float poin[1250][3];

		int k = intarray[i];
		float p[3], ep[3];
		p[0] = op[k].p.x;
		p[1] = op[k].p.y;
		p[2] = op[k].p.z;
		if(p[0] == 0 || p[1] == 0 || p[2] == 0)
		{
			error = true;
			continue;
		}

//		p[0] = Trace[k][0].x;
//		p[1] = Trace[k][0].y;
//		p[2] = Trace[k][0].z;

		fout << p[0] << " " << p[1] << " " << p[2] << endl;
		index[i][0] = i * tau;

		for (int j = 1; j < tau; j++) {
//			generalstreamlineTracing_single(p, bForward, ep, m_x1, m_y1, m_z1,
//					b, d, step, currentDimX, 1);
//			if(ep[0]< b->low.x || ep[0]> b->high.x ||  ep[1]< b->low.y || ep[1]> b->high.y || ep[2]< b->low.z || ep[2]> b->high.z)
//				continue;
			if (Trace[k][j].x == 0){
				error = true;
				Trace[k][j].x = op[k].e.x;
				Trace[k][j].y = op[k].e.y;
				Trace[k][j].z = op[k].e.z;

			}
			poin[j][0] = p[0] = Trace[k][j].x;
			poin[j][1] = p[1] = Trace[k][j].y;
			poin[j][2] = p[2] = Trace[k][j].z;
			index[i][j] = i * tau + j;

			//for (uint32_t row = 0; row < 20; row++)

			//	ASF_vertex v = a[row];
			/*	float p1[3], ep[3];
			 p1[0] = v.e.x;
			 p1[1] = v.e.y;
			 p1[2] = v.e.z;*/
			fout << poin[j][0] << " " << poin[j][1] << " " << poin[j][2]
					<< endl;
			if (poin[j][0] == 0.)
				printf("");
			//	cvertex++;

		}
	}

	fout << endl;

	fout << "VERTICES   " << n << " " << n * 2 << endl;
	for (int k = 0; k < n; k++) {
		fout << 1 << " " << index[k][0];
		fout << endl;

	}
	fout << endl;

	fout << "LINES  " << n << " " << n * (tau + 1) << endl;
//	cvertex = 0;

	for (uint32_t row = 0; row < n; row++) {

		//for (int i = 0; i < 4; i++)

		//	if (a[face.vertexes[i]].getOldRange() == idx)
		{

			fout << tau << " " << index[row][0];
			for (int j = 1; j < tau; j++)
				fout << " " << index[row][j];// << " " << face.vertexes[2] << " " << face.vertexes[3] << endl;
			//cvertex++;
		}
		fout << endl;

	}

	fout << "POINT_DATA " << PARTICLE << endl;
	fout << "SCALARS V_Scalars int " << endl;
	fout << "LOOKUP_TABLE V_Table " << endl;

	for (uint32_t row = 0; row < PARTICLE; row++) {
		fout << 2 << endl;
	}

}

void File_Saver::Save_Streamlines_EndAdvection(ASF_vertex* op, int* intarray,
		bool bForward, string dataName, float* m_x1, float* m_y1, float* m_z1,
		Boundary*b, Dimension* d, Point* step, int currentDimX, int tau,
		int n) {

//==============================================================================
//==============================================================================

	int** index;
	index = new int*[n];
	for (int j = 0; j < n; j++)
		index[j] = new int[1000];
	int PARTICLE = n * tau;
	string dataType = "Lorenz";
	std::ofstream fout;
	string frame_str;
	stringstream ss;
	ss << tau;
	ss >> frame_str;
	string filename = "streamline" + dataName + "--" + frame_str + ".vtk";
	printf("streamline file = %s \n", filename.c_str());
	fout.open(("streamline" + dataName + "--" + frame_str + ".vtk").c_str(),
			ios::out);
	fout << "# vtk DataFile Version 2.0" << endl;
	fout << "hex mesh vtk data - converted from .off" << endl;
	fout << "ASCII" << endl;
	fout << "DATASET POLYDATA " << endl;
	fout << "POINTS " << PARTICLE << " float" << endl;
//	fFace face = fc[idx*3 + 2];

	for (int i = 0; i < n; i++) {

		float poin[1250][3];

		int k = intarray[i];
		float p[3], ep[3];
		p[0] = op[k].e.x;
		p[1] = op[k].e.y;
		p[2] = op[k].e.z;

		fout << p[0] << " " << p[1] << " " << p[2] << endl;
		index[i][0] = i * tau;

		for (int j = 1; j < tau; j++) {
			generalstreamlineTracing_single(p, bForward, ep, m_x1, m_y1, m_z1,
					b, d, step, currentDimX, 1);
//			if(ep[0]< b->low.x || ep[0]> b->high.x ||  ep[1]< b->low.y || ep[1]> b->high.y || ep[2]< b->low.z || ep[2]> b->high.z)
//				continue;
			poin[j][0] = p[0] = ep[0];			//Trace[k][j - 1].x;
			poin[j][1] = p[1] = ep[1];			//Trace[k][j - 1].y;
			poin[j][2] = p[2] = ep[2];			//Trace[k][j - 1].z;
			index[i][j] = i * tau + j;

			//for (uint32_t row = 0; row < 20; row++)

			//	ASF_vertex v = a[row];
			/*	float p1[3], ep[3];
			 p1[0] = v.e.x;
			 p1[1] = v.e.y;
			 p1[2] = v.e.z;*/
			fout << poin[j][0] << " " << poin[j][1] << " " << poin[j][2]
					<< endl;
			if (poin[j][0] == 0.)
				printf("");
			//	cvertex++;

		}
	}

	fout << endl;

	fout << "VERTICES   " << n << " " << n * 2 << endl;
	for (int k = 0; k < n; k++) {
		fout << 1 << " " << index[k][0];
		fout << endl;

	}
	fout << endl;

	fout << "LINES  " << n << " " << n * (tau + 1) << endl;
//	cvertex = 0;

	for (uint32_t row = 0; row < n; row++) {

		//for (int i = 0; i < 4; i++)

		//	if (a[face.vertexes[i]].getOldRange() == idx)
		{

			fout << tau << " " << index[row][0];
			for (int j = 1; j < tau; j++)
				fout << " " << index[row][j];// << " " << face.vertexes[2] << " " << face.vertexes[3] << endl;
			//cvertex++;
		}
		fout << endl;

	}

	fout << "POINT_DATA " << PARTICLE << endl;
	fout << "SCALARS V_Scalars int " << endl;
	fout << "LOOKUP_TABLE V_Table " << endl;

	for (uint32_t row = 0; row < PARTICLE; row++) {
		fout << 2 << endl;
	}

}

/* namespace std */
