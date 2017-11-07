/*
 * load.cpp
 *
 *  Created on: Oct 25, 2016
 *      Author: marzieh
 */

#include <load.h>

void get_Lorenz_Field(Point *p, float& vx, float& vy, float& vz) {

	float sigma = 10.0;
	float ro = 28.0;
	float beta = 8.0 / 3.0;
	vx = sigma * (p->y - p->x);
	vy = (p->x * (ro - p->z)) - p->y;
	vz = p->x * p->y - beta * p->z;
}

void generateData(ASF_vertex* _a, int oxxdim, int oydim, int ozdim,
		Dimension *_size, Boundary* _b, float* vx, float *vy, float*vz,
		Point** v) {

	Point* tempV = new Point[_size->x * _size->y * _size->z];
	Point step = (_b->high - _b->low); // / (*_size);
	step /= (*_size);
	//float step_x = (_b->x_high - _b->x_low) / (_size->x - 1);
	//float step_y = (_b->y_high - _b->y_low) / (_size->y - 1);
	//float step_z = (_b->z_high - _b->z_low) / (_size->z - 1);

	int counter = 0;

	for (int k = 0; k < _size->z; k++) {
		for (int j = 0; j < _size->y; j++) {

			for (int i = 0; i < _size->x; i++) {
				uint32_t index = i + j * _size->x + k * (_size->x * _size->y);

				//_a[counter].e.x =
				_a[counter].p.x = _b->low.x + i * step.x;

				_a[counter].p.y = _b->low.y + j * step.y;
				_a[counter].p.z = _b->low.z + k * step.z;
				_a[counter].e = _a[counter].p;
				//	_a[counter].eb = _a[counter].p;
				_a[counter].setOldRange(index);
				_a[counter].setInBoundary();

				counter++;

			}
			//buffer[counter++] = lowBoundary + j*step;
		}
		//buffer[counter++] = lowBoundary + i*step;
	}

	float vvx, vvy, vvz;

	Point originalstep;
	Point p;
	originalstep.x = (_b->high.x - _b->low.x) / (oxxdim - 1);
	originalstep.y = (_b->high.y - _b->low.y) / (oydim - 1);
	originalstep.z = (_b->high.z - _b->low.z) / (ozdim - 1);

	for (int k = 0; k < oxxdim; k++) {
		for (int j = 0; j < oydim; j++) {

			for (int i = 0; i < ozdim; i++) {
				uint32_t index = i + j * oxxdim + k * (oxxdim * oydim);
				p.x = _b->low.x + i * originalstep.x;
				p.y = _b->low.y + j * originalstep.y;
				p.z = _b->low.z + k * originalstep.z;
				get_Lorenz_Field(&p, vvx, vvy, vvz);
				vx[index] = vvx;
				vy[index] = vvy;
				vz[index] = vvz;
			}
		}
	}

	*v = tempV;

}

int getOceanRawData(float** vx, float** vy, float** vz) {
	int size = NX * NY * NZ;

	string filename1 = "UVEL.1440x720x50.19920102.txt";
	string filename2 = "VVEL.1440x720x50.19920102.txt";

	string filename3 = "WVEL.1440x720x50.19920102.txt";
	float* tvx = new float[size];
	float* tvy = new float[size];
	float* tvz = new float[size];

	std::ifstream myfile1(filename1.c_str());
	std::ifstream myfile2(filename2.c_str());
	std::ifstream myfile3(filename3.c_str());

	if (myfile1.is_open()) {

		for (int i = 0; i < size; i++) {
			myfile1 >> tvx[i];
			myfile2 >> tvy[i];
			myfile3 >> tvz[i];
		}
		myfile1.close();
		myfile2.close();
		myfile3.close();
	}

	*vx = tvx;
	*vy = tvy;
	*vz = tvz;
//	FILE *fp;
//
//	if ((fp = fopen(filename.c_str(), "rb")) == NULL)
//		return 0;
//
////	 fseek(fp, 0, SEEK_END);
////	 long size = ftell(fp);
////	 fseek(fp, 0, SEEK_SET);
//
//	float *uvel = (float *) malloc(sizeof(float) * size);
//	if (uvel == NULL) {
//		fclose(fp);
//		return 0;
//	}
//	int ret = fread(uvel, sizeof(float), size, fp);
//
////	if (fread(uvel, sizeof(float), size, fp) != size) {
////		fclose(fp);
////		return 0;
////	}
//
//	fclose(fp);
//
//	filename = "VVEL.1440x720x50.19920102.raw";
//
//	if ((fp = fopen(filename.c_str(), "rb")) == NULL)
//		return 0;
//
//	float *vvel = (float *) malloc(sizeof(float) * size);
//	if (vvel == NULL) {
//		fclose(fp);
//		return 0;
//	}
//	ret = fread(vvel, sizeof(float), size, fp);
//
////	if (fread(vvel, sizeof(float), size, fp) != size) {
////		fclose(fp);
////		return 0;
////	}
//
//	fclose(fp);
//
//	filename = "VVEL.1440x720x50.19920102.raw";
//
//	if ((fp = fopen(filename.c_str(), "rb")) == NULL)
//		return 0;
//
//	float *wvel = (float *) malloc(sizeof(float) * size);
//	if (wvel == NULL) {
//		fclose(fp);
//		return 0;
//	}
//	ret = fread(wvel, sizeof(float), size, fp);
//
////	if (fread(wvel, sizeof(float), size, fp) != size) {
////		fclose(fp);
////		return 0;
////	}
//
//	fclose(fp);

//	*vx = uvel;
//	*vy = vvel;
//	*vz = vvel;

}

void generateDataOcean(ASF_vertex* _a, Dimension *_size, Boundary* _b,
		int primxDim, float* vx, float* vy, float* vz, Point** v) {

	Point* tempV = new Point[_size->x * _size->y * _size->z];

	Point step = (_b->high - _b->low);	 // / (*_size);
	step /= (*_size);
	//float step_x = (_b->x_high - _b->x_low) / (_size->x - 1);
	//float step_y = (_b->y_high - _b->y_low) / (_size->y - 1);
	//float step_z = (_b->z_high - _b->z_low) / (_size->z - 1);

	int counter = 0;

	for (int k = 0; k < _size->z; k++) {
		for (int j = 0; j < _size->y; j++) {

			for (int i = 0; i < _size->x; i++) {
				uint32_t index = i + j * _size->x + k * (_size->x * _size->y);

				//_a[counter].e.x =
				_a[counter].e.x = _a[counter].p.x = _b->low.x + i * step.x;

				_a[counter].e.y = _a[counter].p.y = _b->low.y + j * step.y;
				_a[counter].e.z = _a[counter].p.z = _b->low.z + k * step.z;

				_a[counter].setOldRange(index);
				_a[counter].setInBoundary();

				if (_size->x < primxDim) {
					tempV[counter].x = vx[index * 2];
					tempV[counter].y = vy[index * 2];
					tempV[counter].z = vz[index * 2];
				} else if (_size->x > primxDim) {
					tempV[counter].x = vx[index / 2];
					tempV[counter].y = vy[index / 2];
					tempV[counter].z = vz[index / 2];
				} else {
					tempV[counter].x = vx[index];
					tempV[counter].y = vy[index];
					tempV[counter].z = vz[index];
				}

				//get_Lorenz_Field(&_a[counter].p, &_a[counter].v);
				counter++;

			}
			//buffer[counter++] = lowBoundary + j*step;
		}
		//buffer[counter++] = lowBoundary + i*step;
	}

	*v = tempV;

}

void generateDataHurricane(ASF_vertex* _a, Dimension *_size, Boundary* _b,
		int primxDim, float* vx, float* vy, float* vz, Point** v) {

	Point* tempV = new Point[_size->x * _size->y * _size->z];

	Point step = (_b->high - _b->low);	 // / (*_size);
	step /= (*_size);
	//float step_x = (_b->x_high - _b->x_low) / (_size->x - 1);
	//float step_y = (_b->y_high - _b->y_low) / (_size->y - 1);
	//float step_z = (_b->z_high - _b->z_low) / (_size->z - 1);

	int counter = 0;

	for (int k = 0; k < _size->z; k++) {
		for (int j = 0; j < _size->y; j++) {

			for (int i = 0; i < _size->x; i++) {
				uint32_t index = i + j * _size->x + k * (_size->x * _size->y);

				//_a[counter].e.x =
				_a[counter].e.x = _a[counter].p.x = _b->low.x + i * step.x;

				_a[counter].e.y = _a[counter].p.y = _b->low.y + j * step.y;
				_a[counter].e.z = _a[counter].p.z = _b->low.z + k * step.z;

				_a[counter].setOldRange(index);
				_a[counter].setInBoundary();

//				if (_size->x < primxDim) {
//					tempV[counter].x = vx[index * 2];
//					tempV[counter].y = vy[index * 2];
//					tempV[counter].z = vz[index * 2];
//				} else if (_size->x > primxDim) {
//					tempV[counter].x = vx[index / 2];
//					tempV[counter].y = vy[index / 2];
//					tempV[counter].z = vz[index / 2];
//				} else {
//					tempV[counter].x = vx[index];
//					tempV[counter].y = vy[index];
//					tempV[counter].z = vz[index];
//				}

				//get_Lorenz_Field(&_a[counter].p, &_a[counter].v);
				counter++;

			}
			//buffer[counter++] = lowBoundary + j*step;
		}
		//buffer[counter++] = lowBoundary + i*step;
	}

	*v = tempV;

}

void getBenardData(float** vx, float** vy, float** vz) {
	FILE * pFile;
	int n = 3;

	int counter = 0;

	int XDIMENSION = 256;
	int YDIMENSION = 64;
	int ZDIMENSION = 128;

	int block[3];
	//FILE* pFile2 = fopen ( "blunt.x" , "rb" );
	FILE* pFile2 = fopen("bernard.raw", "rb");
	if (pFile2 == NULL) {
		fputs("File error", stderr);
		exit(1);
	}
	fseek(pFile2, 0L, SEEK_SET);

	float* m_x1 = new float[XDIMENSION * YDIMENSION * ZDIMENSION];
	float* m_y1 = new float[XDIMENSION * YDIMENSION * ZDIMENSION];
	float* m_z1 = new float[XDIMENSION * YDIMENSION * ZDIMENSION];

	unsigned char* temp_x = new unsigned char[XDIMENSION * YDIMENSION
			* ZDIMENSION];
	unsigned char* temp_y = new unsigned char[XDIMENSION * YDIMENSION
			* ZDIMENSION];
	unsigned char* temp_z = new unsigned char[XDIMENSION * YDIMENSION
			* ZDIMENSION];

	unsigned char t_x;
	unsigned char t_y;
	unsigned char t_z;
	n = 3;
	//int counter = 0;
	//int retVal=static_cast<byte>(fread(temp, sizeof(byte), XDIMENSION*YDIMENSION*ZDIMENSION, pFile2));
	for (int i = 0; i < XDIMENSION * YDIMENSION * ZDIMENSION; i++) {
		int retVal = static_cast<byte>(fread(&t_x, 1, 1, pFile2));
		retVal = static_cast<byte>(fread(&t_y, 1, 1, pFile2));
		retVal = static_cast<byte>(fread(&t_z, 1, 1, pFile2));

		double t_x1 = (double) t_x / 255 - 0.5;
		double t_y1 = (double) t_y / 255 - 0.5;
		double t_z1 = (double) t_z / 255 - 0.5;

		float dist = sqrt((double) t_x1 * t_x1 + t_y1 * t_y1 + t_z1 * t_z1);
		m_x1[i] = (double) (t_x1 / dist);	//+0.5;
		m_y1[i] = (double) (t_y1 / dist);	//+0.5;
		m_z1[i] = (double) (t_z1 / dist);	//+0.5;
		counter++;
	}

	printf("%d \n", counter);
	*vx = m_x1;
	*vy = m_y1;
	*vz = m_z1;
	fclose(pFile2);
}

void generateDataBenard(ASF_vertex* _a, Dimension *_size, Boundary* _b,
		int primxDim, float* vx, float* vy, float* vz, Point** v) {

	Point* tempV = new Point[_size->x * _size->y * _size->z];

	Point step = (_b->high - _b->low);	// / (*_size);
	step /= (*_size);
	//float step_x = (_b->x_high - _b->x_low) / (_size->x - 1);
	//float step_y = (_b->y_high - _b->y_low) / (_size->y - 1);
	//float step_z = (_b->z_high - _b->z_low) / (_size->z - 1);

	int counter = 0;

	for (int k = 0; k < _size->z; k++) {
		for (int j = 0; j < _size->y; j++) {

			for (int i = 0; i < _size->x; i++) {
				uint32_t index = i + j * _size->x + k * (_size->x * _size->y);

				//_a[counter].e.x =
				_a[counter].e.x = _a[counter].p.x = _b->low.x + i * step.x;

				_a[counter].e.y = _a[counter].p.y = _b->low.y + j * step.y;
				_a[counter].e.z = _a[counter].p.z = _b->low.z + k * step.z;

				_a[counter].setOldRange(index);
				_a[counter].setInBoundary();

				if (_size->x < primxDim) {
					tempV[counter].x = vx[index * 2];
					tempV[counter].y = vy[index * 2];
					tempV[counter].z = vz[index * 2];
				} else if (_size->x > primxDim) {
					tempV[counter].x = vx[index / 2];
					tempV[counter].y = vy[index / 2];
					tempV[counter].z = vz[index / 2];
				} else {
					tempV[counter].x = vx[index];
					tempV[counter].y = vy[index];
					tempV[counter].z = vz[index];
				}

				//get_Lorenz_Field(&_a[counter].p, &_a[counter].v);
				counter++;

			}
			//buffer[counter++] = lowBoundary + j*step;
		}
		//buffer[counter++] = lowBoundary + i*step;
	}

	printf("%d \n", counter);
	*v = tempV;

}

void generateDataTornado(ASF_vertex* _a, Dimension *_size, Boundary* _b,
		int primxDim, float* vx, float* vy, float* vz, Point** v) {

	Point* tempV = new Point[_size->x * _size->y * _size->z];
	Point step = (_b->high - _b->low);	// / (*_size);
	step /= (*_size);
	//float step_x = (_b->x_high - _b->x_low) / (_size->x - 1);
	//float step_y = (_b->y_high - _b->y_low) / (_size->y - 1);
	//float step_z = (_b->z_high - _b->z_low) / (_size->z - 1);

	int counter = 0;

	for (int k = 0; k < _size->z; k++) {
		for (int j = 0; j < _size->y; j++) {

			for (int i = 0; i < _size->x; i++) {
				uint32_t index = i + j * _size->x + k * (_size->x * _size->y);

				//_a[counter].e.x =
				_a[counter].e.x = _a[counter].p.x = _b->low.x + i * step.x;

				_a[counter].e.y = _a[counter].p.y = _b->low.y + j * step.y;
				_a[counter].e.z = _a[counter].p.z = _b->low.z + k * step.z;

				_a[counter].setOldRange(index);
				_a[counter].setInBoundary();

				tempV[counter].x = vx[counter];
				tempV[counter].y = vy[counter];
				tempV[counter].z = vz[counter];
				//get_Lorenz_Field(&_a[counter].p, &tempV[counter]);
				counter++;

			}
			//buffer[counter++] = lowBoundary + j*step;
		}
		//buffer[counter++] = lowBoundary + i*step;
	}

	*v = tempV;

}

void getTornadoData(int xs, int ys, int zs, int time, float *m_x1, float* m_y1,
		float* m_z1)

		{

	//m_x1 = new float[xs*ys*zs];
	//m_y1 = new float[xs*ys*zs];
	//m_z1 = new float[xs*ys*zs];
	float x, y, z;
	int ix, iy, iz;
	float r, xc, yc, scale, temp, z0;
	float r2 = 8;
	float SMALL = 0.00000000001;
	float xdelta = 1.0 / (xs - 1.0);
	float ydelta = 1.0 / (ys - 1.0);
	float zdelta = 1.0 / (zs - 1.0);

	//m_x1 = new float[XDIMENSION*YDIMENSION*ZDIMENSION];
	//m_y1 = new float[XDIMENSION*YDIMENSION*ZDIMENSION];
	//m_z1 = new float[XDIMENSION*YDIMENSION*ZDIMENSION];
	int counter = 0;
	for (iz = 0; iz < zs; iz++) {
		z = iz * zdelta;                        // map z to 0->1
		xc = 0.5 + 0.1 * sin(0.04 * time + 10.0 * z); // For each z-slice, determine the spiral circle.
		yc = 0.5 + 0.1 * cos(0.03 * time + 3.0 * z); //    (xc,yc) determine the center of the circle.
		r = 0.1 + 0.4 * z * z + 0.1 * z * sin(8.0 * z); //  The radius also changes at each z-slice.
		r2 = 0.2 + 0.1 * z;      //    r is the center radius, r2 is for damping
		for (iy = 0; iy < ys; iy++) {
			y = iy * ydelta;
			for (ix = 0; ix < xs; ix++) {
				x = ix * xdelta;
				temp = sqrt((y - yc) * (y - yc) + (x - xc) * (x - xc));
				scale = fabs(r - temp);
				/*
				 *  I do not like this next line. It produces a discontinuity
				 *  in the magnitude. Fix it later.
				 *
				 */
				if (scale > r2)
					scale = 0.8 - scale;
				else
					scale = 1.0;
				z0 = 0.1 * (0.1 - temp * z);
				if (z0 < 0.0)
					z0 = 0.0;
				temp = sqrt(temp * temp + z0 * z0);
				scale = (r + r2 - temp) * scale / (temp + SMALL);
				scale = scale / (1 + z);
				/*		m_x1[ix][iy][iz] = scale * (y - yc) + 0.1*(x - xc);
				 m_y1[ix][iy][iz] = scale * -(x - xc) + 0.1*(y - yc);
				 m_z[ix][iy][iz] = scale * z0;*/

				m_x1[((iz * ys + iy) * xs + ix)] = scale * (y - yc)
						+ 0.1 * (x - xc);
				m_y1[((iz * ys + iy) * xs + ix)] = scale * -(x - xc)
						+ 0.1 * (y - yc);
				m_z1[((iz * ys + iy) * xs + ix)] = scale * z0;
			}
		}
	}
}

const char* FindAndJump(const char* buffer, const char* SearchString) {
	const char* FoundLoc = strstr(buffer, SearchString);
	if (FoundLoc)
		return FoundLoc + strlen(SearchString);
	return buffer;
}

int getCylinderData(float* m_x1, float* m_y1, float* m_z1, Boundary b) {
	//const char* FileName = "testscalar.am";
	const char* FileName = "flow_t4048.am";
	//const char* FileName = "testvector3c.am";

	FILE* fp = fopen(FileName, "rb");
	if (!fp) {
		printf("Could not find %s\n", FileName);
		return 1;
	}

	float sum_x = 0.;
	float sum_y = 0.;
	float sum_z = 0.;
	float avg_x = 0.;
	float avg_y = 0.;
	float avg_z = 0.;

	printf("Reading %s\n", FileName);

	//We read the first 2k bytes into memory to parse the header.
	//The fixed buffer size looks a bit like a hack, and it is one, but it gets the job done.
	char buffer[2048];
	fread(buffer, sizeof(char), 2047, fp);
	buffer[2047] = '\0'; //The following string routines prefer null-terminated strings

	if (!strstr(buffer, "# AmiraMesh BINARY-LITTLE-ENDIAN 2.1")) {
		printf("Not a proper AmiraMesh file.\n");
		fclose(fp);
		return 1;
	}

	//Find the Lattice definition, i.e., the dimensions of the uniform grid
	int xDim(0), yDim(0), zDim(0);
	sscanf(FindAndJump(buffer, "define Lattice"), "%d %d %d", &xDim, &yDim,
			&zDim);
	printf("\tGrid Dimensions: %d %d %d\n", xDim, yDim, zDim);

	//Find the BoundingBox
	float xmin(1.0f), ymin(1.0f), zmin(1.0f);
	float xmax(-1.0f), ymax(-1.0f), zmax(-1.0f);
	sscanf(FindAndJump(buffer, "BoundingBox"), "%g %g %g %g %g %g", &xmin,
			&xmax, &ymin, &ymax, &zmin, &zmax);
	printf("\tBoundingBox in x-Direction: [%g ... %g]\n", xmin, xmax);
	printf("\tBoundingBox in y-Direction: [%g ... %g]\n", ymin, ymax);
	printf("\tBoundingBox in z-Direction: [%g ... %g]\n", zmin, zmax);

	//Is it a uniform grid? We need this only for the sanity check below.
	const bool bIsUniform = (strstr(buffer, "CoordType \"uniform\"") != NULL);
	printf("\tGridType: %s\n", bIsUniform ? "uniform" : "UNKNOWN");

	//Type of the field: scalar, vector
	int NumComponents(0);
	if (strstr(buffer, "Lattice { float Data }")) {
		//Scalar field
		NumComponents = 1;
	} else {
		//A field with more than one component, i.e., a vector field
		sscanf(FindAndJump(buffer, "Lattice { float["), "%d", &NumComponents);
	}
	printf("\tNumber of Components: %d\n", NumComponents);

	//Sanity check
	if (xDim <= 0 || yDim <= 0 || zDim <= 0 || xmin > xmax || ymin > ymax
			|| zmin > zmax || !bIsUniform || NumComponents <= 0) {
		printf("Something went wrong\n");
		fclose(fp);
		return 1;
	}

	float step_x = (b.high.x - b.low.x) / xDim;
	float step_y = (b.high.y - b.low.y) / yDim;
	float step_z = (b.high.z - b.low.z) / zDim;
	//Find the beginning of the data section
	const long idxStartData = strstr(buffer, "# Data section follows") - buffer;
	if (idxStartData > 0) {
		//Set the file pointer to the beginning of "# Data section follows"
		fseek(fp, idxStartData, SEEK_SET);
		//Consume this line, which is "# Data section follows"
		fgets(buffer, 2047, fp);
		//Consume the next line, which is "@1"
		fgets(buffer, 2047, fp);

		//Read the data
		// - how much to read
		const size_t NumToRead = xDim * yDim * zDim * NumComponents;
		// - prepare memory; use malloc() if you're using pure C
		float* pData = new float[NumToRead];
		if (pData) {
			// - do it
			const size_t ActRead = fread((void*) pData, sizeof(float),
					NumToRead, fp);
			// - ok?
			if (NumToRead != ActRead) {
				printf(
						"Something went wrong while reading the binary data section.\nPremature end of file?\n");
				delete[] pData;
				fclose(fp);
				return 1;
			}

			//Test: Print all data values
			//Note: Data runs x-fastest, i.e., the loop over the x-axis is the innermost
			printf(
					"\nPrinting all values in the same order in which they are in memory:\n");
			int Idx(0);
			for (int k = 0; k < zDim; k++) {
				for (int j = 0; j < yDim; j++) {
					for (int i = 0; i < xDim; i++) {
						//Note: Random access to the value (of the first component) of the grid point (i,j,k):
						// pData[((k * yDim + j) * xDim + i) * NumComponents]
						assert(
								pData[((k * yDim + j) * xDim + i)
										* NumComponents]
										== pData[Idx * NumComponents]);

						m_x1[i + j * xDim + k * xDim * yDim] = pData[Idx
								* NumComponents];
						m_y1[i + j * xDim + k * xDim * yDim] = pData[Idx
								* NumComponents + 1];
						m_z1[i + j * xDim + k * xDim * yDim] = pData[Idx
								* NumComponents + 2];
						sum_x += pData[Idx * NumComponents];
						sum_y += pData[Idx * NumComponents + 1];
						sum_z += pData[Idx * NumComponents + 2];

						//for(int c=0;c<NumComponents;c++)
						{
							//printf("%g ", pData[Idx * NumComponents + c]);
						}
						//	printf("\n");
						Idx++;
					}
				}
			}

			delete[] pData;
		}
	}
	avg_x = sum_x / (zDim * yDim * xDim);
	avg_y = sum_y / (zDim * yDim * xDim);
	avg_z = sum_z / (zDim * yDim * xDim);

	for (int k = 0; k < zDim; k++) {
		for (int j = 0; j < yDim; j++) {
			for (int i = 0; i < xDim; i++) {
				m_x1[((k * yDim + j) * xDim + i)] = m_x1[((k * yDim + j) * xDim
						+ i)] - avg_x;
				m_y1[((k * yDim + j) * xDim + i)] = m_y1[((k * yDim + j) * xDim
						+ i)] - avg_y;
				m_z1[((k * yDim + j) * xDim + i)] = m_z1[((k * yDim + j) * xDim
						+ i)] - avg_z;

			}
		}
	}
	// 	for (int i = 0;i<XDIMENSION*2;i++)
	// 	{
	// 		for (int j = 0;j < YDIMENSION*2;j++)
	// 		{
	// 			for (int k = 0;k < ZDIMENSION;k++)
	// 			{
	// 			}
	// 		}
	// 	}
	fclose(fp);
}

void generateDataCylinder(ASF_vertex* _a, Dimension *_size, Boundary* _b,
		int primxDim, float* vx, float* vy, float* vz, Point** v) {

	Point* tempV = new Point[_size->x * _size->y * _size->z];

	Point step = (_b->high - _b->low);	// / (*_size);
	step /= (*_size);
	//float step_x = (_b->x_high - _b->x_low) / (_size->x - 1);
	//float step_y = (_b->y_high - _b->y_low) / (_size->y - 1);
	//float step_z = (_b->z_high - _b->z_low) / (_size->z - 1);

	int counter = 0;

	for (int k = 0; k < _size->z; k++) {
		for (int j = 0; j < _size->y; j++) {

			for (int i = 0; i < _size->x; i++) {
				uint32_t index = i + j * _size->x + k * (_size->x * _size->y);

				//_a[counter].e.x =
				_a[counter].e.x = _a[counter].p.x =
						_b->low.x + i * step.x;

				_a[counter].e.y = _a[counter].p.y =
						_b->low.y + j * step.y;
				_a[counter].e.z = _a[counter].p.z =
						_b->low.z + k * step.z;

				_a[counter].setOldRange(index);
				_a[counter].setInBoundary();

				if (_size->x < primxDim) {
					tempV[counter].x = vx[index * 2];
					tempV[counter].y = vy[index * 2];
					tempV[counter].z = vz[index * 2];
				} else if (_size->x > primxDim) {
					tempV[counter].x = vx[index / 2];
					tempV[counter].y = vy[index / 2];
					tempV[counter].z = vz[index / 2];
				} else {
					tempV[counter].x = vx[index];
					tempV[counter].y = vy[index];
					tempV[counter].z = vz[index];
				}

				//get_Lorenz_Field(&_a[counter].p, &_a[counter].v);
				counter++;

			}
			//buffer[counter++] = lowBoundary + j*step;
		}
		//buffer[counter++] = lowBoundary + i*step;
	}

	*v = tempV;

}

void getOceanData(float** uvel, float** vvel, float** wvel) {
	/* This will be the netCDF ID for the file and data variable. */
	int ncid, varid;
//	float data_lat_U[NX][NY][NZ];
//	float data_lat_V[NX][NY][NZ];
//	float data_lat_W[NX][NY][NZ];
//	float data_long_U[NX][NY][NZ];
//	float data_long_V[NX][NY][NZ];
//	float data_long_W[NX][NY][NZ];
//	float data_uvel[NX][NY][NZ];
//	float data_vvel[NX][NY][NZ];
//	float data_wvel[NX][NY][NZ];

	float* lat_u = new float[NX * NY * NZ];
	float* lat_v = new float[NX * NY * NZ];
	float* lat_w = new float[NX * NY * NZ];

	float* long_u = new float[NX * NY * NZ];
	float* long_v = new float[NX * NY * NZ];
	float* long_w = new float[NX * NY * NZ];

	float* test_u = new float[NX * NY * NZ];
	float* test_v = new float[NX * NY * NZ];
	float* test_w = new float[NX * NY * NZ];
	/* Loop indexes, and error handling. */
	int x, y, z, retval;
	/* Open the file. NC_NOWRITE tells netCDF we want read-only access
	 * to the file.*/

	//reading the latitudes
	//========================================
	if ((retval = nc_open(FILE_NAME_U, NC_NOWRITE, &ncid)))
		ERR(retval);
	/* Get the varid of the data variable, based on its name. */
	if ((retval = nc_inq_varid(ncid, "LATITUDE_T", &varid)))
		printf("Error: %s\n", nc_strerror(retval));

	/* Read the data. */
	if ((retval = nc_get_var_float(ncid, varid, &lat_u[0])))	//[0][0])))
		printf("Error: %s\n", nc_strerror(retval));

//	if ((retval = nc_get_var_float(ncid, varid, &data_lat_U[0][0][0])))
//		printf("Error: %s\n", nc_strerror(retval));

	/*float spaceLat[NX*NY*NZ];
	 for (x = 0; x < NX; x++)
	 for (y = 0; y < NY; y++)
	 for (z = 1; z < NZ; z++)
	 spaceLat[x*(NY*NZ) + y*NZ + z] = data_lat[x][y][z] - data_lat[x][y][z - 1];*/
	//reading the longtitude
	//========================================
	if ((retval = nc_inq_varid(ncid, "LONGITUDE_T", &varid)))
		printf("Error: %s\n", nc_strerror(retval));

	/* Read the data. */
	if ((retval = nc_get_var_float(ncid, varid, &long_u[0])))	//[0][0])))
		printf("Error: %s\n", nc_strerror(retval));

//	if ((retval = nc_get_var_float(ncid, varid, &data_long_U[0][0][0])))
//		printf("Error: %s\n", nc_strerror(retval));

	//reading the UVel
	//========================================

	if ((retval = nc_inq_varid(ncid, "UVEL", &varid)))
		printf("Error: %s\n", nc_strerror(retval));

	/* Read the data. */
	/*if ((retval = nc_get_var_float(ncid, varid, &data_uvel[0][0][0])))
	 printf("Error: %s\n", nc_strerror(retval));*/

	if ((retval = nc_get_var_float(ncid, varid, &test_u[0])))
		printf("Error: %s\n", nc_strerror(retval));

	//===============================================================================

	if ((retval = nc_open(FILE_NAME_V, NC_NOWRITE, &ncid)))
		ERR(retval);
	/* Get the varid of the data variable, based on its name. */
	if ((retval = nc_inq_varid(ncid, "LATITUDE_T", &varid)))
		printf("Error: %s\n", nc_strerror(retval));

	/* Read the data. */
	if ((retval = nc_get_var_float(ncid, varid, &lat_v[0])))	//[0][0])))
		printf("Error: %s\n", nc_strerror(retval));

//	if ((retval = nc_get_var_float(ncid, varid, &data_lat_V[0][0][0])))
//		printf("Error: %s\n", nc_strerror(retval));

	if ((retval = nc_inq_varid(ncid, "LONGITUDE_T", &varid)))
		printf("Error: %s\n", nc_strerror(retval));

	/* Read the data. */
	if ((retval = nc_get_var_float(ncid, varid, &long_v[0])))	//[0][0])))
		printf("Error: %s\n", nc_strerror(retval));

//	if ((retval = nc_get_var_float(ncid, varid, &data_long_V[0][0][0])))
//		printf("Error: %s\n", nc_strerror(retval));

	if ((retval = nc_inq_varid(ncid, "VVEL", &varid)))
		printf("Error: %s\n", nc_strerror(retval));

	/* Read the data. */
	if ((retval = nc_get_var_float(ncid, varid, &test_v[0])))	//[0][0])))
		printf("Error: %s\n", nc_strerror(retval));

	/*if ((retval = nc_get_var_float(ncid, varid, &data_vvel[0][0][0])))
	 printf("Error: %s\n", nc_strerror(retval));*/

	//====================================================================
	if ((retval = nc_open(FILE_NAME_W, NC_NOWRITE, &ncid)))
		ERR(retval);
	/* Get the varid of the data variable, based on its name. */
	if ((retval = nc_inq_varid(ncid, "LATITUDE_T", &varid)))
		printf("Error: %s\n", nc_strerror(retval));

	/* Read the data. */
	if ((retval = nc_get_var_float(ncid, varid, &lat_w[0])))	//[0][0])))
		printf("Error: %s\n", nc_strerror(retval));

//	if ((retval = nc_get_var_float(ncid, varid, &data_lat_W[0][0][0])))
//		printf("Error: %s\n", nc_strerror(retval));

	if ((retval = nc_inq_varid(ncid, "LONGITUDE_T", &varid)))
		printf("Error: %s\n", nc_strerror(retval));

	/* Read the data. */
	if ((retval = nc_get_var_float(ncid, varid, &long_w[0])))	//[0][0])))
		printf("Error: %s\n", nc_strerror(retval));

//	if ((retval = nc_get_var_float(ncid, varid, &data_long_W[0][0][0])))
//			printf("Error: %s\n", nc_strerror(retval));

	if ((retval = nc_inq_varid(ncid, "WVEL", &varid)))
		printf("Error: %s\n", nc_strerror(retval));

	/* Read the data. */

	if ((retval = nc_get_var_float(ncid, varid, &test_w[0])))	//[0][0])))
		printf("Error: %s\n", nc_strerror(retval));

	/*if ((retval = nc_get_var_float(ncid, varid, &data_wvel[0][0][0])))
	 printf("Error: %s\n", nc_strerror(retval));*/

	//====================================================================
	/* Check the data. */
//	for (x = 0; x < NX; x++)
//	for (y = 0; y < NY; y++)
//	for (z = 0; z < NZ; z++)
//	{
//		if (test_u[x*NY*NZ + y*NZ + z] != data_uvel[x][y][z])
//			printf("");
//	}
	/*if (data_lat_U[x][y][z] != data_lat_V[x][y][z] && data_lat_V[x][y][z] != data_lat_W[x][y][z])

	 return;*/
	/* Close the file, freeing all resources. */
	if ((retval = nc_close(ncid)))
		ERR(retval);
	*uvel = test_u;	// (float*)data_uvel;
	*vvel = test_v;	// (float*)data_vvel;
	*wvel = test_w;	// (float*)data_wvel;

	//ofstream myfile("UVEL.1440x720x50.19920102.txt");
	//ofstream myfile2("VVEL.1440x720x50.19920102.txt");
	//ofstream myfile3("WVEL.1440x720x50.19920102.txt");

	//if (myfile.is_open())
	//{
	//	for (int k = 0; k < NZ; k++)
	//	for (int j = 0; j < NY;j++)
	//	for (int i = 0; i < NX; i++){
	//		double value = data_uvel[k][j][i];
	//		myfile << value<<" \n";
	//		myfile2 << data_vvel[k][j][i] << " \n";
	//		myfile3 << data_wvel[k][j][i] << " \n";
	//	}
	//	myfile.close();
	//}
	//
	//ofstream binaryFile("UVEL.1440x720x50.19920102.raw", ios::out | ios::binary);
	//binaryFile.write((char*)&data_uvel[x][y][z], sizeof (float)*NX*NY*NZ);
	//binaryFile.close();

	//ofstream binaryFile1("VVEL.1440x720x50.19920102.raw", ios::out | ios::binary);
	//binaryFile1.write((char*)&data_vvel[x][y][z], sizeof (float)*NX*NY*NZ);
	//binaryFile1.close();

	//ofstream binaryFile2("WVEL.1440x720x50.19920102.raw", ios::out | ios::binary);
	//binaryFile2.write((char*)&data_wvel[x][y][z], sizeof (float)*NX*NY*NZ);
	//binaryFile2.close();

	printf("*** SUCCESS reading example file %s!\n", FILE_NAME_U);
}

void Load_AllEdges(Edge * Fc_, uint32_t * Fr_, Edge * Bc_, uint32_t * Br_,
		uint32_t vSize, uint32_t RSize_F, uint32_t RSize_B, Edge ** oFc,
		uint32_t ** oFr, Edge ** oBc, uint32_t ** oBr, uint32_t * oRSize) {
	uint32_t Edges = RSize_F + RSize_B;
	uint32_t Vertices = vSize;
	char tmp[256];
	char tmp_c;
	uint32_t tmp_i, from, to;
	int s = sizeof(uint32_t);

	//std::vector<Edge> v_Fc;
	char* str1 = "edges.txt";
	FILE* f_p1 = fopen(str1, "wb");

	////--- the head of vtk file ---, dataset is DATASET STRUCTURED_POINTS

	uint32_t CSize = Edges;
	uint32_t RSize = Vertices + 2;

	Edge* Fc = new Edge[CSize];
	uint32_t* Fr = new uint32_t[RSize];

	int c = 0;

	for (uint32_t i = 1; i < Vertices; i++) {
		if (c == 30580)
			printf("%d", c);
		if (Fr_[i] - Fr_[i - 1] > 0) {
			//v_Fc.insert(v_Fc.begin(), Fc_[Fr_[i] - 1], Fc_[Fr_[i] ]);
			//memcpy(&Fc[c], &Fc_[Fr_[i] - 1], (Fr_[i] - Fr_[i - 1]));
			for (int j = 0; j < (Fr_[i] - Fr_[i - 1]); j++) {
				Fc[c++] = Fc_[Fr_[i] - 1 + j];
				//	fprintf(f_p1, "%d, %d\n", i, Fc_[Fr_[i] - 1 + j].value);

			}

			//c = c + (Fr_[i] - Fr_[i - 1]);
			//Fc[c] = Fc_[Fr_[i] - 1];
			//Fr[c] = Fr_[i];
			//c++;
		}

		if (Br_[i] - Br_[i - 1] > 0) {
			for (int j = 0; j < (Br_[i] - Br_[i - 1]); j++) {
				Fc[c++] = Bc_[Br_[i] - 1 + j];

				//	fprintf(f_p1, "%d, %d\n", i, Bc_[Fr_[i] - 1 + j].value);

			}
			//v_Fc.insert(v_Fc.begin(), Bc_[Br_[i] - 1], Bc_[Br_[i]]);
			/*memcpy(&Fc[c], &Bc_[Br_[i] - 1], (Br_[i] - Br_[i - 1]));
			 c = c + (Br_[i] - Br_[i - 1]);*/
			//Fc[c] = Bc_[Br_[i] - 1];
			//Fr[c] = Br_[i];
			//c++;
		}
		Fr[i] = Br_[i] + Fr_[i];
	}

	fclose(f_p1);

	//transposition
	//cout<< "Computing the transposition "<< endl;
	uint32_t* Br = new uint32_t[RSize];
	Edge* Bc = new Edge[CSize];

	uint32_t * shift = new uint32_t[RSize];

	uint32_t target_vertex = 0, source_vertex = 0;

	for (unsigned int i = 0; i < RSize; i++) {
		Br[i] = 0;
		shift[i] = 0;
	}

	for (unsigned int i = 0; i < CSize; i++) {
		Br[Fc[i].getValue() + 1]++;
	}

	for (unsigned int i = 0; i < RSize - 1; i++) {
		Br[i + 1] = Br[i] + Br[i + 1];
	}

	for (unsigned int i = 0; i < CSize; i++) {
		while (i >= Fr[target_vertex + 1]) {
			target_vertex++;
		}
		source_vertex = Fc[i].getValue();
		Bc[Br[source_vertex] + shift[source_vertex]].setValue(target_vertex);
		shift[source_vertex]++;
	}
	delete[] shift;

	//*oCSize = CSize;
	*oRSize = CSize;
	*oFc = Fc;
	*oFr = Fr;
	*oBc = Bc;
	*oBr = Br;

	//cout<<"Loading done" <<endl;
}

template<class T>
void endswap(T *objp) {
	unsigned char *memp = reinterpret_cast<unsigned char*>(objp);
	std::reverse(memp, memp + sizeof(T));
}

void getHuricaneData(float** vx, float** vy, float** vz) {

	int tz = 0;

	int sz = 500 * 500 * 100;



	float* datau = new float[sz];

	float* datav = new float[sz];
	float* dataw = new float[sz];

	FILE *file1 = fopen("Uf36.bin", "rb");
	assert(file1);

	int n;
	n = fread(datau, sizeof(float), sz, file1);
	assert(n == sz);

	bool error = false;
	char le[4];
	char be[4];

	for (int i = 0; i < sz; i++) {
		int index = i;	//250 + 500 * (250 + 500 * i);
		//printf("%f----%f \n", datau[0], datau[index]);

		memcpy(le, &datau[index], 4);

		reverse_copy(le, le + 4, be);
		memcpy(&datau[index], be, 4);

		if ((datau[index] < -79.47297 || datau[index] > 85.17703)
				&& datau[index] > 1.00000004e+35) {
			error = true;
		//	printf("%f \n ", datau[index]);
			datau[index] = 0.;
		}
		error = true;
		//printf("%f \n ", datau[index]);
//
//		float te2 = reinterpret_cast<float&>(le);
//		n = 0;
	}
	//====================================================
	FILE *file2 = fopen("Vf36.bin", "rb");
	assert(file2);

	n = fread(datav, sizeof(float), sz, file2);
	assert(n == sz);

	for (int i = 0; i < sz; i++) {
		int index = i;	//250 + 500 * (250 + 500 * i);
		//printf("%f----%f \n", datav[0], datav[index]);

		memcpy(le, &datav[index], 4);

		reverse_copy(le, le + 4, be);
		memcpy(&datav[index], be, 4);

		datav[index] = -datav[index];
		if ((datav[index] < -76.03391 || datav[index] > 82.95293)
				&& datav[index] > 1.00000004e+35) {
			error = true;
		//	printf("%f \n ", datav[index]);
			datav[index] = 0.;
		}

		//printf("%f \n ", datav[index]);
		//
		//		float te2 = reinterpret_cast<float&>(le);
		//		n = 0;
	}
	//========================================================
	FILE *file3 = fopen("Wf36.bin", "rb");
	assert(file3);

	n = fread(dataw, sizeof(float), sz, file3);
	assert(n == sz);

	for (int i = 0; i < sz; i++) {
		int index = i;		//250 + 500 * (250 + 500 * i);
		//	printf("%f----%f \n", dataw[0], dataw[index]);

		memcpy(le, &dataw[index], 4);

		reverse_copy(le, le + 4, be);
		memcpy(&dataw[index], be, 4);

		if ((dataw[index] < -9.06026 || dataw[index] > 28.61434)
				&& dataw[index] > 1.00000004e+35) {
			error = true;
		//	printf("%f \n ", dataw[index]);
			dataw[index] = 0.;
		}
		error = true;
		//	printf("%f \n ", datau[index]);
		//
		//		float te2 = reinterpret_cast<float&>(le);
		//		n = 0;
	}
	//========================================================

	*vx = datav;
	*vy = datau;
	*vz = dataw;
	return;

}

load::load() {
	// TODO Auto-generated constructor stub

}

load::~load() {
	// TODO Auto-generated destructor stub
}

