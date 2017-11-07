#ifndef STREAMLINE_KERNEL_H
#define STREAMLINE_KERNEL_H

#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif

#define BIT_SHIFT ((unsigned)1 << 31)

#define FWD_VISITED_SHIFT ((unsigned)1 << 31)

#define BWD_VISITED_SHIFT ((unsigned)1 << 30)

#define FWD_PROPAGATE_SHIFT ((unsigned)1 << 29)

#define BWD_PROPAGATE_SHIFT ((unsigned)1 << 28)

#define RANGE_SET_SHIFT ((uint32_t)1 << 27)

#ifndef __host__
#define __host__
#endif

#ifndef __device__
#define __device__
#endif

#include <stdint.h>
#include <vector>
#include <queue>
//#include <pthread.h>
#include <assert.h>
#include <cuda_runtime.h>
#include <math.h>
#include <helper_cuda.h>

#include <helper_functions.h>
//#include <helper_timer.h>
//#include "hash_table.h"

#define epsilon exp(-4.0);
using namespace std;

typedef unsigned char byte;

struct Dimension {
	uint32_t x;
	uint32_t y;
	uint32_t z;

};
class Point {
	//private:
public:
	float x;
	float y;
	float z;

	Point() :
			x(0.0), y(0.0), z(0.0) {
	}
	Point(float& _x, float& _y, float& _z) :
			x(_x), y(_y), z(_z) {
	}
	__host__ __device__
	inline float getDist() const {
		return sqrt(x * x + y * y + z * z);
	}
	__host__ __device__
	inline void setValue(float& _x, float& _y, float& _z) {
		x = _x;
		y = _y;
		z = _z;
	}
	__host__ __device__
	inline float dist(Point& other) {

		return sqrtf(
				pow(x - other.x, 2.0) + pow(y - other.y, 2.0)
						+ pow(z - other.z, 2.0));
	}

	__host__ __device__
	inline Dimension& divide(Point _p) {
		Dimension _d;
		if (true) {
			_d.x = (uint32_t) ((x + exp(-4.0)) / _p.x);
			_d.y = (uint32_t) ((y + exp(-4.0)) / _p.y);
			_d.z = ((z + exp(-4.0)) / _p.z);
			return _d;
		}
	}
	__host__ __device__
	inline Point& operator=(const Point& other) {
		x = other.x;
		y = other.y;
		z = other.z;
		return *this;
	}
	__host__ __device__
	inline Point& operator+=(const Point& other) {
		x += other.x;
		y += other.y;
		z += other.z;
		return *this;
	}
	__host__ __device__
	inline Point& operator-=(const Point& other) {
		x -= other.x;
		y -= other.y;
		z -= other.z;
		return *this;
	}
	__host__ __device__
	inline Point& operator*=(const Point& other) {
		x *= other.x;
		y *= other.y;
		z *= other.z;
		return *this;
	}
	__host__ __device__
	inline Point& operator/=(Dimension& other) {
		x = x / (other.x - 1);
		y = y / (other.y - 1);
		z = z / (other.z - 1);
		return *this;
	}
	__host__ __device__
	inline Point& operator/(const float& other) {
		Point p;
		p.x = x / other;
		p.y = y / other;
		p.z = z / other;
		return p;
	}
	__host__ __device__
	inline Point& operator/(const Point& other) {
		x = x / other.x;
		y = y / other.y;
		z = z / other.z;
		return *this;
	}
	__host__ __device__
	inline Point& operator-(Point& other) {
		Point p;
		p.x = x - other.x;
		p.y = y - other.y;
		p.z = z - other.z;
		return p;
	}
	__host__ __device__
	inline bool operator==(const Point& other) {
		if (x == other.x && y == other.y && z == other.z)
			return true;
		return false;
	}
	__host__ __device__
	inline bool operator<(const Point& other) {
		if (x <= other.x || y <= other.y || z <= other.z)
			return true;
		return false;
	}
	__host__ __device__
	inline bool operator>(const Point& other) {
		if (x >= other.x || y >= other.y || z >= other.z)
			return true;
		return false;
	}

	//__host__ __device__ inline uint32_t getLowLink() const { return (low_link & 0x7FFFFFFF); }
	//__host__ __device__ inline void setLowLink(uint32_t n) { low_link = (n | (low_link & BIT_SHIFT)); }
	//__host__ __device__ inline bool isInComponent() const { return (low_link & BIT_SHIFT); }
	//__host__ __device__ inline void setInComponentBit() { low_link = (low_link | BIT_SHIFT); };
	//__host__ __device__ inline void clearInComponentBit() { low_link = (low_link & ~BIT_SHIFT); };
};

struct Boundary {
	Point high;
	Point low;
};

class ASF_vertex {
public:
	uint32_t range;

	uint32_t oldrange;
	Point p;


	//Point v;
	Point e;


	int type;
	uint32_t left;
	uint32_t right;
	uint32_t up;
	uint32_t down;
	//Point eb;

	//Point es;

	uint32_t level;

	//Dimension Fr_xy;	// this array will keep the indices of each face. if we suppose it counter clock wise, we have 1 2 3 4 and accordingly v Fr.x Fr.y and Fr.z
	//Dimension Fr_yz;
	//Dimension Fr_xz;

	ASF_vertex() :
			range(0), oldrange(0), p(), level(1),/*Fr_xy(),Fr_yz(),Fr_xz(),*/e() {
	}
	/*
	 *tag bits distribution in regular (non-col) obf vertex
	 *range=	|Done1	|Visited+Reached_bit	|Propagate+Elim_bit	|real_range
	 *oldrange=	|Done2	|phase-bit#1			|phase-bit#2		|real_oldrange
	 *
	 */

	//general
	__host__ __device__
	inline uint32_t getRange() const {
		return (range & 0x1FFFFFFF);
	}
	//__host__ __device__ inline uint32_t	getRangeBackward() const		{ return (brange & 0x1FFFFFFF); }

	__host__ __device__
	inline void setRange(uint32_t n) {
		range = n;
	}	//((range & 0xE0000000) | n); }
	//__host__ __device__ inline void		setRangeBackward(uint32_t n)	{ brange = ((brange & 0xE0000000) | n); }
	__host__ __device__
	inline uint32_t getOldRange() const {
		return (oldrange & 0x1FFFFFFF);
	}
	__host__ __device__
	inline void setOldRange(uint32_t n) {
		oldrange = ((oldrange & 0xE0000000) | n);
	}
	__host__ __device__
	inline uint32_t getRange(Boundary* _b, Point* _s, Dimension* _d) const {
		int i = (uint32_t) (((e.x - _b->low.x) / _s->x));
		int j = ((uint32_t) ((e.y - _b->low.y) / _s->y));
		int k = ((uint32_t) ((e.z - _b->low.z) / _s->z));
		return (i + j * _d->x + k * _d->x * _d->y);
		if (true)
			return (uint32_t) (((e.x - _b->low.x) / _s->x)
					+ ((uint32_t) ((e.y - _b->low.y) / _s->y) * _d->x)
					+ ((uint32_t) ((e.z - _b->low.z) / _s->z)) * _d->x * _d->y);
		else
			return (uint32_t) (((e.x - _b->low.x) / _s->x)
					+ (((e.y - _b->low.y) / _s->y) * _d->x)
					+ (((e.z - _b->low.z) / _s->z)) * _d->x * _d->y);
	}


	__host__ __device__
		inline uint32_t getRange_tau(Point* ep, Boundary* _b, Point* _s, Dimension* _d,int tau) const {
			int i = (uint32_t) (((ep->x - _b->low.x) / _s->x));
			int j = ((uint32_t) ((ep->y - _b->low.y) / _s->y));
			int k = ((uint32_t) ((ep->z - _b->low.z) / _s->z));
			return (i + j * _d->x + k * _d->x * _d->y);

		}



	__host__ __device__
	inline void getIndex(Boundary* _b, Point* _s, Dimension* _d, int &i, int &j,
			int&k) const {
		i = (uint32_t) (((e.x - _b->low.x) / _s->x));
		j = ((uint32_t) ((e.y - _b->low.y) / _s->y));
		k = ((uint32_t) ((e.z - _b->low.z) / _s->z));
		//else return(uint32_t)(((e.x - _b->low.x) / _s->x) + (((e.y - _b->low.y) / _s->y)*_d->x) + (((e.z - _b->low.z) / _s->z))*_d->x*_d->y);
	}
	__host__ __device__
	inline uint32_t getOldRange(Boundary* _b, Point* _s, Dimension* _d) const {
		return (uint32_t) (((p.x - _b->low.x) / _s->x)
				+ ((uint32_t) ((p.y - _b->low.y) / _s->y) * _d->x)
				+ ((uint32_t) ((p.z - _b->low.z) / _s->z)) * _d->x * _d->y);
	}

	/*__host__ __device__ inline uint32_t getRangeBackward(Boundary* _b, Point* _s, Dimension* _d) const { return (uint32_t)(((uint32_t)((eb.x - _b->low.x) / _s->x)) + ((uint32_t)((eb.y - _b->low.y) / _s->y)*_d->x) + ((uint32_t)((eb.z - _b->low.z) / _s->z)*_d->x*_d->y)); }*/
	__host__ __device__
	inline float getDistance(Point _s) const {
		return (float) sqrt(
				pow((_s.x - e.x), 2.0) + pow((_s.y - e.y), 2.0)
						+ pow((_s.z - e.z), 2.0));
	}

	//__host__ __device__ inline bool checkFaceRange(Dimension d){ if (abs((int)Fr.x - (int)getRange()) <= 1 || abs((int)Fr.x - (int)getRange()) <= d.x || abs((int)Fr.x - (int)getRange()) <= (d.x*d.y))return true; return false; }
	//__host__ __device__ inline void setFaceInd_X(Dimension d){Fr.x = getRange() + 1; Fr.y = getRange() + d.x + 1; Fr.z = getRange() + d.x;	}
	__host__ __device__
	inline bool isInFWD() const {
		return (!(oldrange & 0x40000000) && !(oldrange & 0x20000000));
	}
	__host__ __device__
	inline bool isInSCC() const {
		return (range & 0x80000000);
	}
	__host__ __device__
	inline void setInFWD() {
		oldrange = (oldrange & 0x1FFFFFFF);
		range = (range & 0x1FFFFFFF);
	}
	__host__ __device__
	inline void setInSCC() {
		range = ((range & 0x1FFFFFFF) | 0x80000000);
	}
	__host__ __device__
	inline void setInNextLevel_xy() {
		range = (/*(range & 0x1FFFFFFF)*/range | 0x20000000);
	}
	__host__ __device__
	inline void setInNextLevel_yz() {
		range = (/*(range & 0x1FFFFFFF)*/range | 0x40000000);
	}
	__host__ __device__
	inline void setInNextLevel_xz() {
		range = (/*(range & 0x1FFFFFFF)*/range | 0x80000000);
	}
	__host__ __device__
	inline bool isInOWCTY() const {
		return (!(oldrange & 0x40000000) && (oldrange & 0x20000000));
	}
	__host__ __device__
	inline void setInOWCTY() {
		oldrange = ((oldrange & 0x1FFFFFFF) | 0x20000000);
		range = (range & 0x1FFFFFFF);
	}
	__host__ __device__
	inline bool isInBWD() const {
		return ((oldrange & 0x40000000) && !(oldrange & 0x20000000));
	}
	__host__ __device__
	inline void setInBWD() {
		oldrange = ((oldrange & 0x1FFFFFFF) | 0x40000000);
		range = (range & 0x1FFFFFFF);
	}
	__host__ __device__
	inline bool isInBoundary() const {
		return (oldrange & 0x40000000) && (oldrange & 0x20000000);
	}

	__host__ __device__
	inline bool isInNextLevel_xy() const {
		return (range & 0x20000000);
	}
	__host__ __device__
	inline bool isInNextLevel_yz() const {
		return (range & 0x40000000);
	}
	__host__ __device__
	inline bool isInNextLevel_xz() const {
		return (range & 0x80000000);
	}
	__host__ __device__
	inline void setInBoundary() {
		oldrange = ((oldrange & 0x1FFFFFFF) | 0x60000000);
	}
	__host__ __device__
	inline void unsetInBoundary() {
		oldrange = ((oldrange & 0x1FFFFFFF));
	}
	__host__ __device__
	inline void unsetInFace_xy() {
		range = ((range & 0xDFFFFFFF));
	}
	__host__ __device__
	inline void unsetInFace_yz() {
		range = ((range & 0xBFFFFFFF));
	}
	__host__ __device__
	inline void unsetInFace_xz() {
		range = ((range & 0x7FFFFFFF));
	}
	__host__ __device__
	inline bool checkInBoundary(Boundary* b) {
		if (e < b->low || e > b->high)
			return false;
		return true;
	}

	__host__ __device__
	inline bool checkInBoundary_tau(Point ep, Boundary* b,int tau) {
		if (ep < b->low || ep > b->high)
			return false;
		return true;
	}




	__host__ __device__
	inline bool checkInBoundary_StartPoint(Boundary* b) {
		if (p < b->low || p > b->high)
			return false;
		return true;
	}
	//__host__ __device__ inline bool checkInBoundaryBackward(Boundary* b)			{ if (eb < b->low || eb > b->high)return false; return true; }
	__host__ __device__
	inline bool isInSCC1() const {
		return ((range & 0x40000000) && (range & 0x20000000));
	}
	__host__ __device__
	inline uint32_t getPhase() const {
		return ((oldrange >> 29) & 0x3);
	}
	__host__ __device__
	inline uint32_t getLevel() const {
		return ((oldrange >> 29) & 0x8);
	}

	//pivot only voting bits
	__host__ __device__
	inline bool isDone1() const {
		return (range & 0x80000000);
	}
	__host__ __device__
	inline void setDone1() {
		range = (range | 0x80000000);
	}
	__host__ __device__
	inline void unsetDone1() {
		range = (range & 0x7FFFFFFF);
	}
	__host__ __device__
	inline bool isDone2() const {
		return (oldrange & 0x80000000);
	}
	__host__ __device__
	inline void setDone2() {
		oldrange = (oldrange | 0x80000000);
	}
	__host__ __device__
	inline void unsetDone2() {
		oldrange = (oldrange & 0x7FFFFFFF);
	}

	//phase specific
	__host__ __device__
	inline bool isFWDVisited() const {
		return (range & 0x40000000);
	}
	__host__ __device__
	inline void setFWDVisited() {
		range = (range | 0x40000000);
	}
	__host__ __device__
	inline void setnFace(uint32_t t) {
		range = (range | t);
	}
	__host__ __device__
	inline bool isFWDPropagate() const {
		return (range & 0x20000000);
	}
	__host__ __device__
	inline void setFWDPropagate() {
		range = (range | 0x20000000);
	}

	__host__ __device__
	inline bool isBWDVisited() const {
		return (range & 0x40000000);
	}
	__host__ __device__
	inline void setBWDVisited() {
		range = (range | 0x40000000);
	}
	__host__ __device__
	inline bool isBWDPropagate() const {
		return (range & 0x20000000);
	}
	__host__ __device__
	inline void setBWDPropagate() {
		range = (range | 0x20000000);
	}

	__host__ __device__
	inline bool isReached() const {
		return (range & 0x40000000);
	}
	__host__ __device__
	inline void setReached() {
		range = (range | 0x40000000);
	}
	__host__ __device__
	inline bool isElim() const {
		return (range & 0x20000000);
	}
	__host__ __device__
	inline void setElim() {
		range = (range | 0x20000000);
	}

	//COL
	__host__ __device__
	inline void setInCOL() {
		oldrange = (oldrange | 0xE0000000);
		range = 0x0;
	}
	__host__ __device__
	inline bool isInCOL() const {
		return (oldrange & 0x80000000);
	}
};

class fEdge {
public:
	uint32_t v1;
	uint32_t v2;
	int subedge[2];
	bool bsplit;
	int E2F[4];
	uint32_t E2V;
	uint32_t next;
	uint32_t Prev;
	uint32_t level;
	bool internalEdge;
	bool bconnected;
	int parent;

	//float length;
	__host__ __device__
	inline bool isInBoundary() const {
		return (E2V & 0x20000000);
	}
	__host__ __device__
	inline void setInBoundary() {
		E2V = (E2V | 0x20000000);
	}
	__host__ __device__
	inline void unsetInBoundary() {
		E2V = ((E2V & 0x1FFFFFFF));
	}
	//bool bBoundary;
	//uint32_t E2T[2];
};

class fTriangle {
public:
	uint32_t edge[3];
	uint32_t T2F;

};

class fFace {
public:
	uint32_t edge[4];
	int subface[4];
	uint32_t vertexes[4];
	uint32_t F2V;
	uint32_t parent;
	uint32_t level;
	byte dir;
	bool bsplit;
	bool bconnected;

	uint32_t center;
	fFace() :
			level(0), bsplit(0) {
	}
	//int state;
	//float dx1;
	//float dx2;
	//bool bSplit;
	float Jacobian;__host__ __device__
	inline bool isInBoundary() const {
		return (F2V & 0x20000000);
	}
	__host__ __device__
	inline void setInBoundary() {
		F2V = (F2V | 0x20000000);
	}
	__host__ __device__
	inline void unsetInBoundary() {
		F2V = ((F2V & 0x1FFFFFFF));
	}
	__host__ __device__
	inline uint32_t getFaceVertex() const {
		return (F2V & 0x1FFFFFFF);
	}
};

class fVoxel {
public:
	int id;
	fFace farray[6];
};

//********************************************
#endif // SCC_KERNEL_H
