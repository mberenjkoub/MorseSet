/*
 * scc_kernel.h
 *
 *  Created on: Oct 25, 2016
 *      Author: marzieh
 */




#ifndef SCC_KERNEL_H
#define SCC_KERNEL_H

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

#include <helper_cuda.h>

#include <helper_functions.h>
#include <helper_timer.h>


using namespace std;

const uint32_t MAX_VERTEX = 0x7FFFFFFF; // 2^31 - 1
const uint32_t blockSize = 256;
const uint32_t my_cl = 64;

//!32-bit unsigned integer
 typedef uint32_t ulong_int_t;

enum ordering {MIN, MAX};
enum dataType{Lorenz, Benard, Cylinder, Ocean,Tornado,Hurricane };


class Edge
{
//private:
public:
    uint32_t value;

    Edge() : value(0) {}
	Edge( uint32_t & val ) : value(val) {}
    __host__ __device__ inline uint32_t getValue() const { return (value & 0x7FFFFFFF); }
    __host__ __device__ inline void setValue(uint32_t n) { value = (n | (value & BIT_SHIFT)); }
    __host__ __device__ inline bool isValid() const { return (value & BIT_SHIFT); }
    __host__ __device__ inline void setValidBit() { value = (value | BIT_SHIFT); };
    __host__ __device__ inline void clearValidBit() { value = (value & ~BIT_SHIFT); };
};

class FB_vertex
{
//private:
public:
    uint32_t range;
    uint32_t tags;

    FB_vertex() : range(0), tags(0) {}
    __host__ __device__ inline uint32_t getRange() const { return  range; }
    __host__ __device__ inline void setRange(uint32_t n) { range  = n ; }
    __host__ __device__ inline bool isForwardVisited() const { return ( tags & FWD_VISITED_SHIFT); }
    __host__ __device__ inline bool isForwardPropagate() const { return (tags & FWD_PROPAGATE_SHIFT); }
    __host__ __device__ inline bool isBackwardVisited() const { return (tags & BWD_VISITED_SHIFT); }
    __host__ __device__ inline bool isBackwardPropagate() const { return ( tags & BWD_PROPAGATE_SHIFT); }
    __host__ __device__ inline void setForwardVisitedBit() { tags = ( tags | FWD_VISITED_SHIFT); };
    __host__ __device__ inline void setForwardPropagateBit() { tags = ( tags | FWD_PROPAGATE_SHIFT); };
    __host__ __device__ inline void setBackwardVisitedBit() { tags = ( tags | BWD_VISITED_SHIFT); };
    __host__ __device__ inline void setBackwardPropagateBit() { tags = ( tags | BWD_PROPAGATE_SHIFT); };
    __host__ __device__ inline void clearForwardVisitedBit() { tags = (tags & ~FWD_VISITED_SHIFT); };
    __host__ __device__ inline void clearForwardPropagateBit() { tags = (tags & ~FWD_PROPAGATE_SHIFT); };
    __host__ __device__ inline void clearBackwardVisitedBit() { tags = (tags & ~BWD_VISITED_SHIFT); };
    __host__ __device__ inline void clearBackwardPropagateBit() { tags = (tags & ~BWD_PROPAGATE_SHIFT); };
    __host__ __device__ inline bool isRangeSet() const { return ( tags & RANGE_SET_SHIFT); }
    __host__ __device__ inline void rangeSet() { tags = ( tags | RANGE_SET_SHIFT); };
    __host__ __device__ inline uint32_t getSection() const { return ((bool(tags & FWD_VISITED_SHIFT) << 1) + bool(tags & BWD_VISITED_SHIFT));}
};

class OBF_vertex {
public:
	uint32_t range;
	uint32_t oldrange;
//	uint32_t Fr[3];

	OBF_vertex() : range(0), oldrange(0) {}
/*
 *tag bits distribution in regular (non-col) obf vertex
 *range=	|Done1	|Visited+Reached_bit	|Propagate+Elim_bit	|real_range
 *oldrange=	|Done2	|phase-bit#1			|phase-bit#2		|real_oldrange
 *
 */

//general
	__host__ __device__ inline uint32_t	getRange() const		{ return ( range & 0x1FFFFFFF ); }
	__host__ __device__ inline void		setRange(uint32_t n)	{ range = ( ( range & 0xE0000000 ) | n ); }
	__host__ __device__ inline uint32_t	getOldRange() const		{ return ( oldrange & 0x1FFFFFFF ); }
	__host__ __device__ inline void		setOldRange(uint32_t n)	{ oldrange = ( ( oldrange & 0xE0000000 ) | n ); }

	__host__ __device__ inline bool isInFWD() const		{ return ( !( oldrange & 0x40000000 ) && !( oldrange & 0x20000000 ) ); }
	__host__ __device__ inline void setInFWD()			{ oldrange = ( oldrange & 0x1FFFFFFF ); range = ( range & 0x1FFFFFFF ); }
	__host__ __device__ inline bool isInOWCTY() const	{ return ( !( oldrange & 0x40000000 ) && ( oldrange & 0x20000000 ) ); }
	__host__ __device__ inline void setInOWCTY()		{ oldrange = ( ( oldrange & 0x1FFFFFFF ) | 0x20000000 ); range = ( range & 0x1FFFFFFF ); }
	__host__ __device__ inline bool isInBWD() const		{ return ( ( oldrange & 0x40000000 ) && !( oldrange & 0x20000000 ) ); }
	__host__ __device__ inline void setInBWD()			{ oldrange = ( ( oldrange & 0x1FFFFFFF ) | 0x40000000 ); range = ( range & 0x1FFFFFFF ); }
	__host__ __device__ inline bool isInSCC() const		{ return ((oldrange & 0x40000000) && (oldrange & 0x20000000)); }
	__host__ __device__ inline bool isInSCC1() const		{ return ((range & 0x40000000) && (range & 0x20000000)); }
	__host__ __device__ inline void setInSCC()			{ oldrange = ( ( oldrange & 0x1FFFFFFF ) | 0x60000000 ); }

	__host__ __device__ inline uint32_t getPhase() const { return ( ( oldrange >> 29 ) & 0x3 ); }

//pivot only voting bits
	__host__ __device__ inline bool isDone1() const	{ return ( range & 0x80000000 ); }
	__host__ __device__ inline void setDone1()		{ range = ( range | 0x80000000 ); }
	__host__ __device__ inline void unsetDone1()	{ range = ( range & 0x7FFFFFFF ); }
	__host__ __device__ inline bool isDone2() const	{ return ( oldrange & 0x80000000 ); }
	__host__ __device__ inline void setDone2()		{ oldrange = ( oldrange | 0x80000000 ); }
	__host__ __device__ inline void unsetDone2()	{ oldrange = ( oldrange & 0x7FFFFFFF ); }

//phase specific
	__host__ __device__ inline bool isFWDVisited() const	{ return ( range & 0x40000000 ); }
	__host__ __device__ inline void setFWDVisited()			{ range = ( range | 0x40000000 ); }
	__host__ __device__ inline bool isFWDPropagate() const	{ return ( range & 0x20000000 ); }
	__host__ __device__ inline void setFWDPropagate()		{ range = ( range | 0x20000000 ); }

	__host__ __device__ inline bool isBWDVisited() const	{ return ( range & 0x40000000 ); }
	__host__ __device__ inline void setBWDVisited()			{ range = ( range | 0x40000000 ); }
	__host__ __device__ inline bool isBWDPropagate() const	{ return ( range & 0x20000000 ); }
	__host__ __device__ inline void setBWDPropagate()		{ range = ( range | 0x20000000 ); }

	__host__ __device__ inline bool isReached() const	{ return ( range & 0x40000000 ); }
	__host__ __device__ inline void setReached()		{ range = ( range | 0x40000000 ); }
	__host__ __device__ inline bool isElim() const		{ return ( range & 0x20000000 ); }
	__host__ __device__ inline void setElim()			{ range = ( range | 0x20000000 ); }

//COL
	__host__ __device__ inline void setInCOL()			{ oldrange = ( oldrange | 0xE0000000 ); range = 0x0; }
	__host__ __device__ inline bool isInCOL() const		{ return ( oldrange & 0x80000000 ); }
};

class COL_vertex
{
//private:
public:
	uint32_t map;
	uint32_t oldmap;

	COL_vertex() : map(0), oldmap(0) {}

//COL usage in OBF
	__host__ __device__ inline COL_vertex & operator=(const OBF_vertex & obf)
	{
		if ( obf.isInCOL() ) {
			map = 0;
			oldmap = MAX_VERTEX - obf.getOldRange();
		}
		else {
			map = MAX_VERTEX - obf.getOldRange();
			//oldmap = ( obf.isInSCC() ) ? 0x80000000 : 0;//TODO-test likely unnecessary: vertex not in COL means it's in SCC
			oldmap = 0x80000000;
		}
		return *this;
	}

	__host__ __device__ inline uint32_t getMap() const { return (map & 0x7FFFFFFF); }
	__host__ __device__ inline uint32_t getOldMap() const { return (oldmap & 0x7FFFFFFF); }
	__host__ __device__ inline void setMap(uint32_t n) { map = (n | (map & BIT_SHIFT)); }
	__host__ __device__ inline void setOldMap(uint32_t n) { oldmap = (n | (oldmap & BIT_SHIFT)); }
	__host__ __device__ inline bool isBackwardVisited() const { return (map & BIT_SHIFT); }
	__host__ __device__ inline bool isPropagate() const { return (oldmap & BIT_SHIFT);}
	__host__ __device__ inline void setBackwardVisitedBit() { map = (map | BIT_SHIFT); };
	__host__ __device__ inline void setPropagateBit() { oldmap = (oldmap | BIT_SHIFT); };
	__host__ __device__ inline void clearBackwardVisitedBit() { map = (map & ~BIT_SHIFT); };
	__host__ __device__ inline void clearPropagateBit() { oldmap = (oldmap & ~BIT_SHIFT); };
};

class TR_vertex
{
//private:
public:
	uint32_t visit;
	uint32_t low_link;

	TR_vertex() : visit(0), low_link(0) {}
	__host__ __device__ inline uint32_t getVisited() const { return visit; }
	__host__ __device__ inline void setVisited(uint32_t n) { visit = n; }

	__host__ __device__ inline uint32_t getLowLink() const { return (low_link & 0x7FFFFFFF); }
	__host__ __device__ inline void setLowLink(uint32_t n) { low_link = (n | (low_link & BIT_SHIFT)); }
	__host__ __device__ inline bool isInComponent () const { return (low_link & BIT_SHIFT); }
	__host__ __device__ inline void setInComponentBit() { low_link = (low_link | BIT_SHIFT); };
	__host__ __device__ inline void clearInComponentBit() { low_link = (low_link & ~BIT_SHIFT); };
};

class CPU_OBF_vertex {
//private:
public:
    uint32_t visited;  //visited in the correponding job
    uint32_t inset;    // membership to the corresponding job

    CPU_OBF_vertex() : visited(0), inset(0) {};

    __host__ __device__ inline uint32_t getVisited() const { return (visited & 0x7FFFFFFF); }
    __host__ __device__ inline void setVisited(uint32_t n) { visited = (n | (visited & BIT_SHIFT)); }
    __host__ __device__ inline bool isSeed () const { return (visited & BIT_SHIFT); }
    __host__ __device__ inline void setSeedBit() { visited = (visited | BIT_SHIFT); };
    __host__ __device__ inline void clearSeedBit() { visited = (visited & ~BIT_SHIFT); };

    __host__ __device__ inline uint32_t getInset() const { return (inset & 0x7FFFFFFF); }
    __host__ __device__ inline void setInset(uint32_t n) { inset = (n | (inset & BIT_SHIFT)); }
    __host__ __device__ inline bool isEliminated () const { return (inset & BIT_SHIFT); }
    __host__ __device__ inline void setEliminatedBit() { inset = (inset | BIT_SHIFT); };
    __host__ __device__ inline void cleareliminatedBit() { inset = (inset & ~BIT_SHIFT); };

};

class BFS_vertex {
//private:
public:
    uint32_t distance;

    BFS_vertex() : distance(0) {}

    __host__ __device__ inline uint32_t getDistance() const { return (distance & 0x7FFFFFFF); }
    __host__ __device__ inline void setDistance(uint32_t n) { distance = (n | (distance & BIT_SHIFT)); }
    __host__ __device__ inline bool isVisited () const { return (distance & BIT_SHIFT); }
    __host__ __device__ inline void setVisitedBit() { distance = (distance | BIT_SHIFT); };
    __host__ __device__ inline void clearVisitedBit() { distance = (distance & ~BIT_SHIFT); };
};

class TR_stack_vertex {
//private:
public:
    uint32_t id;
    uint32_t from;

    TR_stack_vertex() : id(0), from(0) {}
    __host__ __device__ inline uint32_t getId() const { return (id & 0x7FFFFFFF); }
    __host__ __device__ inline void setId(uint32_t n) { id = (n | (id & BIT_SHIFT)); }
    __host__ __device__ inline bool isExpanded () const { return (id & BIT_SHIFT); }
    __host__ __device__ inline void setExpandedBit() { id = (id | BIT_SHIFT); };
    __host__ __device__ inline void clearExpandedBit() { id = (id & ~BIT_SHIFT); };

    __host__ __device__ inline uint32_t getFrom() const { return from; }
    __host__ __device__ inline void setFrom(uint32_t n) { from = n ; }

};

//********** for CPU parallel FWD, OBF and parallel generation ****************
struct init_t
{
  uint32_t thread_id;
  uint32_t number_of_threads;
  Edge *Fc;
  uint32_t *Fr;
  uint32_t size;
};


struct OBF_init_t
{
  uint32_t thread_id;
  uint32_t number_of_threads;
  Edge *Fc;
  uint32_t *Fr;
  Edge *Bc;
  uint32_t *Br;
  uint32_t size;
};

struct OBF_info_t
{
 uint32_t inset; //  set of the state in which we performe  the reachibility  (0 for first job we work on all state)
 uint32_t states; // number of visited states in corresponding job

 unsigned short int tag;
 queue<uint32_t>* reached; // corresponding to set Reached in the paper
 uint32_t range; // corresponding to set Range in the paper
 uint32_t executed_by;  // the forward reachibility which executed this job
 uint32_t executed_by_BWD; // the backward reachibility which executed this job
 bool has_root; // has the root ?
 uint32_t root; // root for forward reachibility necessery for next owcty
 queue<uint32_t>* seeds; // corresponding to set Seeds returned by FWDSEED in the paper
 //queue<uint32_t>* SCC_list;
};

struct g_done
{
	uint32_t value1;
	char padding1[ my_cl - sizeof(uint32_t) ];
	bool value2;
	char padding2[ my_cl - sizeof(bool) ];
	bool value3;
	char padding3[ my_cl - sizeof(bool) ];
	bool value4;
	char padding4[ my_cl - sizeof(bool) ];
};

class CFB_vertex {
//private:
public:
/*
 *tag bits distribution in cfb vertex
 *range= | FWDVisited | BWDVisited | Done | real_range
 *
 */
	uint32_t range;

	CFB_vertex() : range(1) {}

	__host__ __device__ inline uint32_t	getRange() const		{ return ( range & 0x1FFFFFFF ); }
	__host__ __device__ inline void		setRange(uint32_t n)	{ range = n; }
	__host__ __device__ inline bool		isDone() const			{ return ( range & 0x20000000 ); }
	__host__ __device__ inline void		setDone()				{ range = ( range | 0x20000000 ); }
    __host__ __device__ inline bool		isFWDVisited() const	{ return ( range & 0x80000000 ); }
    __host__ __device__ inline void		setFWDVisited()			{ range = ( range | 0x80000000 ); }
    __host__ __device__ inline bool		isBWDVisited() const	{ return ( range & 0x40000000 ); }
    __host__ __device__ inline void		setBWDVisited()			{ range = ( range | 0x40000000 ); }
    __host__ __device__ inline bool		isVisited() const		{ return ( !( ( range | 0xC0000000 ) ^ range ) ); }
    __host__ __device__ inline void		setVisited()			{ range = ( range | 0xC0000000 ); }
};

struct _BFS
{
	int thread_id;
	int number_of_threads;
	uint32_t * my_Fr;
	Edge * my_Fc;
	uint32_t * my_Br;
	Edge * my_Bc;
	CFB_vertex * my_m;
	uint32_t max_index;
	uint32_t one_vertices;
	g_done * done;

	uint32_t local_count;
	char padding[ my_cl - sizeof(uint32_t) ];
};

//********************************************
#endif // SCC_KERNEL_H


