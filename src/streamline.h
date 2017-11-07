/*
 * streamline.h
 *
 *  Created on: Oct 25, 2016
 *      Author: marzieh
 */

#ifndef STREAMLINE_H_
#define STREAMLINE_H_

#include <stdint.h>

#include "streamline_kernel .h"
typedef unsigned int uint;

#define FAS(_a, _n)											\
({ __typeof__(_n) __o;										\
	__asm__ __volatile__(									\
		"lock xchg %0,%1"									\
		: "=r" (__o), "=m" (*(volatile unsigned int *)(_a))	\
		:  "0" (_n) );										\
	__o;													\
})

//const int Temp_count = 3;
//const float Time_limit = 50000.0;

//pair <uint32_t, float> FB_Decomposition(uint32_t, uint32_t, Edge *, uint32_t *, Edge *, uint32_t *, FB_vertex *, bool, bool quiet = false);
//
//pair <uint32_t, float> COL_Decomposition(uint32_t, uint32_t, Edge *, uint32_t *, Edge *, uint32_t *, COL_vertex *, ordering,
//	bool qu//iet = false);

//pair <uint32_t, float> Flow_Combinatorialization(ASF_vertex* _a, Point* v, Boundary* _b, Dimension* _d, uint32_t  _tau, Edge ** oFc, uint32_t ** oFr, uint32_t * oRSize, uint32_t whichData, bool bForward);
//pair <uint32_t, float> Flow_Combinatorialization_drawing(ASF_vertex* _a, Point* v, Boundary* _b, Dimension* _d, uint32_t  _tau, Edge ** oFc, uint32_t ** oFr, uint32_t * oRSize, uint32_t whichData, bool bForward);
//pair <uint32_t, float> Adaptive_Sampling(ASF_vertex* _a, ASF_vertex* d_a, Boundary* _b, Dimension* _d,  uint32_t * oRSize, bool bForward);
//
extern "C" {
void main_MGPU(ASF_vertex* _a, Point* _v, Boundary* _b, Dimension* _d,
		uint32_t _tau); //, ASF_vertex** o_a1, ASF_vertex** o_a2, Edge ** oFc, uint32_t ** oFr, uint32_t * oRSize, uint32_t whichData, bool bForward);

}
#endif /* STREAMLINE_H_ */
