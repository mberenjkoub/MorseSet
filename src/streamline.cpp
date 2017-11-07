/*
  cuda-scc: CUDA SCC decomposition

  Author: Milan Ceska, Faculty of Informatics, Masaryk University, Czech Republic
  Mail: xceska@fi.muni.cz
*/

#include "SCC.h"
//#include "scc_kernel.h"
#include "load.h"

#include <stdint.h>
//#include <pthread.h>
#include <set>
#include <map>
#include <string>
#include <queue>
#include <cstdio>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;


void print_help()
{
	printf("========================================================\n");
	printf("CUDA SCC decomposition Tool, usage: \nScc Alg [-w=Num] [Ord] [Trimm] [COL] [COL-limit] [NumT] File \n");
	printf("Input ::= d | g\n");
	printf("Alg ::= t | f | c | o | a | s | p | r | u\n");
	printf("Ord ::= 0 | 1\n");
	printf("Trimm ::= 0 | 1\n");
	printf("Num ::= 1 | ... | 16\n");
	printf("COL ::= 0 | 1\n");
	printf("COL-limit ::= 0 | ...\n");
	printf("--------------------------------------------------------\n");
	printf("Alg: decomposing algorithm\n");
	printf("f - FB algorithm, c - Colouring/Heads-Off algorithm, o - OBF algorithm\n");
	printf("t - Tarjan's algorithm:CPU, u - OBF:CPU\n");
	printf("s - BFS:GPU, p - BFS:CPU(serial), r - BFS:CPU(parallel)\n");
	printf("a - test suite: all available algorithms with most sensible settings\n");
	printf("General parameters:\n");
	printf("Num - number of threads used to generate graph.\n");
	printf("Parameters for FB:\n");
	printf("Trimm - use OWCTY elimination between iterations of FB.\n");
	printf("Parameters for Colouring:\n");
	printf("Ord - use maximal/minimal accepting predecessor.\n");
	printf("Parameters for OBF:\n");
	printf("Trimm - use OWCTY before the actual OBF algorithm,\n");
	printf("COL - use Colouring on small enough SCC-closed subgraphs,\n");
	printf("COL-limit - how small the subgraphs must be.\n");
	printf("Parameters for BFS, OBF:CPU:\n");
	printf("NumT - number of threads used for parallel reachability.\n");
	printf("========================================================\n");
}

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////

//
//int main(int argc, char** argv)
//{
//	//if ( argc < 3 ) {
//	//print_help();
//	//	return 1;
//	//}
//
//	// CSR representation
//	uint32_t CSize; // column arrays size
//	uint32_t RSize; // range arrays size
//	// Forwards arrays
//	Edge * Fc; // forward columns
//	uint32_t * Fr = NULL; // forward ranges
//	// Backwards arrays
//	Edge * Bc; // backward columns
//	uint32_t * Br; // backward ranges
//	GraphGenerator * gen = NULL;
//	_DeviceSet = false;
//	// obtain a CSR graph representation
//	bool dve = true;
//	/*if ( argv[ argc - 1 ][ strlen( argv[ argc - 1 ] ) - 1 ] == 'r' )
//	dve = false;
//	if ( dve ) {
//	int num = 1;
//	for ( int i = 2; i < argc - 1; i++ ) {
//	if ( argv[ i ][ 0 ] == '-' )
//	num = atoi( argv[ i ] + 3 );
//	}
//	gen = new GraphGenerator( argv[ argc - 1 ], num );
//	if ( num > 1 )
//	gen->p_generateGraph( 30000000 );
//	else
//	gen->generateGraph( 30000000 );
//
//	gen->getGraph( &CSize, &RSize, &Fc, &Fr, &Bc, &Br );
//	}
//	else */
//	{
//		//loadFullGraph("scc_all\\Vertices_Forward.txt", &CSize, &RSize, &Fc, &Fr, &Bc, &Br);
//		loadFullGraph("scc_all\\Vertices_Forward_Benard_3000.txt", &CSize, &RSize, &Fc, &Fr, &Bc, &Br);
//	}
//
//	char c = 'o';// argv[1][0];
//	try {
//		switch (c) {
//			/*	case 'f' :
//			{
//			FB_vertex * _m = new FB_vertex[ RSize - 1 ];
//			pair <uint32_t, float> result = FB_Decomposition( CSize, RSize, Fc, Fr, Bc, Br, _m, ( argc > 3 ) ? bool( atoi( argv[ 2 ] )) : false, (argv[ 1 ][ 1 ] != '\0') );
//			if ( argv[ 1 ][ 1 ] != '\0' )
//			printf("%f", result.second);
//			delete [] _m;
//			}
//			break;
//
//			case 'c' :
//			{
//			COL_vertex * _m = new COL_vertex[ RSize - 1 ];
//			pair <uint32_t, float> result = COL_Decomposition( CSize, RSize, Fc, Fr, Bc, Br, _m, ( argc > 3 && !(atoi( argv[ 2 ] )) ) ? MAX : MIN, (argv[ 1 ][ 1 ] != '\0') );
//			if ( argv[ 1 ][ 1 ] != '\0' )
//			printf("%f", result.second);
//			delete [] _m;
//			}
//			break; */
//
//		case 'o':
//		{
//					OBF_vertex * _m = new OBF_vertex[RSize - 1];
//					pair <uint32_t, float> result = OBF_Decomposition(CSize, RSize, Fc, Fr, Bc, Br, _m, /*(argc > 4) ? bool(atoi(argv[3])) :*/ false,
//						/*(argc > 5) ? atoi(argv[4]) :*/ 1000, /*(argc > 3) ? bool(atoi(argv[2])) :*/ true, true/*(argv[1][1] != '\0')*/);
//					//if (argv[1][1] != '\0')
//						printf("%f", result.second);
//					delete[] _m;
//		}
//			break;
//			/*
//			case 't' :
//			{
//			TR_vertex * _m = new TR_vertex[ RSize - 1 ];
//			pair <uint32_t, float> result = TR_Decomposition( CSize, RSize, Fc, Fr, _m, (argv[ 1 ][ 1 ] != '\0') );
//			if ( argv[ 1 ][ 1 ] != '\0' )
//			printf("%f", result.second);
//			delete [] _m;
//			}
//			break;
//
//			case 's' :
//			{
//			FB_vertex * _m = new FB_vertex[ RSize - 1 ];
//			Forward_only( CSize, RSize, Fc, Fr, _m, ( argc > 3 ) ? atoi( argv[ 2 ] ) : 0 );
//			delete [] _m;
//			}
//			break;
//
//			case 'p' :
//			{
//			BFS_vertex * _m = new BFS_vertex[ RSize - 1 ];
//			BFS( CSize, RSize, Fc, Fr, _m );
//			delete [] _m;
//			}
//			break;
//
//			case 'r' :
//			{
//			parallel_FWD( CSize, RSize, Fc, Fr, RSize - 1, ( argc > 3 ) ? atoi( argv[ 2 ] ) : 4 );
//			}
//			break;
//
//			case 'z' :
//			_CPU_BFS( CSize, RSize, Fc, Fr, RSize - 1, ( argc > 3 ) ? atoi( argv[ 2 ] ) : 4 );
//			break;
//
//			case 'y' :
//			{
//			CFB_vertex * _m = new CFB_vertex[ RSize - 1 ];
//			_CPU_FB( CSize, RSize, Fc, Fr, Bc, Br, _m, RSize - 1, ( argc > 3 ) ? atoi( argv[ 2 ] ) : 4 );
//			delete [] _m;
//			break;
//			}
//			case 'a' :
//			{
//			Edge * _Bc = new Edge[ CSize ];
//			uint32_t * _Br = new uint32_t[ RSize ];
//			Edge * _Fc = new Edge[ CSize ];
//			uint32_t * _Fr = new uint32_t[ RSize ];
//
//			printf("SCC decomposition of graph %s ", argv[ argc - 1 ]);
//			printf("of %u Vertices, %u Edges and ", RSize - 2, CSize);
//			pair <uint32_t, float> result;
//
//			memcpy( _Fc, Fc, sizeof(Edge)*CSize );
//			memcpy( _Fr, Fr, sizeof(uint32_t)*RSize );
//			memcpy( _Bc, Bc, sizeof(Edge)*CSize );
//			memcpy( _Br, Br, sizeof(uint32_t)*RSize );
//			TR_vertex * _tm = new TR_vertex[ RSize - 1 ];
//			result = TR_Decomposition( CSize, RSize, _Fc, _Fr, _tm, true );
//			printf("%u Components.\n\n", result.first);
//			printf("\n\t\t\t\tno Trimm\twith Trimm\n");
//			printf("Tarjan's:\t\t\t%f\n", result.second);
//			delete [] _tm;
//
//			memcpy( _Fc, Fc, sizeof(Edge)*CSize );
//			memcpy( _Fr, Fr, sizeof(uint32_t)*RSize );
//			memcpy( _Bc, Bc, sizeof(Edge)*CSize );
//			memcpy( _Br, Br, sizeof(uint32_t)*RSize );
//			BFS_vertex * _bm = new BFS_vertex[ RSize - 1 ];
//			result = BFS( CSize, RSize, _Fc, _Fr, _bm, true );
//			printf("CPU BFS:\t\t\t%f\n", result.second);
//			delete [] _bm;
//
//			/*
//			memcpy( _Fc, Fc, sizeof(Edge)*CSize );
//			memcpy( _Fr, Fr, sizeof(uint32_t)*RSize );
//			memcpy( _Bc, Bc, sizeof(Edge)*CSize );
//			memcpy( _Br, Br, sizeof(uint32_t)*RSize );
//			result = parallel_FWD( CSize, RSize, _Fc, _Fr, RSize - 1, 3, true );
//			printf("CPU BFS, 3 cores\t\t%f, %u vertices found\n", result.second, result.first);
//
//			memcpy( _Fc, Fc, sizeof(Edge)*CSize );
//			memcpy( _Fr, Fr, sizeof(uint32_t)*RSize );
//			memcpy( _Bc, Bc, sizeof(Edge)*CSize );
//			memcpy( _Br, Br, sizeof(uint32_t)*RSize );
//			result = parallel_FWD( CSize, RSize, _Fc, _Fr, RSize - 1, 2, true );
//			printf("CPU BFS, 2 cores\t\t%f, %u vertices found\n", result.second, result.first);
//
//			memcpy( _Fc, Fc, sizeof(Edge)*CSize );
//			memcpy( _Fr, Fr, sizeof(uint32_t)*RSize );
//			memcpy( _Bc, Bc, sizeof(Edge)*CSize );
//			memcpy( _Br, Br, sizeof(uint32_t)*RSize );
//			result = parallel_FWD( CSize, RSize, _Fc, _Fr, RSize - 1, 4, true );
//			printf("CPU BFS, 4 cores\t\t%f, %u vertices found\n", result.second, result.first);
//
//			memcpy( _Fc, Fc, sizeof(Edge)*CSize );
//			memcpy( _Fr, Fr, sizeof(uint32_t)*RSize );
//			memcpy( _Bc, Bc, sizeof(Edge)*CSize );
//			memcpy( _Br, Br, sizeof(uint32_t)*RSize );
//			result = _CPU_BFS( CSize, RSize, _Fc, _Fr, RSize - 1, 1, true );
//			printf("GCPU BFS, 1 cores\t\t%f\n", result.second);
//
//			memcpy( _Fc, Fc, sizeof(Edge)*CSize );
//			memcpy( _Fr, Fr, sizeof(uint32_t)*RSize );
//			memcpy( _Bc, Bc, sizeof(Edge)*CSize );
//			memcpy( _Br, Br, sizeof(uint32_t)*RSize );
//			result = _CPU_BFS( CSize, RSize, _Fc, _Fr, RSize - 1, 2, true );
//			printf("GCPU BFS, 2 cores\t\t%f\n", result.second);
//
//			memcpy( _Fc, Fc, sizeof(Edge)*CSize );
//			memcpy( _Fr, Fr, sizeof(uint32_t)*RSize );
//			memcpy( _Bc, Bc, sizeof(Edge)*CSize );
//			memcpy( _Br, Br, sizeof(uint32_t)*RSize );
//			result = _CPU_BFS( CSize, RSize, _Fc, _Fr, RSize - 1, 3, true );
//			printf("GCPU BFS, 3 cores\t\t%f\n", result.second);
//
//			memcpy( _Fc, Fc, sizeof(Edge)*CSize );
//			memcpy( _Fr, Fr, sizeof(uint32_t)*RSize );
//			memcpy( _Bc, Bc, sizeof(Edge)*CSize );
//			memcpy( _Br, Br, sizeof(uint32_t)*RSize );
//			result = _CPU_BFS( CSize, RSize, _Fc, _Fr, RSize - 1, 4, true );
//			printf("GCPU BFS, 4 cores\t\t%f\n", result.second);
//
//			memcpy( _Fc, Fc, sizeof(Edge)*CSize );
//			memcpy( _Fr, Fr, sizeof(uint32_t)*RSize );
//			memcpy( _Bc, Bc, sizeof(Edge)*CSize );
//			memcpy( _Br, Br, sizeof(uint32_t)*RSize );
//			CFB_vertex * _cfm = new CFB_vertex[ RSize - 1 ];
//			result = _CPU_FB( CSize, RSize, _Fc, _Fr, _Bc, _Br, _cfm, RSize - 1, 4, true );
//			printf("GCPU FB, 4 cores\t\t%f\n", result.second);
//			delete [] _cfm;
//			*/
//			/*
//			memcpy( _Fc, Fc, sizeof(Edge)*CSize );
//			memcpy( _Fr, Fr, sizeof(uint32_t)*RSize );
//			memcpy( _Bc, Bc, sizeof(Edge)*CSize );
//			memcpy( _Br, Br, sizeof(uint32_t)*RSize );
//			FB_vertex * _xm = new FB_vertex[ RSize - 1 ];
//			result = Forward_only( CSize, RSize, _Fc, _Fr, _xm, 0, true );
//			printf("GPU BFS:\t\t\t%f\n", result.second);
//			delete [] _xm;
//
//			// 				memcpy( _Fc, Fc, sizeof(Edge)*CSize );
//			// 				memcpy( _Fr, Fr, sizeof(uint32_t)*RSize );
//			// 				memcpy( _Bc, Bc, sizeof(Edge)*CSize );
//			// 				memcpy( _Br, Br, sizeof(uint32_t)*RSize );
//			// 				FB_vertex * _fm = new FB_vertex[ RSize - 1 ];
//			// 				result = FB_Decomposition( CSize, RSize, _Fc, _Fr, _Bc, _Br, _fm, false, true );
//			// 				printf("%u Components.\n\n", result.first);
//			// 				printf("\n\t\t\t\tno Trimm\twith Trimm\n");
//			// 				if ( result.second < Time_limit )
//			// 					printf("FB:\t\t\t\t%f\t", result.second);
//			// 				else
//			// 					printf("FB:\t\t\t\t>%f\t", Time_limit);
//			// 				delete [] _fm;
//
//
//			memcpy( _Fc, Fc, sizeof(Edge)*CSize );
//			memcpy( _Fr, Fr, sizeof(uint32_t)*RSize );
//			memcpy( _Bc, Bc, sizeof(Edge)*CSize );
//			memcpy( _Br, Br, sizeof(uint32_t)*RSize );
//			FB_vertex * _fm = new FB_vertex[ RSize - 1 ];
//			result = FB_Decomposition( CSize, RSize, _Fc, _Fr, _Bc, _Br, _fm, false, true );
//			if ( result.second < Time_limit )
//			printf("FB:\t\t\t\t%f\t", result.second);
//			else
//			printf("FB:\t\t\t\t>%f\t", Time_limit);
//			delete [] _fm;
//
//
//			memcpy( _Fc, Fc, sizeof(Edge)*CSize );
//			memcpy( _Fr, Fr, sizeof(uint32_t)*RSize );
//			memcpy( _Bc, Bc, sizeof(Edge)*CSize );
//			memcpy( _Br, Br, sizeof(uint32_t)*RSize );
//			FB_vertex * _fm2 = new FB_vertex[ RSize - 1 ];
//			result = FB_Decomposition( CSize, RSize, _Fc, _Fr, _Bc, _Br, _fm2, true, true );
//			if ( result.second < Time_limit )
//			printf("%f\n", result.second);
//			else
//			printf(">%f\n", Time_limit);
//			delete [] _fm2;
//
//			memcpy( _Fc, Fc, sizeof(Edge)*CSize );
//			memcpy( _Fr, Fr, sizeof(uint32_t)*RSize );
//			memcpy( _Bc, Bc, sizeof(Edge)*CSize );
//			memcpy( _Br, Br, sizeof(uint32_t)*RSize );
//			COL_vertex * _cm = new COL_vertex[ RSize - 1 ];
//			result = COL_Decomposition( CSize, RSize, _Fc, _Fr, _Bc, _Br, _cm, MIN, true );
//			if ( result.second < Time_limit )
//			printf("Colouring Max-Min:\t\t%f\n", result.second);
//			else
//			printf("Colouring Max-Min:\t\t>%f\n", Time_limit);
//			delete [] _cm;
//
//
//			memcpy( _Fc, Fc, sizeof(Edge)*CSize );
//			memcpy( _Fr, Fr, sizeof(uint32_t)*RSize );
//			memcpy( _Bc, Bc, sizeof(Edge)*CSize );
//			memcpy( _Br, Br, sizeof(uint32_t)*RSize );
//			OBF_vertex * _om01 = new OBF_vertex[ RSize - 1 ];
//			result = OBF_Decomposition( CSize, RSize, _Fc, _Fr, _Bc, _Br, _om01, false, 10000, false, true );
//			if ( result.second < Time_limit )
//			printf("OBF:\t\t\t\t%f\t", result.second);
//			else
//			printf("OBF:\t\t\t\t>%f\t", Time_limit);
//			delete [] _om01;
//
//			memcpy( _Fc, Fc, sizeof(Edge)*CSize );
//			memcpy( _Fr, Fr, sizeof(uint32_t)*RSize );
//			memcpy( _Bc, Bc, sizeof(Edge)*CSize );
//			memcpy( _Br, Br, sizeof(uint32_t)*RSize );
//			OBF_vertex * _om = new OBF_vertex[ RSize - 1 ];
//			result = OBF_Decomposition( CSize, RSize, _Fc, _Fr, _Bc, _Br, _om, false, 10000, true, true );
//			if ( result.second < Time_limit )
//			printf("%f\n", result.second);
//			else
//			printf(">%f\n", Time_limit);
//			delete [] _om;
//
//			memcpy( _Fc, Fc, sizeof(Edge)*CSize );
//			memcpy( _Fr, Fr, sizeof(uint32_t)*RSize );
//			memcpy( _Bc, Bc, sizeof(Edge)*CSize );
//			memcpy( _Br, Br, sizeof(uint32_t)*RSize );
//			OBF_vertex * _om21 = new OBF_vertex[ RSize - 1 ];
//			uint32_t lim = (10000 < ((RSize - 2) / 100)) ? (RSize - 2) / 100 : 10000 ;
//			result = OBF_Decomposition( CSize, RSize, _Fc, _Fr, _Bc, _Br, _om21, true, lim, false, true );
//			if ( result.second < Time_limit )
//			printf("OBF+COL with %u limit:\t%f\t", lim, result.second);
//			else
//			printf("OBF+COL with %u limit:\t>%f\t", lim, Time_limit);
//			delete [] _om21;
//
//			memcpy( _Fc, Fc, sizeof(Edge)*CSize );
//			memcpy( _Fr, Fr, sizeof(uint32_t)*RSize );
//			memcpy( _Bc, Bc, sizeof(Edge)*CSize );
//			memcpy( _Br, Br, sizeof(uint32_t)*RSize );
//			OBF_vertex * _om2 = new OBF_vertex[ RSize - 1 ];
//			result = OBF_Decomposition( CSize, RSize, _Fc, _Fr, _Bc, _Br, _om2, true, lim, true, true );
//			if ( result.second < Time_limit )
//			printf("%f\n", result.second);
//			else
//			printf(">%f\n", Time_limit);
//			delete [] _om2;
//
//			memcpy( _Fc, Fc, sizeof(Edge)*CSize );
//			memcpy( _Fr, Fr, sizeof(uint32_t)*RSize );
//			memcpy( _Bc, Bc, sizeof(Edge)*CSize );
//			memcpy( _Br, Br, sizeof(uint32_t)*RSize );
//			OBF_vertex * _om31 = new OBF_vertex[ RSize - 1 ];
//			result = OBF_Decomposition( CSize, RSize, _Fc, _Fr, _Bc, _Br, _om31, true, 2 * lim, false, true );
//			if ( result.second < Time_limit )
//			printf("OBF+COL with %u limit:\t%f\t", 2 * lim, result.second);
//			else
//			printf("OBF+COL with %u limit:\t>%f\t", 2 * lim, Time_limit);
//			delete [] _om31;
//
//			memcpy( _Fc, Fc, sizeof(Edge)*CSize );
//			memcpy( _Fr, Fr, sizeof(uint32_t)*RSize );
//			memcpy( _Bc, Bc, sizeof(Edge)*CSize );
//			memcpy( _Br, Br, sizeof(uint32_t)*RSize );
//			OBF_vertex * _om3 = new OBF_vertex[ RSize - 1 ];
//			result = OBF_Decomposition( CSize, RSize, _Fc, _Fr, _Bc, _Br, _om3, true, 2 * lim, true, true );
//			if ( result.second < Time_limit )
//			printf("%f\n", result.second);
//			else
//			printf(">%f\n", Time_limit);
//			delete [] _om3;
//
//			memcpy( _Fc, Fc, sizeof(Edge)*CSize );
//			memcpy( _Fr, Fr, sizeof(uint32_t)*RSize );
//			memcpy( _Bc, Bc, sizeof(Edge)*CSize );
//			memcpy( _Br, Br, sizeof(uint32_t)*RSize );
//			OBF_vertex * _om41 = new OBF_vertex[ RSize - 1 ];
//			result = OBF_Decomposition( CSize, RSize, _Fc, _Fr, _Bc, _Br, _om41, true, 5 * lim, false, true );
//			if ( result.second < Time_limit )
//			printf("OBF+COL with %u limit:\t%f\t", 5 * lim, result.second);
//			else
//			printf("OBF+COL with %u limit:\t>%f\t", 5 * lim, Time_limit);
//			delete [] _om41;
//
//			memcpy( _Fc, Fc, sizeof(Edge)*CSize );
//			memcpy( _Fr, Fr, sizeof(uint32_t)*RSize );
//			memcpy( _Bc, Bc, sizeof(Edge)*CSize );
//			memcpy( _Br, Br, sizeof(uint32_t)*RSize );
//			OBF_vertex * _om4 = new OBF_vertex[ RSize - 1 ];
//			result = OBF_Decomposition( CSize, RSize, _Fc, _Fr, _Bc, _Br, _om4, true, 5 * lim, true, true );
//			if ( result.second < Time_limit )
//			printf("%f\n", result.second);
//			else
//			printf(">%f\n", Time_limit);
//			delete [] _om4;
//			printf("\n=======================================================\n\n");
//			}
//			break;
//			*/
//		}
//	}
//	catch (const char * e)
//	{
//		printf("%s\n", e);
//		return 1;
//	}
//
//	if (gen != NULL)
//		delete gen;
//
//	if (!dve) {
//		delete[] Fr;
//		delete[] Fc;
//		delete[] Br;
//		delete[] Bc;
//	}
//	getchar();
//	return 0;
//}
