/**
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 */
#include <stdio.h>
#include <stdlib.h>

#include <cuda.h>
#include "CPUDriver.h"

//enum data  { Lorenz, Benard, Cylinder, Ocean,Tornado };

int main(int argc, char **argv) {

//	if (argc != 3) {
//		perror("Argument number is wrong!\n");
//		perror("You must have execute file!\n");
//		exit(-1);
//	}
//
//	int currentDataset = atoi(argv[1]);
//	int currenttau = atoi(argv[2]);
	int currentDataset = Lorenz	; // they should be the argument from the input
	int currenttau = 10; //
	printf("Hello");
	CPU_Driver* cp = new CPU_Driver();
	cp->Init(currenttau, currentDataset);

	cp->Run();
	//cp->Run_GPU();

	return 0;
}
