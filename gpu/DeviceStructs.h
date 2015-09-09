/*
 * DeviceStructs.h
 *
 *  Created on: Sep 26, 2012
 *      Author: mwalker
 */

#ifndef DEVICESTRUCTS_H_
#define DEVICESTRUCTS_H_

#include <cuda.h>
#include <curand_kernel.h>
#include "parameters.h"

//Sparse matrix type
//typedef cusp::ell_matrix<int,double,cusp::device_memory> sp_type;

//Domain-specific data structure with device memory pointers
struct DeviceDomainStruct {

	int N_States; //Number of grid states
	int N_Ele; //Number of grid elements
//	sp_type* A; //Sparse matrix array for left-hand side matrices

	double* States; //Grid states
	double* Array_1; //Array of buffer that holds right-hand side terms for each state
	//double* Array_2; //Temporary array used by cusp for sparse matrix multiplication
	double* Boundary_Rates; //Boundary rate of change terms
	int* TropC; //Array of troponin C labels
	int* bJSR; //Array of JSR labels
	double* V0; //Array of cell volumes
	double* TTSurfaceArea; //Array of cell TT surface areas
	double* SRSurfaceArea; //Array of cell SR surface areas
	double Vol_JSR; //JSR volume (um^3)
	double* aij; //Finite volume matrix coefficients (groups of 5)
	int* nij; //Finite volume matrix neighbor indices (groups of 4)

};

//Global variable data structure with device memory pointers
struct DeviceGlobalStruct {

	double* Global_Vars; //Global variable states
	double* GV_src1; //Global variable source terms
	int N_RyR; //Number of RyRs
	int N_LCC; //Number of LCCs
	int N_Channels; //Number of channels
	int* RyR_States; //Pointer to RyR state
	int* RyR_Ele; //RyR grid elements
	int* RyR_JSR_Ele; //RyR JSR grid elements
	int* RyR_Neighb; //RyR neighbors array (4*N_RyR)
	int* LCC_States; //Array of LCC states
	int* LCCV_States; //Array of LCC V-dep states
	int* LCC_Ele; //LCC grid elements
	double* R; //Residual array
	curandState* randStates; //PRNG state array

};

//Function declarations
void InitializeDomainStruct(DeviceDomainStruct* dds, GridStruct* gs, int bSecondaryGradient, ParamStruct* ps);
void InitializeLeastSquares(DeviceDomainStruct* dds, GridStruct* gs, int bSecondaryGradient, ParamStruct* ps);
void InitializeExplicitSparseMatrixFast(GridStruct* gs,int state,double** MLS,int** NLS, int bSecondaryGradient, ParamStruct* ps);
void InitializeDeviceGlobalStruct(DeviceGlobalStruct* dgs, SimData* sd, ParamStruct* ps);

void ExplicitStep(double dt, DeviceDomainStruct* dds, DeviceGlobalStruct* dgs, ParamStruct* ps);
void UpdateChannelGating(double dt, DeviceDomainStruct* dds, DeviceGlobalStruct* dgs);
void InitializeDeviceConstants(ParamStruct* ps);
void FreeDomainStruct(DeviceDomainStruct* dds);
void FreeGlobalStruct(DeviceGlobalStruct* dgs);

#endif /* DEVICESTRUCTS_H_ */
