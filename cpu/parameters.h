/*
 * parameters.h
 *
 *  Created on: Jun 27, 2012
 *      Author: mwalker
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <math.h>
#include <map>

/////////////////////////////
//Input parameters
/////////////////////////////

const int I_OUTPUT_GRID_FLAG = 0; //Flag if grid files should be written
const int I_FLAG_FIDELITY_SIM= 1; //Flag to do fidelity simulations that abort early
const int I_FLAG_CLOCK_SEED = 2; //Flag using clock to seed random number generator
const int I_RAND_SEED_0 = 3; //Random number seed
const int I_RAND_ADD_A = 4; //Random number seed (mulitply by device index)
const int I_RAND_ADD_B = 5; //Random number seed (multiply by ensemble index)
const int I_INITIALIZE_FROM_FILE = 6; //Flag using initial conditions from file
const int I_CUDA_DEVICE = 7; //GPU Device to use
const int I_FLAG_SECONDARY_GRADIENT = 8; //Flag using secondary gradient
const int I_FLAG_NO_FLUX_BOUNDARY = 9; //Flag using no flux at domain boundaries
const int I_N_RYR_FIDELITY = 10; //Number of open RyRs needed to constitute a spark (if FLAG_FIDELITY_SIM is on)
const int I_FLAG_LCC = 11; //Flag for LCC trigger
const int I_LCC_INDEX = 12; //Trigger LCC to open
const int I_FLAG_V_CLAMP = 13; //Flag voltage clamp
const int I_FLAG_NO_GATING = 14; //Flag turning off stochastic gating simulation
const int I_GATING_FFWD = 15; //Flag channel gating w/o diffusion until a channel opens (use for gain simulations at low potentials)
const int I_FLAG_RYR_OUT = 16; //Flag to output RyR states
const int I_N_SIMS = 17; //Number of simulations to run
const int I_INIT_RYR = 18; //Index of an initial RyR to open (set to -1 to ignore)
const int I_T_FINAL = 0; //final simulation time (ms)
const int I_T_STEP = 1; //step size (ms)
const int I_OUTPUT_STATES_INTERVAL = 2; //Output time intervals for states
const int I_OUTPUT_GRID_INTERVAL = 3; //Output time intervals for grids
const int I_C_MIN = 4; //Minimum amount of subspace Ca2+ if FLAG_FIDELITY_SIM is on
const int I_D_CA = 5; //Ca2+ diffusion coefficient (um^2/ms)
const int I_D_CA_JSR = 6; //Ca2+ in JSR
const int I_D_ATP = 7; //ATP
const int I_D_CMDN = 8; //Calmodulin
const int I_D_TRPN = 9; //TroponinC
const int I_D_CSQN = 10; //Calsequestrin binding sites
const int I_D_DYE = 11; //Dye (Fluo-3; Hake et al 2012)
const int I_D_EGTA = 12; //EGTA
const int I_C_0_CA = 13; //Initial cytosolic Ca2+ (uM)
const int I_C_0_CA_JSR = 14; //Initial JSR Ca2+ (uM)
const int I_DEFAULT_VM = 15; //Default membrane potential (mV)
const int I_DEFAULT_CA_NSR = 16; //Default NSR Ca2+ concentration (uM)
const int I_CA_0 = 17; //Extracellular Ca2+ (uM)
const int I_B_TOT_ATP = 18; //Total ATP concentration (uM) (Hake et al 2012)
const int I_K_OFF_ATP = 19; //ATP off rate (1/ms)
const int I_K_ON_ATP = 20; //ATP on rate (1/uM-ms)
const int I_B_TOT_CMDN = 21; //Total calmodulin concentration (uM) (Hake et al 2012)
const int I_K_OFF_CMDN = 22; //calmodulin off rate (1/ms)
const int I_K_ON_CMDN = 23; //calmodulin on rate (1/uM-ms)
const int I_B_TOT_TRPN = 24; //Total high-affinity troponin concentration (uM) (Smith et al 1998)
const int I_K_OFF_TRPN = 25; //High-affinity Troponin off rate (1/ms)
const int I_K_ON_TRPN = 26; //High-affinityTroponin on rate (1/uM-ms)
const int I_B_TOT_SL_HI = 27; //Total high-affinity sarcolemma binding site density (umol/um^2) (Peskoff et al 1992)
const int I_K_D_SL_HI = 28; //uM
const int I_B_TOT_SL_LO = 29; //Total low-affinity sarcolemma binding site density (umol/um^2) (Peskoff et al 1992)
const int I_K_D_SL_LO = 30; //uM
const int I_B_TOT_DYE = 31; //Total cytosolic dye concentration (uM) (Smith et al 1998)
const int I_K_OFF_DYE = 32; //Dye off rate (1/ms)
const int I_K_ON_DYE = 33; //Dye on rate (1/uM-ms)
const int I_B_T_JSR = 34; //uM
const int I_K_M_JSR = 35; //uM
const int I_V_REFILL = 36; //ms-1
const int I_V_RYR = 37; //ms-1 (0.245pA)*0.8
const int I_V_CELL = 38; //um^3
const int I_K_ON_CSQN = 39; //Diffusion-limited on rate
const int I_K_OFF_CSQN = 40;
const int I_RYR_ETA = 41;
const int I_RYR_A_STAR = 42;
const int I_RYR_EPS_CC = 43;
const int I_RYR_EPS_OO = 44;
const int I_RYR_K_PLUS = 45;
const int I_RYR_K_MINUS = 46; //ms-1
const int I_RYR_PHI_M = 47;
const int I_RYR_PHI_B = 48;
const int I_K_D_i = 49; //uM
const int I_K_D_SR = 50; //uM
const int I_A_P = 51; //uM
const int I_LCC_DURATION = 52; //Duration to open LCC for (if FLAG_LCC is on)
const int I_T_CLAMP = 53; //Voltage clamp duration (if FLAG_V_CLAMP is on)
const int I_V_CLAMP = 54; //Voltage clamp setting (mV) (if FLAG_V_CLAMP is on)
const int I_INIT_RYR_HOLD = 55; //Hold initial RyR (INIT_RYR) open for at least this long, does not apply to randomly-opened RyRs
const int I_FILE_BASE = 0; //Input mesh filename base
const int I_OUTPUT_STATE = 1; //Output filename for global states
const int I_OUTPUT_GRID = 2; //Output filename for grid
const int I_PARAM_TITLE = 3; //Parameter set title

/////////////////////////////
//Hard-coded parameters
/////////////////////////////

const int N_PARAMETERS_INT = 19;
const int N_PARAMETERS_REAL = 56;
const int N_PARAMETERS_CHAR = 4;
const char DEFAULT_PARAMETER_FILE[255] = "TestSet.xml";

//Face types
const int FACE_DEFAULT = 0;
const int FACE_BOUND = 1;
const int FACE_NOFLUX = 2;

//Input face types
const int IN_FACE_DEFAULT = 0;
const int IN_FACE_BOUNDARY = 1;
const int IN_FACE_TT = 2;
const int IN_FACE_JSR = 3;

//Cell domains
const int DOMAIN_CYTO = 0;
const int DOMAIN_JSR = 1;

//Channel parameters
const int STATE_RYR_OPEN = 1; //Open state
const int STATE_RYR_CLOSED = 0; //Closed state
const int DEFAULT_STATE_RYR = 0; //Default state

const int STATE_LCC_OPEN_1 = 6;
const int STATE_LCC_OPEN_2 = 12;
const int DEFAULT_STATE_LCC = 1;

const int STATE_LCCV_OPEN = 2;
const int STATE_LCCV_CLOSED = 1;
const int DEFAULT_STATE_LCCV = 2;

//Indices of each state variable
const int INDEX_CA = 0;
const int INDEX_ATP = 1;
const int INDEX_CMDN = 2;
const int INDEX_TRPN = 3;
const int INDEX_CSQN = 4;
const int INDEX_DYE = 5;

//State variable arrays
const int N_STATES = 6; //Number of states per grid cell (cytosol)

//Global variable indices
const int N_GLOBAL_VARS = 2;
const int GLOBAL_VM = 0;
const int GLOBAL_CA_NSR= 1;

//Default globals must match order of global variable indices
//const double DEFAULT_GLOBALS[N_GLOBAL_VARS] = {DEFAULT_VM, DEFAULT_CA_NSR};
//const double DIFF[N_STATES] = {D_CA,D_ATP,D_CMDN,D_TRPN,D_CSQN,D_DYE};


///////////////////////////////
// Data structures
///////////////////////////////

//Grid properties structure
struct GridStruct {

	int N_States; //Number of cell states
	double** States; //2D arrays of cell states
	double* States_Buffer; //Finite volume boundary rates

	int* Ele_1; //Array of each face first element
	int* Ele_2; //Array of each face other element
	int* FaceType; //Array of each face type
	double* aij; //Finite volume matrix coefficients (groups of 5)
	double** Boundary_Rates; //Finite volume boundary rates
	double* Face_C; //Array of each face lumped constant
	double* dXi; //Array of each face Xi vector magnitude
	double* Area_f; //Array of each face area
	double* Afx; //Array of each face A_f vector
	double* Afy; //Array of each face A_f vector
	double* Afz; //Array of each face A_f vector
	double* eXix; //Array of each face e_Xi vector
	double* eXiy; //Array of each face e_Xi vector
	double* eXiz; //Array of each face e_Xi vector
	double* r1x; //Array of each face r1 vector
	double* r1y; //Array of each face r1 vector
	double* r1z; //Array of each face r1 vector
	double* r2x; //Array of each face r2 vector
	double* r2y; //Array of each face r2 vector
	double* r2z; //Array of each face r2 vector
	double* V0; //Array of each cell's volume
	double* TTSurfaceArea; //Array of each cell's TT surface area
	double* SRSurfaceArea; //Array of each cell's SR surface area
	int* nij; //Finite volume matrix neighbor indices (groups of 4)
	int* Domain; //Array of each cell's domain type
	int* TropC; //Array of each cell's troponin amount (0 or 1)
	int N_Ele; //Number of triangular cells
	int N_Face; //Number of edges
	double Vol_JSR;

	double* Boundaries; //Dirichlet boundary conditions

	//Mesh node data
	int N_Nodes;
	double* Node_X;
	double* Node_Y;
	double* Node_Z;
	int* Ele_Nodes;

};

//Simulation data structure
struct SimData {

	GridStruct Grid; //Cytoplasm mesh data

	double Global_Vars[N_GLOBAL_VARS]; //Global variables
	double Global_Src[N_GLOBAL_VARS]; //Global variables buffer

	int N_RyR; //Number of RyRs
	int N_LCC; //Number of LCCs
	int N_Channels; //Number of channels
	int* LCC_States; //Array of LCC states
	int* LCCV_States; //Array of LCC V-dep gating states
	int* LCC_Ele; //LCC grid elements
	int* RyR_States; //Array of RyR States
	int* RyR_Ele; //RyR grid elements
	int* RyR_JSR_Ele; //RyR JSR grid elements
	int* RyR_Neighb; //RyR neighbors array (4*N_RyR)
	double* R; //stochastic residual

	int device;
	int ensemble_index;

};

//Parameters data structure
struct ParamStruct {

	int Integers[N_PARAMETERS_INT];
	double Reals[N_PARAMETERS_REAL];
	char Chars[N_PARAMETERS_CHAR][255];
	double Diff[N_STATES];
	double Default_Globals[N_GLOBAL_VARS];
	//char timestamp[255];
	char rand_id[255];
	char state_file[255];
	char grid_file[255];

};

/////////////////////////////
//Function declarations
/////////////////////////////

void ReadParameters(ParamStruct* ps, const char* param_file);
void ReadGrid(SimData* sd, GridStruct* gs, ParamStruct* ps);
void ReadGridStates(GridStruct* gs, ParamStruct* ps);
void InitializeDefaultGrid(SimData* sd, ParamStruct* ps);
void InitializeDefaultChannels(SimData* sd);
void ResetGlobals(SimData* sd, ParamStruct* ps);
void InitializeGridMaps(GridStruct* gs);
void OpenRandomRyR(ParamStruct* ps, SimData* sd);
void InitializeDomainStruct(GridStruct* gs, int bSecondaryGradient, ParamStruct* ps);
void InitializeExplicitSparseMatrixFast(GridStruct* gs,int state,double** MLS,int** NLS, int bSecondaryGradient, ParamStruct* ps);

void ResetDefaultChannels(SimData* sd);
void ResetDefaultGridStates(SimData* sd, ParamStruct* ps);

void RunSimulation(SimData* sd, ParamStruct* ps);
void ExplicitStep(double dt, SimData* sd, ParamStruct* ps);
void UpdateChannelGating(double dt, SimData* sd, ParamStruct* ps);

void CreateOutputFiles(int run_nmbr, int rank, ParamStruct* ps);
void WriteToGridFiles(int n, double t, GridStruct* gs, SimData* sd, ParamStruct* ps);
void WriteToStateFiles(double t, SimData* sd, ParamStruct* ps);

void FreeGridStruct(GridStruct* gs);
void FreeSimData(SimData* sd);

#endif /* PARAMETERS_H_ */
