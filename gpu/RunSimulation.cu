/*
 * RunSimulation.cpp
 *
 *  Created on: Jun 27, 2012
 *      Author: mwalker
 */

#include <stdio.h>
#include <math.h>
#include <time.h>
#include "DeviceStructs.h"
#include <cuda_profiler_api.h>

//Constant memory variables

__device__ __constant__ double PARAMS_DEVICE[N_PARAMETERS_REAL];

//CUDA function declarations

//Sets all elements in array b to zero
__global__ void resetVectorToZero(double* b, int N)
{
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if (idx<N) b[idx] = 0;
}

//Element-wise addition of b*alpha to a
__global__ void addMultVectors(double* a, double* b, double alpha, int N)
{
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if (idx<N) a[idx] += alpha*b[idx];
}

__global__ void ResetResiduals(curandState * randStates, double* R) {

	int idx = blockIdx.x*blockDim.x + threadIdx.x;

	curandState localState = randStates[idx];
	R[idx] = log((double)curand_uniform( &localState ));
	randStates[idx] =localState;

}

//Fast buffering
__global__ void UpdateFastBuffer(int N_Ele, double* grid_ca, double* buffer_ca, double* TTSurfaceArea,double* V0)
{
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if (idx<N_Ele) {
	  double btothi = TTSurfaceArea[idx]*PARAMS_DEVICE[I_B_TOT_SL_HI] / (V0[idx]*(1e-15));
	  double btotlo= TTSurfaceArea[idx]*PARAMS_DEVICE[I_B_TOT_SL_LO] / (V0[idx]*(1e-15));
	  double B = 1 + (btothi*PARAMS_DEVICE[I_K_D_SL_HI]/((PARAMS_DEVICE[I_K_D_SL_HI]+grid_ca[idx])*(PARAMS_DEVICE[I_K_D_SL_HI]+grid_ca[idx]))) +
			  	     (btotlo*PARAMS_DEVICE[I_K_D_SL_LO]/((PARAMS_DEVICE[I_K_D_SL_LO]+grid_ca[idx])*(PARAMS_DEVICE[I_K_D_SL_LO]+grid_ca[idx])));
	  buffer_ca[idx] /= B;
	/*if (TTSurfaceArea[idx] > 0) { //electrodiffusion approximation (DO NOT USE)
		buffer_ca[idx] /= 38.95;
	}*/
  }
}

//Reaction terms
__global__ void UpdateReactionTerms(int N_Ele, double* grid,
									double* buffer,
									double* GV0,
									double* aij,
									int* nij,
									int* TropC, double* V0, double* TTSurfaceArea,
									double* SRSurfaceArea, double* Boundary_Rates,
									int* bJSR, double Vol_JSR) {

	int j = blockIdx.x*blockDim.x + threadIdx.x;

	if (j < N_Ele) {
		double c_ca = grid[INDEX_CA*N_Ele + j];
		int b_idx;
		double r_net, c_b;

		int idx_ca = INDEX_CA*N_Ele;
		int idx_atp = INDEX_ATP*N_Ele;
		int idx_cmdn = INDEX_CMDN*N_Ele;
		int idx_trpn = INDEX_TRPN*N_Ele;
		int idx_csqn = INDEX_CSQN*N_Ele;
		int idx_dye = INDEX_DYE*N_Ele;

		//Reset buffer values for immobile buffers that may not be set
		buffer[idx_trpn + j] = 0;
		buffer[idx_csqn + j] = 0;

		//Update diffusion
		double a0, a1, a2, a3, a4; //Coefficients
		a0 = aij[5*j];
		a1 = aij[5*j+1];
		a2 = aij[5*j+2];
		a3 = aij[5*j+3];
		a4 = aij[5*j+4];
		int n1, n2, n3, n4; //Neighbors
		n1 = nij[4*j];
		n2 = nij[4*j+1];
		n3 = nij[4*j+2];
		n4 = nij[4*j+3];

		int n0 = idx_ca + j;
		buffer[n0] = PARAMS_DEVICE[I_D_CA]*(c_ca*a0 + grid[idx_ca + n1]*a1 + grid[idx_ca + n2]*a2 + grid[idx_ca + n3]*a3 + grid[idx_ca + n4]*a4) + Boundary_Rates[n0];
		/*buffer[n0] += grid[idx_ca + n1]*a1*PARAMS_DEVICE[I_D_CA];
		buffer[n0] += grid[idx_ca + n2]*a2*PARAMS_DEVICE[I_D_CA];
		buffer[n0] += grid[idx_ca + n3]*a3*PARAMS_DEVICE[I_D_CA];
		buffer[n0] += grid[idx_ca + n4]*a4*PARAMS_DEVICE[I_D_CA];
		buffer[n0] += Boundary_Rates[n0];*/

		n0 = idx_dye + j;
		//buffer[n0] = grid[n0]*a0*PARAMS_DEVICE[I_D_DYE];
		buffer[n0] = PARAMS_DEVICE[I_D_DYE]*(grid[n0]*a0+ grid[idx_dye + n1]*a1+ grid[idx_dye + n2]*a2+ grid[idx_dye + n3]*a3+ grid[idx_dye + n4]*a4)+ Boundary_Rates[n0];
		/*buffer[n0] += grid[idx_dye + n1]*a1*PARAMS_DEVICE[I_D_DYE];
		buffer[n0] += grid[idx_dye + n2]*a2*PARAMS_DEVICE[I_D_DYE];
		buffer[n0] += grid[idx_dye + n3]*a3*PARAMS_DEVICE[I_D_DYE];
		buffer[n0] += grid[idx_dye + n4]*a4*PARAMS_DEVICE[I_D_DYE];
		buffer[n0] += Boundary_Rates[n0];*/

		n0 = idx_atp + j;
		//buffer[n0] = grid[n0]*a0*PARAMS_DEVICE[I_D_ATP];
		buffer[n0] = PARAMS_DEVICE[I_D_ATP]*(grid[n0]*a0+ grid[idx_atp + n1]*a1+ grid[idx_atp + n2]*a2+ grid[idx_atp + n3]*a3+ grid[idx_atp + n4]*a4)+ Boundary_Rates[n0];
		/*buffer[n0] += grid[idx_atp + n1]*a1*PARAMS_DEVICE[I_D_ATP];
		buffer[n0] += grid[idx_atp + n2]*a2*PARAMS_DEVICE[I_D_ATP];
		buffer[n0] += grid[idx_atp + n3]*a3*PARAMS_DEVICE[I_D_ATP];
		buffer[n0] += grid[idx_atp + n4]*a4*PARAMS_DEVICE[I_D_ATP];
		buffer[n0] += Boundary_Rates[n0];*/

		n0 = idx_cmdn + j;
		//buffer[n0] = grid[n0]*a0*PARAMS_DEVICE[I_D_CMDN];
		buffer[n0] = PARAMS_DEVICE[I_D_CMDN]*(grid[n0]*a0+ grid[idx_cmdn + n1]*a1+ grid[idx_cmdn + n2]*a2+ grid[idx_cmdn + n3]*a3+ grid[idx_cmdn + n4]*a4)+ Boundary_Rates[n0];
		/*buffer[n0] += grid[idx_cmdn + n1]*a1*PARAMS_DEVICE[I_D_CMDN];
		buffer[n0] += grid[idx_cmdn + n2]*a2*PARAMS_DEVICE[I_D_CMDN];
		buffer[n0] += grid[idx_cmdn + n3]*a3*PARAMS_DEVICE[I_D_CMDN];
		buffer[n0] += grid[idx_cmdn + n4]*a4*PARAMS_DEVICE[I_D_CMDN];
		buffer[n0] += Boundary_Rates[n0];*/

		if (!bJSR[j]) {
			//Troponin C, SR, and SERCA
			if (TropC[j]) {

				b_idx = idx_trpn + j;
				c_b = grid[b_idx];
				r_net = PARAMS_DEVICE[I_K_ON_TRPN]*(PARAMS_DEVICE[I_B_TOT_TRPN]-c_b)*c_ca - PARAMS_DEVICE[I_K_OFF_TRPN]*c_b;
				buffer[idx_ca + j] += -r_net;
				buffer[idx_trpn + j] = r_net;

				//SERCA pump in same zone as troponin C
				double K_i = c_ca / PARAMS_DEVICE[I_K_D_i];
				K_i = K_i*K_i;
				double K_sr = GV0[GLOBAL_CA_NSR] / PARAMS_DEVICE[I_K_D_SR];
				K_sr = K_sr*K_sr;
				double D_cycle = 0.104217 + 17.923*K_sr + K_i*(1.75583e6 + K_sr*7.61673e6) + K_i*K_i*(6.08463e11 + K_sr*4.50544e11);
				double v_cycle = (K_i*K_i*3.24873e12 + K_i*(9.17846e6 - 11478.2*K_sr) - 0.329904*K_sr) / D_cycle;
				double d_ca = (2e-3) * v_cycle * PARAMS_DEVICE[I_A_P]; //uM / ms (Note it is given in uM/s in Williams et al 2011)
				buffer[idx_ca + j] -= d_ca;

			}
			//Dye buffering
			b_idx = idx_dye + j;
			c_b = grid[b_idx];
			r_net = PARAMS_DEVICE[I_K_ON_DYE]*(PARAMS_DEVICE[I_B_TOT_DYE]-c_b)*c_ca - PARAMS_DEVICE[I_K_OFF_DYE]*c_b;
			buffer[idx_ca + j] += -r_net;
			buffer[idx_dye + j] += r_net;

			//ATP buffering
			b_idx = idx_atp + j;
			c_b = grid[b_idx];
			r_net = PARAMS_DEVICE[I_K_ON_ATP]*(PARAMS_DEVICE[I_B_TOT_ATP]-c_b)*c_ca - PARAMS_DEVICE[I_K_OFF_ATP]*c_b;
			buffer[idx_ca + j] += -r_net;
			buffer[idx_atp + j] += r_net;

			//Calmodulin buffering
			b_idx = idx_cmdn + j;
			c_b = grid[b_idx];
			r_net = PARAMS_DEVICE[I_K_ON_CMDN]*(PARAMS_DEVICE[I_B_TOT_CMDN]-c_b)*c_ca - PARAMS_DEVICE[I_K_OFF_CMDN]*c_b;
			buffer[idx_ca + j] += -r_net;
			buffer[idx_cmdn + j] += r_net;

		} else {

			//Calsequestrin buffering
			b_idx = idx_csqn + j;
			c_b = grid[b_idx];
			r_net = PARAMS_DEVICE[I_K_ON_CSQN]*(PARAMS_DEVICE[I_B_T_JSR]-c_b)*c_ca - PARAMS_DEVICE[I_K_OFF_CSQN]*c_b;
			buffer[idx_ca + j] += -r_net;
			buffer[idx_csqn + j] = r_net;

			//JSR refill
			//buffer_ca[j] += (GV0[GLOBAL_CA_NSR] - c_ca) * PARAMS_DEVICE[I_V_REFILL] * PARAMS_DEVICE[I_V_CELL] / Vol_JSR;
			buffer[idx_ca + j] += (GV0[GLOBAL_CA_NSR] - c_ca) * PARAMS_DEVICE[I_V_REFILL]; //Volume-independent refill

		}


	}

}

__global__ void UpdateRyRFlux(double* Grid_Ca, double *Grid_Buffer,
							   double* GV0, double* V0,
							   int N_RyR, int* RyR_States, int* RyR_Ele,
							   int* RyR_JSR_Ele) {

	int i = blockIdx.x*blockDim.x + threadIdx.x;

	if (i < N_RyR) {

		//RyR flux
		if (RyR_States[i]) {
			double JRyR = PARAMS_DEVICE[I_V_RYR] * (Grid_Ca[RyR_JSR_Ele[i]] - Grid_Ca[RyR_Ele[i]])* PARAMS_DEVICE[I_V_CELL];
			//double JRyR = 5.1822 / V_DYAD; //1pA constant flux
			//double JRyR = (0.2*1e-6*1e15)/(1e3*2.0*96485.3365);

			Grid_Buffer[RyR_Ele[i]] += JRyR / V0[RyR_Ele[i]];
			Grid_Buffer[RyR_JSR_Ele[i]] -= JRyR / V0[RyR_JSR_Ele[i]];
		}

	}
}

__global__ void UpdateLCCFlux(double* Grid_Ca, double *Grid_Buffer,
							   double* GV0, double* V0,
							   int N_LCC, int* LCC_States, int* LCCV_States, int* LCC_Ele) {

	int i = blockIdx.x*blockDim.x + threadIdx.x;

	if (i < N_LCC && (LCC_States[i] == STATE_LCC_OPEN_1 || LCC_States[i] == STATE_LCC_OPEN_2) && LCCV_States[i] == STATE_LCCV_OPEN) {

		const double P_CaL = (9.13e-13)*(1e12)*(1e-3); //um^3 / ms
		const double F = 96.5; // C / mmol
		const double RT = 310 * 8.314; // J / mol
		const double F_over_RT = F/RT;
		double V = GV0[GLOBAL_VM];

		//LCC Flux
		double CA_DYAD = Grid_Ca[LCC_Ele[i]];

		if (fabs(V) < 1e-6) {
			V = 1e-6;
		}
		double J_lcc = - P_CaL * V * (2*F_over_RT) * (CA_DYAD*exp(2*V*F_over_RT) - 0.34*PARAMS_DEVICE[I_CA_0]) / (exp(2*V*F_over_RT)-1) / V0[LCC_Ele[i]];
		Grid_Buffer[LCC_Ele[i]] += J_lcc;

	}
}

__global__ void UpdateRyRStates(double dt, int N_RyR, int* RyR_States,
								double* Grid_Ca, int* RyR_Ele,
								curandState* randStates, double* R,
								double* gv0, int* RyR_Neighb,
								int* RyR_JSR_Ele) {

	int idx = blockIdx.x*blockDim.x + threadIdx.x;

	if (idx < N_RyR) {
		int N_Open = 0;
		int N_Closed = 0;
		for (int i = 0; i < 4; i++) {
			if (RyR_Neighb[idx*4 + i] >= 0) {
				if (RyR_States[RyR_Neighb[idx*4 + i]]) {
					N_Open++;
				} else {
					N_Closed++;
				}
			}
		}

		double X, rate;
		double ca_ss = Grid_Ca[RyR_Ele[idx]];

		if (RyR_States[idx]) {
			//Open -> Closed rate
			X = exp( 0.5*PARAMS_DEVICE[I_RYR_A_STAR]*(N_Open*PARAMS_DEVICE[I_RYR_EPS_OO] - N_Closed*PARAMS_DEVICE[I_RYR_EPS_CC]) );
			rate = X*PARAMS_DEVICE[I_RYR_K_MINUS];

		} else {
			//Closed -> Open rate
			//double phi = PARAMS_DEVICE[I_RYR_PHI_M]*Grid_Ca[RyR_JSR_Ele[idx]] + PARAMS_DEVICE[I_RYR_PHI_B];
			//double phi = (pow(PARAMS_DEVICE[I_RYR_PHI_B],PARAMS_DEVICE[I_RYR_PHI_M])+pow(Grid_Ca[RyR_JSR_Ele[idx]],PARAMS_DEVICE[I_RYR_PHI_M]))/(pow(PARAMS_DEVICE[I_RYR_PHI_B],PARAMS_DEVICE[I_RYR_PHI_M])+pow(1000.0,PARAMS_DEVICE[I_RYR_PHI_M]));
			double phi = (1.0 - pow(1000.0/PARAMS_DEVICE[I_RYR_PHI_B],PARAMS_DEVICE[I_RYR_PHI_M])) + pow(Grid_Ca[RyR_JSR_Ele[idx]]/PARAMS_DEVICE[I_RYR_PHI_B],PARAMS_DEVICE[I_RYR_PHI_M]);
			X = exp( 0.5*PARAMS_DEVICE[I_RYR_A_STAR]*(N_Closed*PARAMS_DEVICE[I_RYR_EPS_CC] - N_Open*PARAMS_DEVICE[I_RYR_EPS_OO]) );
			rate = X*phi*PARAMS_DEVICE[I_RYR_K_PLUS]* pow(ca_ss,PARAMS_DEVICE[I_RYR_ETA]);
			//double p_open = Grid_Ca[RyR_JSR_Ele[idx]] / (Grid_Ca[RyR_JSR_Ele[idx]] + 650);
			//rate *= p_open;
			//rate = min(rate,PARAMS_DEVICE[I_RYR_K_MINUS]/10);

		}

		R[idx] += rate*dt;

		if (R[idx] > 0) {
			//Generate new residual
			curandState localState = randStates[idx];
			R[idx] = log((double)curand_uniform( &localState ));
			randStates[idx] = localState;

			if (RyR_States[idx] == 1) {
				RyR_States[idx] = 0;
			} else {
				RyR_States[idx] = 1;
/*
				localState = randStates[idx];
				double runif = (double)curand_uniform( &localState );
				randStates[idx] = localState;
				double p_open = Grid_Ca[RyR_JSR_Ele[idx]] / (Grid_Ca[RyR_JSR_Ele[idx]] + 650);
				if (runif < p_open) {
				}*/

			}
		}

	}

}

__global__ void UpdateLCCStates(double dt, int N_LCC, int N_RyR,
								int* LCC_States, int* LCCV_States,
								int* LCC_Ele, double* Grid_Ca,
								curandState* randStates, double* R,
								double* gv0) {

	int idx = blockIdx.x*blockDim.x + threadIdx.x;

	if (idx < N_LCC) {


		double V = gv0[GLOBAL_VM];

		double ca_ss = Grid_Ca[LCC_Ele[idx]];

		const double fL=0.85; // transition	rate into open state (1/ms)
		const double gL=2.0; //	transition rate	out	of open	state (1/ms)
//		const double gPhosph = 0.049;
		const double fLprime=0.005;	// transition rate into	Ca mode	open state (1/ms)
		const double gLprime=7.0; // transition	rate out of	Ca mode	open state (1/ms)
		const double bL=1.9356;	// mode	transition parameter
		const double bL2=bL*bL;
		const double bL3=bL*bL*bL;
		const double bL4=bL*bL*bL*bL;
		const double aL=2.0; //	mode transition	parameter
		const double aL2=aL*aL;
		const double aL3=aL*aL*aL;
		const double aL4=aL*aL*aL*aL;
		const double omega=0.83*2.0*1.3*0.01;  // mode transition parameter	(1/ms)

		const double alphacf=4.0*1.2*0.416;
		const double betacf=4.0*0.45*0.049;
		const double gammacf=0.83*1.9*1.3*0.31*7.5*0.09233 / 1000.0;	// (ms-1 uM-1)

		const double CCa0_to_C0	= omega;		// = omega
		const double CCa1_to_C1	= omega/bL;	// = omega/bL
		const double CCa2_to_C2	= omega/bL2;	// = omega/bL^2
		const double CCa3_to_C3	= omega/bL3;	// = omega/bL^3
		const double CCa4_to_C4	= omega/bL4;	// = omega/bL^4

		const double yCa_frac=0.4;	// asymptotic value	for	fraction of	LCCs that

		double alpha =	alphacf	* exp(0.012*(V-35.0));
		double beta = betacf *	exp(-0.05*(V-35.0));
		double alpha_prime	= aL*alpha;
		double beta_prime = beta/bL;

		double gamma_rate =	gammacf*ca_ss;

		int i_lcc = idx + N_RyR;
		int lcc_state = LCC_States[idx];

		int target_states[3];
		double trans_rates[3];

		//LCC state
		if (lcc_state == 1) {
			trans_rates[0] = 4.0*alpha;
			trans_rates[1] = gamma_rate;
			trans_rates[2] = 0;
			target_states[0] = 2;
			target_states[1] = 7;
			target_states[2] = 0;
		} else if (lcc_state == 2) {
			trans_rates[0] = beta;
			trans_rates[1] = 3.0*alpha;
			trans_rates[2] = aL*gamma_rate;
			target_states[0] = 1;
			target_states[1] = 3;
			target_states[2] = 8;
		} else if (lcc_state == 3) {
			trans_rates[0] = 2.0*beta;
			trans_rates[1] = 2.0*alpha;
			trans_rates[2] = aL2*gamma_rate;
			target_states[0] = 2;
			target_states[1] = 4;
			target_states[2] = 9;
		} else if (lcc_state == 4) {
			trans_rates[0] = 3.0*beta;
			trans_rates[1] = alpha;
			trans_rates[2] = aL3*gamma_rate;
			target_states[0] = 3;
			target_states[1] = 5;
			target_states[2] = 10;
		} else if (lcc_state == 5) {
			trans_rates[0] = 4.0*beta;
			trans_rates[1] = fL;
			trans_rates[2] = aL4*gamma_rate;
			target_states[0] = 4;
			target_states[1] = 6;
			target_states[2] = 11;
		} else if (lcc_state == 6) {
			trans_rates[0] = gL;
			trans_rates[1] = 0;
			trans_rates[2] = 0;
			target_states[0] = 5;
			target_states[1] = 0;
			target_states[2] = 0;
		} else if (lcc_state == 7) {
			trans_rates[0] = CCa0_to_C0;
			trans_rates[1] = 4.0*alpha_prime;
			trans_rates[2] = 0;
			target_states[0] = 1;
			target_states[1] = 8;
			target_states[2] = 0;
		} else if (lcc_state == 8) {
			trans_rates[0] = beta_prime;
			trans_rates[1] = CCa1_to_C1;
			trans_rates[2] = 3.0*alpha_prime;
			target_states[0] = 7;
			target_states[1] = 2;
			target_states[2] = 9;
		} else if (lcc_state == 9) {
			trans_rates[0] = 2.0*beta_prime;
			trans_rates[1] = CCa2_to_C2;
			trans_rates[2] = 2.0*alpha_prime;
			target_states[0] = 8;
			target_states[1] = 3;
			target_states[2] = 10;
		} else if (lcc_state == 10) {
			trans_rates[0] = 3.0*beta_prime;
			trans_rates[1] = CCa3_to_C3;
			trans_rates[2] = alpha_prime;
			target_states[0] = 9;
			target_states[1] = 4;
			target_states[2] = 11;
		} else if (lcc_state == 11) {
			trans_rates[0] = 4.0*beta_prime;
			trans_rates[1] = CCa4_to_C4;
			trans_rates[2] = fLprime;
			target_states[0] = 10;
			target_states[1] = 5;
			target_states[2] = 12;
		} else if (lcc_state == 12) {
			trans_rates[0] = gLprime;
			trans_rates[1] = 0;
			trans_rates[2] = 0;
			target_states[0] = 11;
			target_states[1] = 0;
			target_states[2] = 0;
		}

		//LCC V-dependent inactivation gate
		double lccv_rate;
		int lccv_state = LCCV_States[idx];
	    double yCa_inf	= yCa_frac/(1.0+exp((V + 12.5)/5.0)) + (1.0-yCa_frac);
	    double tau_yCa	= 60.0 + 340.0/(1.0	+ exp((V+30.0)/12.0));
		if (lccv_state == 1) { //Inactivated state
			lccv_rate = yCa_inf/tau_yCa;
		} else { //Activated state
			lccv_rate = (1.0-yCa_inf)/tau_yCa;
		}

		R[i_lcc] += (trans_rates[0] + trans_rates[1] + trans_rates[2] + lccv_rate)*dt;

		if (R[i_lcc] > 0) {

			//Generate new residual
			curandState localState = randStates[i_lcc];
			R[i_lcc] = log((double)curand_uniform( &localState ));
			randStates[i_lcc] = localState;

			//Generate random uniform
			localState = randStates[i_lcc];
			double runif = (double)curand_uniform(&localState) * (trans_rates[0] + trans_rates[1] + trans_rates[2] + lccv_rate);
			randStates[i_lcc] = localState;

			//Update state
			if ( runif < trans_rates[0]) {
				LCC_States[idx] = target_states[0];
			} else if ( runif < trans_rates[0]+trans_rates[1]) {
				LCC_States[idx] = target_states[1];
			} else if ( runif < trans_rates[0] + trans_rates[1] + trans_rates[2]) {
				LCC_States[idx] = target_states[2];
			} else {
				LCCV_States[idx] = (lccv_state == 1) ? 2 : 1;
			}
		}

	}
}

void RunSimulation(SimData* sd, ParamStruct* ps) {

	fprintf(stdout,"Initializing simulation...\n");
	
	//Enumerate GPU devices
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, sd->device);

	printf("Using device %d with compute capability %d.%d.\n", sd->device, deviceProp.major, deviceProp.minor);
	cudaSetDevice(sd->device);

	//Declare device data structures
	DeviceDomainStruct dds;
	DeviceGlobalStruct dgs;

	if ( cudaSuccess != cudaGetLastError() )
		fprintf(stderr,"*****************Cuda Error setting up device!\n");
	//Initialize data structures
	fprintf(stdout,"Initializing cytosol domain (%d)...\n",sd->device);
	InitializeDomainStruct(&dds, &sd->Grid, ps->Integers[I_FLAG_SECONDARY_GRADIENT], ps);

	if ( cudaSuccess != cudaGetLastError() )
		fprintf(stderr,"*****************Cuda Error in InitializeDomainStruct()!\n");

	fprintf(stdout,"Initializing global domain (%d)...\n",sd->device);
	InitializeDeviceGlobalStruct(&dgs, sd, ps);

	if ( cudaSuccess != cudaGetLastError() )
		fprintf(stderr,"*****************Cuda Error in InitializeDeviceGlobalStruct()!\n");

	//Initialize device constants
	fprintf(stdout,"Initializing device constants (%d)...\n",sd->device);
	InitializeDeviceConstants(ps);

	if ( cudaSuccess != cudaGetLastError() )
		fprintf(stderr,"*****************Cuda Error in InitializeDeviceConstants!\n");

	//Initialize output files
	CreateOutputFiles(sd->ensemble_index,sd->device,ps);

	//Simulation variables
	double t = 0; //time (ms)
	int output_counter = 0; //index for grid outputs

	//Time step size
	double dt = ps->Reals[I_T_STEP];
	double dt_last = 0;

	fprintf(stdout,"Starting simulation (%d)...\n",sd->device);

	//Output initial state
	printf("Writing grid files (t = %f, #%d).\n",t,output_counter);
	WriteToGridFiles(output_counter, t, &sd->Grid, sd, ps);
	output_counter++;
	printf("Writing state files (t = %f).\n",t);
	WriteToStateFiles(t, sd, ps);

	//Timer
	clock_t clock_start = clock();

	//Keep track of number of open RyRs, terminate when none are open
	int N_RyR_Open = 0;
	for (int i = 0; i < sd->N_RyR; i++) {
		N_RyR_Open += sd->RyR_States[i];
	}

	//Start CUDA profiler
	cudaProfilerStart();

	//Pull in some parameters for readability
	int OUTPUT_GRID_FLAG = ps->Integers[I_OUTPUT_GRID_FLAG];
	int FLAG_FIDELITY_SIM = ps->Integers[I_FLAG_FIDELITY_SIM];
	double T_FINAL = ps->Reals[I_T_FINAL];
	double C_MIN = ps->Reals[I_C_MIN];
	int N_OPEN_MAX = ps->Integers[I_N_RYR_FIDELITY];
	double OUTPUT_GRID_INTERVAL = ps->Reals[I_OUTPUT_GRID_INTERVAL];
	double OUTPUT_STATES_INTERVAL = ps->Reals[I_OUTPUT_STATES_INTERVAL];
	double ca_ss_max = 0;

	//LCC protocol
	int bLCC_Open = ps->Integers[I_FLAG_LCC];
	int seed = 93843 + ps->Integers[I_RAND_SEED_0] + sd->device*ps->Integers[I_RAND_ADD_A] + sd->ensemble_index*ps->Integers[I_RAND_ADD_B];
	if (ps->Integers[I_FLAG_CLOCK_SEED]) {
		seed += (int)time(NULL);
	}
	srand(seed);
	double T_LCC_Close = -ps->Reals[I_LCC_DURATION]*log(((double)rand())/((double)RAND_MAX));
	if (bLCC_Open) {
		fprintf(stdout,"LCC will be closed at %g ms.\n",T_LCC_Close);
	}

	//Channel gating fast-forward option
	//This assumes the system is at steady-state so that Ca2+ concentrations are not changing
	//Simulates only channel gating until a channel opens or T_FINAL is reached
	//Useful for accelerating gain simulations at -40mV where openings are rare
	if (ps->Integers[I_GATING_FFWD]) {
		fprintf(stdout,"Using accelerated steady-state gating...\n");
		int bDone = 0;
		while (!bDone) {

			UpdateChannelGating( dt, &dds, &dgs);
			t += dt;
			dt_last = dt;

			cudaMemcpy(sd->RyR_States, dgs.RyR_States, sizeof(int)*dgs.N_RyR, cudaMemcpyDeviceToHost);
			cudaMemcpy(sd->LCC_States, dgs.LCC_States, sizeof(int)*dgs.N_LCC, cudaMemcpyDeviceToHost);
			cudaMemcpy(sd->LCCV_States, dgs.LCCV_States, sizeof(int)*dgs.N_LCC, cudaMemcpyDeviceToHost);
			cudaDeviceSynchronize();
			N_RyR_Open = 0;
			for (int i = 0; i < sd->N_RyR; i++) {
				N_RyR_Open += sd->RyR_States[i];
			}
			int n_open_lcc = 0;
			for (int i = 0; i < sd->N_LCC; i++) {
				if (sd->LCCV_States[i] == 2 && (sd->LCC_States[i] == 6 || sd->LCC_States[i] == 12)) {
					n_open_lcc++;
				}
			}

			//Check if a channel opened or we've reached T_FINAL
			//fprintf(stdout,"RyR open = %d, LCC_open = %d, t = %f\n",N_RyR_Open,n_open_lcc,t);
			if (N_RyR_Open > 0 || n_open_lcc > 0 || t > T_FINAL) {
				bDone = 1;
			}

			//Output states if not finished
			if (!bDone) {
				//Output grid and states
				if (OUTPUT_GRID_FLAG && floor(t/OUTPUT_GRID_INTERVAL) > floor((t-dt_last)/OUTPUT_GRID_INTERVAL)) {

					//Copy data from device
					cudaDeviceSynchronize();
					for (int i = 0; i < dds.N_States; i++) {
						cudaMemcpy(sd->Grid.States[i],dds.States + (i*dds.N_Ele), sizeof(double)*dds.N_Ele, cudaMemcpyDeviceToHost);
					}
					cudaDeviceSynchronize();
					printf("Writing grid file (gating accelerated)  (t = %f, #%d).\n",t,output_counter);
					WriteToGridFiles(output_counter, t, &sd->Grid, sd, ps);
					output_counter++;
				}
				if (floor(t/OUTPUT_STATES_INTERVAL) > floor((t-dt_last)/OUTPUT_STATES_INTERVAL)) {
					//Copy data from device
					cudaDeviceSynchronize();
					for (int i = 0; i < dds.N_States; i++) { //Note: grid values only needed if calculating total Ca2+
						cudaMemcpy(sd->Grid.States[i],dds.States + (i*dds.N_Ele), sizeof(double)*dds.N_Ele, cudaMemcpyDeviceToHost);
					}
					cudaMemcpy(sd->Global_Vars, dgs.Global_Vars, sizeof(double)*N_GLOBAL_VARS, cudaMemcpyDeviceToHost);
					cudaDeviceSynchronize();
					printf("Writing to states file (gating accelerated) (t = %f, LCC_Open = %d, RyR_open = %d).\n",t,n_open_lcc,N_RyR_Open);
					WriteToStateFiles(t, sd, ps);
				}
			}

			//Voltage clamp protocol
			if (ps->Integers[I_FLAG_V_CLAMP]) {
				if (t <= ps->Reals[I_T_CLAMP]) {
					sd->Global_Vars[GLOBAL_VM] = ps->Reals[I_V_CLAMP];
				} else {
					sd->Global_Vars[GLOBAL_VM] = ps->Reals[I_DEFAULT_VM];
				}
				cudaMemcpy(&dgs.Global_Vars[GLOBAL_VM], &sd->Global_Vars[GLOBAL_VM], sizeof(double), cudaMemcpyHostToDevice);
			}
		}
	}

	//Main time loop
	while (t <= T_FINAL && (!FLAG_FIDELITY_SIM || ( (bLCC_Open || N_RyR_Open > 0 || ca_ss_max > C_MIN) && N_RyR_Open < N_OPEN_MAX ))) {

		//LCC protocol
		if (ps->Integers[I_FLAG_LCC] && bLCC_Open && t >= T_LCC_Close ) {
			bLCC_Open = 0;
			sd->LCC_States[ps->Integers[I_LCC_INDEX]] = DEFAULT_STATE_LCC;
			cudaMemcpy(&dgs.LCC_States[ps->Integers[I_LCC_INDEX]], &sd->LCC_States[ps->Integers[I_LCC_INDEX]], sizeof(int), cudaMemcpyHostToDevice);
		}

		//Voltage clamp protocol
		if (ps->Integers[I_FLAG_V_CLAMP]) {
			if (t <= ps->Reals[I_T_CLAMP]) {
				sd->Global_Vars[GLOBAL_VM] = ps->Reals[I_V_CLAMP];
			} else {
				sd->Global_Vars[GLOBAL_VM] = ps->Reals[I_DEFAULT_VM];
			}
			cudaMemcpy(&dgs.Global_Vars[GLOBAL_VM], &sd->Global_Vars[GLOBAL_VM], sizeof(double), cudaMemcpyHostToDevice);
		}

		//RyR Hold Option

		//RyR Hold Option
		if (ps->Integers[I_INIT_RYR] != -1 && t < ps->Reals[I_INIT_RYR_HOLD]) {
			cudaMemcpy(&sd->RyR_States[ps->Integers[I_INIT_RYR]], &dgs.RyR_States[ps->Integers[I_INIT_RYR]], sizeof(int), cudaMemcpyDeviceToHost);
			if (!sd->RyR_States[ps->Integers[I_INIT_RYR]]) {
			    sd->RyR_States[ps->Integers[I_INIT_RYR]] = 1;
				cudaMemcpy(&dgs.RyR_States[ps->Integers[I_INIT_RYR]], &sd->RyR_States[ps->Integers[I_INIT_RYR]], sizeof(int), cudaMemcpyHostToDevice);
			}
		}

		/*if (t > 2.0) {
			for (int i = 0; i < sd->N_RyR; i++) sd->RyR_States[i] = 0;
			cudaMemcpy(&dgs.RyR_States[ps->Integers[I_INIT_RYR]], &sd->RyR_States[ps->Integers[I_INIT_RYR]], sizeof(int), cudaMemcpyHostToDevice);
		}*/

		//Output states
		if (OUTPUT_GRID_FLAG && floor(t/OUTPUT_GRID_INTERVAL) > floor((t-dt_last)/OUTPUT_GRID_INTERVAL)) {
			
			//Copy data from device
			cudaDeviceSynchronize();
			for (int i = 0; i < dds.N_States; i++) {
				cudaMemcpy(sd->Grid.States[i],dds.States + (i*dds.N_Ele), sizeof(double)*dds.N_Ele, cudaMemcpyDeviceToHost);
			}
			cudaDeviceSynchronize();
			printf("Writing grid file (t = %f, #%d).\n",t,output_counter);
			WriteToGridFiles(output_counter, t, &sd->Grid, sd, ps);
			output_counter++;
		}

		if (floor(t/OUTPUT_STATES_INTERVAL) > floor((t-dt_last)/OUTPUT_STATES_INTERVAL)) {

			//Copy data from device
			cudaDeviceSynchronize();
			for (int i = 0; i < dds.N_States; i++) { //Note: grid values only needed if calculating total Ca2+
				cudaMemcpy(sd->Grid.States[i],dds.States + (i*dds.N_Ele), sizeof(double)*dds.N_Ele, cudaMemcpyDeviceToHost);
			}
			cudaMemcpy(sd->RyR_States, dgs.RyR_States, sizeof(int)*dgs.N_RyR, cudaMemcpyDeviceToHost);
			cudaMemcpy(sd->LCC_States, dgs.LCC_States, sizeof(int)*dgs.N_LCC, cudaMemcpyDeviceToHost);
			cudaMemcpy(sd->LCCV_States, dgs.LCCV_States, sizeof(int)*dgs.N_LCC, cudaMemcpyDeviceToHost);
			cudaMemcpy(sd->Global_Vars, dgs.Global_Vars, sizeof(double)*N_GLOBAL_VARS, cudaMemcpyDeviceToHost);
			cudaDeviceSynchronize();

			N_RyR_Open = 0;
			ca_ss_max = 0;
			for (int i = 0; i < sd->N_RyR; i++) {
				N_RyR_Open += sd->RyR_States[i];
				double c = sd->Grid.States[INDEX_CA][sd->RyR_Ele[i]];
				if (c > ca_ss_max) {
					ca_ss_max = c;
				}
			}

			printf("Writing to states file (t = %f, N_open = %d, Ca_SS_Max = %g).\n",t,N_RyR_Open,ca_ss_max);
			WriteToStateFiles(t, sd, ps);
		}

		//Fast explicit step routine
		ExplicitStep(dt, &dds, &dgs, ps);

		//Increment time step
		t += dt;
		dt_last = dt;

	}

	//Stop CUDA profiler
	cudaProfilerStop();

	//Timer
	clock_t clock_end = clock();

	//Output final state
	cudaDeviceSynchronize();
	for (int i = 0; i < dds.N_States; i++) {
		cudaMemcpy(sd->Grid.States[i],dds.States + (i*dds.N_Ele), sizeof(double)*dds.N_Ele, cudaMemcpyDeviceToHost);
	}
	cudaMemcpy(sd->RyR_States, dgs.RyR_States, sizeof(int)*dgs.N_RyR, cudaMemcpyDeviceToHost);
	cudaMemcpy(sd->LCC_States, dgs.LCC_States, sizeof(int)*dgs.N_LCC, cudaMemcpyDeviceToHost);
	cudaMemcpy(sd->LCCV_States, dgs.LCCV_States, sizeof(int)*dgs.N_LCC, cudaMemcpyDeviceToHost);
	cudaMemcpy(sd->Global_Vars, dgs.Global_Vars, sizeof(double)*N_GLOBAL_VARS, cudaMemcpyDeviceToHost);
	cudaMemcpy(sd->R, dgs.R, sizeof(double)*dgs.N_Channels, cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
	printf("Writing grid files (t = %f, #%d).\n",t,output_counter);
	WriteToGridFiles(output_counter, t, &sd->Grid, sd, ps);
	printf("Writing state files (t = %f).\n",t);
	WriteToStateFiles(t, sd, ps);

	//Write time for simulation
	printf("Simulation completed in %f secs.\n",((float)clock_end-clock_start)/CLOCKS_PER_SEC);

	//Free memory
	FreeDomainStruct(&dds);
	FreeGlobalStruct(&dgs);
	cudaDeviceSynchronize();


}

void ExplicitStep(double dt, DeviceDomainStruct* dds, DeviceGlobalStruct* dgs, ParamStruct* ps) {

	int N_Ele = dds->N_Ele;
	int block_size = 256;
	int n_blocks_ele = N_Ele/block_size + (N_Ele%block_size == 0?0:1);
	int n_blocks_addmult = ((N_Ele*dds->N_States)/block_size) + ((N_Ele*dds->N_States)%block_size == 0?0:1);
	int block_size_gv = 8;
	int n_blocks_gv = N_GLOBAL_VARS/block_size_gv + (N_GLOBAL_VARS%block_size_gv == 0?0:1);

	//Initialize global variable delta vector to zero
	resetVectorToZero<<< n_blocks_gv, block_size_gv >>> (dgs->GV_src1, N_GLOBAL_VARS);

	//Diffusion step
	/*for (int i = 0; i < dds->N_States; i++) {
		if (ps->Diff[i] > 0) { //Diffusible states

			//Know this is synchronous as long as the block size is >64KB
			cudaMemcpy(thrust::raw_pointer_cast(&dds->Array_2[0]), &dds->States[i*N_Ele], sizeof(double)*N_Ele, cudaMemcpyDeviceToDevice);

			cusp::multiply(dds->A[i],dds->Array_2,dds->Array_1[i]);

		} else { //Non-diffusible states
			resetVectorToZero<<< n_blocks_ele, block_size >>> (thrust::raw_pointer_cast(&dds->Array_1[i][0]), N_Ele);
		}

	}*/

	/*for (int i = 0; i < dds->N_States; i++) {
		resetVectorToZero<<< n_blocks_ele, block_size >>> (thrust::raw_pointer_cast(&dds->Array_1[i][0]), N_Ele);
	}*/

	//Buffering reactions
	UpdateReactionTerms<<< n_blocks_ele, block_size >>>(N_Ele,
													dds->States,
													dds->Array_1,
													dgs->Global_Vars,
													dds->aij,
													dds->nij,
													dds->TropC, dds->V0, dds->TTSurfaceArea, dds->SRSurfaceArea, dds->Boundary_Rates,
													dds->bJSR, dds->Vol_JSR);

	//Update Channel Ca fluxes
	UpdateRyRFlux<<< 1, dgs->N_RyR >>>(&dds->States[INDEX_CA*N_Ele],
						&dds->Array_1[INDEX_CA*N_Ele],
						dgs->Global_Vars,
						dds->V0,
						dgs->N_RyR, dgs->RyR_States, dgs->RyR_Ele,
						dgs->RyR_JSR_Ele);
	UpdateLCCFlux<<< 1, dgs->N_LCC >>>(&dds->States[INDEX_CA*N_Ele],
						&dds->Array_1[INDEX_CA*N_Ele],
						dgs->Global_Vars, dds->V0,
						dgs->N_LCC, dgs->LCC_States, dgs->LCCV_States, dgs->LCC_Ele);

	//Fast buffering (must be performed last)
	UpdateFastBuffer<<< n_blocks_ele, block_size >>>(N_Ele,&dds->States[INDEX_CA*N_Ele], &dds->Array_1[INDEX_CA*N_Ele], dds->TTSurfaceArea,dds->V0);

	//Update new states
	/*for (int k = 0; k < dds->N_States; k++) {
		addMultVectors<<< n_blocks_ele, block_size >>> (&dds->States[k*N_Ele], &dds->Array_1[k*N_Ele], dt, N_Ele);
	}*/
	addMultVectors<<< n_blocks_addmult, block_size >>> (dds->States, dds->Array_1, dt, N_Ele*dds->N_States);
	addMultVectors<<< n_blocks_gv, block_size_gv >>> (dgs->Global_Vars, dgs->GV_src1, dt, N_GLOBAL_VARS);

	//Do channel gating
	if (!ps->Integers[I_FLAG_NO_GATING]) {
		UpdateChannelGating( dt, dds, dgs);
	}

	if ( cudaSuccess != cudaGetLastError() )
		fprintf(stderr,"*****************Cuda Error during step computations!\n");

}

void UpdateChannelGating(double dt, DeviceDomainStruct* dds, DeviceGlobalStruct* dgs) {

	//Do channel gating
	int N_Ele = dds->N_Ele;
	UpdateRyRStates <<< 1, dgs->N_RyR >>> ( dt, dgs->N_RyR, dgs->RyR_States,
											&dds->States[INDEX_CA*N_Ele], dgs->RyR_Ele,
											dgs->randStates, dgs->R, dgs->Global_Vars, dgs->RyR_Neighb,
											dgs->RyR_JSR_Ele);
	UpdateLCCStates <<< 1, dgs->N_LCC >>> ( dt, dgs->N_LCC, dgs->N_RyR, dgs->LCC_States, dgs->LCCV_States,
											dgs->LCC_Ele, &dds->States[INDEX_CA*N_Ele],
											dgs->randStates, dgs->R, dgs->Global_Vars);

}


void InitializeDeviceConstants(ParamStruct* ps) {

	cudaMemcpyToSymbol(PARAMS_DEVICE,(void*)ps->Reals,sizeof(double)*N_PARAMETERS_REAL,0,cudaMemcpyHostToDevice);
}

void FreeDomainStruct(DeviceDomainStruct* dds) {

	//delete[] dds->A;
	cudaFree(dds->States);
	cudaFree(dds->Boundary_Rates);
	cudaFree(dds->TropC);
	cudaFree(dds->bJSR);
	cudaFree(dds->V0);
	cudaFree(dds->TTSurfaceArea);
	cudaFree(dds->SRSurfaceArea);

}

void FreeGlobalStruct(DeviceGlobalStruct* dgs) {

	cudaFree(dgs->Global_Vars);
	cudaFree(dgs->GV_src1);
	cudaFree(dgs->RyR_States);
	cudaFree(dgs->RyR_Ele);
	cudaFree(dgs->RyR_JSR_Ele);
	cudaFree(dgs->RyR_Neighb);
	cudaFree(dgs->LCC_States);
	cudaFree(dgs->LCCV_States);
	cudaFree(dgs->LCC_Ele);
	cudaFree(dgs->R);
	cudaFree(dgs->randStates);

}
