/*
 * OutputFuncs.cpp
 *
 *  Created on: Jul 5, 2012
 *      Author: mwalker
 */

#include <iostream>
#include <stdio.h>
#include "parameters.h"

//Creates global variable output file
void CreateOutputFiles(int run_nmbr, int rank, ParamStruct* ps) {

	sprintf(ps->state_file,"%s/states_0_%d.out",ps->Chars[I_OUTPUT_STATE],run_nmbr);
	sprintf(ps->grid_file,"%s/grid",ps->Chars[I_OUTPUT_GRID]);

}

void WriteToGridFiles(int n, double t, GridStruct* gs, SimData* sd, ParamStruct* ps) {

	const char* filename = ps->Chars[I_OUTPUT_GRID];
	int run_nmbr = sd->ensemble_index;
	int rank = sd->device;
	const int buffer_size = (int)10e6; //10MB buffer

	FILE * pFile;
	if (ps->Integers[I_OUTPUT_GRID_FLAG] || (rank == 0 && run_nmbr == 0 && n == 0)) {
		char grid_file[255];
		sprintf(grid_file,"%s_%d_%d_%d.vtk",ps->grid_file,rank,run_nmbr,n);

		//Delete old grid file if it exists
		//remove(grid_file);

		pFile = fopen (grid_file,"wb");
		if (pFile != NULL) {
			
			//Buffer
			char* buffer = new char[buffer_size];
			int ptr = 0;

			//VTK format
			ptr += sprintf(buffer,"# vtk DataFile Version 3.0\n");
			ptr += sprintf(&buffer[ptr],"Grid data for run %d-%d-%d at t = %f\n",run_nmbr,rank,n,t);
			ptr += sprintf(&buffer[ptr],"ASCII\n\n");
			ptr += sprintf(&buffer[ptr],"DATASET UNSTRUCTURED_GRID\n");
			ptr += sprintf(&buffer[ptr],"POINTS %d float\n",gs->N_Nodes);
			for (int j = 0; j < gs->N_Nodes; j++) {
				ptr += sprintf(&buffer[ptr],"%f %f %f\n",gs->Node_X[j],gs->Node_Y[j],gs->Node_Z[j]);
			}
			ptr += sprintf(&buffer[ptr],"\nCELLS %d %d\n",gs->N_Ele,5*gs->N_Ele);
			for (int j = 0; j < gs->N_Ele; j++) {
				int idx = j*4;
				ptr += sprintf(&buffer[ptr],"4 %d %d %d %d\n",gs->Ele_Nodes[idx],gs->Ele_Nodes[idx+1],gs->Ele_Nodes[idx+2],gs->Ele_Nodes[idx+3]);
			}
			ptr += sprintf(&buffer[ptr],"\nCELL_TYPES %d\n",gs->N_Ele);
			for (int j = 0; j < gs->N_Ele; j++) {
				ptr += sprintf(&buffer[ptr],"10\n");
			}

			ptr += sprintf(&buffer[ptr],"\nCELL_DATA %d\n",gs->N_Ele);

			ptr += sprintf(&buffer[ptr],"SCALARS C_CA float\n");
			ptr += sprintf(&buffer[ptr],"LOOKUP_TABLE default\n");
			for (int j = 0; j < gs->N_Ele; j++) {
				ptr += sprintf(&buffer[ptr],"%f\n",gs->States[INDEX_CA][j]);
			}
			/*
			ptr += sprintf(&buffer[ptr],"SCALARS C_ATP float\n");
			ptr += sprintf(&buffer[ptr],"LOOKUP_TABLE default\n");
			for (int j = 0; j < gs->N_Ele; j++) {
				ptr += sprintf(&buffer[ptr],"%f\n",gs->States[INDEX_ATP][j]);
			}
			ptr += sprintf(&buffer[ptr],"SCALARS C_CMDN float\n");
			ptr += sprintf(&buffer[ptr],"LOOKUP_TABLE default\n");
			for (int j = 0; j < gs->N_Ele; j++) {
				ptr += sprintf(&buffer[ptr],"%f\n",gs->States[INDEX_CMDN][j]);
			}
			ptr += sprintf(&buffer[ptr],"SCALARS C_TRPN float\n");
			ptr += sprintf(&buffer[ptr],"LOOKUP_TABLE default\n");
			for (int j = 0; j < gs->N_Ele; j++) {
				ptr += sprintf(&buffer[ptr],"%f\n",gs->States[INDEX_TRPN][j]);
			}
			ptr += sprintf(&buffer[ptr],"SCALARS C_CSQN float\n");
			ptr += sprintf(&buffer[ptr],"LOOKUP_TABLE default\n");
			for (int j = 0; j < gs->N_Ele; j++) {
				ptr += sprintf(&buffer[ptr],"%f\n",gs->States[INDEX_CSQN][j]);
			}*/
			ptr += sprintf(&buffer[ptr],"SCALARS C_DYE float\n");
			ptr += sprintf(&buffer[ptr],"LOOKUP_TABLE default\n");
			for (int j = 0; j < gs->N_Ele; j++) {
				ptr += sprintf(&buffer[ptr],"%f\n",gs->States[INDEX_DYE][j]);
			}
			/*ptr += sprintf(&buffer[ptr],"SCALARS C_SL_LO float\n");
			ptr += sprintf(&buffer[ptr],"LOOKUP_TABLE default\n");
			for (int j = 0; j < gs->N_Ele; j++) {
				double c = ps->Reals[I_B_TOT_SL_LO]*gs->TTSurfaceArea[j]*gs->States[INDEX_CA][j]/(ps->Reals[I_K_D_SL_LO]+gs->States[INDEX_CA][j]);
				ptr += sprintf(&buffer[ptr],"%g\n",c);
			}
			ptr += sprintf(&buffer[ptr],"SCALARS C_SL_HI float\n");
			ptr += sprintf(&buffer[ptr],"LOOKUP_TABLE default\n");
			for (int j = 0; j < gs->N_Ele; j++) {
				double c = ps->Reals[I_B_TOT_SL_HI]*gs->TTSurfaceArea[j]*gs->States[INDEX_CA][j]/(ps->Reals[I_K_D_SL_HI]+gs->States[INDEX_CA][j]);
				ptr += sprintf(&buffer[ptr],"%g\n",c);
			}
*/
			ptr += sprintf(&buffer[ptr],"SCALARS B_RYR int\n");
			ptr += sprintf(&buffer[ptr],"LOOKUP_TABLE default\n");
			int* bRyR = new int[gs->N_Ele];
			for (int i = 0; i < gs->N_Ele; i++) {
				bRyR[i] = 0;
			}
			for (int i = 0; i < sd->N_RyR; i++) {
				bRyR[sd->RyR_Ele[i]] = 1;
				if (sd->RyR_States[i]) {
					bRyR[sd->RyR_Ele[i]] = 2;
				}
			}
			for (int j = 0; j < gs->N_Ele; j++) {
				ptr += sprintf(&buffer[ptr],"%d\n",bRyR[j]);
			}
			delete[] bRyR;

			ptr += sprintf(&buffer[ptr],"SCALARS B_JSR int\n");
			ptr += sprintf(&buffer[ptr],"LOOKUP_TABLE default\n");
			for (int j = 0; j < gs->N_Ele; j++) {
				ptr += sprintf(&buffer[ptr],"%d\n",gs->Domain[j]==DOMAIN_JSR ? 1 : 0);
			}
/*
			ptr += sprintf(&buffer[ptr],"SCALARS B_RYR_JSR int\n");
			ptr += sprintf(&buffer[ptr],"LOOKUP_TABLE default\n");
			int* bRyRJSR = new int[gs->N_Ele];
			for (int i = 0; i < gs->N_Ele; i++) {
				bRyRJSR[i] = 0;
			}
			for (int i = 0; i < sd->N_RyR; i++) {
				bRyRJSR[sd->RyR_JSR_Ele[i]] = 1;
			}
			for (int j = 0; j < gs->N_Ele; j++) {
				ptr += sprintf(&buffer[ptr],"%d\n",bRyRJSR[j]);
			}
			delete[] bRyRJSR;
*/
			ptr += sprintf(&buffer[ptr],"SCALARS B_LCC int\n");
			ptr += sprintf(&buffer[ptr],"LOOKUP_TABLE default\n");
			int* bLCC = new int[gs->N_Ele];
			for (int i = 0; i < gs->N_Ele; i++) {
				bLCC[i] = 0;
			}
			for (int i = 0; i < sd->N_LCC; i++) {
				bLCC[sd->LCC_Ele[i]] = 1;
			}
			for (int j = 0; j < gs->N_Ele; j++) {
				ptr += sprintf(&buffer[ptr],"%d\n",bLCC[j]);
			}
			delete[] bLCC;

			if (ptr > buffer_size) {
				fprintf(stderr,"Error: buffer overrun writing grid file.\n");
				exit(1);
			}

		    fwrite(buffer, sizeof(char), ptr, pFile);

			delete[] buffer;

			fclose(pFile);
		} else {
			fprintf(stderr,"Warning: could not open grid file %s.\n", filename);
		}
	}

}


void WriteToStateFiles(double t, SimData* sd, ParamStruct* ps) {

	int run_nmbr = sd->ensemble_index;
	int rank = sd->device;
	const int buffer_size = 1024*64; //64kb buffer
	const int buffer_real_size = 255;

	FILE * pFile;
	char state_file[255];
	pFile = fopen (ps->state_file,"ab");
	if (pFile != NULL) {

		

		//int ptr = 0;
		//char* buffer = new char[buffer_size];
		int ptr_real = 0;
		float* buf_real = new float[buffer_real_size];

		//Output time
		//ptr += sprintf(&buffer[ptr],"%g ",t);
		buf_real[ptr_real] = t;
		ptr_real++;

		//Compartment Ca2+ concentrations
		//ptr += sprintf(&buffer[ptr],"%g ",sd->Global_Vars[GLOBAL_CA_NSR]);
		//buf_real[ptr_real] = (float)sd->Global_Vars[GLOBAL_CA_NSR];
		//ptr_real++;

		//Total RyR current in pA
		double i_ryr = 0;
		for (int i = 0; i < sd->N_RyR; i++) {
			double CaJSR = sd->Grid.States[INDEX_CA][sd->RyR_JSR_Ele[i]];
			double CaSS = sd->Grid.States[INDEX_CA][sd->RyR_Ele[i]];
			i_ryr += sd->RyR_States[i] * ps->Reals[I_V_RYR]*(CaJSR-CaSS)*(ps->Reals[I_V_CELL])*96485*2*(1e-6); //pA (Using Williams model)
		}
		//ptr += sprintf(&buffer[ptr],"%f ",i_ryr);
		buf_real[ptr_real] = (float)i_ryr;
		ptr_real++;


		//Number of open RyRs
		int N_ryr_open = 0;
		for (int i = 0; i < sd->N_RyR; i++) {
			N_ryr_open += sd->RyR_States[i];
		}
		//ptr += sprintf(&buffer[ptr],"%d ",N_ryr_open);
		buf_real[ptr_real] = (float)N_ryr_open;
		ptr_real++;


		//Number of open LCCs
		int n_open_lcc = 0;
		for (int i = 0; i < sd->N_LCC; i++) {
			if (sd->LCCV_States[i] == 2 && (sd->LCC_States[i] == 6 || sd->LCC_States[i] == 12)) {
				n_open_lcc++;
			}
		}
		//ptr += sprintf(&buffer[ptr],"%d ",n_open_lcc);
		buf_real[ptr_real] = (float)n_open_lcc;
		ptr_real++;

		//Total LCC current
		const double P_CaL = (9.13e-13)*(1e12)*(1e-3); //um^3 / ms
		const double F = 96.5; // C / mmol
		const double RT = 310 * 8.314; // J / mol
		const double F_over_RT = F/RT;
		double V = sd->Global_Vars[GLOBAL_VM];
		if (fabs(V) < 1e-6) {
			V = 1e-6;
		}
		double i_lcc = 0;
		for (int i = 0; i < sd->N_LCC; i++) {
			double CA_DYAD = sd->Grid.States[INDEX_CA][sd->LCC_Ele[i]];

			if ((sd->LCC_States[i]==STATE_LCC_OPEN_1||sd->LCC_States[i]==STATE_LCC_OPEN_2)&&(sd->LCCV_States[i]==STATE_LCCV_OPEN)) {
				i_lcc += - P_CaL * V * (2*F_over_RT) * (CA_DYAD*exp(2*V*F_over_RT) - 0.34*ps->Reals[I_CA_0]) / (exp(2*V*F_over_RT)-1)*96485*2*(1e-6);
			}
		}
		//ptr += sprintf(&buffer[ptr],"%f ",i_lcc);
		buf_real[ptr_real] = (float)i_lcc;
		ptr_real++;

		//Mean RyR Ca_SR
		double mean_ca_sr = 0;
		for (int i = 0; i < sd->N_RyR; i++) {
			mean_ca_sr += sd->Grid.States[INDEX_CA][sd->RyR_JSR_Ele[i]];
		}
		mean_ca_sr /= sd->N_RyR;
		buf_real[ptr_real] = (float)mean_ca_sr;
		ptr_real++;

		//Mean RyR Ca_cyto
		double mean_ca_i = 0;
		for (int i = 0; i < sd->N_RyR; i++) {
			mean_ca_i += sd->Grid.States[INDEX_CA][sd->RyR_Ele[i]];
		}
		mean_ca_i /= sd->N_RyR;
		buf_real[ptr_real] = (float)mean_ca_i;
		ptr_real++;

		
		//RyR Ca2+ levels
		for (int i = 0; i < sd->N_RyR; i++) {
			//ptr += sprintf(&buffer[ptr],"%g ",sd->Grid.States[INDEX_CA][sd->RyR_Ele[i]]);
			buf_real[ptr_real] = (float)sd->Grid.States[INDEX_CA][sd->RyR_Ele[i]];
			ptr_real++;
		}
		
		/*
		//RyR JSR Ca2+ levels
		for (int i = 0; i < sd->N_RyR; i++) {
			//ptr += sprintf(&buffer[ptr],"%g ",sd->Grid.States[INDEX_CA][sd->RyR_JSR_Ele[i]]);
			buf_real[ptr_real] = (float)sd->Grid.States[INDEX_CA][sd->RyR_JSR_Ele[i]];
			ptr_real++;
		}
		 */

		if (ps->Integers[I_FLAG_RYR_OUT]) {
			//RyR states
			for (int i = 0; i < sd->N_RyR; i++) {
				//ptr += sprintf(&buffer[ptr],"%d ",sd->RyR_States[i]);
				buf_real[ptr_real] = (float)sd->RyR_States[i];
				ptr_real++;
			}
		}


		/*
		//RyR rates
		double* rate_Open = new double[sd->N_RyR];
		double* rate_Close = new double[sd->N_RyR];
		double* rate_X_Open = new double[sd->N_RyR];
		double* rate_X_Close = new double[sd->N_RyR];
		double* rate_phi = new double[sd->N_RyR];
		double* rate_Ca = new double[sd->N_RyR];
		for (int i = 0; i < sd->N_RyR; i++) {
			int N_Open = 0;
			int N_Closed = 0;
			for (int j = 0; j < 4; j++) {
				if (sd->RyR_Neighb[i*4 + j] >= 0) {
					if (sd->RyR_States[sd->RyR_Neighb[i*4 + j]]) {
						N_Open++;
					} else {
						N_Closed++;
					}
				}
			}

			double X, rate;
			double ca_ss = sd->Grid.States[INDEX_CA][sd->RyR_Ele[i]];

			//Open -> Closed rate
			X = exp( 0.5*ps->Reals[I_RYR_A_STAR]*(N_Open*ps->Reals[I_RYR_EPS_OO] - N_Closed*ps->Reals[I_RYR_EPS_CC]) );
			rate = X*ps->Reals[I_RYR_K_MINUS];
			rate_Close[i] = rate;
			rate_X_Close[i] = X;

			//Closed -> Open rate
			double phi = ps->Reals[I_RYR_PHI_M]*sd->Grid.States[INDEX_CA][sd->RyR_JSR_Ele[i]] + ps->Reals[I_RYR_PHI_B];
			X = exp( 0.5*ps->Reals[I_RYR_A_STAR]*(N_Closed*ps->Reals[I_RYR_EPS_CC] - N_Open*ps->Reals[I_RYR_EPS_OO]) );
			rate = X*phi*ps->Reals[I_RYR_K_PLUS]* pow(ca_ss,ps->Reals[I_RYR_ETA]);
			rate_Open[i] = rate;
			rate_X_Open[i] = X;
			rate_phi[i] = phi;
			rate_Ca[i] = pow(ca_ss,ps->Reals[I_RYR_ETA]);
		}
		for (int i = 0; i < sd->N_RyR; i++) {
			ptr += sprintf(&buffer[ptr],"%f ",rate_Open[i]);
		}
		for (int i = 0; i < sd->N_RyR; i++) {
			ptr += sprintf(&buffer[ptr],"%f ",rate_Close[i]);
		}
		for (int i = 0; i < sd->N_RyR; i++) {
			ptr += sprintf(&buffer[ptr],"%f ",rate_X_Open[i]);
		}
		for (int i = 0; i < sd->N_RyR; i++) {
			ptr += sprintf(&buffer[ptr],"%f ",rate_X_Close[i]);
		}
		for (int i = 0; i < sd->N_RyR; i++) {
			ptr += sprintf(&buffer[ptr],"%f ",rate_phi[i]);
		}
		for (int i = 0; i < sd->N_RyR; i++) {
			ptr += sprintf(&buffer[ptr],"%f ",rate_Ca[i]);
		}
		delete[] rate_Open;
		delete[] rate_Close;
		delete[] rate_X_Open;
		delete[] rate_X_Close;
		delete[] rate_phi;
		delete[] rate_Ca;
*/
/*
		//LCC states
		for (int i = 0; i < sd->N_LCC; i++) {
			ptr += sprintf(&buffer[ptr],"%d ",sd->LCC_States[i]);
			buf_real[ptr_real] = (double)sd->LCC_States[i];
			ptr_real++;
		}

		//LCCV states
		for (int i = 0; i < sd->N_LCC; i++) {
			ptr += sprintf(&buffer[ptr],"%d ",sd->LCCV_States[i]);
			buf_real[ptr_real] = (double)sd->LCCV_States[i];
			ptr_real++;
		}

		//LCC Ca2+ levels
		for (int i = 0; i < sd->N_LCC; i++) {
			ptr += sprintf(&buffer[ptr],"%g ",sd->Grid.States[INDEX_CA][sd->LCC_Ele[i]]);
			buf_real[ptr_real] = sd->Grid.States[INDEX_CA][sd->LCC_Ele[i]];
			ptr_real++;
		}
*/
		//ptr += sprintf(&buffer[ptr],"\n");

		if (ptr_real > buffer_real_size) {
			fprintf(stderr,"Error: buffer overrun writing states file.\n");
			exit(1);
		}

		fwrite(buf_real, sizeof(float), ptr_real, pFile);

		//delete[] buffer;
		delete[] buf_real;

		fclose(pFile);
	} else {
		fprintf(stderr,"Warning: could not open global state file.\n");
	}
}
