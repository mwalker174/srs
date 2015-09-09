/*
 * main.cpp
 *
 *  Created on: Jun 27, 2012
 *      Author: mwalker
 */

#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <string.h>
#include <sys/stat.h>
#include "parameters.h"

using namespace std;

int main(int argc, char **argv) {

	SimData sd;
	ParamStruct ps;
	char* parameters_file = (char*)DEFAULT_PARAMETER_FILE;
	int myrank = 0; //Ported from MPI code
	
	if (myrank == 0) {
	
		time_t timer;
		time(&timer);
		
  		struct tm y2k;
  		y2k.tm_hour = 0;   y2k.tm_min = 0; y2k.tm_sec = 0;
  		y2k.tm_year = 100; y2k.tm_mon = 0; y2k.tm_mday = 1;
  		
  		int seconds = (int)difftime(timer,mktime(&y2k));
		srand(seconds);
		int id = rand() % 1000;
		
		int id_length = sprintf(ps.rand_id, "%d_%d", seconds, id);
		
	}

	//Set default simulation identifiers
	int idx0 = 0;
	char* output_dir = 0;
	char* mesh_base = 0;

	//Read command line options
	for (int i = 1; i < argc; i += 2) {
		if (!strcmp(argv[i],"-rank")) {
			myrank = atoi(argv[i+1]); //Override rank
		} else if (!strcmp(argv[i],"-start")) {
			idx0 = atoi(argv[i+1]); //Starting simulation index (use for resuming ensemble simulations)
		} else if (!strcmp(argv[i],"-param")) {
			parameters_file = argv[i+1]; //Specify parameter file name
		} else if (!strcmp(argv[i],"-out")) {
			output_dir = argv[i+1]; //Specify output directory
		} else if (!strcmp(argv[i],"-mesh")) {
			mesh_base = argv[i+1]; //Specify mesh location
		} else {
			fprintf(stderr,"Error: unrecognized command line option: %s\n",argv[i]);
			exit(1);
		}
	}

	//Read in parameters from file
	ReadParameters(&ps,parameters_file);
	
	//Set device
	sd.device = ps.Integers[I_CUDA_DEVICE];

	if (output_dir != 0) {
		sprintf(ps.Chars[I_OUTPUT_STATE],"%s",output_dir);
		sprintf(ps.Chars[I_OUTPUT_GRID],"%s",output_dir);
	} else {
		//If output directory not specified, create time-stamped directory
		sprintf(ps.Chars[I_OUTPUT_STATE],"%s/%s_%s",ps.Chars[I_OUTPUT_STATE],ps.Chars[I_PARAM_TITLE],ps.rand_id);
		sprintf(ps.Chars[I_OUTPUT_GRID],"%s/%s_%s",ps.Chars[I_OUTPUT_GRID],ps.Chars[I_PARAM_TITLE],ps.rand_id);
	}
	fprintf(stdout,"Output will be written to %s\n",ps.Chars[I_OUTPUT_STATE]);

	if (mesh_base != 0) {
		sprintf(ps.Chars[I_FILE_BASE],"%s",mesh_base);
		fprintf(stdout,"Reading mesh from %s\n",ps.Chars[I_FILE_BASE]);
	}

	//Let master process handle error handling and output directory setup
	if (myrank == 0) {

		//Create directory
	    if (output_dir == 0) {
			char newdir[255];
			sprintf(newdir,"%s",ps.Chars[I_OUTPUT_STATE]);
			mkdir(newdir,S_IRWXU);
	    }

	    //Get parameters file name
	    char* parameters_filename = strrchr(parameters_file,'/');
	    if (parameters_filename == NULL) {
	    	parameters_filename = strrchr(parameters_file,'\\');
	    }
	    if (parameters_filename == NULL) {
	    	parameters_filename = parameters_file;
	    }

		//Copy parameters file
	    char buf[BUFSIZ];
	    char dst_paramfile[1024];
	    if (output_dir == 0) { //If output directory not specified, add time-stamped directory to path
	    	sprintf(dst_paramfile,"%s/%s",ps.Chars[I_OUTPUT_GRID],parameters_filename);
	    } else {
	    	sprintf(dst_paramfile,"%s/%s",ps.Chars[I_OUTPUT_GRID],parameters_filename);
	    }
	    size_t size;

	    FILE* source = fopen(parameters_file, "rb");
	    FILE* dest = fopen(dst_paramfile, "wb");

	    if (source == NULL) {
	    	fprintf(stderr,"Error: could not open parameters file %s for copying.\n",parameters_file);
	    	exit(1);
	    }
	    if (dest == NULL) {
	    	fprintf(stderr,"Error: could not open destination file %s for copying.\n",dst_paramfile);
	    	exit(1);
	    }

	    while (size = fread(buf, 1, BUFSIZ, source)) {
	        fwrite(buf, 1, size, dest);
	    }

	    fclose(source);
	    fclose(dest);
	}

	//Read in grid properties from files
	ReadGrid(&sd,&sd.Grid,&ps);

	//Initialize global state variables
	fprintf(stdout,"Initializing default global variables (%d)...\n",myrank);
	ResetGlobals(&sd,&ps);

	if (ps.Integers[I_INITIALIZE_FROM_FILE]) {
		fprintf(stdout,"Initializing grid states from files (%d)...\n",myrank);
		ReadGridStates(&sd.Grid,&ps);
	} else {
		fprintf(stdout,"Initializing default grid states (%d)...\n",myrank);
		InitializeDefaultGrid(&sd,&ps);
	}

	//Initialize default channel states and allocate memory for residuals
	fprintf(stdout,"Initializing default channel states (%d)...\n",myrank);
	InitializeDefaultChannels(&sd);

	for (int i = idx0; i < ps.Integers[I_N_SIMS]; i++) {

		fprintf(stdout,"Beginning ensemble simulation %d of %d (rank %d).\n",i+1,ps.Integers[I_N_SIMS],myrank);

		//Set ensemble index
		sd.ensemble_index = i;

		//Channel state initialization
		if (ps.Integers[I_FLAG_LCC]) {

			//Open LCC
			fprintf(stdout,"Opening LCC %d...\n",ps.Integers[I_LCC_INDEX]);
			sd.LCC_States[ps.Integers[I_LCC_INDEX]] = STATE_LCC_OPEN_1;

		} else if (ps.Integers[I_FLAG_V_CLAMP]) {

			//Do not open any channels for voltage clamp

		} else {

			//Open RyRs initially
			if (ps.Integers[I_INIT_RYR] != -1) {
			
				if (ps.Integers[I_INIT_RYR] >= 0 && ps.Integers[I_INIT_RYR] < sd.N_RyR) {

					//Open single specified RyR
					sd.RyR_States[ps.Integers[I_INIT_RYR]] = 1;
					
				} else {
					fprintf(stderr,"Error: Invalid initial RyR index (%d), must be -1 to ignore, or between 0 and the number of RyRs (%d).",ps.Integers[I_INIT_RYR],sd.N_RyR);
					exit(1);
				}

			} else {

				//Default: Open random RyR
				OpenRandomRyR(&ps, &sd);

			}
		}

		//Run simulation
		RunSimulation(&sd,&ps);
		//Reset states
		ResetDefaultGridStates(&sd,&ps);
		ResetDefaultChannels(&sd);
		ResetGlobals(&sd,&ps);
	}

	//Free memory
	FreeSimData(&sd);

	return 0;
}

//Cleanup function
void FreeSimData(SimData* sd) {

	FreeGridStruct(&sd->Grid);

	delete[] sd->LCC_States;
	delete[] sd->LCCV_States;
	delete[] sd->LCC_Ele;
	delete[] sd->RyR_States;
	delete[] sd->RyR_Ele;
	delete[] sd->RyR_JSR_Ele;
	delete[] sd->RyR_Neighb;
	delete[] sd->R;
}

//Helper cleanup function
void FreeGridStruct(GridStruct* gs) {

	for (int i = 0; i < gs->N_States; i++) {
		delete[] gs->States[i];
	}
	delete[] gs->States;

	delete[] gs->Boundaries;

	delete[] gs->Domain;
	delete[] gs->V0;
	delete[] gs->TTSurfaceArea;
	delete[] gs->SRSurfaceArea;
	delete[] gs->TropC;
	delete[] gs->Ele_1;
	delete[] gs->Ele_2;
	delete[] gs->FaceType;
	delete[] gs->Face_C;
	delete[] gs->dXi;
	delete[] gs->Area_f;
	delete[] gs->Afx;
	delete[] gs->Afy;
	delete[] gs->Afz;
	delete[] gs->eXix;
	delete[] gs->eXiy;
	delete[] gs->eXiz;
	delete[] gs->r1x;
	delete[] gs->r1y;
	delete[] gs->r1z;
	delete[] gs->r2x;
	delete[] gs->r2y;
	delete[] gs->r2z;

	delete[] gs->Node_X;
	delete[] gs->Node_Y;
	delete[] gs->Node_Z;
	delete[] gs->Ele_Nodes;
}
