/*
 * InitFuncs.cpp
 *
 *  Created on: Jun 27, 2012
 *      Author: mwalker
 */

#include <stdio.h>
#include <iostream>
#include <map>
#include <vector>
#include <stdexcept>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/sax/HandlerBase.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include "DeviceStructs.h"

using namespace std;
using namespace xercesc;

//Reads parameters
void ReadParameters(ParamStruct* ps, const char* param_file) {

	//Assign default parameters
	ps->Integers[I_OUTPUT_GRID_FLAG] = 0;
	ps->Integers[I_FLAG_FIDELITY_SIM] = 0;
	ps->Integers[I_FLAG_CLOCK_SEED] = 1;
	ps->Integers[I_RAND_SEED_0] = 38376427;
	ps->Integers[I_RAND_ADD_A] = 12783;
	ps->Integers[I_RAND_ADD_B] = 83;
	ps->Integers[I_INITIALIZE_FROM_FILE] = 0;
	ps->Integers[I_CUDA_DEVICE] = 0;
	ps->Integers[I_FLAG_SECONDARY_GRADIENT] = 0;
	ps->Integers[I_FLAG_NO_FLUX_BOUNDARY] = 1;
	ps->Integers[I_N_RYR_FIDELITY] = 1000;
	ps->Integers[I_FLAG_LCC] = 0;
	ps->Integers[I_LCC_INDEX] = 0;
	ps->Integers[I_FLAG_V_CLAMP] = 0;
	ps->Integers[I_FLAG_NO_GATING] = 0;
	ps->Integers[I_GATING_FFWD] = 0;
	ps->Integers[I_FLAG_RYR_OUT] = 0;
	ps->Integers[I_N_SIMS] = 1;
	ps->Integers[I_INIT_RYR] = -1;

	ps->Reals[I_T_FINAL] = 500;
	ps->Reals[I_T_STEP] = 12e-6;
	ps->Reals[I_OUTPUT_STATES_INTERVAL] = 0.01;
	ps->Reals[I_OUTPUT_GRID_INTERVAL] = 1;
	ps->Reals[I_C_MIN] = 1.0;
	ps->Reals[I_D_CA] = 0.25;
	ps->Reals[I_D_CA_JSR] = 0.25;
	ps->Reals[I_D_ATP] = 0.14;
	ps->Reals[I_D_CMDN] = 0.025;
	ps->Reals[I_D_TRPN] = 0;
	ps->Reals[I_D_CSQN] = 0;
	ps->Reals[I_D_DYE] = 0.042;
	ps->Reals[I_D_EGTA] = 0;
	ps->Reals[I_C_0_CA] = 0.1;
	ps->Reals[I_C_0_CA_JSR] = 1000;
	ps->Reals[I_DEFAULT_VM] = -80;
	ps->Reals[I_DEFAULT_CA_NSR] = 1000;
	ps->Reals[I_CA_0] = 2000;
	ps->Reals[I_B_TOT_ATP] = 455;
	ps->Reals[I_K_OFF_ATP] = 45;
	ps->Reals[I_K_ON_ATP] = 0.225;
	ps->Reals[I_B_TOT_CMDN] = 24;
	ps->Reals[I_K_OFF_CMDN] = 0.238;
	ps->Reals[I_K_ON_CMDN] = 0.034;
	ps->Reals[I_B_TOT_TRPN] = 70;
	ps->Reals[I_K_OFF_TRPN] = 0.02;
	ps->Reals[I_K_ON_TRPN] = 0.039;
	ps->Reals[I_B_TOT_SL_HI] = 1.6e-13;
	ps->Reals[I_K_D_SL_HI] = 13;
	ps->Reals[I_B_TOT_SL_LO] = 0;
	ps->Reals[I_K_D_SL_LO] = 1100;
	ps->Reals[I_B_TOT_DYE] = 50;
	ps->Reals[I_K_OFF_DYE] = 0.110;
	ps->Reals[I_K_ON_DYE] = 0.10;
	ps->Reals[I_B_T_JSR] = 30000;
	ps->Reals[I_K_M_JSR] = 638;
	ps->Reals[I_V_REFILL] = 0.095;
	ps->Reals[I_V_RYR] = 3e-8;
	ps->Reals[I_V_CELL] = 25.84e3;
	ps->Reals[I_K_ON_CSQN] = 0.10;
	ps->Reals[I_K_OFF_CSQN] = 63.8;
	ps->Reals[I_RYR_ETA] = 2.1;
	ps->Reals[I_RYR_A_STAR] = 0;
	ps->Reals[I_RYR_EPS_CC] = -0.92;
	ps->Reals[I_RYR_EPS_OO] = -0.85;
	ps->Reals[I_RYR_K_PLUS] = 0.1107e-3;
	ps->Reals[I_RYR_K_MINUS] = 0.500;
	ps->Reals[I_RYR_PHI_M] = 4.0;
	ps->Reals[I_RYR_PHI_B] = 1500;
	ps->Reals[I_K_D_i] = 910;
	ps->Reals[I_K_D_SR] = 2240;
	ps->Reals[I_A_P] = 150;
	ps->Reals[I_LCC_DURATION] = 0.5;
	ps->Reals[I_T_CLAMP] = 200;
	ps->Reals[I_V_CLAMP] = 0;
	ps->Reals[I_INIT_RYR_HOLD] = 0;

	sprintf(ps->Chars[I_FILE_BASE],"mesh");
	sprintf(ps->Chars[I_OUTPUT_STATE],"output");
	sprintf(ps->Chars[I_OUTPUT_GRID],"output");
	sprintf(ps->Chars[I_PARAM_TITLE],"Untitled");

	//Initialize XML parser
	try {
        XMLPlatformUtils::Initialize();
    }
    catch (const XMLException& toCatch) {
        char* message = XMLString::transcode(toCatch.getMessage());
        cout << "Error during initialization! :\n"
             << message << "\n";
        XMLString::release(&message);
        return;
    }

    XercesDOMParser* parser = new XercesDOMParser();
    parser->setValidationScheme(XercesDOMParser::Val_Always);
    parser->setDoNamespaces(true);    // optional

    ErrorHandler* errHandler = (ErrorHandler*) new HandlerBase();
    parser->setErrorHandler(errHandler);

	//Parse XML file
	fprintf(stdout,"Parsing parameter file %s...\n",param_file);
    char* xmlFile = (char*)param_file;

    try {
        parser->parse(xmlFile);
    }
    catch (const XMLException& toCatch) {
        char* message = XMLString::transcode(toCatch.getMessage());
        cout << "Exception message is: \n"
             << message << "\n";
        XMLString::release(&message);
        return;
    }
    catch (const DOMException& toCatch) {
        char* message = XMLString::transcode(toCatch.msg);
        cout << "Exception message is: \n"
             << message << "\n";
        XMLString::release(&message);
        return;
    }
    catch (...) {
        cout << "Unexpected Exception \n" ;
        return;
    }

    DOMDocument* xmlDoc = parser->getDocument();

    XMLCh* TAG_PARAMETER = XMLString::transcode("parameter");
    XMLCh* TAG_SYMBOL = XMLString::transcode("symbol");
    XMLCh* TAG_VALUE = XMLString::transcode("value");

    DOMNodeList* paramElements = xmlDoc->getElementsByTagName(TAG_PARAMETER);
    XMLSize_t nodeCount = paramElements->getLength();

    if (nodeCount == 0) {
    	throw(std::runtime_error("No parameters were found!"));
    }

    for (XMLSize_t xx = 0; xx < nodeCount; ++xx) {
    	DOMNode* currentNode = paramElements->item(xx);
    	if (currentNode->getNodeType() && currentNode->getNodeType() == DOMNode::ELEMENT_NODE) {
    		DOMElement* currentElement = dynamic_cast< xercesc::DOMElement* >(currentNode);
			DOMNodeList* children_param = currentElement->getChildNodes();
			XMLSize_t nodeCount_param = children_param->getLength();
			char* symbol;
			char* value;

			for (XMLSize_t i = 0; i < nodeCount_param; ++i) {
				DOMNode* currentNode_param = children_param->item(i);
				if (currentNode_param->getNodeType() && currentNode_param->getNodeType() == DOMNode::ELEMENT_NODE) {
					DOMElement* currentElement_param = dynamic_cast< xercesc::DOMElement* >(currentNode_param);
					if (XMLString::equals(currentElement_param->getTagName(), TAG_SYMBOL)) {
						symbol = XMLString::transcode(currentElement_param->getTextContent());
					} else if (XMLString::equals(currentElement_param->getTagName(), TAG_VALUE)) {
						value = XMLString::transcode(currentElement_param->getTextContent());
					}
				}
			}

			cout << "symbol = " << symbol << ", value = " << value << endl;

			if (!strcmp(symbol,"flag_noflux")) {
				if (!strcmp(value,"true")) {
					ps->Integers[I_FLAG_NO_FLUX_BOUNDARY] = 1;
				} else {
					ps->Integers[I_FLAG_NO_FLUX_BOUNDARY] = 0;
				}
			} else if (!strcmp(symbol,"flag_nogating")) {
				if (!strcmp(value,"true")) {
					ps->Integers[I_FLAG_NO_GATING] = 1;
				} else {
					ps->Integers[I_FLAG_NO_GATING] = 0;
				}
			} else if (!strcmp(symbol,"flag_ffwd")) {
				if (!strcmp(value,"true")) {
					ps->Integers[I_GATING_FFWD] = 1;
				} else {
					ps->Integers[I_GATING_FFWD] = 0;
				}
			} else if (!strcmp(symbol,"flag_use_clock")) {
				if (!strcmp(value,"true")) {
					ps->Integers[I_FLAG_CLOCK_SEED] = 1;
				} else {
					ps->Integers[I_FLAG_CLOCK_SEED] = 0;
				}
			} else if (!strcmp(symbol,"flag_output_grid")) {
				if (!strcmp(value,"true")) {
					ps->Integers[I_OUTPUT_GRID_FLAG] = 1;
				} else {
					ps->Integers[I_OUTPUT_GRID_FLAG] = 0;
				}
			} else if (!strcmp(symbol,"protocol")) {
				if (!strcmp(value,"gain")) {
					ps->Integers[I_FLAG_V_CLAMP] = 1;
					ps->Integers[I_FLAG_FIDELITY_SIM] = 0;
				} else if (!strcmp(value,"fidelity")) {
					ps->Integers[I_FLAG_V_CLAMP] = 0;
					ps->Integers[I_FLAG_FIDELITY_SIM] = 1;
				} else {
					ps->Integers[I_FLAG_V_CLAMP] = 0;
					ps->Integers[I_FLAG_FIDELITY_SIM] = 0;
				}
			} else if (!strcmp(symbol,"flag_output_ryr_states")) {
				if (!strcmp(value,"true")) {
					ps->Integers[I_FLAG_RYR_OUT] = 1;
				} else {
					ps->Integers[I_FLAG_RYR_OUT] = 0;
				}
			} else if (!strcmp(symbol,"sims_per_proc")) {
				ps->Integers[I_N_SIMS] = atoi(value);
			} else if (!strcmp(symbol,"ryr_open_init")) {
				ps->Integers[I_INIT_RYR] = atoi(value);
			} else if (!strcmp(symbol,"c_min")) {
				ps->Reals[I_C_MIN] = atof(value);
			} else if (!strcmp(symbol,"ryr_open_time")) {
				ps->Reals[I_INIT_RYR_HOLD] = atof(value);
			} else if (!strcmp(symbol,"ryr_open_max")) {
				ps->Integers[I_N_RYR_FIDELITY] = atoi(value);
			} else if (!strcmp(symbol,"v_clamp")) {
				ps->Reals[I_V_CLAMP] = atof(value);
			} else if (!strcmp(symbol,"t_clamp")) {
				ps->Reals[I_T_CLAMP] = atof(value);
			} else if (!strcmp(symbol,"t_final")) {
				ps->Reals[I_T_FINAL] = atof(value);
			} else if (!strcmp(symbol,"t_step")) {
				ps->Reals[I_T_STEP] = atof(value)*(1e-6); //convert nanosecond to milliseconds
			} else if (!strcmp(symbol,"states_interval")) {
				ps->Reals[I_OUTPUT_STATES_INTERVAL] = atof(value);
			} else if (!strcmp(symbol,"grid_interval")) {
				ps->Reals[I_OUTPUT_GRID_INTERVAL] = atof(value);
			} else if (!strcmp(symbol,"seed_base")) {
				ps->Integers[I_RAND_SEED_0] = atoi(value);
			} else if (!strcmp(symbol,"seed_A")) {
				ps->Integers[I_RAND_ADD_A] = atoi(value);
			} else if (!strcmp(symbol,"seed_B")) {
				ps->Integers[I_RAND_ADD_B] = atoi(value);
			} else if (!strcmp(symbol,"CA_I")) {
				ps->Reals[I_C_0_CA] = atof(value);
			} else if (!strcmp(symbol,"CA_SR")) {
				ps->Reals[I_C_0_CA_JSR] = atof(value);
				ps->Reals[I_DEFAULT_CA_NSR] = atof(value);
			} else if (!strcmp(symbol,"CA_O")) {
				ps->Reals[I_CA_0] = atof(value);
			} else if (!strcmp(symbol,"V_M")) {
				ps->Reals[I_DEFAULT_VM] = atof(value);
			} else if (!strcmp(symbol,"D_CA")) {
				ps->Reals[I_D_CA] = atof(value);
			} else if (!strcmp(symbol,"D_CAJSR")) {
				ps->Reals[I_D_CA_JSR] = atof(value);
			} else if (!strcmp(symbol,"D_ATP")) {
				ps->Reals[I_D_ATP] = atof(value);
			} else if (!strcmp(symbol,"D_CMDN")) {
				ps->Reals[I_D_CMDN] = atof(value);
			} else if (!strcmp(symbol,"D_DYE")) {
				ps->Reals[I_D_DYE] = atof(value);
			} else if (!strcmp(symbol,"B_TOT_ATP")) {
				ps->Reals[I_B_TOT_ATP] = atof(value);
			} else if (!strcmp(symbol,"K_OFF_ATP")) {
				ps->Reals[I_K_OFF_ATP] = atof(value);
			} else if (!strcmp(symbol,"K_ON_ATP")) {
				ps->Reals[I_K_ON_ATP] = atof(value);
			} else if (!strcmp(symbol,"B_TOT_CMDN")) {
				ps->Reals[I_B_TOT_CMDN] = atof(value);
			} else if (!strcmp(symbol,"K_OFF_CMDN")) {
				ps->Reals[I_K_OFF_CMDN] = atof(value);
			} else if (!strcmp(symbol,"K_ON_CMDN")) {
				ps->Reals[I_K_ON_CMDN] = atof(value);
			} else if (!strcmp(symbol,"B_TOT_TRPN")) {
				ps->Reals[I_B_TOT_TRPN] = atof(value);
			} else if (!strcmp(symbol,"K_OFF_TRPN")) {
				ps->Reals[I_K_OFF_TRPN] = atof(value);
			} else if (!strcmp(symbol,"K_ON_TRPN")) {
				ps->Reals[I_K_ON_TRPN] = atof(value);
			} else if (!strcmp(symbol,"B_TOT_SL")) {
				ps->Reals[I_B_TOT_SL_HI] = atof(value);
			} else if (!strcmp(symbol,"K_D_SL")) {
				ps->Reals[I_K_D_SL_HI] = atof(value);
			} else if (!strcmp(symbol,"B_TOT_DYE")) {
				ps->Reals[I_B_TOT_DYE] = atof(value);
			} else if (!strcmp(symbol,"K_OFF_DYE")) {
				ps->Reals[I_K_OFF_DYE] = atof(value);
			} else if (!strcmp(symbol,"K_ON_DYE")) {
				ps->Reals[I_K_ON_DYE] = atof(value);
			} else if (!strcmp(symbol,"B_TOT_CSQN")) {
				ps->Reals[I_B_T_JSR] = atof(value);
			} else if (!strcmp(symbol,"K_OFF_CSQN")) {
				ps->Reals[I_K_OFF_CSQN] = atof(value);
			} else if (!strcmp(symbol,"K_ON_CSQN")) {
				ps->Reals[I_K_ON_CSQN] = atof(value);
			} else if (!strcmp(symbol,"V_REFILL")) {
				ps->Reals[I_V_REFILL] = atof(value);
			} else if (!strcmp(symbol,"V_RYR")) {
				ps->Reals[I_V_RYR] = atof(value);
			} else if (!strcmp(symbol,"RYR_ETA")) {
				ps->Reals[I_RYR_ETA] = atof(value);
			} else if (!strcmp(symbol,"RYR_A_STAR")) {
				ps->Reals[I_RYR_A_STAR] = atof(value);
			} else if (!strcmp(symbol,"RYR_EPS_CC")) {
				ps->Reals[I_RYR_EPS_CC] = atof(value);
			} else if (!strcmp(symbol,"RYR_EPS_OO")) {
				ps->Reals[I_RYR_EPS_OO] = atof(value);
			} else if (!strcmp(symbol,"RYR_K_MINUS")) {
				ps->Reals[I_RYR_K_MINUS] = atof(value);
			} else if (!strcmp(symbol,"RYR_K_PLUS")) {
				ps->Reals[I_RYR_K_PLUS] = atof(value);
			} else if (!strcmp(symbol,"RYR_PHI_N")) {
				ps->Reals[I_RYR_PHI_M] = atof(value);
			} else if (!strcmp(symbol,"RYR_PHI_KD")) {
				ps->Reals[I_RYR_PHI_B] = atof(value);
			} else if (!strcmp(symbol,"K_D_I")) {
				ps->Reals[I_K_D_i] = atof(value);
			} else if (!strcmp(symbol,"K_D_SR")) {
				ps->Reals[I_K_D_SR] = atof(value);
			} else if (!strcmp(symbol,"A_P")) {
				ps->Reals[I_A_P] = atof(value);
			} else if (!strcmp(symbol,"paramset_title")) {
				sprintf(ps->Chars[I_PARAM_TITLE],"%s",value);
    		}
    	}
    }

    delete parser;
    delete errHandler;

	//Initialize parameter arrays

	//Diffusion coefficients
	ps->Diff[0] = ps->Reals[I_D_CA];
	ps->Diff[1] = ps->Reals[I_D_ATP];
	ps->Diff[2] = ps->Reals[I_D_CMDN];
	ps->Diff[3] = ps->Reals[I_D_TRPN];
	ps->Diff[4] = ps->Reals[I_D_CSQN];
	ps->Diff[5] = ps->Reals[I_D_DYE];

	//Default dynamic global variable values
	ps->Default_Globals[0] = ps->Reals[I_DEFAULT_VM];
	ps->Default_Globals[1] = ps->Reals[I_DEFAULT_CA_NSR];

	for (int i = 0; i < N_PARAMETERS_INT; i++) {
		fprintf(stdout,"Integer Parameter[%d] = %d\n",i,ps->Integers[i]);
	}
	for (int i = 0; i < N_PARAMETERS_REAL; i++) {
		fprintf(stdout,"Real Parameter[%d] = %g\n",i,ps->Reals[i]);
	}
	for (int i = 0; i < N_PARAMETERS_CHAR; i++) {
		fprintf(stdout,"Char Parameter[%d] = %s\n",i,ps->Chars[i]);
	}


}

//Reads grid structure from files
void ReadGrid(SimData* sd, GridStruct* gs, ParamStruct* ps) {

	char filename[255];
	sprintf(filename,"%s_tet_properties.txt",ps->Chars[I_FILE_BASE]);

	FILE* pFile;
	pFile = fopen(filename,"r");
	char line[255];

	fprintf(stdout,"Reading element property file %s...\n",filename);
	if (pFile != NULL) {

		//Read number of elements
		if (!fscanf(pFile,"%d\n",&(gs->N_Ele) )) {
			fprintf(stderr,"Error reading element property : number of elements\n");
		}

		//Allocate space for element property arrays
		try {
			gs->Domain = new int[gs->N_Ele];
			gs->V0 = new double[gs->N_Ele];
			gs->TTSurfaceArea = new double[gs->N_Ele];
			gs->SRSurfaceArea = new double[gs->N_Ele];
			gs->TropC = new int[gs->N_Ele];
		} catch (std::bad_alloc& ba)
		{
			fprintf(stderr,"bad_alloc caught: %s\n",ba.what());
		}

		//Skip column labels
		if (fgets(line,255,pFile) == NULL) {
			fprintf(stderr,"Error skipping first header line in element property file.\n");
		}

		//Read properties
		for (int i = 0; i < gs->N_Ele; i++) {
			double j1,j2,j3;
			if (!fscanf(pFile,"%lg %lg %lg %lg %lg %lg %d %d\n",&j1,&j2,&j3,
						&(gs->V0[i]),&(gs->TTSurfaceArea[i]),&(gs->SRSurfaceArea[i]),&(gs->TropC[i]),&(gs->Domain[i]) )) {
				fprintf(stderr,"Error reading element property data at line %d.\n",i+3);
			}
		}

		fclose(pFile);
	} else {
		fprintf(stderr,"Error opening element properties file: %s\n",filename);
	}

	//Read face file
	sprintf(filename,"%s_face_properties.txt",ps->Chars[I_FILE_BASE]);
	fprintf(stdout,"Reading face property file %s...\n",filename);
	pFile = fopen(filename,"r");

	if (pFile != NULL) {

		//Read number of faces
		if (!fscanf(pFile,"%d\n",&(gs->N_Face) )) {
			fprintf(stderr,"Error reading face property : number of faces.\n");
		}

		//Allocate space for triangle property arrays
		try {
			gs->Ele_1 = new int[gs->N_Face];
			gs->Ele_2 = new int[gs->N_Face];
			gs->FaceType = new int[gs->N_Face];
			gs->Face_C = new double[gs->N_Face];
			gs->dXi = new double[gs->N_Face];
			gs->Area_f = new double[gs->N_Face];
			gs->Afx = new double[gs->N_Face];
			gs->Afy = new double[gs->N_Face];
			gs->Afz = new double[gs->N_Face];
			gs->eXix = new double[gs->N_Face];
			gs->eXiy = new double[gs->N_Face];
			gs->eXiz = new double[gs->N_Face];
			gs->r1x = new double[gs->N_Face];
			gs->r1y = new double[gs->N_Face];
			gs->r1z = new double[gs->N_Face];
			gs->r2x = new double[gs->N_Face];
			gs->r2y = new double[gs->N_Face];
			gs->r2z = new double[gs->N_Face];
		} catch (std::bad_alloc& ba)
		{
			fprintf(stderr,"bad_alloc caught: %s\n",ba.what());
		}

		//Skip column labels
		if (fgets(line,255,pFile) == NULL) {
			fprintf(stderr,"Error skipping second header line in face property file.\n");
		}

		//Read properties
		for (int i = 0; i < gs->N_Face; i++) {

			double j1,j2,j3,j4,j5,j6;
			if ( !fscanf(pFile,"%d %d %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %d\n",
					&(gs->Ele_1[i]),&(gs->Ele_2[i]),
					&j1,&j2,&j3, //do not load cell normals
					&j4,&j5,&j6, //do not load cell centers
					&(gs->Area_f[i]),
					&(gs->Afx[i]),&(gs->Afy[i]),&(gs->Afz[i]),
					&(gs->r1x[i]),&(gs->r1y[i]),&(gs->r1z[i]),
					&(gs->r2x[i]),&(gs->r2y[i]),&(gs->r2z[i]),
					&(gs->eXix[i]),&(gs->eXiy[i]),&(gs->eXiz[i]),
					&(gs->dXi[i]),
					&(gs->Face_C[i]),
					&(gs->FaceType[i])) ) {
				fprintf(stderr,"Error reading face data at line %d.\n",i+3);
			}

			//Convert from geometry interface type to simulation interface type
			switch (gs->FaceType[i]) {
			case IN_FACE_DEFAULT:
				gs->FaceType[i] = FACE_DEFAULT;
				break;
			case IN_FACE_BOUNDARY:
				if (ps->Integers[I_FLAG_NO_FLUX_BOUNDARY]) {
					gs->FaceType[i] = FACE_NOFLUX;
				} else {
					gs->FaceType[i] = FACE_BOUND;
				}
				break;
			case IN_FACE_TT:
				gs->FaceType[i] = FACE_NOFLUX;
				break;
			case IN_FACE_JSR:
				gs->FaceType[i] = FACE_NOFLUX;
				break;
			default:
				fprintf(stderr,"Warning: undefined input face type %d.\n",gs->FaceType[i]);
				break;
			}
		}
		fclose(pFile);
	} else {
		fprintf(stderr,"Error opening face properties file: %s\n",filename);
	}

	sprintf(filename,"%s.node",ps->Chars[I_FILE_BASE]);
	fprintf(stdout,"Reading node file %s...\n",filename);
	pFile = fopen(filename,"r");
	if (pFile != NULL) {

		//Read number of nodes
		if (!fscanf(pFile,"%d",&(gs->N_Nodes) )) {
			fprintf(stderr,"Error node file - number of nodes\n");
		}

		//Allocate space for element property arrays
		try {
			gs->Node_X = new double[gs->N_Nodes];
			gs->Node_Y = new double[gs->N_Nodes];
			gs->Node_Z = new double[gs->N_Nodes];
		} catch (std::bad_alloc& ba)
		{
			fprintf(stderr,"bad_alloc caught: %s\n",ba.what());
		}

		//Skip column labels
		if (fgets(line,255,pFile) == NULL) {
			fprintf(stderr,"Error skipping first header line in element node file.\n");
		}

		//Read properties
		for (int i = 0; i < gs->N_Nodes; i++) {
			int j1;
			if (!fscanf(pFile,"%d %lg %lg %lg\n",&j1,&(gs->Node_X[i]),&(gs->Node_Y[i]),&(gs->Node_Z[i])) ) {
				fprintf(stderr,"Error reading element property data at line %d.\n",i+3);
			}
		}

		fclose(pFile);
	} else {
		fprintf(stderr,"Error opening element properties file: %s\n",filename);
	}

	sprintf(filename,"%s.ele",ps->Chars[I_FILE_BASE]);
	fprintf(stdout,"Reading cell-node file %s...\n",filename);
	pFile = fopen(filename,"r");
	if (pFile != NULL) {

		//Skip column labels
		if (fgets(line,255,pFile) == NULL) {
			fprintf(stderr,"Error skipping first header line in element node file.\n");
		}

		//Allocate space for element property arrays
		try {
			gs->Ele_Nodes = new int[gs->N_Ele*4];
		} catch (std::bad_alloc& ba)
		{
			fprintf(stderr,"bad_alloc caught: %s\n",ba.what());
		}

		//Read properties
		for (int i = 0; i < gs->N_Ele; i++) {
			int j1,j2;
			int idx = 4*i;
			if (!fscanf(pFile,"%d %d %d %d %d %d\n",&j1,&(gs->Ele_Nodes[idx]),&(gs->Ele_Nodes[idx+1]),&(gs->Ele_Nodes[idx+2]),&(gs->Ele_Nodes[idx+3]),&j2) ) {
				fprintf(stderr,"Error reading element property data at line %d.\n",i+3);
			}

			//Convert to 0-based indexes
			gs->Ele_Nodes[idx]--;
			gs->Ele_Nodes[idx+1]--;
			gs->Ele_Nodes[idx+2]--;
			gs->Ele_Nodes[idx+3]--;
		}

		fclose(pFile);
	} else {
		fprintf(stderr,"Error opening element properties file: %s\n",filename);
	}

	sprintf(filename,"%s_lcc.txt",ps->Chars[I_FILE_BASE]);
	fprintf(stdout,"Reading LCC file %s...\n",filename);
	pFile = fopen(filename,"r");
	if (pFile != NULL) {

		//Read number of channels
		if (!fscanf(pFile,"%d\n",&(sd->N_LCC) )) {
			fprintf(stderr,"Error node file - number of channels\n");
		}

		//Allocate space for element property arrays
		try {
			sd->LCC_Ele = new int[sd->N_LCC];
		} catch (std::bad_alloc& ba)
		{
			fprintf(stderr,"bad_alloc caught: %s\n",ba.what());
		}

		//Read properties
		for (int i = 0; i < sd->N_LCC; i++) {
			if (!fscanf(pFile,"%d\n",&(sd->LCC_Ele[i]) )) {
				fprintf(stderr,"Error reading LCC property data at line %d.\n",i+3);
			}
		}
		fclose(pFile);
	} else {
		fprintf(stderr,"Error opening file: %s\n",filename);
	}

	sprintf(filename,"%s_ryr.txt",ps->Chars[I_FILE_BASE]);
	fprintf(stdout,"Reading RyR file %s...\n",filename);
	pFile = fopen(filename,"r");
	if (pFile != NULL) {

		//Read number of channels
		if (!fscanf(pFile,"%d\n",&(sd->N_RyR) )) {
			fprintf(stderr,"Error node file - number of channels\n");
		}

		//Allocate space for element property arrays
		try {
			sd->RyR_Ele = new int[sd->N_RyR];
			sd->RyR_JSR_Ele = new int[sd->N_RyR];
			sd->RyR_Neighb = new int[sd->N_RyR*4];
		} catch (std::bad_alloc& ba)
		{
			fprintf(stderr,"bad_alloc caught: %s\n",ba.what());
		}

		//Read properties
		for (int i = 0; i < sd->N_RyR; i++) {
			if (!fscanf(pFile,"%d %d %d %d %d %d\n",&(sd->RyR_Ele[i]),&(sd->RyR_JSR_Ele[i]),&(sd->RyR_Neighb[i*4+0]),&(sd->RyR_Neighb[i*4+1]),&(sd->RyR_Neighb[i*4+2]),&(sd->RyR_Neighb[i*4+3]) )) {
				fprintf(stderr,"Error reading RyR property data at line %d.\n",i+3);
			}
		}
		fclose(pFile);
	} else {
		fprintf(stderr,"Error opening file: %s\n",filename);
	}

	for (int i = 0; i < sd->N_RyR; i++) {
		for (int j = i+1; j < sd->N_RyR; j++) {
			if (sd->RyR_Ele[i] == sd->RyR_Ele[j] || sd->RyR_JSR_Ele[i] == sd->RyR_JSR_Ele[j]) {
				fprintf(stderr,"**Warning: RyRs %d and %d located in same element. Parallel race conditions may exist.\n",i,j);
			}
		}
	}

	for (int i = 0; i < sd->N_LCC; i++) {
		for (int j = i+1; j < sd->N_LCC; j++) {
			if (sd->LCC_Ele[i] == sd->LCC_Ele[j]) {
				fprintf(stderr,"**Warning: LCCs %d and %d located in same element. Parallel race conditions may exist.\n",i,j);
			}
		}
	}

	//Add total number of channels
	sd->N_Channels = sd->N_RyR + sd->N_LCC;

	fprintf(stdout,"Done reading input mesh. # Faces = %d, # Elements = %d\n",gs->N_Face,gs->N_Ele);

	fprintf(stdout,"Optimizing grid numbering structure...\n");
	vector<int> ele_1; //lists of adjacent elements
	vector<int> ele_2;
	for (int i = 0; i < sd->Grid.N_Face; i++) {
		if (gs->Ele_2[i] >= 0) {
			ele_1.push_back(gs->Ele_1[i]);
			ele_2.push_back(gs->Ele_2[i]);
		}
	}

	//The next section re-orders elements to improve memory access times
	//Simulation gets up to a ~10-50% speed increase depending on mesh size

	//Modified version of Reverse Cuthill-McKee (RCM) algorithm
	vector<int> vertices;
	vector<int> visited;
	vertices.push_back(0); //Starting node
	while (!vertices.empty()) {
		int v = vertices[0];
		vertices.erase(vertices.begin());
		visited.push_back(v);

		//Find unvisited neighbors of v
		for (int i = 0; i < ele_1.size(); i++) {
			if (ele_1[i] == v || ele_2[i] == v) {
				int neighb;
				if (ele_1[i] == v) {
					neighb = ele_2[i];
				} else {
					neighb = ele_1[i];
				}
				int bSkip = 0;
				for (int j = 0; j < vertices.size(); j++) {
					if (vertices[j] == neighb) {
						bSkip = 1;
						break;
					}
				}
				if (!bSkip) {
					for (int j = 0; j < visited.size(); j++) {
						if (visited[j] == neighb) {
							bSkip = 1;
							break;
						}
					}
					if (!bSkip) {
						vertices.push_back(neighb);
						//ele_1.erase(ele_1.begin() + i);
						//ele_2.erase(ele_2.begin() + i);
					}
				}
			}
		}
	}

	//Remap elements
	fprintf(stdout,"Remapping elements...\n");
	int* Domain = new int[gs->N_Ele];
	double* V0 = new double[gs->N_Ele];
	double* TTSurfaceArea = new double[gs->N_Ele];
	double* SRSurfaceArea = new double[gs->N_Ele];
	int* TropC = new int[gs->N_Ele];
	int* Ele_Nodes = new int[gs->N_Ele*4];
	int* Ele_1 = new int[gs->N_Face];
	int* Ele_2 = new int[gs->N_Face];
	int* RyR_Ele = new int[sd->N_RyR];
	int* RyR_JSR_Ele = new int[sd->N_RyR];
	int* LCC_Ele = new int[sd->N_LCC];

	memcpy(Domain,gs->Domain,gs->N_Ele*sizeof(int));
	memcpy(V0,gs->V0,gs->N_Ele*sizeof(double));
	memcpy(TTSurfaceArea,gs->TTSurfaceArea,gs->N_Ele*sizeof(double));
	memcpy(SRSurfaceArea,gs->SRSurfaceArea,gs->N_Ele*sizeof(double));
	memcpy(TropC,gs->TropC,gs->N_Ele*sizeof(int));
	memcpy(Ele_Nodes,gs->Ele_Nodes,4*gs->N_Ele*sizeof(int));
	memcpy(Ele_1,gs->Ele_1,gs->N_Face*sizeof(int));
	memcpy(Ele_2,gs->Ele_2,gs->N_Face*sizeof(int));
	memcpy(RyR_Ele,sd->RyR_Ele,sd->N_RyR*sizeof(int));
	memcpy(RyR_JSR_Ele,sd->RyR_JSR_Ele,sd->N_RyR*sizeof(int));
	memcpy(LCC_Ele,sd->LCC_Ele,sd->N_LCC*sizeof(int));

	for (int i = 0; i < gs->N_Ele; i++) {
		gs->Domain[i] = Domain[visited[i]];
		gs->V0[i] = V0[visited[i]];
		gs->TTSurfaceArea[i] = TTSurfaceArea[visited[i]];
		gs->SRSurfaceArea[i] = SRSurfaceArea[visited[i]];
		gs->TropC[i] = TropC[visited[i]];
		gs->Ele_Nodes[4*i+0] = Ele_Nodes[4*visited[i]+0];
		gs->Ele_Nodes[4*i+1] = Ele_Nodes[4*visited[i]+1];
		gs->Ele_Nodes[4*i+2] = Ele_Nodes[4*visited[i]+2];
		gs->Ele_Nodes[4*i+3] = Ele_Nodes[4*visited[i]+3];

		for (int j = 0; j < gs->N_Face; j++) {
			if (Ele_1[j] == visited[i]) {
				gs->Ele_1[j] = i;
			}
			if (Ele_2[j] == visited[i]) {
				gs->Ele_2[j] = i;
			}
		}
		for (int j = 0; j < sd->N_RyR; j++) {
			if (RyR_Ele[j] == visited[i]) {
				sd->RyR_Ele[j] = i;
			}
		}
		for (int j = 0; j < sd->N_RyR; j++) {
			if (RyR_JSR_Ele[j] == visited[i]) {
				sd->RyR_JSR_Ele[j] = i;
			}
		}
		for (int j = 0; j < sd->N_LCC; j++) {
			if (LCC_Ele[j] == visited[i]) {
				sd->LCC_Ele[j] = i;
			}
		}
	}

	fprintf(stdout,"Computing matrix properties...\n");
	//Compute maximum and mean distance from a neighbor
	int max_dist = 0;
	double mean_dist = 0;
	for (int i = 0; i < gs->N_Face; i++) {
		if (gs->Ele_2[i] >= 0) {
			int dist = abs(gs->Ele_1[i]-gs->Ele_2[i]);
			mean_dist += (double)dist;
			if (dist > max_dist) {
				max_dist = dist;
			}
		}
	}
	mean_dist /= gs->N_Face;
	fprintf(stdout,"\tMax / Mean neighbor distance: %d / %g\n",max_dist,mean_dist);

	delete[] Domain, V0, TTSurfaceArea, SRSurfaceArea, TropC, Ele_Nodes, Ele_1, Ele_2, RyR_Ele, LCC_Ele;

}

void ReadGridStates(GridStruct* gs, ParamStruct* ps) {

	char filename[255];
	sprintf(filename,"%s_init_states.txt",ps->Chars[I_FILE_BASE]);

	//Read Grid File
	fprintf(stdout,"Reading grid initial conditions file %s...\n",filename);
	FILE* pFile;
	pFile = fopen(filename,"r");
	char line[255];

	if (pFile != NULL) {

		//Skip header
		fgets(line,255,pFile);
		fgets(line,255,pFile);

		//Get number of states
		int temp1, temp2;
		if (!fscanf(pFile,"%lf %d %d\n",&temp1, &temp2, &gs->N_States)) {
			fprintf(stderr,"Error reading grid header, could not get number of states.\n");
			exit(0);
		}

		fgets(line,255,pFile);
		fgets(line,255,pFile);

		//Allocate memory
		gs->States = new double*[gs->N_States];
		for (int i = 0; i < gs->N_States; i++) {
			gs->States[i] = new double[gs->N_Ele];
		}

		//Read grid data
		for (int k = 0; k < gs->N_States; k++) {
			for (int i = 0; i < gs->N_Ele; i++) {
				if (!fscanf(pFile,"%lf ",&(gs->States[k][i]))) {
					fprintf(stderr,"Error reading grid initial conditions data for state %d, cell %d.\n",k, i);
				}
			}
			fgets(line,255,pFile);
		}

		fclose(pFile);
	} else {
		fprintf(stderr,"Error: Could not open grid initial conditions file %s.\n",filename);
		exit(0);
	}


}

void InitializeDefaultGrid(SimData* sd, ParamStruct* ps) {

	////////////////////////
	//Cytosolic grid
	////////////////////////

	GridStruct* gs = &sd->Grid;

	//Use states specified in parameters.h
	gs->N_States = N_STATES;

	//Allocate grid memory
	gs->States = new double*[gs->N_States];
	for (int i = 0; i < gs->N_States; i++) {
		gs->States[i] = new double[gs->N_Ele];
	}
	gs->Boundaries = new double[gs->N_Ele];

	//Initialize default values
	ResetDefaultGridStates(sd,ps);

}

void ResetDefaultGridStates(SimData* sd, ParamStruct* ps) {

	GridStruct* gs = &sd->Grid;

	//Calculate equilibrium concentrations given fixed Ca2+
	double ca0 = ps->Reals[I_C_0_CA];
	double ca0_jsr = ps->Reals[I_C_0_CA_JSR];
	double Kd, btot;

	//TroponinC
	Kd = ps->Reals[I_K_OFF_TRPN]/ps->Reals[I_K_ON_TRPN];
	double c_trpn  = ca0 * ps->Reals[I_B_TOT_TRPN] / (Kd + ca0);

	//Calsequestrin
	Kd = ps->Reals[I_K_OFF_CSQN]/ps->Reals[I_K_ON_CSQN];
	double c_csqn  = ca0_jsr * ps->Reals[I_B_T_JSR] / (Kd + ca0_jsr);

	//ATP
	Kd = ps->Reals[I_K_OFF_ATP]/ps->Reals[I_K_ON_ATP];
	double c_atp  = ca0 * ps->Reals[I_B_TOT_ATP] / (Kd + ca0);

	//Calmodulin
	Kd = ps->Reals[I_K_OFF_CMDN]/ps->Reals[I_K_ON_CMDN];
	double c_cmdn  = ca0 * ps->Reals[I_B_TOT_CMDN] / (Kd + ca0);

	//Dye
	Kd = ps->Reals[I_K_OFF_DYE]/ps->Reals[I_K_ON_DYE];
	double c_dye  = ca0 * ps->Reals[I_B_TOT_DYE] / (Kd + ca0);

	//Set boundary conditions
	gs->Boundaries[INDEX_CA] = ca0;
	gs->Boundaries[INDEX_TRPN] = c_trpn;
	gs->Boundaries[INDEX_CSQN] = 0;
	gs->Boundaries[INDEX_ATP] = c_atp;
	gs->Boundaries[INDEX_CMDN] = c_cmdn;
	gs->Boundaries[INDEX_DYE] = c_dye;

	//Initialize grid values
	for (int j = 0; j < gs->N_Ele; j++) {

		if (gs->Domain[j] == DOMAIN_CYTO) {
			gs->States[INDEX_CA][j] = ca0;
			gs->States[INDEX_TRPN][j] = c_trpn*gs->TropC[j];
			gs->States[INDEX_ATP][j] = c_atp;
			gs->States[INDEX_CMDN][j] = c_cmdn;

			/*
			//Sarcolemmal binding sites must be calculated for each cell
			btot = gs->TTSurfaceArea[j]*B_TOT_SL / (gs->V0[j]*(1e-15));
			Kd = K_OFF_SL/K_ON_SL;
			double c_sl  = ca0 * btot / (Kd + ca0);
			gs->States[INDEX_SL][j] = c_sl;
			*/

			gs->States[INDEX_DYE][j] = c_dye;

		} else if (gs->Domain[j] == DOMAIN_JSR) {
			gs->States[INDEX_CA][j] = ca0_jsr;
			gs->States[INDEX_TRPN][j] = 0;
			gs->States[INDEX_CSQN][j] = c_csqn;
			gs->States[INDEX_ATP][j] = 0;
			gs->States[INDEX_CMDN][j] = 0;
			gs->States[INDEX_DYE][j] = 0;
		}
	}
}

void InitializeDefaultChannels(SimData* sd) {

	sd->LCC_States = new int[sd->N_LCC];
	sd->LCCV_States = new int[sd->N_LCC];

	//Initialize residuals
	sd->R = new double[sd->N_Channels];

	//Reset to default values
	ResetDefaultChannels(sd);

}

void ResetDefaultChannels(SimData* sd) {

	//Initialize channel states
	for (int i = 0; i < sd->N_LCC; i++) {
		sd->LCC_States[i] = DEFAULT_STATE_LCC;
		sd->LCCV_States[i] = DEFAULT_STATE_LCCV;
	}

	sd->RyR_States = new int[sd->N_RyR];
	for (int i = 0; i < sd->N_RyR; i++) {
		sd->RyR_States[i] = DEFAULT_STATE_RYR;
	}
}

void ResetGlobals(SimData* sd, ParamStruct* ps) {

	//Initialize default global variables
	for (int i = 0; i < N_GLOBAL_VARS; i++) {
		sd->Global_Vars[i] = ps->Default_Globals[i];
	}
}

__global__ void FindDiagonals(int N_Ele, int* diag, int n_cols, int* col_indices) {

	//Note the columns indices are stored in column-major order
	for (int j = 0; j < n_cols; j++) {
		for (int k = 0; k < N_Ele; k++) {
			if (col_indices[j*N_Ele + k] == k) {
				diag[k] = j*N_Ele + k;
			}
		}
	}
}

void InitializeDomainStruct(DeviceDomainStruct* dds, GridStruct* gs, int bSecondaryGradient, ParamStruct* ps) {

	//Create sparse matrices
	//fprintf(stdout,"Initializing domain structure...\n");
	//dds->A = new sp_type[gs->N_States];
	InitializeLeastSquares(dds,gs,bSecondaryGradient,ps);

	if ( cudaSuccess != cudaGetLastError() )
		fprintf(stderr,"*****************Cuda Error during InitializeLeastSquares()!\n");

	//fprintf(stdout,"Allocating memory for dense vectors...\n");
	//Allocate dense vector arrays
	//dds->Array_2 = cusp::array1d<double,cusp::device_memory>(gs->N_Ele);


	if ( cudaSuccess != cudaGetLastError() )
		fprintf(stderr,"*****************Cuda Error!\n");
	//Grid buffer
	//fprintf(stdout,"Allocating memory for buffer vectors...\n");
	/*dds->Array_1 = new cusp::array1d<double,cusp::device_memory>[gs->N_States];
	for (int i = 0; i < gs->N_States; i++) {
		dds->Array_1[i] = cusp::array1d<double,cusp::device_memory>(gs->N_Ele);
	}*/

	if ( cudaSuccess != cudaGetLastError() )
		fprintf(stderr,"*****************Cuda Error!\n");

	//fprintf(stdout,"Computing boundary rates...\n");
	//Boundary fluxes for diffusible states
	double** boundary_rates = new double*[gs->N_States];
	for (int k = 0; k < gs->N_States; k++) {
		boundary_rates[k] = new double[gs->N_Ele];
		for (int i = 0; i < gs->N_Ele; i++) {
			boundary_rates[k][i] = 0;
		}
		for (int i = 0; i < gs->N_Face; i++) {
			int e = i;
			int tri1 = gs->Ele_1[e];
			if (gs->FaceType[e] == FACE_BOUND) {
				boundary_rates[k][tri1] += ps->Diff[k] * (gs->Face_C[e]) * gs->Boundaries[k] / gs->V0[tri1];
			}
		}
	}

	//JSR labels and volume
	int* bJSR = new int[gs->N_Ele];
	dds->Vol_JSR = 0;
	for (int i = 0; i < gs->N_Ele; i++) {
		bJSR[i] = gs->Domain[i] == DOMAIN_JSR ? 1 : 0;
		if (gs->Domain[i] == DOMAIN_JSR) {
			dds->Vol_JSR += gs->V0[i];
		}
	}

	//Number of grid states
	dds->N_States = gs->N_States;
	//Number of grid elements
	dds->N_Ele = gs->N_Ele;


	if ( cudaSuccess != cudaGetLastError() )
		fprintf(stderr,"*****************Cuda Error initializing boundaries and JSR!\n");
	//fprintf(stdout,"Allocating device memory...\n");
	//Allocate arrays on device
	cudaMalloc((void**) &dds->Array_1, sizeof(double)*(dds->N_States)*dds->N_Ele);
	cudaMalloc((void**) &dds->States, sizeof(double)*(dds->N_States)*dds->N_Ele);
	cudaMalloc((void**) &dds->Boundary_Rates, sizeof(double)*(dds->N_States)*dds->N_Ele);
	cudaMalloc((void**) &dds->TropC, sizeof(int)*dds->N_Ele);
	cudaMalloc((void**) &dds->bJSR, sizeof(int)*dds->N_Ele);
	cudaMalloc((void**) &dds->V0, sizeof(double)*dds->N_Ele);
	cudaMalloc((void**) &dds->TTSurfaceArea, sizeof(double)*dds->N_Ele);
	cudaMalloc((void**) &dds->SRSurfaceArea, sizeof(double)*dds->N_Ele);
	cudaMalloc((void**) &dds->aij, sizeof(double)*5*dds->N_Ele);
	cudaMalloc((void**) &dds->nij, sizeof(int)*4*dds->N_Ele);


	if ( cudaSuccess != cudaGetLastError() )
		fprintf(stderr,"*****************Cuda Error cudaMalloc'ing space!\n");

	//Copy data to device
	//fprintf(stdout,"Copying to device memory...\n");
	for (int i = 0; i < dds->N_States; i++) {
		cudaMemcpy(dds->States + (i*gs->N_Ele), gs->States[i], sizeof(double)*dds->N_Ele, cudaMemcpyHostToDevice);
		cudaMemcpy(dds->Boundary_Rates + (i*gs->N_Ele), boundary_rates[i], sizeof(double)*dds->N_Ele, cudaMemcpyHostToDevice);
	}
	cudaMemcpy(dds->TropC, gs->TropC, sizeof(int)*dds->N_Ele, cudaMemcpyHostToDevice);
	cudaMemcpy(dds->bJSR, bJSR, sizeof(int)*dds->N_Ele, cudaMemcpyHostToDevice);
	cudaMemcpy(dds->V0, gs->V0, sizeof(double)*dds->N_Ele, cudaMemcpyHostToDevice);
	cudaMemcpy(dds->TTSurfaceArea, gs->TTSurfaceArea, sizeof(double)*dds->N_Ele, cudaMemcpyHostToDevice);
	cudaMemcpy(dds->SRSurfaceArea, gs->SRSurfaceArea, sizeof(double)*dds->N_Ele, cudaMemcpyHostToDevice);
	cudaMemcpy(dds->aij, gs->aij, sizeof(double)*5*dds->N_Ele, cudaMemcpyHostToDevice);
	cudaMemcpy(dds->nij, gs->nij, sizeof(int)*4*dds->N_Ele, cudaMemcpyHostToDevice);

	if ( cudaSuccess != cudaGetLastError() )
		fprintf(stderr,"*****************Cuda Error copying data to device!\n");

	//Free memory
	//fprintf(stdout,"Freeing memory...\n");
	for (int k = 0; k < gs->N_States; k++) {
		delete[] boundary_rates[k];
	}
	delete[] boundary_rates;
	delete[] bJSR;

}

void InitializeLeastSquares(DeviceDomainStruct* dds, GridStruct* gs, int bSecondaryGradient, ParamStruct* ps) {

	//Sparse matrix for gradient calculation

	if (!bSecondaryGradient) {

		//Create sparse matrix
		fprintf(stdout,"Computing sparse matrix elements (primary gradient only)...\n");
		/*for (int i = 0; i < gs->N_States; i++) {
			dds->A[i] = InitializeExplicitSparseMatrixFast(gs,i,0,0,bSecondaryGradient,ps);
		}*/
		InitializeExplicitSparseMatrixFast(gs,0,0,0,bSecondaryGradient,ps);

	} else {
		fprintf(stdout,"Computing sparse matrix elements (secondary gradient)...\n");
		//Least-squares matrices (estimate = MLS*DLS)
		double** MLS = new double*[gs->N_Ele];
		for (int i = 0; i < gs->N_Ele; i++) {
			MLS[i] = new double[12];
		}

		//2D array of faces of each cell
		int** NLS = new int*[gs->N_Ele];
		for (int i = 0; i < gs->N_Ele; i++) {
			NLS[i] = new int[4];
			for (int j = 0; j < 4; j++) {
				NLS[i][j] = -1;
			}
		}

		//Calculate least-squares matrices
		int* numSet = new int[gs->N_Ele];
		for (int i = 0; i < gs->N_Ele; i++) {
			numSet[i] = 0;
		}

		double** M0 = new double*[gs->N_Ele];
		for (int i = 0; i < gs->N_Ele; i++) {
			M0[i] = new double[12];
		}

		double** M1 = new double*[gs->N_Ele];
		for (int i = 0; i < gs->N_Ele; i++) {
			M1[i] = new double[9];
		}

		double** W = new double*[gs->N_Ele];
		for (int i = 0; i < gs->N_Ele; i++) {
			W[i] = new double[4];
		}

		//Calculate M0=M matrix
		int tri1,tri2,idx;
		double x1,x2,x3,w;
		for (int i = 0; i < gs->N_Face; i++) {
			tri1 = gs->Ele_1[i];
			tri2 = gs->Ele_2[i];

			if (gs->FaceType[i] == FACE_DEFAULT) {
				x1 = gs->eXix[i]*gs->dXi[i];
				x2 = gs->eXiy[i]*gs->dXi[i];
				x3 = gs->eXiz[i]*gs->dXi[i];
			} else {
				x1 = gs->r1x[i];
				x2 = gs->r1y[i];
				x3 = gs->r1z[i];
			}

			idx = 3*numSet[tri1];

			//Weighted regression
			w = (x1*x1 + x2*x2 + x3*x3);
			//w = 1;
			if (w == 0) {
				fprintf(stderr,"Warning: infinite weight, assigning to 1e6.\n");
				w = 1e6;
			}
			W[tri1][numSet[tri1]] = 1/w;

			M0[tri1][idx] = x1/w;
			M0[tri1][idx+1] = x2/w;
			M0[tri1][idx+2] = x3/w;
			NLS[tri1][numSet[tri1]] = i;
			numSet[tri1]++;


			if (gs->FaceType[i] != FACE_BOUND) {

				if (tri2 >= 0) {
					if (gs->FaceType[i] == FACE_DEFAULT) {
						x1 = -gs->eXix[i]*gs->dXi[i];
						x2 = -gs->eXiy[i]*gs->dXi[i];
						x3 = -gs->eXiz[i]*gs->dXi[i];
					} else { //else interface
						x1 = gs->r2x[i];
						x2 = gs->r2y[i];
						x3 = gs->r2z[i];
					}

					idx = 3*numSet[tri2];

					//Weighted regression
					w = (x1*x1 + x2*x2 + x3*x3);
					//w = 1;
					if (w == 0) {
						fprintf(stderr,"Warning: infinite weight, assigning to 1e6.\n");
						w = 1e6;
					}
					W[tri2][numSet[tri2]] = 1/w;

					M0[tri2][idx] = x1/w;
					M0[tri2][idx+1] = x2/w;
					M0[tri2][idx+2] = x3/w;
					NLS[tri2][numSet[tri2]] = i;
					numSet[tri2]++;
				}
			}
		}

		//Calculate M1 = M0^T * M0
		for (int i = 0; i < gs->N_Ele; i++) {
			M1[i][0] = M0[i][0]*M0[i][0] + M0[i][3]*M0[i][3] + M0[i][6]*M0[i][6] + M0[i][9]*M0[i][9];
			M1[i][1] = M0[i][0]*M0[i][1] + M0[i][3]*M0[i][4] + M0[i][6]*M0[i][7] + M0[i][9]*M0[i][10];
			M1[i][2] = M0[i][0]*M0[i][2] + M0[i][3]*M0[i][5] + M0[i][6]*M0[i][8] + M0[i][9]*M0[i][11];
			M1[i][3] = M0[i][1]*M0[i][0] + M0[i][4]*M0[i][3] + M0[i][7]*M0[i][6] + M0[i][10]*M0[i][9];
			M1[i][4] = M0[i][1]*M0[i][1] + M0[i][4]*M0[i][4] + M0[i][7]*M0[i][7] + M0[i][10]*M0[i][10];
			M1[i][5] = M0[i][1]*M0[i][2] + M0[i][4]*M0[i][5] + M0[i][7]*M0[i][8] + M0[i][10]*M0[i][11];
			M1[i][6] = M0[i][2]*M0[i][0] + M0[i][5]*M0[i][3] + M0[i][8]*M0[i][6] + M0[i][11]*M0[i][9];
			M1[i][7] = M0[i][2]*M0[i][1] + M0[i][5]*M0[i][4] + M0[i][8]*M0[i][7] + M0[i][11]*M0[i][10];
			M1[i][8] = M0[i][2]*M0[i][2] + M0[i][5]*M0[i][5] + M0[i][8]*M0[i][8] + M0[i][11]*M0[i][11];
		}

		//Calculate inverse of M1
		double a,b,c,d,e,f,g,h,k, det;
		double A, B, C, D, E, F, G, H, K;
		for (int i = 0; i < gs->N_Ele; i++) {
			a = M1[i][0];
			b = M1[i][1];
			c = M1[i][2];
			d = M1[i][3];
			e = M1[i][4];
			f = M1[i][5];
			g = M1[i][6];
			h = M1[i][7];
			k = M1[i][8];
			det = a*e*k + b*f*g + c*d*h - c*e*g - b*d*k - a*f*h;

			if (det==0) {
				fprintf(stderr,"Warning: singular least-square matrix (%d), setting det=1e-6.\n",i);
				det = 1e-6;
			}

			A = e*k - f*h;
			B = f*g - d*k;
			C = d*h - e*g;
			D = c*h - b*k;
			E = a*k - c*g;
			F = g*b - a*h;
			G = b*f - c*e;
			H = c*d - a*f;
			K = a*e - b*d;

			M1[i][0] = A/det;
			M1[i][1] = D/det;
			M1[i][2] = G/det;
			M1[i][3] = B/det;
			M1[i][4] = E/det;
			M1[i][5] = H/det;
			M1[i][6] = C/det;
			M1[i][7] = F/det;
			M1[i][8] = K/det;
		}

		//Calculate MLS = M1^-1 * M0^T
		for (int i = 0; i < gs->N_Ele; i++) {
			MLS[i][0] = M1[i][0]*M0[i][0] + M1[i][1]*M0[i][1] + M1[i][2]*M0[i][2];
			MLS[i][1] = M1[i][0]*M0[i][3] + M1[i][1]*M0[i][4] + M1[i][2]*M0[i][5];
			MLS[i][2] = M1[i][0]*M0[i][6] + M1[i][1]*M0[i][7] + M1[i][2]*M0[i][8];
			MLS[i][3] = M1[i][0]*M0[i][9] + M1[i][1]*M0[i][10] + M1[i][2]*M0[i][11];
			MLS[i][4] = M1[i][3]*M0[i][0] + M1[i][4]*M0[i][1] + M1[i][5]*M0[i][2];
			MLS[i][5] = M1[i][3]*M0[i][3] + M1[i][4]*M0[i][4] + M1[i][5]*M0[i][5];
			MLS[i][6] = M1[i][3]*M0[i][6] + M1[i][4]*M0[i][7] + M1[i][5]*M0[i][8];
			MLS[i][7] = M1[i][3]*M0[i][9] + M1[i][4]*M0[i][10] + M1[i][5]*M0[i][11];
			MLS[i][8] = M1[i][6]*M0[i][0] + M1[i][7]*M0[i][1] + M1[i][8]*M0[i][2];
			MLS[i][9] = M1[i][6]*M0[i][3] + M1[i][7]*M0[i][4] + M1[i][8]*M0[i][5];
			MLS[i][10] = M1[i][6]*M0[i][6] + M1[i][7]*M0[i][7] + M1[i][8]*M0[i][8];
			MLS[i][11] = M1[i][6]*M0[i][9] + M1[i][7]*M0[i][10] + M1[i][8]*M0[i][11];
			//printf("%d %f %f %f %f %f %f\n",i,MLS[i][0],MLS[i][1],MLS[i][2],MLS[i][3],MLS[i][4],MLS[i][5]);
		}

		//Factor in weights
		for (int i = 0; i < gs->N_Ele; i++) {
			MLS[i][0] *= W[i][0];
			MLS[i][1] *= W[i][1];
			MLS[i][2] *= W[i][2];
			MLS[i][3] *= W[i][3];
			MLS[i][4] *= W[i][0];
			MLS[i][5] *= W[i][1];
			MLS[i][6] *= W[i][2];
			MLS[i][7] *= W[i][3];
			MLS[i][8] *= W[i][0];
			MLS[i][9] *= W[i][1];
			MLS[i][10] *= W[i][2];
			MLS[i][11] *= W[i][3];
		}

		//Create sparse matrix
		/*for (int i = 0; i < gs->N_States; i++) {
			dds->A[i] = InitializeExplicitSparseMatrixFast(gs,i,MLS,NLS,bSecondaryGradient,ps);
		}*/
		InitializeExplicitSparseMatrixFast(gs,0,MLS,NLS,bSecondaryGradient,ps);

		//Free memory
		for (int i = 0; i < gs->N_Ele; i++) {
			delete[] MLS[i];
		}
		delete[] MLS;
		for (int i = 0; i < gs->N_Ele; i++) {
			delete[] NLS[i];
		}
		delete[] NLS;
		delete[] numSet;
		for (int i = 0; i < gs->N_Ele; i++) {
			delete[] M0[i];
		}
		delete[] M0;
		for (int i = 0; i < gs->N_Ele; i++) {
			delete[] M1[i];
		}
		delete[] M1;
		for (int i = 0; i < gs->N_Ele; i++) {
			delete[] W[i];
		}
		delete[] W;
	}
}

//Class used for constructing sparse matrices with hash maps
class Entry2D {
public:
    int i;
    int j;
    int N_Cols;
    bool operator<(const Entry2D &other) const {
        return (((long int)this->i)*(this->N_Cols) + (this->j) < ((long int)other.i)*(this->N_Cols) + (other.j));
    }
    bool operator==(const Entry2D &other) const {
        return (this->i == other.i && this->j == other.j);
    }
};

//int Entry2D::N_Rows = 0; //Should be set to number of cells in code


void InitializeExplicitSparseMatrixFast(GridStruct* gs,int state,double** MLS,int** NLS, int bSecondaryGradient, ParamStruct* ps) {

	//Creates sparse matrix from using a map object, should work better for larger meshes
	//This matrix should be used for explicit time stepping

	int N_Face = gs->N_Face; //Number of edges
	int N_Cells = gs->N_Ele; //Number of control volume elements

	map< Entry2D ,double> valuemap;

	//Primary gradient
	double k1,k2;
	int tri1,tri2;
	Entry2D key;
	key.N_Cols = N_Cells;
	for (int i = 0; i < N_Face; i++) {
		tri1 = gs->Ele_1[i];
		tri2 = gs->Ele_2[i];

		/*double diff = ps->Diff[state];
		if (gs->Domain[tri1] == DOMAIN_JSR) {
			diff = ps->Reals[I_D_CA_JSR];
		}*/
		double diff = 1; //TODO

		k1 = diff * (gs->Face_C[i]) / gs->V0[tri1];

		if (gs->FaceType[i] == FACE_DEFAULT) {

			key.i = tri1;
			key.j = tri2;
			valuemap[key] += k1;

			key.i = tri1;
			key.j = tri1;
			valuemap[key] += -k1;

			k2 = diff * (gs->Face_C[i]) / gs->V0[tri2];

			key.i = tri2;
			key.j = tri1;
			valuemap[key] += k2;

			key.i = tri2;
			key.j = tri2;
			valuemap[key] += -k2;

		} else if (gs->FaceType[i] == FACE_BOUND) {

			key.i = tri1;
			key.j = tri1;
			//valuemap[key] += -k1; //Uncomment to enable diffusion out of domain boundaries

		}
	}

	//Secondary gradient
	if (bSecondaryGradient) {
		map< Entry2D ,double> secondarymap;
		int c0, c1, f, n;
		double g_c1x[5],g_c0x[5];
		double g_c1y[5],g_c0y[5];
		double g_c1z[5],g_c0z[5];
		for (int i = 0; i < N_Face; i++) {
			c0 = gs->Ele_1[i];
			c1 = gs->Ele_2[i];

			double diff = ps->Diff[state];
			if (gs->Domain[c0] == DOMAIN_JSR) {
				diff = ps->Reals[I_D_CA_JSR];
			}

			if (gs->FaceType[i] == FACE_DEFAULT) {

				//X direction
				g_c1x[1] = 0.5 * MLS[c1][0] * diff * (-gs->Afx[i] + (gs->Face_C[i]*gs->dXi[i]*gs->eXix[i]));
				g_c1x[2] = 0.5 * MLS[c1][1] * diff * (-gs->Afx[i] + (gs->Face_C[i]*gs->dXi[i]*gs->eXix[i]));
				g_c1x[3] = 0.5 * MLS[c1][2] * diff * (-gs->Afx[i] + (gs->Face_C[i]*gs->dXi[i]*gs->eXix[i]));
				g_c1x[4] = 0.5 * MLS[c1][3] * diff * (-gs->Afx[i] + (gs->Face_C[i]*gs->dXi[i]*gs->eXix[i]));
				g_c1x[0] = -(g_c1x[1] + g_c1x[2] + g_c1x[3] + g_c1x[4]);

				g_c0x[1] = 0.5 * MLS[c0][0] * diff * (-gs->Afx[i] + (gs->Face_C[i]*gs->dXi[i]*gs->eXix[i]));
				g_c0x[2] = 0.5 * MLS[c0][1] * diff * (-gs->Afx[i] + (gs->Face_C[i]*gs->dXi[i]*gs->eXix[i]));
				g_c0x[3] = 0.5 * MLS[c0][2] * diff * (-gs->Afx[i] + (gs->Face_C[i]*gs->dXi[i]*gs->eXix[i]));
				g_c0x[4] = 0.5 * MLS[c0][3] * diff * (-gs->Afx[i] + (gs->Face_C[i]*gs->dXi[i]*gs->eXix[i]));
				g_c0x[0] = -(g_c0x[1] + g_c0x[2] + g_c0x[3] + g_c0x[4]);

				//Y direction
				g_c1y[1] = 0.5 * MLS[c1][4] * diff * (-gs->Afy[i] + (gs->Face_C[i]*gs->dXi[i]*gs->eXiy[i]));
				g_c1y[2] = 0.5 * MLS[c1][5] * diff * (-gs->Afy[i] + (gs->Face_C[i]*gs->dXi[i]*gs->eXiy[i]));
				g_c1y[3] = 0.5 * MLS[c1][6] * diff * (-gs->Afy[i] + (gs->Face_C[i]*gs->dXi[i]*gs->eXiy[i]));
				g_c1y[4] = 0.5 * MLS[c1][7] * diff * (-gs->Afy[i] + (gs->Face_C[i]*gs->dXi[i]*gs->eXiy[i]));
				g_c1y[0] = -(g_c1y[1] + g_c1y[2] + g_c1y[3] + g_c1y[4]);

				g_c0y[1] = 0.5 * MLS[c0][4] * diff * (-gs->Afy[i] + (gs->Face_C[i]*gs->dXi[i]*gs->eXiy[i]));
				g_c0y[2] = 0.5 * MLS[c0][5] * diff * (-gs->Afy[i] + (gs->Face_C[i]*gs->dXi[i]*gs->eXiy[i]));
				g_c0y[3] = 0.5 * MLS[c0][6] * diff * (-gs->Afy[i] + (gs->Face_C[i]*gs->dXi[i]*gs->eXiy[i]));
				g_c0y[4] = 0.5 * MLS[c0][7] * diff * (-gs->Afy[i] + (gs->Face_C[i]*gs->dXi[i]*gs->eXiy[i]));
				g_c0y[0] = -(g_c0y[1] + g_c0y[2] + g_c0y[3] + g_c0y[4]);

				//Z direction
				g_c1z[1] = 0.5 * MLS[c1][8] * diff * (-gs->Afz[i] + (gs->Face_C[i]*gs->dXi[i]*gs->eXiz[i]));
				g_c1z[2] = 0.5 * MLS[c1][9] * diff * (-gs->Afz[i] + (gs->Face_C[i]*gs->dXi[i]*gs->eXiz[i]));
				g_c1z[3] = 0.5 * MLS[c1][10] * diff * (-gs->Afz[i] + (gs->Face_C[i]*gs->dXi[i]*gs->eXiz[i]));
				g_c1z[4] = 0.5 * MLS[c1][11] * diff * (-gs->Afz[i] + (gs->Face_C[i]*gs->dXi[i]*gs->eXiz[i]));
				g_c1z[0] = -(g_c1z[1] + g_c1z[2] + g_c1z[3] + g_c1z[4]);

				g_c0z[1] = 0.5 * MLS[c0][8] * diff * (-gs->Afz[i] + (gs->Face_C[i]*gs->dXi[i]*gs->eXiz[i]));
				g_c0z[2] = 0.5 * MLS[c0][9] * diff * (-gs->Afz[i] + (gs->Face_C[i]*gs->dXi[i]*gs->eXiz[i]));
				g_c0z[3] = 0.5 * MLS[c0][10] * diff * (-gs->Afz[i] + (gs->Face_C[i]*gs->dXi[i]*gs->eXiz[i]));
				g_c0z[4] = 0.5 * MLS[c0][11] * diff * (-gs->Afz[i] + (gs->Face_C[i]*gs->dXi[i]*gs->eXiz[i]));
				g_c0z[0] = -(g_c0z[1] + g_c0z[2] + g_c0z[3] + g_c0z[4]);


				//Visit C0's neighbors
				//fprintf(stdout,"%g\n",-g_c0z[0] / gs->V0[c1]);
				for (int j = 0; j < 4; j++) {
					f = NLS[c0][j];
					if (gs->FaceType[f] == FACE_DEFAULT) {
						n = gs->Ele_1[f];
						if (n == c0) {
							n = gs->Ele_2[f];
						}

						key.i = c0;
						key.j = c0;
						valuemap[key] += -(g_c0x[j+1] + g_c0y[j+1] + g_c0z[j+1]) / gs->V0[c0];
						secondarymap[key] += -(g_c0x[j+1] + g_c0y[j+1] + g_c0z[j+1]) / gs->V0[c0];

						key.i = c0;
						key.j = n;
						valuemap[key] += g_c0x[j+1] / gs->V0[c0];
						valuemap[key] += g_c0y[j+1] / gs->V0[c0];
						valuemap[key] += g_c0z[j+1] / gs->V0[c0];
						secondarymap[key] += g_c0x[j+1] / gs->V0[c0];
						secondarymap[key] += g_c0y[j+1] / gs->V0[c0];
						secondarymap[key] += g_c0z[j+1] / gs->V0[c0];

						key.i = c1;
						key.j = c1;
						valuemap[key] -= -(g_c0x[j+1] + g_c0y[j+1] + g_c0z[j+1]) / gs->V0[c1];
						secondarymap[key] -= -(g_c0x[j+1] + g_c0y[j+1] + g_c0z[j+1]) / gs->V0[c1];

						key.i = c1;
						key.j = n;
						valuemap[key] -= g_c0x[j+1] / gs->V0[c1];
						valuemap[key] -= g_c0y[j+1] / gs->V0[c1];
						valuemap[key] -= g_c0z[j+1] / gs->V0[c1];
						secondarymap[key] -= g_c0x[j+1] / gs->V0[c1];
						secondarymap[key] -= g_c0y[j+1] / gs->V0[c1];
						secondarymap[key] -= g_c0z[j+1] / gs->V0[c1];
					}
				}

				//Visit C1's neighbors
				for (int j = 0; j < 4; j++) {
					f = NLS[c1][j];
					if (gs->FaceType[f] == FACE_DEFAULT) {
						n = gs->Ele_1[f];
						if (n == c1) {
							n = gs->Ele_2[f];
						}

						key.i = c0;
						key.j = c0;
						valuemap[key] += -(g_c1x[j+1] + g_c1y[j+1] + g_c1z[j+1]) / gs->V0[c0];
						secondarymap[key] += -(g_c1x[j+1] + g_c1y[j+1] + g_c1z[j+1]) / gs->V0[c0];

						key.i = c0;
						key.j = n;
						valuemap[key] += g_c1x[j+1] / gs->V0[c0];
						valuemap[key] += g_c1y[j+1] / gs->V0[c0];
						valuemap[key] += g_c1z[j+1] / gs->V0[c0];
						secondarymap[key] += g_c1x[j+1] / gs->V0[c0];
						secondarymap[key] += g_c1y[j+1] / gs->V0[c0];
						secondarymap[key] += g_c1z[j+1] / gs->V0[c0];

						key.i = c1;
						key.j = c1;
						valuemap[key] -= -(g_c1x[j+1] + g_c1y[j+1] + g_c1z[j+1]) / gs->V0[c1];
						secondarymap[key] -= -(g_c1x[j+1] + g_c1y[j+1] + g_c1z[j+1]) / gs->V0[c1];

						key.i = c1;
						key.j = n;
						valuemap[key] -= g_c1x[j+1] / gs->V0[c1];
						valuemap[key] -= g_c1y[j+1] / gs->V0[c1];
						valuemap[key] -= g_c1z[j+1] / gs->V0[c1];
						secondarymap[key] -= g_c1x[j+1] / gs->V0[c1];
						secondarymap[key] -= g_c1y[j+1] / gs->V0[c1];
						secondarymap[key] -= g_c1z[j+1] / gs->V0[c1];
					}
				}

			}
		}
/*
		FILE * pFile;
		map<Entry2D,double>::iterator it;
		char file[255];
		sprintf(file,"secondary_matrix_%d.txt",state);
		pFile = fopen (file,"w");
		if (pFile != NULL) {
			for (it = secondarymap.begin(); it != secondarymap.end(); it++) {
				if ((*it).second != 0) {
					fprintf(pFile,"%d %d %g\n",(*it).first.i,(*it).first.j,(*it).second);
				}
			}
			fclose(pFile);
		} else {
			fprintf(stderr,"Warning: could not open matrix file %s.\n", file);
		}
*/
	}


	//Compute coordinate-based data
	int N_Nonzero = 0;
	map<Entry2D,double>::iterator it;
	for (it = valuemap.begin(); it != valuemap.end(); it++) {
		if ((*it).second != 0) {
			N_Nonzero++;
		}
	}
	double* vals = new double[N_Nonzero];
	int* r = new int[N_Nonzero];
	int* c = new int[N_Nonzero];
	int idx = 0;
	for (it = valuemap.begin(); it != valuemap.end(); it++) {
		if ((*it).second != 0) {
			vals[idx] = (*it).second;
			r[idx] = (*it).first.i;
			c[idx] = (*it).first.j;
			//fprintf(stdout,"val[%d] = %g, r[%d] = %d, c[%d] = %d\n",idx,vals[idx],idx,r[idx],idx,c[idx]);
			idx++;
		}
	}

	//Construct custom sparse matrix
	//If an element does not have 4 neighbors, make coefficient 0 and set neighbor to self
	gs->aij = new double[5*gs->N_Ele];
	gs->nij = new int[4*gs->N_Ele];
	int* n_faces = new int[gs->N_Ele];
	for (int i = 0; i < 5*gs->N_Ele; i++) {
		gs->aij[i] = 0;
	}
	for (int i = 0; i < gs->N_Ele; i++) {
		gs->nij[i*4] = i;
		gs->nij[i*4+1] = i;
		gs->nij[i*4+2] = i;
		gs->nij[i*4+3] = i;
	}
	for (int i = 0; i < gs->N_Ele; i++) {
		n_faces[i] = 0;
	}
	for (int i = 0; i < N_Nonzero; i++) {
		int row = r[i];
		int col = c[i];
		if (row == col) {
			gs->aij[5*row] = vals[i];
		} else {
			gs->aij[5*row + 1 + n_faces[row]] = vals[i];
			gs->nij[4*row + n_faces[row]] = col;
			n_faces[row]++;
		}
	}

	delete[] n_faces;

/*
	//First construct coordinate-based matrix since it's easier
	cusp::coo_matrix<int,double,cusp::host_memory> Ac(N_Cells,N_Cells,N_Nonzero);

	for (int i = 0; i < N_Nonzero; i++) {
		Ac.row_indices[i] = r[i];
		Ac.column_indices[i] = c[i];
		Ac.values[i] = vals[i];
	}
*/


	//Check matrix validity
	/*fprintf(stdout,"Verifying matrix validity...\n");
	if (N_Nonzero > 0) {

		//Check validity of the row and column array order
		for (int i = 1; i < N_Nonzero; i++) {
			if (r[i] < r[i-1]) {
				fprintf(stderr,"Error: rows out of order.\n");
				exit(0);
			}
		}

		int current_row = r[0];
		int last_col = c[0];
		for (int i = 1; i < N_Nonzero; i++) {
			if (r[i] > current_row) {
				current_row = r[i];
			} else {
				if (c[i] < last_col) {
					fprintf(stderr,"Error: columns out of order.\n");
					exit(0);
				}
			}
		}

		//Check validity of row and column values
		for (int i = 0; i < N_Nonzero; i++) {
			if (r[i] < 0 || r[i] >= N_Cells) {
				fprintf(stderr,"Error: invalid row r[%d] = %d\n",i,r[i]);
				exit(0);
			}
			if (c[i] < 0 || c[i] >= N_Cells) {
				fprintf(stderr,"Error: invalid column c[%d] = %d\n",i,c[i]);
				exit(0);
			}
		}

		//Check for NaN values
		for (int i = 0; i < N_Nonzero; i++) {
			if (vals[i] != vals[i]) {
				fprintf(stderr,"Error: sparse matrix value vals[%d] = %d\n",i,vals[i]);
				exit(0);
			}
		}
	}
*/

/*
	FILE * pFile;
	if (OUTPUT_GRID_FLAG) {
		char file[255];
		sprintf(file,"sparsematrix_%d.txt",state);
		pFile = fopen (file,"w");
		if (pFile != NULL) {
			//Output grid states.
			//Each state variable's grid is written one at a time
			for (int j = 0; j < N_Nonzero; j++) {
				fprintf(pFile,"%d %d %g\n",r[j],c[j],vals[j]);
			}
			fclose(pFile);
		} else {
			fprintf(stderr,"Warning: could not open matrix file %s.\n", file);
		}
	}
*/

	//fprintf(stdout,"Constructing ELL sparse matrix... ");
	//sp_type Cc(Ac);
	//fprintf(stdout,"done.\n");

	delete[] vals;
	delete[] r;
	delete[] c;

	//Return ELL sparse matrix
	//return Cc;


}

//Device function for initializing PRNG and residuals
__global__ void InitializeRandom(curandState * randStates, int seed, double* R, int n_channels) {

	int idx = blockIdx.x*blockDim.x + threadIdx.x;

	if (idx < n_channels) {
		curand_init ( seed, idx, 0, &randStates[idx] );

		curandState localState = randStates[idx];
		R[idx] = log((double)curand_uniform( &localState ));
		randStates[idx] =localState;
	}

}

void InitializeDeviceGlobalStruct(DeviceGlobalStruct* dgs, SimData* sd, ParamStruct* ps) {

	//Initialize variables
	dgs->N_Channels = sd->N_Channels;
	dgs->N_RyR = sd->N_RyR;
	dgs->N_LCC = sd->N_LCC;

	//Allocate space on device
	cudaMalloc((void**) &dgs->Global_Vars, sizeof(double)*N_GLOBAL_VARS);
	cudaMalloc((void**) &dgs->GV_src1, sizeof(double)*N_GLOBAL_VARS);
	cudaMalloc((void**) &dgs->RyR_States, dgs->N_RyR*sizeof(int));
	cudaMalloc((void**) &dgs->RyR_Neighb, 4*dgs->N_RyR*sizeof(int));
	cudaMalloc((void**) &dgs->RyR_Ele, dgs->N_RyR*sizeof(int));
	cudaMalloc((void**) &dgs->RyR_JSR_Ele, dgs->N_RyR*sizeof(int));
	cudaMalloc((void**) &dgs->LCC_States, dgs->N_LCC*sizeof(int));
	cudaMalloc((void**) &dgs->LCCV_States, dgs->N_LCC*sizeof(int));
	cudaMalloc((void**) &dgs->LCC_Ele, dgs->N_LCC*sizeof(int));
	cudaMalloc((void**) &dgs->R, sizeof(double)*dgs->N_Channels);

	//Copy data to device
	cudaMemcpy(dgs->Global_Vars, sd->Global_Vars, sizeof(double)*N_GLOBAL_VARS, cudaMemcpyHostToDevice);
	cudaMemcpy(dgs->RyR_States, sd->RyR_States, sizeof(int)*dgs->N_RyR, cudaMemcpyHostToDevice);
	cudaMemcpy(dgs->RyR_Ele, sd->RyR_Ele, sizeof(int)*dgs->N_RyR, cudaMemcpyHostToDevice);
	cudaMemcpy(dgs->RyR_JSR_Ele, sd->RyR_JSR_Ele, sizeof(int)*dgs->N_RyR, cudaMemcpyHostToDevice);
	cudaMemcpy(dgs->RyR_Neighb, sd->RyR_Neighb, 4*sizeof(int)*dgs->N_RyR, cudaMemcpyHostToDevice);
	cudaMemcpy(dgs->LCC_States, sd->LCC_States, sizeof(int)*dgs->N_LCC, cudaMemcpyHostToDevice);
	cudaMemcpy(dgs->LCCV_States, sd->LCCV_States, sizeof(int)*dgs->N_LCC, cudaMemcpyHostToDevice);
	cudaMemcpy(dgs->LCC_Ele, sd->LCC_Ele, sizeof(int)*dgs->N_LCC, cudaMemcpyHostToDevice);
	cudaMemcpy(dgs->R, sd->R, sizeof(double)*dgs->N_Channels, cudaMemcpyHostToDevice);

	//Random numbers
    cudaMalloc ( &dgs->randStates, dgs->N_Channels*sizeof( curandState ) );

    // setup seeds
	int block_size = 32;
	int n_blocks = dgs->N_Channels/block_size + (dgs->N_Channels%block_size == 0?0:1);
	int seed = ps->Integers[I_RAND_SEED_0] + sd->device*ps->Integers[I_RAND_ADD_A] + sd->ensemble_index*ps->Integers[I_RAND_ADD_B];
	if (ps->Integers[I_FLAG_CLOCK_SEED]) {
		seed += (int)time(NULL);
	}
    InitializeRandom <<< n_blocks, block_size >>> ( dgs->randStates, seed, dgs->R, dgs->N_Channels );
    cudaDeviceSynchronize();

	cudaMemcpy(sd->R, dgs->R, sizeof(double)*dgs->N_Channels, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();

}

void OpenRandomRyR(ParamStruct* ps, SimData* sd) {

	double* RyR_Rates = new double[sd->N_RyR];
	double Rate_Total = 0;

	for (int idx = 0; idx < sd->N_RyR; idx++) {
		int N_Open = 0;
		int N_Closed = 0;
		for (int i = 0; i < 4; i++) {
			if (sd->RyR_Neighb[idx*4 + i] >= 0) {
				if (sd->RyR_States[sd->RyR_Neighb[idx*4 + i]]) {
					N_Open++;
				} else {
					N_Closed++;
				}
			}
		}

		double X, rate;
		double ca_ss = sd->Grid.States[INDEX_CA][sd->RyR_Ele[idx]];

		if (sd->RyR_States[idx]) {
			//Open -> Closed rate
			X = exp( 0.5*ps->Reals[I_RYR_A_STAR]*(N_Open*ps->Reals[I_RYR_EPS_OO] - N_Closed*ps->Reals[I_RYR_EPS_CC]) );
			rate = X*ps->Reals[I_RYR_K_MINUS];

		} else {
			//Closed -> Open rate
			double phi = ps->Reals[I_RYR_PHI_M]*sd->Grid.States[INDEX_CA][sd->RyR_JSR_Ele[idx]] + ps->Reals[I_RYR_PHI_B];
			X = exp( 0.5*ps->Reals[I_RYR_A_STAR]*(N_Closed*ps->Reals[I_RYR_EPS_CC] - N_Open*ps->Reals[I_RYR_EPS_OO]) );
			rate = X*phi*ps->Reals[I_RYR_K_PLUS]* pow(ca_ss,ps->Reals[I_RYR_ETA]);
		}
		RyR_Rates[idx] = rate;
		Rate_Total += rate;
	}

	//Randomly choose RyR to open
	int seed = 93843 + ps->Integers[I_RAND_SEED_0] + sd->device*ps->Integers[I_RAND_ADD_A] + sd->ensemble_index*ps->Integers[I_RAND_ADD_B];
	if (ps->Integers[I_FLAG_CLOCK_SEED]) {
		seed += (int)time(NULL);
	}
	srand(seed);
	double r = Rate_Total*((double)rand())/((double)RAND_MAX);
	int RyR_Open = 0;
	double sum = 0;
	for (int i = 0; i < sd->N_RyR; i++) {
		sum += RyR_Rates[i];
		if (sum >= r) {
			RyR_Open = i;
			break;
		}
	}
	fprintf(stdout,"Opening RyR #%d...\n",RyR_Open);
	sd->RyR_States[RyR_Open] = 1;
	delete[] RyR_Rates;
}
