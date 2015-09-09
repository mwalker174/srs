/*
 * main.cpp
 *
 *  Created on: Jul 13, 2012
 *      Author: mwalker
 */

#include <iostream>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <string.h>
#include <stdexcept>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/sax/HandlerBase.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include "parameters.h"

using namespace std;
using namespace xercesc;

ParamStruct ps;
MeshStruct ms;

int main(int argc, char **argv) {

	//Parameters file location
	char* parameters_file = (char*)DEFAULT_PARAM_FILE;

	//Number of scan lines
	int scanlines_option = -1;
	int bSetFilename = 0;
	char* filename;
	int index_option = -1;
	int bSetOutputname = 0;
	char* outname;

	//Read command line options
	for (int i = 1; i < argc; i += 2) {
		if (!strcmp(argv[i],"-param")) {
			parameters_file = argv[i+1];
		} else if (!strcmp(argv[i],"-n")) {
			scanlines_option = atoi(argv[i+1]);
		} else if (!strcmp(argv[i],"-f")) {
			bSetFilename = 1;
			filename = argv[i+1];
		} else if (!strcmp(argv[i],"-i")) {
			index_option = atoi(argv[i+1]);
		} else if (!strcmp(argv[i],"-out")) {
			bSetOutputname = 1;
			outname = argv[i+1];
		} else {
			fprintf(stderr,"Error: unrecognized command line option: %s\n",argv[i]);
			exit(1);
		}
	}

	//Read in parameters
	ReadParameters(parameters_file);

	//Override number of scan lines if specified at the command line
	if (scanlines_option > 0) {
		ps.Integers[I_N_SCANS] = scanlines_option;
	}
	if (bSetFilename) {
		sprintf(ps.Chars[I_FILE_BASE],"%s",filename);
	}
	if (bSetOutputname) {
		sprintf(ps.Chars[I_FILE_CA],"%s",outname);
		sprintf(ps.Chars[I_FILE_CAF],"%s",outname);
		sprintf(ps.Chars[I_FILE_CAFB],"%s",outname);
	}
	if (index_option >= 0) {
		sprintf(ps.Chars[I_FILE_BASE],"%s_%d_0",ps.Chars[I_FILE_BASE],index_option);
		sprintf(ps.Chars[I_FILE_CA],"%s_%d.txt",ps.Chars[I_FILE_CA],index_option);
		sprintf(ps.Chars[I_FILE_CAF],"%s_%d.txt",ps.Chars[I_FILE_CAF],index_option);
		sprintf(ps.Chars[I_FILE_CAFB],"%s_%d.txt",ps.Chars[I_FILE_CAFB],index_option);
	}

	//Read in geometry
	int N_Points, N_Cells;
	Vector3D *points;
	int** cells;
	ReadGridGeometry(0, N_Points, points, N_Cells, cells);
	int N_Cells_Mesh = N_Cells; //Size of cells array

	//Compute cell volumes and centers, copy states
	fprintf(stdout,"Adding computing cell volumes and centers...\n");
	ms.CellCenters = vector<Vector3D>(N_Cells);
	ms.V0 = vector<double>(N_Cells);
	ms.C_Ca = vector<double>(N_Cells);
	ms.C_CaF = vector<double>(N_Cells);
	for (int i = 0; i < N_Cells; i++) {
		Vector3D c1, c2, c3, c4;
		c1 = points[cells[i][0]];
		c2 = points[cells[i][1]];
		c3 = points[cells[i][2]];
		c4 = points[cells[i][3]];
		ms.CellCenters[i].x = (c1.x + c2.x + c3.x + c4.x)/4.0;
		ms.CellCenters[i].y = (c1.y + c2.y + c3.y + c4.y)/4.0;
		ms.CellCenters[i].z = (c1.z + c2.z + c3.z + c4.z)/4.0;

		Vector3D v1, v2, v3;
		v1 = vectorMinus(c1,c4);
		v2 = vectorMinus(c2,c4);
		v3 = vectorMinus(c3,c4);
		ms.V0[i] = fabs( dot(v1,cross(v2.x,v2.y,v2.z,v3.x,v3.y,v3.z))/6.0 );
	}

	//Add padding to the domain
	fprintf(stdout,"Adding domain padding...\n");

	double pad_min = ps.Reals[I_PADDING_MIN];
	double pad_max = ps.Reals[I_PADDING_MAX];
	double domain_x_min = -ps.Reals[I_DOMAIN_X]/2;
	double domain_x_max = ps.Reals[I_DOMAIN_X]/2;
	double domain_y_min = -ps.Reals[I_DOMAIN_Y]/2;
	double domain_y_max = ps.Reals[I_DOMAIN_Y]/2;
	double domain_z_min = -ps.Reals[I_DOMAIN_Z]/2;
	double domain_z_max = ps.Reals[I_DOMAIN_Z]/2;

	AddPointGrid(pad_min,pad_max, pad_min,pad_max, domain_z_max,pad_max);
	AddPointGrid(pad_min,pad_max, pad_min,pad_max, domain_z_min,pad_min);

	AddPointGrid(domain_x_min,pad_min, pad_min,pad_max, domain_z_min,domain_z_max);
	AddPointGrid(domain_x_max,pad_max, pad_min,pad_max, domain_z_min,domain_z_max);

	AddPointGrid(domain_x_min,domain_x_max, domain_y_min,pad_min, domain_z_min,domain_z_max);
	AddPointGrid(domain_x_min,domain_x_max, domain_y_max,pad_max, domain_z_min,domain_z_max);

	//Compute linescan coordinates
	int N_X_F = (int)((ps.Reals[I_X_MAX]-ps.Reals[I_X_MIN])/ps.Reals[I_RES_X_F]);
	int N_X_Ca = (int)((ps.Reals[I_X_MAX]-ps.Reals[I_X_MIN])/ps.Reals[I_RES_X_Ca]);
	double* X_F = new double[N_X_F];
	double* X_Ca = new double[N_X_Ca];
	for (int i = 0; i < N_X_F; i++) {
		X_F[i] = ps.Reals[I_X_MIN] + i*ps.Reals[I_RES_X_F];
	}
	for (int i = 0; i < N_X_Ca; i++) {
		X_Ca[i] = ps.Reals[I_X_MIN] + i*ps.Reals[I_RES_X_Ca];
	}
	int N_Scans = ps.Integers[I_N_SCANS];
	double** I_Ca = new double*[N_Scans];
	double** I_CaF = new double*[N_Scans];
	double** I_CaFB = new double*[N_Scans];
	for (int i = 0; i < N_Scans; i++) {
		I_Ca[i] = new double[N_X_Ca];
		I_CaF[i] = new double[N_X_F];
		I_CaFB[i] = new double[N_X_F];
	}

	//Y and Z offsets
	double Y0 = ps.Reals[I_Y_OFFSET];
	double Z0 = ps.Reals[I_Z_OFFSET];

	//PSF parameters
	double sigma_x, sigma_y, sigma_z;

	//Main computation
	fprintf(stdout,"Computing linescans...\n");
	Vector3D center_0;
	N_Cells = ms.C_CaF.size();

	int axis; //AXIS_X, AXIS_Y, OR AXIS_Z
	if (!strcmp("X",ps.Chars[I_SCAN_AXIS])) {
		center_0.y = Y0;
		center_0.z = Z0;
		sigma_x = sqrt(ps.Reals[I_VAR_X]);
		sigma_y = sqrt(ps.Reals[I_VAR_Y]);
		sigma_z = sqrt(ps.Reals[I_VAR_Z]);
		axis = AXIS_X;
	} else if (!strcmp("Y",ps.Chars[I_SCAN_AXIS])) {
		center_0.x = Y0;
		center_0.z = Z0;
		sigma_x = sqrt(ps.Reals[I_VAR_X]);
		sigma_y = sqrt(ps.Reals[I_VAR_Y]);
		sigma_z = sqrt(ps.Reals[I_VAR_Z]);
		axis = AXIS_Y;
	} else if (!strcmp("Z",ps.Chars[I_SCAN_AXIS])) {
		center_0.x = Y0;
		center_0.y = Z0;
		sigma_x = sqrt(ps.Reals[I_VAR_Z]);
		sigma_y = sqrt(ps.Reals[I_VAR_X]);
		sigma_z = sqrt(ps.Reals[I_VAR_Y]);
		axis = AXIS_Z;
	} else {
		fprintf(stderr,"Error: invalid scan axis: %s\n",ps.Chars[I_SCAN_AXIS]);
		exit(0);
	}
	double r_thresh = 1; //9*max(ps.Reals[I_VAR_X],max(ps.Reals[I_VAR_Y],ps.Reals[I_VAR_Z]));
	for (int k = 0; k < N_Scans; k++) {

		ReadGridGStates(k*ps.Integers[I_N_SCANS_SKIP], N_Cells_Mesh, ms.C_Ca, ms.C_CaF);

		for (int i = 0; i < N_X_F; i++) {

			if (axis == AXIS_X) {
				center_0.x = X_F[i];
			} else if (axis == AXIS_Y) {
				center_0.y = X_F[i];
			} else {
				center_0.z = X_F[i];
			}

			double sum = 0;
			Vector3D r = vectorMinus(ms.CellCenters[0],center_0);
			double r_min = r.x*r.x + r.y*r.y + r.z*r.z;
			int cell_min = 0;
			for (int j = 0; j < N_Cells; j++) {
				r = vectorMinus(ms.CellCenters[j],center_0);
				double dist = r.x*r.x + r.y*r.y + r.z*r.z;
				if (dist < r_thresh) {
					double w = normpdf3d(sigma_x, sigma_y, sigma_z, r);
					sum += ms.V0[j] * ms.C_CaF[j] * w;
				}
				if (dist < r_min) {
					r_min = dist;
					cell_min = j;
				}
			}
			I_CaFB[k][i] = sum;
			I_CaF[k][i] = ms.C_CaF[cell_min];
		}

		for (int i = 0; i < N_X_Ca; i++) {

			if (axis == AXIS_X) {
				center_0.x = X_Ca[i];
			} else if (axis == AXIS_Y) {
				center_0.y = X_Ca[i];
			} else {
				center_0.z = X_Ca[i];
			}
			Vector3D r = vectorMinus(ms.CellCenters[0],center_0);
			double r_min = r.x*r.x + r.y*r.y + r.z*r.z;
			int cell_min = 0;
			for (int j = 0; j < N_Cells; j++) {
				r = vectorMinus(ms.CellCenters[j],center_0);
				double dist = r.x*r.x + r.y*r.y + r.z*r.z;
				if (dist < r_min) {
					r_min = dist;
					cell_min = j;
				}
			}
			I_Ca[k][i] = ms.C_Ca[cell_min];
		}
	}

	fprintf(stdout,"Saving to disk...\n");
	WriteOutputFile(I_Ca,N_X_Ca,N_Scans,ps.Chars[I_FILE_CA]);
	WriteOutputFile(I_CaF,N_X_F,N_Scans,ps.Chars[I_FILE_CAF]);
	WriteOutputFile(I_CaFB,N_X_F,N_Scans,ps.Chars[I_FILE_CAFB]);

	delete[] X_F;
	delete[] X_Ca;
	for (int i = 0; i < N_Scans; i++) {
		delete[] I_Ca[i];
		delete[] I_CaF[i];
		delete[] I_CaFB[i];
	}
	delete[] I_Ca;
	delete[] I_CaF;
	delete[] I_CaFB;
	delete[] points;
	for (int i = 0; i < N_Cells_Mesh; i++) {
		delete[] cells[i];
	}
	delete[] cells;

	return 0;
}

void AddPointGrid( double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) {

	//Grid spacing
	double dx = ps.Reals[I_PADDING_RES];
	double dy = ps.Reals[I_PADDING_RES];
	double dz = ps.Reals[I_PADDING_RES];

	//Space the grid in correct directions if min and max are flipped (always start at min and work outwards)
	if (xmin > xmax)
		dx = -dx;
	if (ymin > ymax)
		dy = -dy;
	if (zmin > zmax)
		dz = -dz;

	//Number of elements in each dimension
	int N_X = (int)((xmax-xmin)/dx);
	int N_Y = (int)((ymax-ymin)/dy);
	int N_Z = (int)((zmax-zmin)/dz);

	fprintf(stdout,"Adding %d padding elements...\n",N_X*N_Y*N_Z);

	//Create grid
	Vector3D v;
	for (int i = 0; i < N_X; i++) {
		for (int j = 0; j < N_Y; j++) {
			for (int k = 0; k < N_Z; k++) {
				v.x = xmin + i*dx + 0.5*dx;
				v.y = ymin + j*dy + 0.5*dy;
				v.z = zmin + k*dz + 0.5*dz;
				ms.CellCenters.push_back(v);
				ms.V0.push_back(fabs(dx*dy*dz));
				ms.C_Ca.push_back(ps.Reals[I_C0]);
				ms.C_CaF.push_back(ps.Reals[I_F0]);
			}
		}
	}

}


void WriteOutputFile(double** I, int N_I, int N_Scans, char* filename) {

	//Write linescan to file
	FILE * pFile;
	fprintf(stdout,"Writing output file %s\n",filename);
	pFile = fopen (filename,"w");
	if (pFile != NULL) {

		for (int j = 0; j < N_Scans; j++) {
			for (int i = 0; i < N_I; i++) {
				fprintf(pFile,"%g ",I[j][i]);
			}
			fprintf(pFile,"\n");
		}

		fclose(pFile);

	} else {
		fprintf(stderr,"Warning: could not open output file.\n");
	}

}

//Reads parameters
void ReadParameters(const char* param_file) {

	//Assign default parameters
	ps.Integers[I_N_SCANS] = 26;
	ps.Integers[I_N_SCANS_SKIP] = 2;
	ps.Reals[I_RES_X_F] = 0.15;
	ps.Reals[I_RES_X_Ca] = 0.01;
	ps.Reals[I_RES_T] = 2.0;
	ps.Reals[I_X_MIN] = -4.0;
	ps.Reals[I_X_MAX] = 4.0;
	ps.Reals[I_Y_OFFSET] = 0.1075;
	ps.Reals[I_Z_OFFSET] = 0;
	ps.Reals[I_DOMAIN_X] = 4;
	ps.Reals[I_DOMAIN_Y] = 4;
	ps.Reals[I_DOMAIN_Z] = 4;
	ps.Reals[I_PADDING_MIN] = -5.0;
	ps.Reals[I_PADDING_MAX] = 5.0;
	ps.Reals[I_PADDING_RES] = 0.3;
	ps.Reals[I_VAR_X] = 0.0289;
	ps.Reals[I_VAR_Y] = 0.1154;
	ps.Reals[I_VAR_Z] = 0.0289;
	ps.Reals[I_F0] = 4.1667;
	ps.Reals[I_C0] = 0.10;
	sprintf(ps.Chars[I_FILE_BASE],"../output/cru3d_grid_0_0");
	sprintf(ps.Chars[I_FILE_CA],"scan_output_Ca.txt");
	sprintf(ps.Chars[I_FILE_CAF],"scan_output_CaF.txt");
	sprintf(ps.Chars[I_FILE_CAFB],"scan_output.txt");
	sprintf(ps.Chars[I_SCAN_AXIS],"Z");

	//Initialize XML parser
	try {
        XMLPlatformUtils::Initialize();
    }
    catch (const XMLException& toCatch) {
        char* message = XMLString::transcode(toCatch.getMessage());
        cout << "Error during initialization! :\n"
             << message << "\n";
        XMLString::release(&message);
        exit(1);
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
        exit(1);
    }
    catch (const DOMException& toCatch) {
        char* message = XMLString::transcode(toCatch.msg);
        cout << "Exception message is: \n"
             << message << "\n";
        XMLString::release(&message);
        exit(1);
    }
    catch (...) {
        cout << "Unexpected Exception \n" ;
        exit(1);
    }
    
    DOMDocument* xmlDoc = parser->getDocument();
    DOMElement* elementRoot = xmlDoc->getDocumentElement();
    if (!elementRoot) {
    	throw(std::runtime_error("empty XML document"));
        exit(1);
    }
    
    DOMNodeList* children = elementRoot->getChildNodes();
    XMLSize_t nodeCount = children->getLength();
    
    XMLCh* TAG_ROOT = XMLString::transcode("analysis");
    
    if (!XMLString::equals(elementRoot->getTagName(), TAG_ROOT)) {
    	throw(std::runtime_error("Invalid XML root"));
        exit(1);
    }
    children = elementRoot->getChildNodes();
    nodeCount = children->getLength();
    
    XMLCh* TAG_PARAMETER = XMLString::transcode("parameter");
    XMLCh* TAG_SYMBOL = XMLString::transcode("symbol");
    XMLCh* TAG_VALUE = XMLString::transcode("value");
    
    char RyR_Cluster[1024];
    char LCC_Cluster[1024];
    
    for (XMLSize_t xx = 0; xx < nodeCount; ++xx) {
    	DOMNode* currentNode = children->item(xx);
    	if (currentNode->getNodeType() && currentNode->getNodeType() == DOMNode::ELEMENT_NODE) {
    		DOMElement* currentElement = dynamic_cast< xercesc::DOMElement* >(currentNode);
    		if (XMLString::equals(currentElement->getTagName(), TAG_PARAMETER)) {
    		
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
    			
    			if (!strcmp(symbol,"scan_n")) {
    				ps.Integers[I_N_SCANS] = atoi(value);
    			} else if (!strcmp(symbol,"scan_n_per_line")) {
    				ps.Integers[I_N_SCANS_SKIP] = atoi(value);
    			} else if (!strcmp(symbol,"scan_x_res_f")) {
    				ps.Reals[I_RES_X_F] = atof(value);
    			} else if (!strcmp(symbol,"flag_x_res_ca")) {
    				ps.Integers[I_RES_X_Ca] = atoi(value);
    			} else if (!strcmp(symbol,"scan_x_min")) {
    				ps.Reals[I_X_MIN] = atof(value);
    			} else if (!strcmp(symbol,"scan_x_max")) {
    				ps.Reals[I_X_MAX] = atof(value);
    			} else if (!strcmp(symbol,"scan_offset_1")) {
    				ps.Reals[I_Y_OFFSET] = atof(value);
    			} else if (!strcmp(symbol,"scan_offset_2")) {
    				ps.Reals[I_Z_OFFSET] = atof(value);
    			} else if (!strcmp(symbol,"scan_domain_x")) {
    				ps.Reals[I_DOMAIN_X] = atof(value);
    			} else if (!strcmp(symbol,"scan_domain_y")) {
    				ps.Reals[I_DOMAIN_Y] = atof(value);
    			} else if (!strcmp(symbol,"scan_domain_z")) {
    				ps.Reals[I_DOMAIN_Z] = atof(value);
    			} else if (!strcmp(symbol,"scan_padding_min")) {
    				ps.Reals[I_PADDING_MIN] = atof(value);
    			} else if (!strcmp(symbol,"scan_padding_max")) {
    				ps.Reals[I_PADDING_MAX] = atof(value);
    			} else if (!strcmp(symbol,"scan_padding_res")) {
    				ps.Reals[I_PADDING_RES] = atof(value);
    			} else if (!strcmp(symbol,"scan_var_x")) {
    				ps.Reals[I_VAR_X] = atof(value);
    			} else if (!strcmp(symbol,"scan_var_y")) {
    				ps.Reals[I_VAR_Y] = atof(value);
    			} else if (!strcmp(symbol,"scan_var_z")) {
    				ps.Reals[I_VAR_Z] = atof(value);
    			} else if (!strcmp(symbol,"scan_f0")) {
    				ps.Reals[I_F0] = atof(value);
    			} else if (!strcmp(symbol,"scan_c0")) {
    				ps.Reals[I_C0] = atof(value);
    			} else if (!strcmp(symbol,"scan_gridfile_base")) {
    				sprintf(ps.Chars[I_FILE_BASE],"%s",value);
    			} else if (!strcmp(symbol,"outdir")) {
    				sprintf(ps.Chars[I_FILE_CA],"%s/scan_output_Ca.txt",value);
    				sprintf(ps.Chars[I_FILE_CAF],"%s/scan_output_CaF.txt",value);
    				sprintf(ps.Chars[I_FILE_CAFB],"%s/scan_output.txt",value);
    			} else if (!strcmp(symbol,"scan_axis")) {
    				sprintf(ps.Chars[I_SCAN_AXIS],"%s",value);
    			}
    			
    		}
    	}
    }
    
    delete parser;
    delete errHandler;

	for (int i = 0; i < N_PARAMETERS_INT; i++) {
		fprintf(stdout,"Integer Parameter[%d] = %d\n",i,ps.Integers[i]);
	}
	for (int i = 0; i < N_PARAMETERS_REAL; i++) {
		fprintf(stdout,"Real Parameter[%d] = %g\n",i,ps.Reals[i]);
	}
	for (int i = 0; i < N_PARAMETERS_CHAR; i++) {
		fprintf(stdout,"Char Parameter[%d] = %s\n",i,ps.Chars[i]);
	}
}

//Reads grid file
void ReadGridGeometry(int i, int &N_Points, Vector3D *&points, int &N_Cells, int** &cells) {

	FILE* pFile;
	char filename[255];
	sprintf(filename,"%s_%d.vtk",ps.Chars[I_FILE_BASE],i);
	pFile = fopen(filename,"r");
	char line[255];

	fprintf(stdout,"Reading grid file %s...\n",filename);
	if (pFile != NULL) {

		//Skip to points section
		fprintf(stdout,"Skipping to points section...\n");
		while (strcmp(fgets(line,255,pFile),"DATASET UNSTRUCTURED_GRID\n")) {
		}

		fprintf(stdout,"Reading number of points...\n");
		if (!fscanf(pFile,"POINTS %d float\n",&N_Points)) {
			fprintf(stderr,"Error reading number of points.\n");
			exit(1);
		}

		fprintf(stdout,"Reading points data...\n");
		points = new Vector3D[N_Points];
		for (int i = 0; i < N_Points; i++) {
			if (!fscanf(pFile,"%lg %lg %lg",&(points[i].x),&(points[i].y),&(points[i].z))) {
				fprintf(stderr,"Error reading points data.\n");
				exit(1);
			}
		}

		//Skip line
		fgets(line,255,pFile);
		fgets(line,255,pFile);

		//Get cells data
		if (!fscanf(pFile,"CELLS %d %*d\n",&N_Cells)) {
			fprintf(stderr,"Error reading number of cells.\n");
			exit(1);
		}

		fprintf(stdout,"Reading cells data...\n");
		cells = new int*[N_Cells];
		for (int i = 0; i < N_Cells; i++) {
			cells[i] = new int[4];
			if (!fscanf(pFile,"4 %d %d %d %d\n",cells[i],cells[i]+1,cells[i]+2,cells[i]+3)) {
				fprintf(stderr,"Error reading cells data.\n");
				exit(1);
			}
		}

		fclose(pFile);
	} else {
		fprintf(stderr,"Error opening file: %s\n",filename);
		exit(1);
	}

}

//Reads grid file
void ReadGridGStates(int i, int N_Cells, vector<double> &states_ca, vector<double> &states_caf) {

	FILE* pFile;
	char filename[255];
	sprintf(filename,"%s_%d.vtk",ps.Chars[I_FILE_BASE],i);
	pFile = fopen(filename,"r");
	char line[255];

	fprintf(stdout,"Reading grid file %s...\n",filename);
	if (pFile != NULL) {

		//Skip to Ca section
		while (strcmp(fgets(line,255,pFile),"SCALARS C_CA float\n")) {
		}

		//Skip line
		fgets(line,255,pFile);

		//Get cell states
		for (int i = 0; i < N_Cells; i++) {
			if (!fscanf(pFile,"%lg\n",&states_ca[i])) {
				fprintf(stderr,"Error reading cells data.\n");
				exit(1);
			}
		}

		//Skip to CaF section
		while (strcmp(fgets(line,255,pFile),"SCALARS C_DYE float\n")) {
		}

		//Skip line
		fgets(line,255,pFile);

		//Get cell states
		for (int i = 0; i < N_Cells; i++) {
			if (!fscanf(pFile,"%lg\n",&states_caf[i])) {
				fprintf(stderr,"Error reading cells data.\n");
				exit(1);
			}
		}

		fclose(pFile);
	} else {
		fprintf(stderr,"Error opening file: %s\n",filename);
		exit(1);
	}

}
