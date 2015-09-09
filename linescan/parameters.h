/*
 * parameters.h
 *
 *  Created on: Feb 6, 2013
 *      Author: mwalker
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <vector>
#include <math.h>
using namespace std;

//Number of input parameters of each type
const int N_PARAMETERS_INT = 2;
const int N_PARAMETERS_REAL = 18;
const int N_PARAMETERS_CHAR = 5;

const int I_N_SCANS = 0; //Number of scan lines
const int I_N_SCANS_SKIP = 1; //Number of files to skip between scans
const int I_RES_X_F = 0; //Spatial resolution for fluorescence
const int I_RES_X_Ca = 1; //Spatial resolution for Ca concentration
const int I_RES_T = 2; //Time resolution
const int I_X_MIN = 3; //Starting X location
const int I_X_MAX = 4; //Ending X location
const int I_Y_OFFSET = 5; //Y offset from origin
const int I_Z_OFFSET = 6; //Z offset from origin
const int I_DOMAIN_X = 7; //Domain size in X
const int I_DOMAIN_Y = 8; //Domain size in Y
const int I_DOMAIN_Z = 9; //Domain size in Z
const int I_PADDING_MIN = 10; //Domain padding min coordinate
const int I_PADDING_MAX = 11; //Domain padding max coordinate
const int I_PADDING_RES = 12; //Domain padding spacing
const int I_VAR_X = 13; //PSF variance in X direction
const int I_VAR_Y = 14; //PSF variance in Y direction
const int I_VAR_Z = 15; //PSF variance in Z direction
const int I_F0 = 16; //Baseline signal
const int I_C0 = 17; //Baseline Ca signal
const int I_FILE_BASE = 0; //Input file base
const int I_FILE_CA = 1; //Output file name for Ca
const int I_FILE_CAF = 2; //Output file name for CaF
const int I_FILE_CAFB = 3; //Output file name for CaF (blurred)
const int I_SCAN_AXIS = 4; //Axis along which to scan (X, Y, or Z)

//Default parameters file
const char DEFAULT_PARAM_FILE[255] = "DefaultScan.xml";

struct Vector3D {
	double x, y, z;
};

const double PI = 3.14159;
const double PDF3D_ALPHA = 1.0 / pow(2*PI,1.5);
const int AXIS_X = 0;
const int AXIS_Y = 1;
const int AXIS_Z = 2;

inline double norm(Vector3D v) { return sqrt(v.x*v.x + v.y*v.y + v.z*v.z); };
inline void normalize(Vector3D &v) { double n = norm(v); v.x /= n; v.y /= n; v.z /= n; };
inline double dot(Vector3D v, Vector3D w) { return v.x*w.x + v.y*w.y + v.z*w.z; };
inline double distance3D(Vector3D v, Vector3D w) { return sqrt((v.x-w.x)*(v.x-w.x) + (v.y-w.y)*(v.y-w.y) + (v.z-w.z)*(v.z-w.z)); }
inline double normpdf3d(double sigma_x, double sigma_y, double sigma_z, Vector3D r) {
	return PDF3D_ALPHA * exp( -0.5*( (r.x*r.x/(sigma_x*sigma_x)) + (r.y*r.y/(sigma_y*sigma_y))  + (r.z*r.z/(sigma_z*sigma_z)) ) ) / (sigma_x*sigma_y*sigma_z);
};

Vector3D negative(Vector3D v) {
	Vector3D c;
	c.x = -v.x;
	c.y = -v.y;
	c.z = -v.z;
	return c;
}

Vector3D vectorMinus(Vector3D v1, Vector3D v2) {
	Vector3D c;
	c.x = v1.x - v2.x;
	c.y = v1.y - v2.y;
	c.z = v1.z - v2.z;
	return c;
}

Vector3D cross(double x1, double y1, double z1, double x2, double y2, double z2) {
	Vector3D c;
	c.x = y1*z2 - z1*y2; //Cross product of two vectors in the face's plane
	c.y = z1*x2 - x1*z2;
	c.z = x1*y2 - y1*x2;
	return c;
}

//Parameters data structure
struct ParamStruct {

	int Integers[N_PARAMETERS_INT];
	double Reals[N_PARAMETERS_REAL];
	char Chars[N_PARAMETERS_CHAR][255];

};

struct MeshStruct {

	//Cell properties
	vector<Vector3D> CellCenters; //Cell centers
	vector<double> V0; //Cell volumes
	vector<double> C_Ca; //Cell states (Ca)
	vector<double> C_CaF; //Cell states (Ca-Dye)

};

//Function declarations
void AddPointGrid(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);

void ReadParameters(const char* param_file);
void ReadGridGeometry(int i, int &N_Points, Vector3D *&points, int &N_Cells, int** &cells);
void ReadGridGStates(int i, int N_Cells, vector<double> &states_ca, vector<double> &states_caf);
void WriteOutputFile(double** I, int N_I, int N_Scans, char* filename);

#endif /* PARAMETERS_H_ */
