/*
 * parameters.h
 *
 *  Created on: Jul 13, 2012
 *      Author: mwalker
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <vector>
#include <math.h>

using namespace std;

//Number of input parameters of each type
const int N_PARAMETERS_INT = 7;
const int N_PARAMETERS_REAL = 12;
const int N_PARAMETERS_CHAR = 3;

const int I_N_TT_SEGMENTS = 0; //Number of sides on TT, each of length DYAD_RES
const int I_N_TT_SEGMENTS_Z = 1; //Number of vertical segments on TT, each of length DYAD_RES
const int I_N_JSR_SEGMENTS = 2; //Number of segments that make the width of the JSR
const int I_N_JSR_SEGMENTS_Z = 3; //Number of segments that make the height of the JSR
const int I_FLAG_TT = 4; //Flag to include TT
const int I_FLAG_DYAD_GUIDE = 5; //Flag to guidepoints 15nm above SR
const int I_FLAG_TWO_CRUS = 6; //Flag to add second identical CRU
const int I_VMAX = 0; //Maximum element volume
const int I_QMAX = 1; //Maximum Q ratio (quality metric, should be between 1.414 and 2)
const int I_GRID_X = 2; //Domain size in X
const int I_GRID_Y = 3; //Domain size in Y
const int I_GRID_Z = 4; //Domain size in Z
const int I_DYAD_R = 5; //Dyad height
const int I_JSR_R = 6; //JSR depth
const int I_DYAD_RES = 7; //Dyad resolution, recommended 30nm
const int I_R_TRPN = 8; //Radius beyond which Troponin and SERCA are present
const int I_RYR_SPACING = 9; //Spacing between adjacent RyRs
const int I_JSR_TRANS = 10; //X translation for JSR
const int I_CRU_SEPARATION = 11; //If FLAG_TWO_CRUS, specifies distance between them
const int I_OUTPUT_NAME = 0; //Output file base
const int I_RYR_FILE = 1; //RyR placement file
const int I_LCC_FILE = 2; //LCC placement file

/*
int N_RYR_X = 7; //Number of RyRs in horizontal direction
int N_RYR_Y = 7; //Number of RyRs in vertical direction (along the T-Tubule)
int N_LCC_X = 3; //Number of LCCs in horizontal direction
int N_LCC_Y = 3; //Number of LCCs in vertical direction
const int N_TT_SEGMENTS = 21; //Number of sides on TT, each of length DYAD_RES
const int N_TT_SEGMENTS_Z = 33; //Number of vertical segments on TT, each of length DYAD_RES
const int N_JSR_SEGMENTS = 17; //Number of segments that make the width of the JSR
const int N_JSR_SEGMENTS_Z = 13; //Number of segments that make the height of the JSR
const double VMAX = 0.005; //Maximum element volume
const double QMAX = 2; //Maximum Q ratio (quality metric, should be between 1.414 and 2)
const double GRID_X = 4; //Domain size in X
const double GRID_Y = 4; //Domain size in Y
const double GRID_Z = 4; //Domain size in Z
const double DYAD_R = 15e-3; //Dyad height
const double JSR_R = 30e-3; //JSR depth
const double DYAD_RES = 30e-3; //Dyad resolution, recommended 30nm
const double R_TRPN = 0.2; //Radius beyond which Troponin and SERCA are present
const double RYR_SPACING = 30e-3; //Spacing between adjacent RyRs
const char OUTPUT_NAME[255] = "../mesh"; //Output file base
*/

//Other constants
const double PI = 3.14159;

//Face attributes
const int ATTR_NONE = 0;
const int ATTR_BOUNDARY = 1;
const int ATTR_TT = 2;
const int ATTR_JSR = 3;

//Region attributes
const int REGION_CYTO = 0;
const int REGION_JSR = 1;

//Default parameters file
const char DEFAULT_PARAM_FILE[255] = "DefaultMesh.xml";

struct Vector3D {
	double x, y, z;
};

struct IndexPair {
	int i,j;
};

inline double norm(Vector3D v) { return sqrt(v.x*v.x + v.y*v.y + v.z*v.z); };
inline void normalize(Vector3D &v) { double n = norm(v); v.x /= n; v.y /= n; v.z /= n; };
inline double dot(Vector3D v, Vector3D w) { return v.x*w.x + v.y*w.y + v.z*w.z; };
inline double distance3D(Vector3D v, Vector3D w) { return sqrt((v.x-w.x)*(v.x-w.x) + (v.y-w.y)*(v.y-w.y) + (v.z-w.z)*(v.z-w.z)); }

Vector3D negative(Vector3D v) {
	Vector3D c;
	c.x = -v.x;
	c.y = -v.y;
	c.z = -v.z;
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
	double JSR_Z;
	double TT_RADIUS; //Circumradius of TT
	double CRU_1_Z; //Z position of first CRU
	double CRU_2_Z; //Z position of second CRU

	int N_RyR_X;
	int N_RyR_Y;
	int N_RyR;
	int** bRyR; //2D array of RyR placement
	int N_LCC_X;
	int N_LCC_Y;
	int N_LCC;
	int** bLCC; //2D array of LCC placement

};

struct GeometryStruct {

	vector<double> VerticesX; //Vertex coordinates
	vector<double> VerticesY;
	vector<double> VerticesZ;

	vector<double> HolesX; //3D hole coordinates
	vector<double> HolesY;
	vector<double> HolesZ;

	vector< vector<int> > Facets;

	vector<vector<int> > Polys; //Polygon vertex list
	vector<int> PolyLabels; //Polygon attribute labels
	vector<int> PolyCount; //Number of polygons
	vector<int> PolyHolesCount; //Number of polygons with holes
	vector<double> PolyHolesX; //Polygon hole coordinates
	vector<double> PolyHolesY;
	vector<double> PolyHolesZ;

	vector<double> RegionsX; //Region coordinates
	vector<double> RegionsY;
	vector<double> RegionsZ;
	vector<int> RegionAttr; //Region attributes

};

struct MeshStruct {

	//Cell properties
	vector<Vector3D> CellCenters;
	vector<double> V0;
	vector<int> bTroponinC;
	vector<double> TTSurfaceArea;
	vector<double> SRSurfaceArea;
	vector<int> CellAttrList;

	//Face properties
	vector<IndexPair> FaceTet;
	vector<Vector3D> FaceNormals;
	vector<Vector3D> FaceCenters;
	vector<Vector3D> Af;
	vector<double> FaceArea;
	vector<Vector3D> r1;
	vector<Vector3D> r2;
	vector<Vector3D> eXi;
	vector<double> dXi;
	vector<double> FaceCoeff;
	vector<int> FaceAttrList;
	vector<int> bElectroDiff;
};

struct InterfaceStruct {

	//Cell properties
	vector<IndexPair> cell_pairs;

	//Face properties
	vector<IndexPair> face_pairs;
	vector<Vector3D> FaceNormals;
	vector<Vector3D> FaceCenters;
	vector<Vector3D> Af;
	vector<double> FaceArea;
	vector<Vector3D> r1;
	vector<Vector3D> r2;
	vector<Vector3D> eXi;
	vector<double> dXi;
	vector<double> FaceCoeff;
};


//Function declarations
void ReadParameters(ParamStruct* ps, const char* param_file);
void ReadChannelGeometry(char* filename, int** &bChannel, int &N_Channels, int &N_Dyad_X, int &N_Dyad_Y);
void WritePropertiesFile(MeshStruct* ms, const char* filebase);
void WritePolyFile(vector<double> &VX, vector<double> &VY, vector<int> &C1, vector<int> &C2, vector<double> &HX, vector<double> &HY, vector<int> &A);
void WriteChannelFiles(vector<Vector3D> &CellCenters,vector<int> &CellAttrList);

void InitializeVertices(GeometryStruct* gs);

void AddFace(GeometryStruct* gs, int attribute,
			 double x1, double y1, double z1,
			 double x2, double y2, double z2,
			 double x3, double y3, double z3,
			 double x4, double y4, double z4);

void AddBox(GeometryStruct* gs, int attribute, double x_center, double y_center, double z_center, double dx, double dy, double dz);
void AddSegmentedTT(GeometryStruct* gs, int Poly_Attr, double radius, int bCapped);
void AddJSR(GeometryStruct* gs, double z0);
void JSRHelper(GeometryStruct* gs, vector<int> v_stitch_in, vector<int> v_stitch_out,
	  	   vector<int> &stitch_next_in, vector<int> &stitch_next_out, double z0, int bCapTop);
void AddDyadGuide(GeometryStruct* gs, double z_center);
void ChannelHelper(vector<Vector3D> &CellCenters, vector<int> &CellAttrList,
				   vector<int> &RyR_Ele,vector<int> &RyR_JSR_Ele,vector<int> &LCC_Ele,
				   vector<int> &RYR_i, double z_center);

#endif /* PARAMETERS_H_ */
