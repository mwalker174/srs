/*
 * main.cpp
 *
 *  Created on: Jul 13, 2012
 *      Author: mwalker
 */

#include <iostream>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <stdexcept>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/sax/HandlerBase.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include "tetgen.h"
#include "parameters.h"

using namespace std;
using namespace xercesc;

ParamStruct ps;

int main(int argc, char **argv) {

	//Parameters file location
	char* parameters_file = (char*)DEFAULT_PARAM_FILE;
	char* output_file = 0;

	//Read command line options
	for (int i = 1; i < argc; i += 2) {
		if (!strcmp(argv[i],"-param")) {
			parameters_file = argv[i+1];
		} else if (!strcmp(argv[i],"-out")) {
			output_file = argv[i+1];
		} else {
			fprintf(stderr,"Error: unrecognized command line option: %s\n",argv[i]);
			exit(1);
		}
	}

	//Read in parameters
	ReadParameters(&ps, parameters_file);
	ps.JSR_Z = ps.Integers[I_N_JSR_SEGMENTS_Z]*ps.Reals[I_DYAD_RES];
	double tt_side = 27e-3;
	ps.TT_RADIUS = tt_side / (2*sin(PI/ps.Integers[I_N_TT_SEGMENTS])); //Circumradius of n-sided regular polygon
	if (ps.Integers[I_FLAG_TWO_CRUS]) {
		ps.CRU_1_Z = -ps.Reals[I_CRU_SEPARATION]/2;
		ps.CRU_2_Z = ps.Reals[I_CRU_SEPARATION]/2;
	} else {
		ps.CRU_1_Z = 0;
		ps.CRU_2_Z = 0;
	}
	
	if (output_file != 0) {
		sprintf(ps.Chars[I_OUTPUT_NAME],"%s",output_file);
		fprintf(stdout,"Mesh will be written to %s\n",ps.Chars[I_OUTPUT_NAME]);
	}

	fprintf(stdout,"Creating geometry...\n");

	vector<double> BoundaryX; //Vertex coordinates for outer layer of the boundary layer
	vector<double> BoundaryY;
	vector<double> BoundaryZ;

	GeometryStruct gs;

	fprintf(stdout,"Initializing cytosolic mesh...\n");
	InitializeVertices(&gs);

	//Check lists
	for (int i = 0; i < gs.VerticesX.size(); i++) {
		double x = gs.VerticesX[i];
	}
	for (int i = 0; i < gs.VerticesX.size(); i++) {
		double x = gs.VerticesY[i];
	}
	for (int i = 0; i < gs.VerticesX.size(); i++) {
		double x = gs.VerticesZ[i];
	}
	for (int i = 0; i < gs.HolesX.size(); i++) {
		double x = gs.HolesX[i];
	}
	for (int i = 0; i < gs.HolesX.size(); i++) {
		double x = gs.HolesY[i];
	}
	for (int i = 0; i < gs.HolesX.size(); i++) {
		double x = gs.HolesZ[i];
	}
	for (int i = 0; i < gs.Facets.size(); i++) {
		for (int j = 0; j < gs.Facets[i].size(); j++) {
			int x = gs.Facets[i][j];
		}
	}
	for (int i = 0; i < gs.Polys.size(); i++) {
		for (int j = 0; j < gs.Polys[i].size(); j++) {
			int x = gs.Polys[i][j];
		}
	}
	for (int i = 0; i < gs.PolyLabels.size(); i++) {
		int x = gs.PolyLabels[i];
	}
	for (int i = 0; i < gs.PolyCount.size(); i++) {
		int x = gs.PolyCount[i];
	}
	for (int i = 0; i < gs.PolyHolesCount.size(); i++) {
		int x = gs.PolyHolesCount[i];
	}
	for (int i = 0; i < gs.PolyHolesX.size(); i++) {
		double x = gs.PolyHolesX[i];
	}
	for (int i = 0; i < gs.PolyHolesX.size(); i++) {
		double x = gs.PolyHolesY[i];
	}
	for (int i = 0; i < gs.PolyHolesX.size(); i++) {
		double x = gs.PolyHolesZ[i];
	}
	for (int i = 0; i < gs.RegionsX.size(); i++) {
		double x = gs.RegionsX[i];
	}
	for (int i = 0; i < gs.RegionsX.size(); i++) {
		double x = gs.RegionsY[i];
	}
	for (int i = 0; i < gs.RegionsX.size(); i++) {
		double x = gs.RegionsZ[i];
	}
	for (int i = 0; i < gs.RegionAttr.size(); i++) {
		int x = gs.RegionAttr[i];
	}

	fprintf(stdout,"Creating tetgen object...\n");
	tetgenio in;
	tetgenio out;
	in.firstnumber = 1;

	//Nodes
	in.numberofpoints = gs.VerticesX.size();
	in.pointlist = new REAL[in.numberofpoints * 3];
	for (int i = 0; i < gs.VerticesX.size(); i++) {
		int i3 = i*3;
		in.pointlist[i3] = gs.VerticesX[i];
		in.pointlist[i3+1] = gs.VerticesY[i];
		in.pointlist[i3+2] = gs.VerticesZ[i];
	}

	//Holes
	in.numberofholes = gs.HolesX.size();
	if (in.numberofholes > 0) {
		in.holelist = new double[in.numberofholes*3];
		for (int i = 0; i < in.numberofholes; i++) {
			in.holelist[i*3] = gs.HolesX[i];
			in.holelist[i*3+1] = gs.HolesY[i];
			in.holelist[i*3+2] = gs.HolesZ[i];
		}
	} else {
		in.holelist = NULL;
	}

	//Regions
	in.numberofregions = gs.RegionAttr.size();
	if (in.numberofregions > 0) {
		in.regionlist = new double[in.numberofregions*5];
		for (int i = 0; i < in.numberofregions; i++) {
			in.regionlist[i*5] = gs.RegionsX[i];
			in.regionlist[i*5+1] = gs.RegionsY[i];
			in.regionlist[i*5+2] = gs.RegionsZ[i];
			in.regionlist[i*5+3] = gs.RegionAttr[i];
			in.regionlist[i*5+4] = 0;
		}
	} else {
		in.regionlist = NULL;
	}

	//Surfaces
	in.numberoffacets = gs.Facets.size();

	/*
	if (!TEST_FLAG) {
		in.numberoffacets -= 2;
	}
	*/

	in.facetlist = new tetgenio::facet[in.numberoffacets];
	in.facetmarkerlist = new int[in.numberoffacets];
	tetgenio::facet *f;
	tetgenio::polygon *p;
	int holeidx = 0;
	for (int i = 0; i < gs.Facets.size(); i++) {
		f = &in.facetlist[i];
		f->numberofpolygons = gs.Facets[i].size();
		f->polygonlist = new tetgenio::polygon[(f->numberofpolygons)];
		f->numberofholes = gs.PolyHolesCount[i];
		if (f->numberofholes > 0) {
			f->numberofholes = gs.PolyHolesCount[i];
			f->holelist = new double[(f->numberofholes)*3];
			for (int j = 0; j < f->numberofholes; j++) {
				f->holelist[j*3] = gs.PolyHolesX[holeidx];
				f->holelist[j*3+1] = gs.PolyHolesY[holeidx];
				f->holelist[j*3+2] = gs.PolyHolesZ[holeidx];
				holeidx++;
			}
		} else {
			f->holelist = NULL;
		}

		in.facetmarkerlist[i] = gs.PolyLabels[i];
		for (int k = 0; k < gs.Facets[i].size(); k++) {
			int polyidx = gs.Facets[i][k];
			p = &f->polygonlist[k];
			p->numberofvertices = gs.Polys[polyidx].size();
			p->vertexlist = new int[p->numberofvertices];
			for (int j = 0; j < gs.Polys[polyidx].size(); j++) {
				p->vertexlist[j] = gs.Polys[polyidx][j];
			}
		}
		int z = 0;
	}

	fprintf(stdout,"Invoking tetgen tetrahedralizer...\n");
	tetgenbehavior b;
	char switches[255]; //a0.001
	sprintf(switches,"nnfApq%fa%f",ps.Reals[I_QMAX],ps.Reals[I_VMAX]);
	//char switches[255] = "nnfp";
	b.parse_commandline(switches);
	tetrahedralize(&b,&in,&out);

	//Calculate tetrahedra volumes, centers, and troponin labels
	fprintf(stdout,"Calculating tetrahedra properties...\n");
	MeshStruct ms;
	ms.CellCenters = vector<Vector3D> (out.numberoftetrahedra);
	ms.V0 = vector<double>(out.numberoftetrahedra);
	ms.bTroponinC = vector<int>(out.numberoftetrahedra);
	ms.TTSurfaceArea = vector<double>(out.numberoftetrahedra);
	ms.SRSurfaceArea = vector<double>(out.numberoftetrahedra);
	ms.CellAttrList = vector<int> (out.numberoftetrahedra);
	for (int i = 0; i < out.numberoftetrahedra; i++) {

		//Atribute list
		if (out.tetrahedronattributelist != NULL) {
			ms.CellAttrList[i] = (int)out.tetrahedronattributelist[i];
		} else {
			ms.CellAttrList[i] = REGION_CYTO;
		}

		int j1 = out.tetrahedronlist[i*4] - out.firstnumber;
		int j2 = out.tetrahedronlist[i*4+1] - out.firstnumber;
		int j3 = out.tetrahedronlist[i*4+2] - out.firstnumber;
		int j4 = out.tetrahedronlist[i*4+3] - out.firstnumber;

		Vector3D v1, v2, v3, v4;

		v1.x = out.pointlist[j1*3];
		v1.y = out.pointlist[j1*3+1];
		v1.z = out.pointlist[j1*3+2];

		v2.x = out.pointlist[j2*3];
		v2.y = out.pointlist[j2*3+1];
		v2.z = out.pointlist[j2*3+2];

		v3.x = out.pointlist[j3*3];
		v3.y = out.pointlist[j3*3+1];
		v3.z = out.pointlist[j3*3+2];

		v4.x = out.pointlist[j4*3];
		v4.y = out.pointlist[j4*3+1];
		v4.z = out.pointlist[j4*3+2];

		//Simplex center
		Vector3D w;
		w.x = (v1.x+v2.x+v3.x+v4.x)/4;
		w.y = (v1.y+v2.y+v3.y+v4.y)/4;
		w.z = (v1.z+v2.z+v3.z+v4.z)/4;
		ms.CellCenters[i] = w;

		/*
		//TroponinC located everywhere in cytosolic region
		if (ms.CellAttrList[i] == REGION_CYTO) {
			ms.bTroponinC[i] = 1;
		} else {
			ms.bTroponinC[i] = 0;
		}
		*/
		//TODO
		//TroponinC located beyond a given radius
		if (sqrt(w.x*w.x + w.y*w.y) > ps.Reals[I_R_TRPN] && ms.CellAttrList[i] != REGION_JSR) {
			ms.bTroponinC[i] = 1;
		} else {
			ms.bTroponinC[i] = 0;
		}

		//Volume calculation
		v1.x -= v4.x;
		v1.y -= v4.y;
		v1.z -= v4.z;

		v2.x -= v4.x;
		v2.y -= v4.y;
		v2.z -= v4.z;

		v3.x -= v4.x;
		v3.y -= v4.y;
		v3.z -= v4.z;

		double T = dot(v1,cross(v2.x,v2.y,v2.z,v3.x,v3.y,v3.z));
		ms.V0[i] = fabs(T) / 6;
	}

	//Count number of valid faces
	int N_Faces = 0;
	for (int i = 0; i < out.numberoftrifaces; i++) {
		if (out.adjtetlist[i*2] > 0 && out.adjtetlist[i*2] <= out.numberoftetrahedra) {
			N_Faces++;
		} else {
			break;
		}
	}

	//Calculate face properties
	fprintf(stdout,"Calculating face properties...\n");
	ms.FaceTet = vector<IndexPair>(N_Faces);
	ms.FaceNormals = vector<Vector3D>(N_Faces);
	ms.FaceCenters = vector<Vector3D>(N_Faces);
	ms.Af = vector<Vector3D>(N_Faces);
	ms.FaceArea = vector<double>(N_Faces);
	ms.r1 = vector<Vector3D>(N_Faces);
	ms.r2 = vector<Vector3D>(N_Faces);
	ms.eXi = vector<Vector3D>(N_Faces);
	ms.dXi = vector<double>(N_Faces);
	ms.FaceCoeff = vector<double>(N_Faces);
	ms.FaceAttrList = vector<int>(N_Faces);

	for (int i = 0; i < N_Faces; i++) {

		//Face attributes
		ms.FaceAttrList[i] = out.trifacemarkerlist[i];

		//Face neighbors
		IndexPair p;
		p.i = out.adjtetlist[i*2]-out.firstnumber;
		p.j = out.adjtetlist[i*2+1]-out.firstnumber;
		ms.FaceTet[i] = p;

		Vector3D v;

		int v1 = out.trifacelist[i*3]-out.firstnumber;
		int v2 = out.trifacelist[i*3+1]-out.firstnumber;
		int v3 = out.trifacelist[i*3+2]-out.firstnumber;

		double x1 = out.pointlist[v1*3];
		double y1 = out.pointlist[v1*3+1];
		double z1 = out.pointlist[v1*3+2];

		double x2 = out.pointlist[v2*3];
		double y2 = out.pointlist[v2*3+1];
		double z2 = out.pointlist[v2*3+2];

		double x3 = out.pointlist[v3*3];
		double y3 = out.pointlist[v3*3+1];
		double z3 = out.pointlist[v3*3+2];

		//Face center
		double cx = (x1+x2+x3)/3;
		double cy = (y1+y2+y3)/3;
		double cz = (z1+z2+z3)/3;

		v.x = cx;
		v.y = cy;
		v.z = cz;

		ms.FaceCenters[i] = v;

		//Calculate a normal vector
		double xw = x2-x1;
		double yw = y2-y1;
		double zw = z2-z1;

		double xu = x3-x1;
		double yu = y3-y1;
		double zu = z3-z1;

		Vector3D n = cross(xw,yw,zw,xu,yu,zu);

		//Calculate face area
		ms.FaceArea[i] = norm(n)/2; //Cross product magnitude is twice the area of the triangle

		//Calculate area vector
		v.x = n.x/2; //Cross product magnitude is twice the area of the triangle
		v.y = n.y/2;
		v.z = n.z/2;

		//Check if it is Af points in correct direction (away from face neighbor i)
		int ni = p.i;
		double xs = ms.CellCenters[ni].x;
		double ys = ms.CellCenters[ni].y;
		double zs = ms.CellCenters[ni].z;

		Vector3D d1, d2;

		d1.x = cx + v.x - xs;
		d1.y = cy + v.y - ys;
		d1.z = cz + v.z - zs;

		d2.x = cx - v.x - xs;
		d2.y = cy - v.y - ys;
		d2.z = cz - v.z - zs;

		if (norm(d1) < norm(d2)) {
			v.x = -v.x;
			v.y = -v.y;
			v.z = -v.z;
		}

		ms.Af[i] = v;

		//Unit normal
		normalize(v);
		ms.FaceNormals[i] = v;

		//r1 vector from simplex i's center to face center
		xs = ms.CellCenters[ni].x;
		ys = ms.CellCenters[ni].y;
		zs = ms.CellCenters[ni].z;
		v.x = cx - xs;
		v.y = cy - ys;
		v.z = cz - zs;
		ms.r1[i] = v;

		//r2 vector from simplex j's center to face center
		int nj = p.j;
		double xt, yt, zt;
		if (nj >= 0) {
			xt = ms.CellCenters[nj].x;
			yt = ms.CellCenters[nj].y;
			zt = ms.CellCenters[nj].z;
			v.x = cx - xt;
			v.y = cy - yt;
			v.z = cz - zt;
		} /*else if (ms.FaceAttrList[i] == ATTR_DYAD) {
			//If a dyad face, r2 points from center of dyad to face center
			if (!FLAG_TEST) {
				xt = (TT_RADIUS + (DYAD_R/2)) * cos(DYAD_THETA);
				yt = (TT_RADIUS + (DYAD_R/2)) * sin(DYAD_THETA);
				zt = 0;
			} else if (FLAG_TEST_POINT_SOURCE) {
				xt = 0;
				yt = 0;
				zt = 0;
			} else {
				xt = cx;
				yt = cy;
				zt = cz;
			}
			v.x = cx - xt;
			v.y = cy - yt;
			v.z = cz - zt;
		}*/ else {
			v.x = 0;
			v.y = 0;
			v.z = 0;
		}
		ms.r2[i] = v;

		//Xi vector and magnitude between simplex centers
		if (nj >= 0) {// || ms.FaceAttrList[i] == ATTR_DYAD) {
			v.x = xt - xs;
			v.y = yt - ys;
			v.z = zt - zs;
		} else {
			v.x = cx - xs;
			v.y = cy - ys;
			v.z = cz - zs;
		}
		ms.dXi[i] = norm(v);
		normalize(v);
		ms.eXi[i] = v;

		//Numerical coefficient for the face
		ms.FaceCoeff[i] = dot(ms.Af[i],ms.Af[i]) / (ms.dXi[i] * dot(ms.Af[i],ms.eXi[i]));

	}

	//Add TT and SR surface area adjacent to each cell
	for (int i = 0; i < ms.TTSurfaceArea.size(); i++) {
		ms.TTSurfaceArea[i] = 0;
		ms.SRSurfaceArea[i] = 0;
	}

	for (int i = 0; i < ms.FaceAttrList.size(); i++) {
		if (ms.FaceAttrList[i] == ATTR_TT && ps.Integers[I_FLAG_TT]) {
			ms.TTSurfaceArea[ms.FaceTet[i].i] += ms.FaceArea[i];
		} else if (ms.FaceAttrList[i] == ATTR_JSR) {
			ms.SRSurfaceArea[ms.FaceTet[i].i] += ms.FaceArea[i];
		}
	}

	//Save basic information to standard tetgen files
	out.save_nodes((char*)ps.Chars[I_OUTPUT_NAME]);
	out.save_faces((char*)ps.Chars[I_OUTPUT_NAME]);
	out.save_neighbors((char*)ps.Chars[I_OUTPUT_NAME]);
	out.save_elements((char*)ps.Chars[I_OUTPUT_NAME]);

	WriteChannelFiles(ms.CellCenters,ms.CellAttrList);
	WritePropertiesFile(&ms,ps.Chars[I_OUTPUT_NAME]);

	return 0;
}


void WriteChannelFiles(vector<Vector3D> &CellCenters, vector<int> &CellAttrList) {

	int N_LCC = ps.N_LCC;

	vector<int> RyR_Ele;
	vector<int> RyR_JSR_Ele;
	vector<int> LCC_Ele;
	vector<int> RYR_i;

	ChannelHelper(CellCenters,CellAttrList,
			   RyR_Ele,RyR_JSR_Ele,LCC_Ele,
			   RYR_i, ps.CRU_1_Z);
	if (ps.Integers[I_FLAG_TWO_CRUS]) {
		ChannelHelper(CellCenters,CellAttrList,
				   RyR_Ele,RyR_JSR_Ele,LCC_Ele,
				   RYR_i, ps.CRU_2_Z);
	}

	//Write RyR to file
	FILE * pFile;
	char outfile[255];
	sprintf(outfile,"%s_ryr.txt",ps.Chars[I_OUTPUT_NAME]);
	fprintf(stdout,"Writing RyR file %s\n",outfile);
	pFile = fopen (outfile,"w");
	if (pFile != NULL) {

		fprintf(pFile,"%d\n",RyR_Ele.size()); //header line is # of ryr's
		for (int i = 0; i < RyR_Ele.size(); i++) {
			fprintf(pFile,"%d %d %d %d %d %d\n",RyR_Ele[i],RyR_JSR_Ele[i],RYR_i[4*i+0],RYR_i[4*i+1],RYR_i[4*i+2],RYR_i[4*i+3]);
		}

		fclose(pFile);

	} else {
		fprintf(stderr,"Warning: could not open RyR output file.\n");
		exit(1);
	}

	//Write LCC to file
	sprintf(outfile,"%s_lcc.txt",ps.Chars[I_OUTPUT_NAME]);
	fprintf(stdout,"Writing LCC file %s\n",outfile);
	pFile = fopen (outfile,"w");
	if (pFile != NULL) {

		fprintf(pFile,"%d\n",LCC_Ele.size()); //header line is # of lcc's
		for (int i = 0; i < LCC_Ele.size(); i++) {
			fprintf(pFile,"%d\n",LCC_Ele[i]);
		}

		fclose(pFile);

	} else {
		fprintf(stderr,"Warning: could not open LCC output file.\n");
		exit(1);
	}

}

void ChannelHelper(vector<Vector3D> &CellCenters, vector<int> &CellAttrList,
				   vector<int> &RyR_Ele,vector<int> &RyR_JSR_Ele,vector<int> &LCC_Ele,
				   vector<int> &RYR_i, double z_center) {

	fprintf(stdout,"Adding channels z = %f...\n",z_center);

	int N_LCC = ps.N_LCC;

	vector<double> LCC_X(N_LCC);
	vector<double> LCC_Y(N_LCC);
	vector<double> LCC_Z(N_LCC);

	vector<double> RYR_X(ps.N_RyR);
	vector<double> RYR_Y(ps.N_RyR);
	vector<double> RYR_Z(ps.N_RyR);

	vector<double> RYR_JSR_X(ps.N_RyR);
	vector<double> RYR_JSR_Y(ps.N_RyR);
	vector<double> RYR_JSR_Z(ps.N_RyR);

	double x0,y0,z0,angle;
	double ryr_y = ps.Reals[I_RYR_SPACING]*ps.N_RyR_Y;
	if (ps.N_RyR_Y % 2 == 1) {
	    ryr_y = ps.Reals[I_RYR_SPACING]*ps.N_RyR_Y;
	} else {
	    ryr_y = ps.Reals[I_RYR_SPACING]*(ps.N_RyR_Y+1);
	}
	double lcc_y;
	if (ps.N_LCC_Y % 2 == 1) {
	    lcc_y = ps.Reals[I_RYR_SPACING]*ps.N_LCC_Y;
	} else {
	    lcc_y = ps.Reals[I_RYR_SPACING]*(ps.N_LCC_Y+1);
	}
	int ryr_offset = RyR_Ele.size();

	//double dtheta_ryr = (ps.Reals[I_RYR_SPACING]/ps.Reals[I_DYAD_RES])*2*PI/ps.Integers[I_N_TT_SEGMENTS];
	double dtheta_ryr = ps.Reals[I_DYAD_RES]/(ps.Reals[I_DYAD_R]+ps.TT_RADIUS);
	double dtheta_lcc = 2*PI/ps.Integers[I_N_TT_SEGMENTS];
	double dyad_angle = (ps.N_RyR_X)*dtheta_ryr;

	//LCC placement
	int n_segments = ps.Integers[I_N_JSR_SEGMENTS];
	double angle0;
	if (ps.N_LCC_X % 2 == 1) {
	    angle0 = (0.5 - (double)((ps.N_LCC_X-1)/2))*dtheta_lcc;
	} else {
	    angle0 = (0.5 -(double)((ps.N_LCC_X)/2))*dtheta_lcc;
	}
	double r0 = ps.TT_RADIUS;
	int idx = 0;
	for (int i = 0; i < ps.N_LCC_Y; i++) {
		for (int j = 0; j < ps.N_LCC_X; j++) {
			if (ps.bLCC[i][j]) {
				angle = angle0 + (j*dtheta_lcc);
				x0 = r0*cos(angle);
				y0 = r0*sin(angle);
				z0 = z_center+(lcc_y/2) - i*ps.Reals[I_RYR_SPACING] - (ps.Reals[I_RYR_SPACING]/2);
				LCC_X[idx] = x0;
				LCC_Y[idx] = y0;
				LCC_Z[idx] = z0;
				idx++;
			}
		}
	}

	//RyR placement
	n_segments = ps.Integers[I_N_JSR_SEGMENTS];
	angle0 = (0.5 - (double)((ps.N_RyR_X-1)/2))*dtheta_ryr;
	if (ps.N_RyR_X % 2 == 1) {
	    angle0 = (0.5 - (double)((ps.N_RyR_X-1)/2))*dtheta_ryr;
	} else {
	    angle0 = (0.5 - (double)((ps.N_RyR_X)/2))*dtheta_ryr;
	}
	r0 = ps.TT_RADIUS + ps.Reals[I_DYAD_R] - (10e-3);
	idx = 0;
	for (int i = 0; i < ps.N_RyR_Y; i++) {
		for (int j = 0; j < ps.N_RyR_X; j++) {
			if (ps.bRyR[i][j]) {
				angle = angle0 + (j*dtheta_ryr);
				x0 = r0*cos(angle) + ps.Reals[I_JSR_TRANS];
				y0 = r0*sin(angle);
				z0 = z_center+(ryr_y/2) - i*ps.Reals[I_RYR_SPACING] - (ps.Reals[I_RYR_SPACING]/2);
				RYR_X[idx] = x0;
				RYR_Y[idx] = y0;
				RYR_Z[idx] = z0;
				idx++;
			}
		}
	}

	//RyR placement (JSR)
	r0 = ps.TT_RADIUS + ps.Reals[I_DYAD_R];
	idx = 0;
	for (int i = 0; i < ps.N_RyR_Y; i++) {
		for (int j = 0; j < ps.N_RyR_X; j++) {
			if (ps.bRyR[i][j]) {
				angle = angle0 + (j*dtheta_ryr);
				x0 = r0*cos(angle) + ps.Reals[I_JSR_TRANS];
				y0 = r0*sin(angle);
				z0 = z_center+(ryr_y/2) - i*ps.Reals[I_RYR_SPACING] - (ps.Reals[I_RYR_SPACING]/2);
				RYR_JSR_X[idx] = x0;
				RYR_JSR_Y[idx] = y0;
				RYR_JSR_Z[idx] = z0;
				idx++;
			}
		}
	}

	//Determine which element centers each channel is closest to

	double r2,r2min;
	double dx,dy,dz;
	int emin;
	for (int i = 0; i < N_LCC; i++) {
		r2min = ps.Reals[I_GRID_X]*ps.Reals[I_GRID_X] + ps.Reals[I_GRID_Y]*ps.Reals[I_GRID_Y] + ps.Reals[I_GRID_Z]*ps.Reals[I_GRID_Z];
		emin = -1;
		for (int j = 0; j < CellCenters.size(); j++) {
			dx = CellCenters[j].x-LCC_X[i];
			dy = CellCenters[j].y-LCC_Y[i];
			dz = CellCenters[j].z-LCC_Z[i];
			r2 = dx*dx + dy*dy + dz*dz;
			if (r2 < r2min) {
				r2min = r2;
				emin = j;
			}
		}
		LCC_Ele.push_back(emin);
	}

	for (int i = 0; i < ps.N_RyR; i++) {
		r2min = ps.Reals[I_GRID_X]*ps.Reals[I_GRID_X] + ps.Reals[I_GRID_Y]*ps.Reals[I_GRID_Y] + ps.Reals[I_GRID_Z]*ps.Reals[I_GRID_Z];
		emin = -1;
		for (int j = 0; j < CellCenters.size(); j++) {
			if (CellAttrList[j] == REGION_CYTO) {
				dx = CellCenters[j].x-RYR_X[i];
				dy = CellCenters[j].y-RYR_Y[i];
				dz = CellCenters[j].z-RYR_Z[i];
				r2 = dx*dx + dy*dy + dz*dz;
				if (r2 < r2min) {
					r2min = r2;
					emin = j;
				}
			}
		}
		RyR_Ele.push_back(emin);
	}

	//JSR elements
	for (int i = 0; i < ps.N_RyR; i++) {
		r2min = ps.Reals[I_GRID_X]*ps.Reals[I_GRID_X] + ps.Reals[I_GRID_Y]*ps.Reals[I_GRID_Y] + ps.Reals[I_GRID_Z]*ps.Reals[I_GRID_Z];
		emin = -1;
		for (int j = 0; j < CellCenters.size(); j++) {
			if (CellAttrList[j] == REGION_JSR) {
				dx = CellCenters[j].x-RYR_JSR_X[i];
				dy = CellCenters[j].y-RYR_JSR_Y[i];
				dz = CellCenters[j].z-RYR_JSR_Z[i];
				r2 = dx*dx + dy*dy + dz*dz;
				if (r2 < r2min) {
					r2min = r2;
					emin = j;
				}
			}
		}
		RyR_JSR_Ele.push_back(emin);
	}

	//For RyR connectivity, need to re-label the arrays according to RyR index (+1)
	idx = 1;
	for (int i = 0; i < ps.N_RyR_Y; i++) {
		for (int j = 0; j < ps.N_RyR_X; j++) {
			if (ps.bRyR[i][j]) {
				ps.bRyR[i][j] = idx;
				idx++;
			}
		}
	}

	//Compute RyR connectivity
	for (int i = 0; i < ps.N_RyR_Y; i++) {
		for (int j = 0; j < ps.N_RyR_X; j++) {
			if (ps.bRyR[i][j]) {
				if (i > 0 && ps.bRyR[i-1][j])
					RYR_i.push_back(ryr_offset+ps.bRyR[i-1][j]-1);
				else
					RYR_i.push_back(-1);

				if (i < ps.N_RyR_Y-1 && ps.bRyR[i+1][j])
					RYR_i.push_back(ryr_offset+ps.bRyR[i+1][j]-1);
				else
					RYR_i.push_back(-1);

				if (j > 0 && ps.bRyR[i][j-1])
					RYR_i.push_back(ryr_offset+ps.bRyR[i][j-1]-1);
				else
					RYR_i.push_back(-1);

				if (j < ps.N_RyR_X-1 && ps.bRyR[i][j+1])
					RYR_i.push_back(ryr_offset+ps.bRyR[i][j+1]-1);
				else
					RYR_i.push_back(-1);
			}
		}
	}

}

void InitializeVertices(GeometryStruct* gs) {

	double x0,y0,z0;
	vector <int> v;

	//Outer boundaries
	int f0_bound = gs->Facets.size();
	AddBox(gs, ATTR_BOUNDARY, 0, 0, 0, ps.Reals[I_GRID_X], ps.Reals[I_GRID_Y], ps.Reals[I_GRID_Z]);

	//JSR
	AddJSR(gs,ps.CRU_1_Z);
	if (ps.Integers[I_FLAG_DYAD_GUIDE]) {
		AddDyadGuide(gs,ps.CRU_1_Z); //Adds points to help keep RyR elements similar across varying dyad height values
		//AddSegmentedTT(gs, ATTR_NONE, ps.TT_RADIUS + ps.Reals[I_DYAD_R] - (15e-3), 0);
	}

	if (ps.Integers[I_FLAG_TWO_CRUS]) {
		AddJSR(gs,ps.CRU_2_Z);
		if (ps.Integers[I_FLAG_DYAD_GUIDE]) {
			AddDyadGuide(gs,ps.CRU_2_Z); //Adds points to help keep RyR elements similar across varying dyad height values
			//AddSegmentedTT(gs, ATTR_NONE, ps.TT_RADIUS + ps.Reals[I_DYAD_R] - (15e-3), 0);
		}
	}

	//TT
	if (ps.Integers[I_FLAG_TT]) {
		AddSegmentedTT(gs, ATTR_TT, ps.TT_RADIUS, 1);
	}

	//Volume holes

	if (ps.Integers[I_FLAG_TT]) {
		//TT
		gs->HolesX.push_back(0);
		gs->HolesY.push_back(0);
		gs->HolesZ.push_back(0);
	}

	//Regions

	//Cytosol
	x0 = ps.TT_RADIUS+ps.Reals[I_DYAD_R]/2;
	y0 = 0;
	z0 = 0;
	gs->RegionsX.push_back(x0);
	gs->RegionsY.push_back(y0);
	gs->RegionsZ.push_back(z0);
	gs->RegionAttr.push_back(REGION_CYTO);

	//JSR
	double jsr_theta_center = 0;
	x0 = (ps.TT_RADIUS+ps.Reals[I_DYAD_R]+(ps.Reals[I_JSR_R]/2))*cos(jsr_theta_center) + ps.Reals[I_JSR_TRANS];
	y0 = (ps.TT_RADIUS+ps.Reals[I_DYAD_R]+(ps.Reals[I_JSR_R]/2))*sin(jsr_theta_center);
	z0 = ps.CRU_1_Z;
	gs->RegionsX.push_back(x0);
	gs->RegionsY.push_back(y0);
	gs->RegionsZ.push_back(z0);
	gs->RegionAttr.push_back(REGION_JSR);

	if (ps.Integers[I_FLAG_TWO_CRUS]) {
		z0 = ps.CRU_2_Z;
		gs->RegionsX.push_back(x0);
		gs->RegionsY.push_back(y0);
		gs->RegionsZ.push_back(z0);
		gs->RegionAttr.push_back(REGION_JSR);
	}

}

void AddJSR(GeometryStruct* gs, double z_center) {

	vector<int> stitch_in;
	vector<int> stitch_out;
	vector<int> stitch_next_in;
	vector<int> stitch_next_out;

	double dz = ps.Reals[I_DYAD_RES];
	
	double z0;
	if (ps.Integers[I_N_JSR_SEGMENTS_Z] % 2 == 1) {
	    z0 = -(((double)ps.Integers[I_N_JSR_SEGMENTS_Z])*ps.Reals[I_DYAD_RES])/2.0;
	} else {
	    z0 = -(((double)ps.Integers[I_N_JSR_SEGMENTS_Z]+1)*ps.Reals[I_DYAD_RES])/2.0;
	}
	for (int i = 0; i <= ps.Integers[I_N_JSR_SEGMENTS_Z]; i++) {

		int bCapTop = (i==0||i==ps.Integers[I_N_JSR_SEGMENTS_Z]) ? 1 : 0;

		JSRHelper(gs,stitch_in,stitch_out,stitch_next_in,stitch_next_out,z_center+z0+i*dz,bCapTop);
		stitch_in = stitch_next_in;
		stitch_out = stitch_next_out;
	}

}

void JSRHelper(GeometryStruct* gs, vector<int> v_stitch_in, vector<int> v_stitch_out,
		  	   vector<int> &stitch_next_in, vector<int> &stitch_next_out, double z0, int bCapTop) {

	vector<int> v;
	vector<int> f;
	stitch_next_in.clear();
	stitch_next_out.clear();
	int v0 = gs->VerticesX.size() + 1;
	double x0, y0;
	int bStitch = (v_stitch_in.size()>0)&&(v_stitch_out.size()>0);

	//double d_theta = 2*PI*ps.TT_RADIUS / ((ps.Reals[I_DYAD_R]+ps.TT_RADIUS)*(double)ps.Integers[I_N_TT_SEGMENTS]);
	//double d_theta = 2*PI / (ps.Integers[I_N_TT_SEGMENTS]);
	double d_theta = ps.Reals[I_DYAD_RES]/(ps.Reals[I_DYAD_R]+ps.TT_RADIUS);
	//fprintf(stdout,"d_theta = %f, d_theta0 = %f\n",d_theta,d_theta0);

	int n_segments = ps.Integers[I_N_JSR_SEGMENTS];
	int jsr_start;
	if (ps.Integers[I_N_JSR_SEGMENTS] % 2 == 1) {
	    jsr_start = -(n_segments-1)/2;
	} else {
	    jsr_start = -(n_segments)/2;
	}

	//Inner side
	int idx = 0;
	for (int i = jsr_start; i <= jsr_start+n_segments; i++) {
		double angle = d_theta*i;
		x0 = (ps.TT_RADIUS+ps.Reals[I_DYAD_R])*cos(angle) + ps.Reals[I_JSR_TRANS];
		y0 = (ps.TT_RADIUS+ps.Reals[I_DYAD_R])*sin(angle);

		gs->VerticesX.push_back(x0);
		gs->VerticesY.push_back(y0);
		gs->VerticesZ.push_back(z0);

		stitch_next_in.push_back(v0 + idx);

		if (bStitch && idx > 0) {
			v.clear();
			v.push_back(v_stitch_in[idx-1]);
			v.push_back(v0 + idx - 1);
			v.push_back(v0 + idx);
			v.push_back(v_stitch_in[idx]);
			gs->Polys.push_back(v);
			gs->PolyLabels.push_back(ATTR_JSR);
			gs->PolyCount.push_back(1);
			gs->PolyHolesCount.push_back(0);
			f.clear();
			f.push_back(gs->Polys.size()-1);
			gs->Facets.push_back(f);
		}

		idx++;
	}


	//Outer side
	v0 = gs->VerticesX.size() + 1;
	idx = 0;
	for (int i = jsr_start; i <= jsr_start+n_segments; i++) {

		double angle = d_theta*i;
		x0 = (ps.TT_RADIUS+ps.Reals[I_DYAD_R]+ps.Reals[I_JSR_R])*cos(angle) + ps.Reals[I_JSR_TRANS];
		y0 = (ps.TT_RADIUS+ps.Reals[I_DYAD_R]+ps.Reals[I_JSR_R])*sin(angle);

		gs->VerticesX.push_back(x0);
		gs->VerticesY.push_back(y0);
		gs->VerticesZ.push_back(z0);

		stitch_next_out.push_back(v0 + idx);

		if (bStitch && idx > 0) {
			v.clear();
			v.push_back(v_stitch_out[idx-1]);
			v.push_back(v0 + idx - 1);
			v.push_back(v0 + idx);
			v.push_back(v_stitch_out[idx]);
			gs->Polys.push_back(v);
			gs->PolyLabels.push_back(ATTR_JSR);
			gs->PolyCount.push_back(1);
			gs->PolyHolesCount.push_back(0);
			f.clear();
			f.push_back(gs->Polys.size()-1);
			gs->Facets.push_back(f);
		}
		idx++;

	}
	if (bCapTop) {
		//Join inner/outer JSR
		for (int i = 1; i < stitch_next_out.size(); i++) {
			v.clear();
			v.push_back(stitch_next_in[i]);
			v.push_back(stitch_next_out[i]);
			v.push_back(stitch_next_out[i-1]);
			v.push_back(stitch_next_in[i-1]);
			gs->Polys.push_back(v);
			gs->PolyLabels.push_back(ATTR_JSR);
			gs->PolyCount.push_back(1);
			gs->PolyHolesCount.push_back(0);
			f.clear();
			f.push_back(gs->Polys.size()-1);
			gs->Facets.push_back(f);
		}
	}
	//Cap sides
	if (bStitch) {
		v.clear();
		v.push_back(v_stitch_in[0]);
		v.push_back(v_stitch_out[0]);
		v.push_back(stitch_next_out[0]);
		v.push_back(stitch_next_in[0]);
		gs->Polys.push_back(v);
		gs->PolyLabels.push_back(ATTR_JSR);
		gs->PolyCount.push_back(1);
		gs->PolyHolesCount.push_back(0);
		f.clear();
		f.push_back(gs->Polys.size()-1);
		gs->Facets.push_back(f);

		v.clear();
		v.push_back(v_stitch_in[v_stitch_in.size()-1]);
		v.push_back(v_stitch_out[v_stitch_out.size()-1]);
		v.push_back(stitch_next_out[stitch_next_out.size()-1]);
		v.push_back(stitch_next_in[stitch_next_in.size()-1]);
		gs->Polys.push_back(v);
		gs->PolyLabels.push_back(ATTR_JSR);
		gs->PolyCount.push_back(1);
		gs->PolyHolesCount.push_back(0);
		f.clear();
		f.push_back(gs->Polys.size()-1);
		gs->Facets.push_back(f);
	}
}

void AddDyadGuide(GeometryStruct* gs, double z_center) {

	double x0, y0, z0;

	double d_theta = ps.Reals[I_DYAD_RES]/(ps.Reals[I_DYAD_R]+ps.TT_RADIUS);

	int n_segments = ps.Integers[I_N_JSR_SEGMENTS];
	int jsr_start;
	if (ps.Integers[I_N_JSR_SEGMENTS] % 2 == 1) {
	    jsr_start = -(n_segments-1)/2;
	} else {
	    jsr_start = -(n_segments)/2;
	}
	double dz = ps.Reals[I_DYAD_RES];
	
	double z_base;
	if (ps.Integers[I_N_JSR_SEGMENTS_Z] % 2 == 1) {
	    z_base = -(((double)ps.Integers[I_N_JSR_SEGMENTS_Z])*ps.Reals[I_DYAD_RES])/2.0;
	} else {
	    z_base = -(((double)ps.Integers[I_N_JSR_SEGMENTS_Z]+1)*ps.Reals[I_DYAD_RES])/2.0;
	}

	//jSR side always added
	for (int k = 0; k < ps.Integers[I_N_JSR_SEGMENTS_Z]; k++) {
		z0 = z_center + z_base + k*dz;
		for (int i = jsr_start; i <= jsr_start+n_segments; i++) {
			double angle = d_theta*i;
			x0 = (ps.TT_RADIUS+ps.Reals[I_DYAD_R]-(15e-3))*cos(angle) + ps.Reals[I_JSR_TRANS];
			y0 = (ps.TT_RADIUS+ps.Reals[I_DYAD_R]-(15e-3))*sin(angle);

			gs->VerticesX.push_back(x0);
			gs->VerticesY.push_back(y0);
			gs->VerticesZ.push_back(z0);
		}
	}
	
	//If dyad is wide enough, add TT side
	if (ps.Reals[I_DYAD_R] >= 45e-3) {
	
	for (int k = 0; k < ps.Integers[I_N_JSR_SEGMENTS_Z]; k++) {
		z0 = z_center + z_base +k*dz;
		for (int i = jsr_start; i <= jsr_start+n_segments; i++) {
			double angle = d_theta*i;
			x0 = (ps.TT_RADIUS+(15e-3))*cos(angle) + ps.Reals[I_JSR_TRANS];
			y0 = (ps.TT_RADIUS+(15e-3))*sin(angle);

			gs->VerticesX.push_back(x0);
			gs->VerticesY.push_back(y0);
			gs->VerticesZ.push_back(z0);
		}
		
	}
	
	}
}

void AddFace(GeometryStruct* gs, int attribute,
			 double x1, double y1, double z1,
			 double x2, double y2, double z2,
			 double x3, double y3, double z3,
			 double x4, double y4, double z4) {

	int v0 = gs->VerticesX.size();

	gs->VerticesX.push_back(x1);
	gs->VerticesY.push_back(y1);
	gs->VerticesZ.push_back(z1);

	gs->VerticesX.push_back(x2);
	gs->VerticesY.push_back(y2);
	gs->VerticesZ.push_back(z2);

	gs->VerticesX.push_back(x3);
	gs->VerticesY.push_back(y3);
	gs->VerticesZ.push_back(z3);

	gs->VerticesX.push_back(x4);
	gs->VerticesY.push_back(y4);
	gs->VerticesZ.push_back(z4);

	vector <int> v;
	v.push_back(v0+0);
	v.push_back(v0+1);
	v.push_back(v0+2);
	v.push_back(v0+3);
	gs->Polys.push_back(v);
	gs->PolyLabels.push_back(attribute);
	gs->PolyCount.push_back(1);
	gs->PolyHolesCount.push_back(0);
	vector <int> f;
	f.push_back(gs->Polys.size()-1);
	gs->Facets.push_back(f);


}

void AddBox(GeometryStruct* gs, int attribute, double x_center, double y_center, double z_center, double dx, double dy, double dz) {

	fprintf(stdout,"Adding box with dimensions %g x %g x %g\n",dx,dy,dz);
	int v0 = gs->VerticesX.size();

	double x0 = (dx/2.0);
	double y0 = (dy/2.0);
	double z0 = (dz/2.0);

	//1
	gs->VerticesX.push_back(x_center + x0);
	gs->VerticesY.push_back(y_center + y0);
	gs->VerticesZ.push_back(z_center + z0);

	//2
	gs->VerticesX.push_back(x_center - x0);
	gs->VerticesY.push_back(y_center + y0);
	gs->VerticesZ.push_back(z_center + z0);

	//3
	gs->VerticesX.push_back(x_center + x0);
	gs->VerticesY.push_back(y_center - y0);
	gs->VerticesZ.push_back(z_center + z0);

	//4
	gs->VerticesX.push_back(x_center + x0);
	gs->VerticesY.push_back(y_center + y0);
	gs->VerticesZ.push_back(z_center - z0);

	//5
	gs->VerticesX.push_back(x_center - x0);
	gs->VerticesY.push_back(y_center - y0);
	gs->VerticesZ.push_back(z_center + z0);

	//6
	gs->VerticesX.push_back(x_center - x0);
	gs->VerticesY.push_back(y_center + y0);
	gs->VerticesZ.push_back(z_center - z0);

	//7
	gs->VerticesX.push_back(x_center + x0);
	gs->VerticesY.push_back(y_center - y0);
	gs->VerticesZ.push_back(z_center - z0);

	//8
	gs->VerticesX.push_back(x_center - x0);
	gs->VerticesY.push_back(y_center - y0);
	gs->VerticesZ.push_back(z_center - z0);

	vector <int> v;
	vector <int> f;

	//XY
	v.clear();
	v.push_back(v0+1);
	v.push_back(v0+2);
	v.push_back(v0+5);
	v.push_back(v0+3);
	gs->Polys.push_back(v);
	gs->PolyLabels.push_back(attribute);
	gs->PolyCount.push_back(1);
	gs->PolyHolesCount.push_back(0);
	f.clear();
	f.push_back(gs->Polys.size()-1);
	gs->Facets.push_back(f);

	v.clear();
	v.push_back(v0+4);
	v.push_back(v0+6);
	v.push_back(v0+8);
	v.push_back(v0+7);
	gs->Polys.push_back(v);
	gs->PolyLabels.push_back(attribute);
	gs->PolyCount.push_back(1);
	gs->PolyHolesCount.push_back(0);
	f.clear();
	f.push_back(gs->Polys.size()-1);
	gs->Facets.push_back(f);

	//XZ
	v.clear();
	v.push_back(v0+1);
	v.push_back(v0+2);
	v.push_back(v0+6);
	v.push_back(v0+4);
	gs->Polys.push_back(v);
	gs->PolyLabels.push_back(attribute);
	gs->PolyCount.push_back(1);
	gs->PolyHolesCount.push_back(0);
	f.clear();
	f.push_back(gs->Polys.size()-1);
	gs->Facets.push_back(f);

	v.clear();
	v.push_back(v0+3);
	v.push_back(v0+5);
	v.push_back(v0+8);
	v.push_back(v0+7);
	gs->Polys.push_back(v);
	gs->PolyLabels.push_back(attribute);
	gs->PolyCount.push_back(1);
	gs->PolyHolesCount.push_back(0);
	f.clear();
	f.push_back(gs->Polys.size()-1);
	gs->Facets.push_back(f);

	//YZ
	v.clear();
	v.push_back(v0+1);
	v.push_back(v0+3);
	v.push_back(v0+7);
	v.push_back(v0+4);
	gs->Polys.push_back(v);
	gs->PolyLabels.push_back(attribute);
	gs->PolyCount.push_back(1);
	gs->PolyHolesCount.push_back(0);
	f.clear();
	f.push_back(gs->Polys.size()-1);
	gs->Facets.push_back(f);

	v.clear();
	v.push_back(v0+2);
	v.push_back(v0+5);
	v.push_back(v0+8);
	v.push_back(v0+6);
	gs->Polys.push_back(v);
	gs->PolyLabels.push_back(attribute);
	gs->PolyCount.push_back(1);
	gs->PolyHolesCount.push_back(0);
	f.clear();
	f.push_back(gs->Polys.size()-1);
	gs->Facets.push_back(f);

}

void AddSegmentedTT(GeometryStruct* gs, int Poly_Attr, double radius, int bCapped) {

	//Angle between points
	double theta = 2.0*PI / (double)ps.Integers[I_N_TT_SEGMENTS];
	double theta0 = 0;
	double z0;
	if (ps.Integers[I_N_TT_SEGMENTS_Z] % 2 == 1) {
	    z0 = -ps.Integers[I_N_TT_SEGMENTS_Z]*ps.Reals[I_DYAD_RES]/2;
	} else {
	    z0 = -(ps.Integers[I_N_TT_SEGMENTS_Z]+1)*ps.Reals[I_DYAD_RES]/2;
	}
	double dz = ps.Reals[I_DYAD_RES];

	for (int i = 0; i <= ps.Integers[I_N_TT_SEGMENTS_Z]; i++) {

		double z1 = z0 + i*dz;

		vector<int> v_cap;

		for (int j = 0; j < ps.Integers[I_N_TT_SEGMENTS]; j++) {

			double angle = theta0 + j*theta;
			double x1 = radius * cos(angle);
			double y1 = radius * sin(angle);

			gs->VerticesX.push_back(x1);
			gs->VerticesY.push_back(y1);
			gs->VerticesZ.push_back(z1);

			v_cap.push_back(gs->VerticesX.size());

			if (i > 0) {

				int v0 = gs->VerticesX.size();

				int dj = (j==0) ? ps.Integers[I_N_TT_SEGMENTS]-1 : -1;

				vector <int> v;
				v.push_back(v0);
				v.push_back(v0+dj);
				v.push_back(v0+dj-ps.Integers[I_N_TT_SEGMENTS]);
				v.push_back(v0-ps.Integers[I_N_TT_SEGMENTS]);
				gs->Polys.push_back(v);
				gs->PolyLabels.push_back(Poly_Attr);
				gs->PolyCount.push_back(1);
				gs->PolyHolesCount.push_back(0);
				vector<int> f;
				f.push_back(gs->Polys.size()-1);
				gs->Facets.push_back(f);

			}

		}

		if (bCapped && (i == 0 || i == ps.Integers[I_N_TT_SEGMENTS_Z])) {
			gs->Polys.push_back(v_cap);
			gs->PolyLabels.push_back(Poly_Attr);
			gs->PolyCount.push_back(1);
			gs->PolyHolesCount.push_back(0);
			vector<int> f;
			f.push_back(gs->Polys.size()-1);
			gs->Facets.push_back(f);
		}

	}

}


void WritePropertiesFile(MeshStruct* ms, const char* filebase) {

	//Write face properties to file
	FILE * pFile;
	char outfile[255];
	sprintf(outfile,"%s_face_properties.txt",filebase);
	fprintf(stdout,"Writing face properties file %s\n",outfile);
	pFile = fopen (outfile,"w");
	if (pFile != NULL) {
		//fprintf(pFile,"#NFaces\n");
		//fprintf(pFile,"%d\n",N_Faces);
		fprintf(pFile,"%d\n",ms->FaceTet.size());
		fprintf(pFile,"Neighb1 Neighb2 NormX NormY NormZ CentX CentY CentZ Area AfX AfY AfZ r1X r1Y r1Z r2X r2Y r2Z eXiX eXiY eXiZ dXi FaceCoeff FaceType bElectroDiff\n");

		for (int i = 0; i < ms->FaceTet.size(); i++) {
			fprintf(pFile,"%d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %d\n",
					ms->FaceTet[i].i,ms->FaceTet[i].j,
					ms->FaceNormals[i].x,ms->FaceNormals[i].y,ms->FaceNormals[i].z,
					ms->FaceCenters[i].x,ms->FaceCenters[i].y,ms->FaceCenters[i].z,
					ms->FaceArea[i],
					ms->Af[i].x,ms->Af[i].y,ms->Af[i].z,
					ms->r1[i].x,ms->r1[i].y,ms->r1[i].z,
					ms->r2[i].x,ms->r2[i].y,ms->r2[i].z,
					ms->eXi[i].x,ms->eXi[i].y,ms->eXi[i].z,
					ms->dXi[i],
					ms->FaceCoeff[i],
					ms->FaceAttrList[i]);
		}

		fclose(pFile);

	} else {
		fprintf(stderr,"Warning: could not open face properties output file.\n");
		exit(1);
	}

	//Write tetrahedra properties file
	sprintf(outfile,"%s_tet_properties.txt",filebase);
	fprintf(stdout,"Writing tetrahedra properties file %s\n",outfile);
	pFile = fopen (outfile,"w");
	if (pFile != NULL) {
		//fprintf(pFile,"#NTetrahedra\n");
		//fprintf(pFile,"%d\n",out.numberoftetrahedra);

		fprintf(pFile,"%d\n",ms->CellCenters.size());
		fprintf(pFile,"CentX CentY CentZ V0 TTSurfaceArea SRSurfaceArea bTroponinC Type\n");

		for (int i = 0; i < ms->CellCenters.size(); i++) {
			fprintf(pFile,"%g %g %g %g %g %g %d %d",ms->CellCenters[i].x,ms->CellCenters[i].y,ms->CellCenters[i].z,ms->V0[i],
					ms->TTSurfaceArea[i],ms->SRSurfaceArea[i],ms->bTroponinC[i],ms->CellAttrList[i]);
			fprintf(pFile,"\n");
		}

		fclose(pFile);

	} else {
		fprintf(stderr,"Warning: could not open tetrahedra properties output file.\n");
		exit(1);
	}
}

void WritePolyFile(vector<double> &VX, vector<double> &VY, vector<int> &C1, vector<int> &C2, vector<double> &HX, vector<double> &HY, vector<int> &A) {

	FILE * pFile;
	char outfile[255];
	sprintf(outfile,"%s.poly",ps.Chars[I_OUTPUT_NAME]);
	pFile = fopen (outfile,"w");
	if (pFile != NULL) {

		//First line
		fprintf(pFile,"%d 2 0 0\n",VX.size());

		//Output vertex coordinates
		for (int i = 0; i < VX.size(); i++) {
			fprintf(pFile,"%d %f %f 0\n",i,VX[i],VY[i]);
		}

		//Segments header
		fprintf(pFile,"%d 1\n",C1.size());

		//Segments list
		for (int i = 0; i < C1.size(); i++) {
			fprintf(pFile,"%d %d %d %d\n",i,C1[i],C2[i],A[i]);
		}

		//Holes header
		fprintf(pFile,"%d\n",HX.size());

		//Holes list
		for (int i = 0; i < HX.size(); i++) {
			fprintf(pFile,"%d %f %f\n",i,HX[i],HY[i]);
		}

		fclose(pFile);

	} else {
		fprintf(stderr,"Warning: could not open vertices output file.\n");
		exit(1);
	}

}


//Reads parameters
void ReadParameters(ParamStruct* ps, const char* param_file) {

	//Assign default parameters
	ps->Integers[I_N_TT_SEGMENTS] = 23;
	ps->Integers[I_N_TT_SEGMENTS_Z] = 67;
	ps->Integers[I_N_JSR_SEGMENTS] = 15;
	ps->Integers[I_N_JSR_SEGMENTS_Z] = 15;
	ps->Integers[I_FLAG_TT] = 1;
	ps->Integers[I_FLAG_DYAD_GUIDE] = 0;
	ps->Integers[I_FLAG_TWO_CRUS] = 0;
	ps->Reals[I_VMAX] = 1.0;
	ps->Reals[I_QMAX] = 2.0;
	ps->Reals[I_GRID_X] = 4.0;
	ps->Reals[I_GRID_Y] = 4.0;
	ps->Reals[I_GRID_Z] = 4.0;
	ps->Reals[I_DYAD_R] = 15e-3;
	ps->Reals[I_JSR_R] = 40e-3;
	ps->Reals[I_DYAD_RES] = 31e-3;
	ps->Reals[I_R_TRPN] = 0.20;
	ps->Reals[I_RYR_SPACING] = 31e-3;
	ps->Reals[I_JSR_TRANS] = 0;
	ps->Reals[I_CRU_SEPARATION] = 0.806;
	sprintf(ps->Chars[I_OUTPUT_NAME],"../mesh");
	sprintf(ps->Chars[I_RYR_FILE],"");
	sprintf(ps->Chars[I_LCC_FILE],"");

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
    
    XMLCh* TAG_ROOT = XMLString::transcode("mesh");
    
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
    			
    			if (!strcmp(symbol,"flag_tt")) {
    				if (!strcmp(value,"true")) {
    					ps->Integers[I_FLAG_TT] = 1;
    				} else {
    					ps->Integers[I_FLAG_TT] = 0;
    				}
    			} else if (!strcmp(symbol,"n_tt_seg")) {
    				ps->Integers[I_N_TT_SEGMENTS] = atoi(value);
    			} else if (!strcmp(symbol,"n_tt_seg_z")) {
    				ps->Integers[I_N_TT_SEGMENTS_Z] = atoi(value);
    			} else if (!strcmp(symbol,"n_jsr_seg")) {
    				ps->Integers[I_N_JSR_SEGMENTS] = atoi(value);
    			} else if (!strcmp(symbol,"n_jsr_seg_z")) {
    				ps->Integers[I_N_JSR_SEGMENTS_Z] = atoi(value);
    			} else if (!strcmp(symbol,"dyad_r")) {
    				ps->Reals[I_DYAD_R] = atof(value);
    			} else if (!strcmp(symbol,"flag_dyad_guide")) {
    				if (!strcmp(value,"true")) {
    					ps->Integers[I_FLAG_DYAD_GUIDE] = 1;
    				} else {
    					ps->Integers[I_FLAG_DYAD_GUIDE] = 0;
    				}
    			} else if (!strcmp(symbol,"grid_x")) {
    				ps->Reals[I_GRID_X] = atof(value);
    			} else if (!strcmp(symbol,"grid_y")) {
    				ps->Reals[I_GRID_Y] = atof(value);
    			} else if (!strcmp(symbol,"grid_z")) {
    				ps->Reals[I_GRID_Z] = atof(value);
    			} else if (!strcmp(symbol,"jsr_r")) {
    				ps->Reals[I_JSR_R] = atof(value);
    			} else if (!strcmp(symbol,"r_trpn")) {
    				ps->Reals[I_R_TRPN] = atof(value);
    			} else if (!strcmp(symbol,"cluster_ryr")) {
    				memcpy(RyR_Cluster,value,strlen(value)*sizeof(char));
    			} else if (!strcmp(symbol,"cluster_lcc")) {
    				memcpy(LCC_Cluster,value,strlen(value)*sizeof(char));
    			}
    			
    		}
    	}
    }
    
    delete parser;
    delete errHandler;
    
    //Generate channel geometry
	ReadChannelGeometry(RyR_Cluster, ps->bRyR, ps->N_RyR, ps->N_RyR_X, ps->N_RyR_Y);
	ReadChannelGeometry(LCC_Cluster, ps->bLCC, ps->N_LCC, ps->N_LCC_X, ps->N_LCC_Y);
    
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


//Reads parameters
void ReadChannelGeometry(char* geometry, int** &bChannel, int &N_Channels, int &N_Dyad_X, int &N_Dyad_Y) {

	int N = strlen(geometry);
	N_Channels = 0;
	
	int idx = 0;
	while (idx < N && (geometry[idx] == '0' || geometry[idx] == '1')) {
		idx++;
	}
	int N_First_Line = idx;
	//fprintf(stdout,"N_First_Line = %d\n",N_First_Line);
	int N_This_Line = 0;
	int Line_Started = 0;
	int N_Lines = 1;
	while (idx < N) {
		if (geometry[idx] == '0' || geometry[idx] == '1') {
			if (!Line_Started) {
				Line_Started = 1;
				N_Lines++;
			}
			N_This_Line++;
		} else if (Line_Started) {
			//fprintf(stdout,"N_This_Line = %d\n",N_This_Line);
			if (N_This_Line != N_First_Line) {
				throw(std::runtime_error("Channel geometry non-rectangular."));
				exit(1);
			}
			N_This_Line = 0;
			Line_Started = 0;
		}
		idx++;
	}
	
	N_Dyad_X = N_First_Line;
	N_Dyad_Y = N_Lines;
	
	for (int i = 0; i < N; i++) {
		if (geometry[i] == '1') {
			N_Channels++;
		}
	}
	
	bChannel = new int*[N_Dyad_Y];
	for (int i = 0; i < N_Dyad_Y; i++) {
		bChannel[i] = new int[N_Dyad_X];
	}
	
	idx = 0;
	Line_Started = 0;
	int X = 0;
	int Y = -1;
	while (idx < N) {
		if (geometry[idx] == '0' || geometry[idx] == '1') {
			//fprintf(stdout,"X = %d, Y = %d\n",X,Y);
			if (!Line_Started) {
				Line_Started = 1;
				Y++;
				X = 0;
			} else {
				X++;
			}
			if (geometry[idx] == '0') {
				bChannel[Y][X] = 0;
			} else {
				bChannel[Y][X] = 1;
			}
		} else {
			Line_Started = 0;
		}
		idx++;
	}
	
	fprintf(stdout,"Array dimensions: %d x %d\n",N_Dyad_X,N_Dyad_Y);
	fprintf(stdout,"# Channels : %d\n",N_Channels);
	for (int i = 0; i < N_Dyad_Y; i++) {
		for (int j = 0; j < N_Dyad_X; j++) {
			fprintf(stdout,"%d",bChannel[i][j]);
		}
		fprintf(stdout,"\n");
	}

}
