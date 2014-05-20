#include<fstream>
#include<iostream>
#include<string>
#include<sstream>
#include<iomanip>
#include<cmath>
#include<ctime>
#include<ctype.h>
#include<vector>
extern "C" {

#include<cpgplot.h>

}

#ifndef RJJ_ObjGen_Plots
#define RJJ_ObjGen_Plots

using namespace std;

#define SPEC_PLOT_LIMIT 1.0e5

extern void HeapSort(int n, float ra[]);

// functions using floats

extern void CreateObjPlots(int m, int plot_mode, std::string output_code, int * xyz_order, int NOx, int NOy, int NOf, vector<object_props *> & detections, int NOobj, int obj_limit, int ctype3, float crpix3, float crval3, float cdelt3, float restfreq);
extern void CreateObjPlots(int m, int plot_mode, std::string output_code, int * xyz_order, int NOx, int NOy, int NOf, vector<object_props *> & detections, long int NOobj, int obj_limit, int ctype3, float crpix3, float crval3, float cdelt3, float restfreq);

// functions using doubles

extern void CreateObjPlots(int m, int plot_mode, std::string output_code, int * xyz_order, int NOx, int NOy, int NOf, vector<object_props_dbl *> & detections, int NOobj, int obj_limit, int ctype3, float crpix3, float crval3, float cdelt3, float restfreq);
extern void CreateObjPlots(int m, int plot_mode, std::string output_code, int * xyz_order, int NOx, int NOy, int NOf, vector<object_props_dbl *> & detections, long int NOobj, int obj_limit, int ctype3, float crpix3, float crval3, float cdelt3, float restfreq);

#endif

