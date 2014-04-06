#include<fstream>
#include<iostream>
#include<string>
#include<sstream>
#include<iomanip>
#include<cmath>
#include<ctime>
#include<ctype.h>
extern "C" {

#include<cpgplot.h>

}
#define SPEC_PLOT_LIMIT 1.0e5

extern void HeapSort(int n, float ra[]);

extern void CreateObjPlots(int m, int plot_mode, std::string output_code, int * xyz_order, int NOx, int NOy, int NOf, object_props ** detections, int NOobj, int obj_limit, int obj_batch_limit, int ctype3, float crpix3, float crval3, float cdelt3, float restfreq);


