#include<iostream>
#include<RJJ_ObjGen.h>
#include<vector>
extern "C" {

#include<fitsio.h>
#include<wcslib.h>

}

#ifndef RJJ_ObjGen_WCS
#define RJJ_ObjGen_WCS

using namespace std;

void PrintCatalogueHeader_WCS(std::fstream& output_file, int cat_mode);

// functions using floats

int CreateCatalogue_WCS(std::fstream& output_file, vector<object_props *> & detections, int NOobj, int obj_limit, int cat_mode, fitsfile * fits_input);
long int CreateCatalogue_WCS(std::fstream& output_file, vector<object_props *> & detections, long int NOobj, int obj_limit, int cat_mode, fitsfile * fits_input);

// functions using doubles

int CreateCatalogue_WCS(std::fstream& output_file, vector<object_props_dbl *> & detections, int NOobj, int obj_limit, int cat_mode, fitsfile * fits_input);
long int CreateCatalogue_WCS(std::fstream& output_file, vector<object_props_dbl *> & detections, long int NOobj, int obj_limit, int cat_mode, fitsfile * fits_input);

#endif

