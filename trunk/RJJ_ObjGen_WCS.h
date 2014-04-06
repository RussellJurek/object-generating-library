#include<iostream>
extern "C" {

#include<fitsio.h>
#include<wcslib.h>

}

void PrintCatalogueHeader_WCS(std::fstream& output_file, int cat_mode);

int CreateCatalogue_WCS(std::fstream& output_file, object_props ** detections, int NOobj, int obj_limit, int cat_mode, fitsfile * fits_input);

