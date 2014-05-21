#include<RJJ_ObjGen.h>
#include<RJJ_ObjGen_WCS.h>

using namespace std;

void PrintCatalogueHeader_WCS(std::fstream& output_file, int cat_mode){

  if(cat_mode > 0){
    
    output_file << "# ID NO_voxels RA Dec Freq RA_L Dec_L Freq_L RA_pix Dec_pix Freq_pix RA_L_pix Dec_L_pix Freq_L_pix RA_L_Ppix Dec_L_Ppix Freq_L_Ppix RA_L_Npix Dec_L_Npix Freq_L_Npix Total_L Total_L_Ppix Total_L_Npix Avg_L Std_Dev_L RMS_L Min_L Max_L RA_min RA_max DEC_min DEC_max Freq_min Freq_max Ispec_max_chan Ispec_w50_min Ispec_w50_max Ispec_w20_min Ispec_w20_max Ispec_cmax_chan Ispec_cw50_min Ispec_cw50_max Ispec_cw20_min Ispec_cw20_max : #strings grid : strings range " << std::endl;

  } else {

    output_file << "# ID NO_voxels RA Dec Freq RA_L Dec_L Freq_L RA_pix Dec_pix Freq_pix RA_L_pix Dec_L_pix Freq_L_pix RA_L_Ppix Dec_L_Ppix Freq_L_Ppix RA_L_Npix Dec_L_Npix Freq_L_Npix Total_L Total_L_Ppix Total_L_Npix Avg_L Std_Dev_L RMS_L Min_L Max_L RA_min RA_max DEC_min DEC_max Freq_min Freq_max Ispec_max_chan Ispec_w50_min Ispec_w50_max Ispec_w20_min Ispec_w20_max Ispec_cmax_chan Ispec_cw50_min Ispec_cw50_max Ispec_cw20_min Ispec_cw20_max " << std::endl;

  }


}

// functions using floats

int CreateCatalogue_WCS(std::fstream& output_file, vector<object_props *> & detections, int NOobj, int obj_limit, int cat_mode, fitsfile * fits_input){

  int k,m,obj_batch;
  float progress;

  char * header;
  int status,nkeys,nreject,nwcs,stat[NWCSFIX];
  struct wcsprm *wcs;
  double wcs_vals[6],world[1][2],phi[1],theta[1],imgcrd[1][2],pixcrd[1][2];

  // extract wcs from fits_input
  status = 0;
  fits_hdr2str(fits_input,1,NULL,0,&header,&nkeys,&status);
  wcspih(header,nkeys,WCSHDR_all,2,&nreject,&nwcs,&wcs);
  fits_read_wcstab(fits_input,wcs->nwtb,(wtbarr *)wcs->wtb,&status);
  wcsfix(7,0,wcs,stat);  
  wcsset(wcs);

  progress = 0.0;
  std::cout << "0 | |:| | : | |:| | 100% complete" << std::endl;
  m = 1;
  for(k = 0; k < NOobj; ++k){
      
    // calculate obj_batch
    obj_batch = (int) floorf(((float) k / (float) obj_limit));
    
    // move on if this object has been re-initialised
    if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ 
      
      while(progress <= (((float) (k + 1)) / ((float) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
      continue; 
      
    }
    
    // calculate properties
    detections[obj_batch][(k - (obj_batch * obj_limit))].CalcProps();
    
    // calculate RA,Dec values
    pixcrd[0][0] = (double) detections[obj_batch][(k - (obj_batch * obj_limit))].GetRA();
    pixcrd[0][1] = (double) detections[obj_batch][(k - (obj_batch * obj_limit))].GetDEC();
    wcsp2s(wcs,1,2,pixcrd[0],imgcrd[0],phi,theta,world[0],stat);
    wcs_vals[0] = world[0][0];
    wcs_vals[1] = world[0][1];
    wcs_vals[2] = -99;
    pixcrd[0][0] = (double) detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAi();
    pixcrd[0][1] = (double) detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECi();
    wcsp2s(wcs,1,2,pixcrd[0],imgcrd[0],phi,theta,world[0],stat);
    wcs_vals[3] = world[0][0];
    wcs_vals[4] = world[0][1];
    wcs_vals[5] = -99;
   
    // write object properties to output file/terminal
    detections[obj_batch][(k - (obj_batch * obj_limit))].ShowAll_file_WCS(m,output_file,cat_mode,wcs_vals);
    ++m;

    while(progress <= (((float) (k + 1)) / ((float) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }

    //for(k = 0; k < NOobj; k++) 
  }
  std::cout << "* done." << std::endl;
  std::cout << (m - 1) << " objects written to output catalogue." << std::endl;

  // free up the memory used by the wcs struct
  wcsvfree(&nwcs,&wcs);

  return m;

}

long int CreateCatalogue_WCS(std::fstream& output_file, vector<object_props *> & detections, long int NOobj, int obj_limit, int cat_mode, fitsfile * fits_input){

  long int k,m,obj_batch;
  float progress;

  char * header;
  int status,nkeys,nreject,nwcs,stat[NWCSFIX];
  struct wcsprm *wcs;
  double wcs_vals[6],world[1][2],phi[1],theta[1],imgcrd[1][2],pixcrd[1][2];

  // extract wcs from fits_input
  status = 0;
  fits_hdr2str(fits_input,1,NULL,0,&header,&nkeys,&status);
  wcspih(header,nkeys,WCSHDR_all,2,&nreject,&nwcs,&wcs);
  fits_read_wcstab(fits_input,wcs->nwtb,(wtbarr *)wcs->wtb,&status);
  wcsfix(7,0,wcs,stat);  
  wcsset(wcs);

  progress = 0.0;
  std::cout << "0 | |:| | : | |:| | 100% complete" << std::endl;
  m = 1;
  for(k = 0; k < NOobj; ++k){
      
    // calculate obj_batch
    obj_batch = (long int) floor(((double) k / (double) obj_limit));
    
    // move on if this object has been re-initialised
    if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ 
      
      while(progress <= (((double) (k + 1)) / ((double) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
      continue; 
      
    }
    
    // calculate properties
    detections[obj_batch][(k - (obj_batch * obj_limit))].CalcProps();
    
    // calculate RA,Dec values
    pixcrd[0][0] = (double) detections[obj_batch][(k - (obj_batch * obj_limit))].GetRA();
    pixcrd[0][1] = (double) detections[obj_batch][(k - (obj_batch * obj_limit))].GetDEC();
    wcsp2s(wcs,1,2,pixcrd[0],imgcrd[0],phi,theta,world[0],stat);
    wcs_vals[0] = world[0][0];
    wcs_vals[1] = world[0][1];
    wcs_vals[2] = -99;
    pixcrd[0][0] = (double) detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAi();
    pixcrd[0][1] = (double) detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECi();
    wcsp2s(wcs,1,2,pixcrd[0],imgcrd[0],phi,theta,world[0],stat);
    wcs_vals[3] = world[0][0];
    wcs_vals[4] = world[0][1];
    wcs_vals[5] = -99;
   
    // write object properties to output file/terminal
    detections[obj_batch][(k - (obj_batch * obj_limit))].ShowAll_file_WCS(m,output_file,cat_mode,wcs_vals);
    ++m;

    while(progress <= (((double) (k + 1)) / ((double) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }

    //for(k = 0; k < NOobj; k++) 
  }
  std::cout << "* done." << std::endl;
  std::cout << (m - 1) << " objects written to output catalogue." << std::endl;

  // free up the memory used by the wcs struct
  wcsvfree(&nwcs,&wcs);

  return m;

}

// functions using doubles

int CreateCatalogue_WCS(std::fstream& output_file, vector<object_props_dbl *> & detections, int NOobj, int obj_limit, int cat_mode, fitsfile * fits_input){

  int k,m,obj_batch;
  float progress;

  char * header;
  int status,nkeys,nreject,nwcs,stat[NWCSFIX];
  struct wcsprm *wcs;
  double wcs_vals[6],world[1][2],phi[1],theta[1],imgcrd[1][2],pixcrd[1][2];

  // extract wcs from fits_input
  status = 0;
  fits_hdr2str(fits_input,1,NULL,0,&header,&nkeys,&status);
  wcspih(header,nkeys,WCSHDR_all,2,&nreject,&nwcs,&wcs);
  fits_read_wcstab(fits_input,wcs->nwtb,(wtbarr *)wcs->wtb,&status);
  wcsfix(7,0,wcs,stat);  
  wcsset(wcs);

  progress = 0.0;
  std::cout << "0 | |:| | : | |:| | 100% complete" << std::endl;
  m = 1;
  for(k = 0; k < NOobj; ++k){
      
    // calculate obj_batch
    obj_batch = (int) floorf(((float) k / (float) obj_limit));
    
    // move on if this object has been re-initialised
    if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ 
      
      while(progress <= (((float) (k + 1)) / ((float) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
      continue; 
      
    }
    
    // calculate properties
    detections[obj_batch][(k - (obj_batch * obj_limit))].CalcProps();
    
    // calculate RA,Dec values
    pixcrd[0][0] = (double) detections[obj_batch][(k - (obj_batch * obj_limit))].GetRA();
    pixcrd[0][1] = (double) detections[obj_batch][(k - (obj_batch * obj_limit))].GetDEC();
    wcsp2s(wcs,1,2,pixcrd[0],imgcrd[0],phi,theta,world[0],stat);
    wcs_vals[0] = world[0][0];
    wcs_vals[1] = world[0][1];
    wcs_vals[2] = -99;
    pixcrd[0][0] = (double) detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAi();
    pixcrd[0][1] = (double) detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECi();
    wcsp2s(wcs,1,2,pixcrd[0],imgcrd[0],phi,theta,world[0],stat);
    wcs_vals[3] = world[0][0];
    wcs_vals[4] = world[0][1];
    wcs_vals[5] = -99;
   
    // write object properties to output file/terminal
    detections[obj_batch][(k - (obj_batch * obj_limit))].ShowAll_file_WCS(m,output_file,cat_mode,wcs_vals);
    ++m;

    while(progress <= (((float) (k + 1)) / ((float) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }

    //for(k = 0; k < NOobj; k++) 
  }
  std::cout << "* done." << std::endl;
  std::cout << (m - 1) << " objects written to output catalogue." << std::endl;

  // free up the memory used by the wcs struct
  wcsvfree(&nwcs,&wcs);

  return m;

}

long int CreateCatalogue_WCS(std::fstream& output_file, vector<object_props_dbl *> & detections, long int NOobj, int obj_limit, int cat_mode, fitsfile * fits_input){

  long int k,m,obj_batch;
  float progress;

  char * header;
  int status,nkeys,nreject,nwcs,stat[NWCSFIX];
  struct wcsprm *wcs;
  double wcs_vals[6],world[1][2],phi[1],theta[1],imgcrd[1][2],pixcrd[1][2];

  // extract wcs from fits_input
  status = 0;
  fits_hdr2str(fits_input,1,NULL,0,&header,&nkeys,&status);
  wcspih(header,nkeys,WCSHDR_all,2,&nreject,&nwcs,&wcs);
  fits_read_wcstab(fits_input,wcs->nwtb,(wtbarr *)wcs->wtb,&status);
  wcsfix(7,0,wcs,stat);  
  wcsset(wcs);

  progress = 0.0;
  std::cout << "0 | |:| | : | |:| | 100% complete" << std::endl;
  m = 1;
  for(k = 0; k < NOobj; ++k){
      
    // calculate obj_batch
    obj_batch = (long int) floor(((double) k / (double) obj_limit));
    
    // move on if this object has been re-initialised
    if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ 
      
      while(progress <= (((double) (k + 1)) / ((double) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
      continue; 
      
    }
    
    // calculate properties
    detections[obj_batch][(k - (obj_batch * obj_limit))].CalcProps();
    
    // calculate RA,Dec values
    pixcrd[0][0] = (double) detections[obj_batch][(k - (obj_batch * obj_limit))].GetRA();
    pixcrd[0][1] = (double) detections[obj_batch][(k - (obj_batch * obj_limit))].GetDEC();
    wcsp2s(wcs,1,2,pixcrd[0],imgcrd[0],phi,theta,world[0],stat);
    wcs_vals[0] = world[0][0];
    wcs_vals[1] = world[0][1];
    wcs_vals[2] = -99;
    pixcrd[0][0] = (double) detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAi();
    pixcrd[0][1] = (double) detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECi();
    wcsp2s(wcs,1,2,pixcrd[0],imgcrd[0],phi,theta,world[0],stat);
    wcs_vals[3] = world[0][0];
    wcs_vals[4] = world[0][1];
    wcs_vals[5] = -99;
   
    // write object properties to output file/terminal
    detections[obj_batch][(k - (obj_batch * obj_limit))].ShowAll_file_WCS(m,output_file,cat_mode,wcs_vals);
    ++m;

    while(progress <= (((double) (k + 1)) / ((double) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }

    //for(k = 0; k < NOobj; k++) 
  }
  std::cout << "* done." << std::endl;
  std::cout << (m - 1) << " objects written to output catalogue." << std::endl;

  // free up the memory used by the wcs struct
  wcsvfree(&nwcs,&wcs);

  return m;

}

