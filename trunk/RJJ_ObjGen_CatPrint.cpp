#include<iostream>
#include<RJJ_ObjGen.h>

void PrintCatalogueHeader(std::fstream& output_file, int cat_mode){

  if(cat_mode > 0){
    
    output_file << "# ID NO_voxels RA Dec Freq RA_L Dec_L Freq_L Total_L Avg_L Std_Dev_L RMS_L Min_L Max_L RA_min RA_max DEC_min DEC_max Freq_min Freq_max Ispec_max_chan Ispec_w50_min Ispec_w50_max Ispec_w20_min Ispec_w20_max Ispec_cmax_chan Ispec_cw50_min Ispec_cw50_max Ispec_cw20_min Ispec_cw20_max : #strings grid : strings range " << std::endl;

  } else {

    output_file << "# ID NO_voxels RA Dec Freq RA_L Dec_L Freq_L Total_L Avg_L Std_Dev_L RMS_L Min_L Max_L RA_min RA_max DEC_min DEC_max Freq_min Freq_max Ispec_max_chan Ispec_w50_min Ispec_w50_max Ispec_w20_min Ispec_w20_max Ispec_cmax_chan Ispec_cw50_min Ispec_cw50_max Ispec_cw20_min Ispec_cw20_max " << std::endl;

  }


}

int CreateCatalogue(std::fstream& output_file, object_props ** detections, int NOobj, int obj_limit, int cat_mode){

  int k,m,obj_batch;
  float progress;

  progress = 0.0;
  std::cout << "0 | |:| | : | |:| | 100% complete" << std::endl;
  m = 1;
  for(k = 0; k < NOobj; k++){
      
    // calculate obj_batch
    obj_batch = floorf(((float) k / (float) obj_limit));
    
    // move on if this object has been re-initialised
    if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ 
      
      while(progress <= (((float) (k + 1)) / ((float) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
      continue; 
      
    }
    
    // calculate properties
    detections[obj_batch][(k - (obj_batch * obj_limit))].CalcProps();
    
    // write object properties to output file/terminal
    detections[obj_batch][(k - (obj_batch * obj_limit))].ShowAll_file(m,output_file,cat_mode);
    m++;

    while(progress <= (((float) (k + 1)) / ((float) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }

    //for(k = 0; k < NOobj; k++) 
  }
  std::cout << "* done." << std::endl;
  std::cout << (m - 1) << " objects written to output catalogue." << std::endl;

  return m;

}

