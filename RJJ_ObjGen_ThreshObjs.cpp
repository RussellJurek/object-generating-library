#include<iostream>
#include<RJJ_ObjGen.h>

void ThresholdObjs(object_props ** detections, int NOobj, int obj_limit, int min_x_size, int min_y_size, int min_z_size, int min_v_size, float intens_thresh_min, float intens_thresh_max, int min_LoS_count){
  
  float progress;
  int i, j, k, obj_batch;

  progress = 0.0;
  std::cout << "0 | |:| | : | |:| | 100% complete" << std::endl;
  for(k = 0; k < NOobj; k++){
	
    // calculate the obj_batch value for the existing object
    obj_batch = (int) floorf(((float) k / (float) obj_limit));
	
    // move to the next object if this one has been re-initialised
    if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ 

      while(progress <= (((float) (k + 1)) / ((float) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
      continue; 

    }
	
    // count the number of LOSs through this object that contain an object section
    i = 0;
    for(j = 0; j < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); j++){

      if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_grid((j + 1)) > detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_grid(j)){ i++; }

    }

    // apply the size threshold, and if it fails re-initialise the object and pop it's id to the list
    // of available obj_ids
    if(((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) < min_x_size)  || ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1) < min_y_size) || ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1) < min_z_size) || (detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < min_v_size) || (detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() < intens_thresh_min) || (detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() > intens_thresh_max) || (i < min_LoS_count)){

      // re-initialise object
      detections[obj_batch][(k - (obj_batch * obj_limit))].ReInit();
      if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_update() != 0) && (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) >= 0)){
	
	detections[obj_batch][(k - (obj_batch * obj_limit))].ReInit_srep();
	detections[obj_batch][(k - (obj_batch * obj_limit))].ReInit_mini();
	
      }
      detections[obj_batch][(k - (obj_batch * obj_limit))].ReInit_size();
      detections[obj_batch][(k - (obj_batch * obj_limit))].Set_srep_update(0);

    }
    
    while(progress <= (((float) (k + 1)) / ((float) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
      
    // for(i = 0; i < NO_check_obj_ids; i++)
  }
  std::cout << "* done." << std::endl;

}

