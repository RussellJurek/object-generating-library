#include<iostream>
#include<RJJ_ObjGen.h>

using namespace std;

// functions using floats

int AddObjsToChunk(int * flag_vals, vector<object_props *> & detections, int NOobj, int obj_limit, int chunk_x_start, int chunk_y_start, int chunk_z_start, int chunk_x_size, int chunk_y_size, int chunk_z_size, vector<int> & check_obj_ids, int * data_metric, int * xyz_order){

  float progress;
  int obj, obj_batch, x, y, g, f;
  int temp_x[2], temp_y[2], temp_z[2];

  // reorder the datacube and subcube limits to be in x,y,z order
  temp_x[0] = chunk_x_start;
  temp_y[0] = chunk_y_start;
  temp_z[0] = chunk_z_start;
  temp_x[1] = chunk_x_size;
  temp_y[1] = chunk_y_size;
  temp_z[1] = chunk_z_size;
  switch(xyz_order[0]){
  case 1:
    chunk_x_start = temp_x[0];
    chunk_x_size = temp_x[1];
    break;
  case 2:
    chunk_x_start = temp_y[0];
    chunk_x_size = temp_y[1];
    break;
  case 3:
    chunk_x_start = temp_z[0];
    chunk_x_size = temp_z[1];
    break;
  default:
    chunk_x_start = temp_x[0];
    chunk_x_size = temp_x[1];
    break;
  }
  switch(xyz_order[1]){
  case 1:
    chunk_y_start = temp_x[0];
    chunk_y_size = temp_x[1];
    break;
  case 2:
    chunk_y_start = temp_y[0];
    chunk_y_size = temp_y[1];
    break;
  case 3:
    chunk_y_start = temp_z[0];
    chunk_y_size = temp_z[1];
    break;
  default:
    chunk_y_start = temp_y[0];
    chunk_y_size = temp_y[1];
    break;
  }
  switch(xyz_order[2]){
  case 1:
    chunk_z_start = temp_x[0];
    chunk_z_size = temp_x[1];
    break;
  case 2:
    chunk_z_start = temp_y[0];
    chunk_z_size = temp_y[1];
    break;
  case 3:
    chunk_z_start = temp_z[0];
    chunk_z_size = temp_z[1];
    break;
  default:
    chunk_z_start = temp_z[0];
    chunk_z_size = temp_z[1];
    break;
  }

  // add sources into the flag_vals array
  check_obj_ids.resize(0);
  if(NOobj > 0){
    
    std::cout << "Updating flag_vals array for existing objects . . . " << std::endl;
    progress = 0.0;
    std::cout << "0 | |:| | : | |:| | 100% complete" << std::endl;
    
    for(obj = 0; obj < NOobj; ++obj){
      
      // calculate obj_batch number for this object
      obj_batch = (int) floorf(((float) obj / (float) obj_limit));
      
      // move to the next object if it's been re-initialised
      if(detections[obj_batch][(obj - (obj_batch * obj_limit))].ShowVoxels() < 1){ 
	
	while(progress <= (((float) (obj + 1)) / ((float) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
	continue; 
	
      }
      
      // move to the next object if the bounding box is outside of the current chunk
      if((detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmax() < chunk_x_start) || (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmax() < chunk_y_start) || (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmin() >= (chunk_x_start + chunk_x_size)) || (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmin() >= (chunk_y_start + chunk_y_size))){ 
	
	while(progress <= (((float) (obj + 1)) / ((float) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
	continue; 
	
      }
      
      // add the object to the initial list of object ids to be checked
      check_obj_ids.push_back(obj);

      // reconstruct the object from the sparse representation and update the flag_vals array
      for(x = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0); x <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(1); ++x){
	
	for(y = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(2); y <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(3); ++y){
	  
	  // x,y position of object lies inside boundaries of chunk, so flag the flag_vals array with the obj number for every object string
	  if((x >= chunk_x_start) && (x < (chunk_x_start + chunk_x_size)) && (y >= chunk_y_start) && (y < (chunk_y_start + chunk_y_size))){
	    
	    for(g = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((((y - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + x - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0)))); g < detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((((y - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + x - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0) + 1))); ++g){
	      
	      for(f = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings((2 * g)); ((f <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings(((2 * g) + 1))) && (f < (chunk_z_size + chunk_z_start))); ++f){
		
		//flag_vals[(((f - chunk_z_start) * (data_x_size - chunk_x_start) * (data_y_size - chunk_y_start)) + ((y - chunk_y_start) * (data_x_size - chunk_x_start)) + x - chunk_x_start)] = obj;
		flag_vals[(((f - chunk_z_start) * data_metric[2]) + ((y - chunk_y_start) * data_metric[1]) + ((x - chunk_x_start) * data_metric[0]))] = obj;

		// for(f . . . )
	      }
	      
	      // for(g . . . )
	    }
	    
	    // if((x >= chunk_x_start) && (x < (chunk_x_start + chunk_x_size)) && (y >= chunk_y_start) && (y < (chunk_y_start + chunk_y_size)))
	  }
	  
	  // for(y = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmin(); y <= detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmax(); y++)
	}
	
	// for(x = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmin(); x <= detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmax(); x++)
      }
      
      while(progress <= (((float) (obj + 1)) / ((float) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
      
      // for(obj = 0; obj < NOobj; obj++)
    }
    
    std::cout << "* done." << std::endl; 
    
    //if(NOobj > 0)
  } 
  
  return check_obj_ids.size();
  
}

int AddObjsToChunk(long int * flag_vals, vector<object_props *> & detections, long int NOobj, int obj_limit, int chunk_x_start, int chunk_y_start, int chunk_z_start, int chunk_x_size, int chunk_y_size, int chunk_z_size, vector<long int> & check_obj_ids, int * data_metric, int * xyz_order){

  float progress;
  long int obj, obj_batch;
  int x, y, g, f;
  int temp_x[2], temp_y[2], temp_z[2];

  // reorder the datacube and subcube limits to be in x,y,z order
  temp_x[0] = chunk_x_start;
  temp_y[0] = chunk_y_start;
  temp_z[0] = chunk_z_start;
  temp_x[1] = chunk_x_size;
  temp_y[1] = chunk_y_size;
  temp_z[1] = chunk_z_size;
  switch(xyz_order[0]){
  case 1:
    chunk_x_start = temp_x[0];
    chunk_x_size = temp_x[1];
    break;
  case 2:
    chunk_x_start = temp_y[0];
    chunk_x_size = temp_y[1];
    break;
  case 3:
    chunk_x_start = temp_z[0];
    chunk_x_size = temp_z[1];
    break;
  default:
    chunk_x_start = temp_x[0];
    chunk_x_size = temp_x[1];
    break;
  }
  switch(xyz_order[1]){
  case 1:
    chunk_y_start = temp_x[0];
    chunk_y_size = temp_x[1];
    break;
  case 2:
    chunk_y_start = temp_y[0];
    chunk_y_size = temp_y[1];
    break;
  case 3:
    chunk_y_start = temp_z[0];
    chunk_y_size = temp_z[1];
    break;
  default:
    chunk_y_start = temp_y[0];
    chunk_y_size = temp_y[1];
    break;
  }
  switch(xyz_order[2]){
  case 1:
    chunk_z_start = temp_x[0];
    chunk_z_size = temp_x[1];
    break;
  case 2:
    chunk_z_start = temp_y[0];
    chunk_z_size = temp_y[1];
    break;
  case 3:
    chunk_z_start = temp_z[0];
    chunk_z_size = temp_z[1];
    break;
  default:
    chunk_z_start = temp_z[0];
    chunk_z_size = temp_z[1];
    break;
  }

  // add sources into the flag_vals array
  check_obj_ids.resize(0);
  if(NOobj > 0){
    
    std::cout << "Updating flag_vals array for existing objects . . . " << std::endl;
    progress = 0.0;
    std::cout << "0 | |:| | : | |:| | 100% complete" << std::endl;
    
    for(obj = 0; obj < NOobj; ++obj){
      
      // calculate obj_batch number for this object
      obj_batch = (long int) floor(((double) obj / (double) obj_limit));
      
      // move to the next object if it's been re-initialised
      if(detections[obj_batch][(obj - (obj_batch * obj_limit))].ShowVoxels() < 1){ 
	
	while(progress <= (((double) (obj + 1)) / ((double) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
	continue; 
	
      }
      
      // move to the next object if the bounding box is outside of the current chunk
      if((detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmax() < chunk_x_start) || (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmax() < chunk_y_start) || (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmin() >= (chunk_x_start + chunk_x_size)) || (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmin() >= (chunk_y_start + chunk_y_size))){ 
	
	while(progress <= (((double) (obj + 1)) / ((double) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
	continue; 
	
      }
      
      // add the object to the initial list of object ids to be checked
      check_obj_ids.push_back(obj);

      // reconstruct the object from the sparse representation and update the flag_vals array
      for(x = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0); x <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(1); ++x){
	
	for(y = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(2); y <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(3); ++y){
	  
	  // x,y position of object lies inside boundaries of chunk, so flag the flag_vals array with the obj number for every object string
	  if((x >= chunk_x_start) && (x < (chunk_x_start + chunk_x_size)) && (y >= chunk_y_start) && (y < (chunk_y_start + chunk_y_size))){
	    
	    for(g = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((((y - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + x - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0)))); g < detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((((y - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + x - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0) + 1))); ++g){
	      
	      for(f = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings((2 * g)); ((f <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings(((2 * g) + 1))) && (f < (chunk_z_size + chunk_z_start))); ++f){
		
		//flag_vals[(((f - chunk_z_start) * (data_x_size - chunk_x_start) * (data_y_size - chunk_y_start)) + ((y - chunk_y_start) * (data_x_size - chunk_x_start)) + x - chunk_x_start)] = obj;
		flag_vals[(((f - chunk_z_start) * data_metric[2]) + ((y - chunk_y_start) * data_metric[1]) + ((x - chunk_x_start) * data_metric[0]))] = obj;

		// for(f . . . )
	      }
	      
	      // for(g . . . )
	    }
	    
	    // if((x >= chunk_x_start) && (x < (chunk_x_start + chunk_x_size)) && (y >= chunk_y_start) && (y < (chunk_y_start + chunk_y_size)))
	  }
	  
	  // for(y = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmin(); y <= detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmax(); y++)
	}
	
	// for(x = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmin(); x <= detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmax(); x++)
      }
      
      while(progress <= (((double) (obj + 1)) / ((double) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
      
      // for(obj = 0; obj < NOobj; obj++)
    }
    
    std::cout << "* done." << std::endl; 
    
    //if(NOobj > 0)
  } 
  
  return check_obj_ids.size();
  
}

// functions using doubles

int AddObjsToChunk(int * flag_vals, vector<object_props_dbl *> & detections, int NOobj, int obj_limit, int chunk_x_start, int chunk_y_start, int chunk_z_start, int chunk_x_size, int chunk_y_size, int chunk_z_size, vector<int> & check_obj_ids, int * data_metric, int * xyz_order){

  float progress;
  int obj, obj_batch, x, y, g, f;
  int temp_x[2], temp_y[2], temp_z[2];

  // reorder the datacube and subcube limits to be in x,y,z order
  temp_x[0] = chunk_x_start;
  temp_y[0] = chunk_y_start;
  temp_z[0] = chunk_z_start;
  temp_x[1] = chunk_x_size;
  temp_y[1] = chunk_y_size;
  temp_z[1] = chunk_z_size;
  switch(xyz_order[0]){
  case 1:
    chunk_x_start = temp_x[0];
    chunk_x_size = temp_x[1];
    break;
  case 2:
    chunk_x_start = temp_y[0];
    chunk_x_size = temp_y[1];
    break;
  case 3:
    chunk_x_start = temp_z[0];
    chunk_x_size = temp_z[1];
    break;
  default:
    chunk_x_start = temp_x[0];
    chunk_x_size = temp_x[1];
    break;
  }
  switch(xyz_order[1]){
  case 1:
    chunk_y_start = temp_x[0];
    chunk_y_size = temp_x[1];
    break;
  case 2:
    chunk_y_start = temp_y[0];
    chunk_y_size = temp_y[1];
    break;
  case 3:
    chunk_y_start = temp_z[0];
    chunk_y_size = temp_z[1];
    break;
  default:
    chunk_y_start = temp_y[0];
    chunk_y_size = temp_y[1];
    break;
  }
  switch(xyz_order[2]){
  case 1:
    chunk_z_start = temp_x[0];
    chunk_z_size = temp_x[1];
    break;
  case 2:
    chunk_z_start = temp_y[0];
    chunk_z_size = temp_y[1];
    break;
  case 3:
    chunk_z_start = temp_z[0];
    chunk_z_size = temp_z[1];
    break;
  default:
    chunk_z_start = temp_z[0];
    chunk_z_size = temp_z[1];
    break;
  }

  // add sources into the flag_vals array
  check_obj_ids.resize(0);
  if(NOobj > 0){
    
    std::cout << "Updating flag_vals array for existing objects . . . " << std::endl;
    progress = 0.0;
    std::cout << "0 | |:| | : | |:| | 100% complete" << std::endl;
    
    for(obj = 0; obj < NOobj; ++obj){
      
      // calculate obj_batch number for this object
      obj_batch = (int) floorf(((float) obj / (float) obj_limit));
      
      // move to the next object if it's been re-initialised
      if(detections[obj_batch][(obj - (obj_batch * obj_limit))].ShowVoxels() < 1){ 
	
	while(progress <= (((float) (obj + 1)) / ((float) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
	continue; 
	
      }
      
      // move to the next object if the bounding box is outside of the current chunk
      if((detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmax() < chunk_x_start) || (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmax() < chunk_y_start) || (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmin() >= (chunk_x_start + chunk_x_size)) || (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmin() >= (chunk_y_start + chunk_y_size))){ 
	
	while(progress <= (((float) (obj + 1)) / ((float) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
	continue; 
	
      }
      
      // add the object to the initial list of object ids to be checked
      check_obj_ids.push_back(obj);

      // reconstruct the object from the sparse representation and update the flag_vals array
      for(x = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0); x <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(1); ++x){
	
	for(y = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(2); y <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(3); ++y){
	  
	  // x,y position of object lies inside boundaries of chunk, so flag the flag_vals array with the obj number for every object string
	  if((x >= chunk_x_start) && (x < (chunk_x_start + chunk_x_size)) && (y >= chunk_y_start) && (y < (chunk_y_start + chunk_y_size))){
	    
	    for(g = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((((y - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + x - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0)))); g < detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((((y - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + x - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0) + 1))); ++g){
	      
	      for(f = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings((2 * g)); ((f <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings(((2 * g) + 1))) && (f < (chunk_z_size + chunk_z_start))); ++f){
		
		//flag_vals[(((f - chunk_z_start) * (data_x_size - chunk_x_start) * (data_y_size - chunk_y_start)) + ((y - chunk_y_start) * (data_x_size - chunk_x_start)) + x - chunk_x_start)] = obj;
		flag_vals[(((f - chunk_z_start) * data_metric[2]) + ((y - chunk_y_start) * data_metric[1]) + ((x - chunk_x_start) * data_metric[0]))] = obj;

		// for(f . . . )
	      }
	      
	      // for(g . . . )
	    }
	    
	    // if((x >= chunk_x_start) && (x < (chunk_x_start + chunk_x_size)) && (y >= chunk_y_start) && (y < (chunk_y_start + chunk_y_size)))
	  }
	  
	  // for(y = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmin(); y <= detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmax(); y++)
	}
	
	// for(x = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmin(); x <= detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmax(); x++)
      }
      
      while(progress <= (((float) (obj + 1)) / ((float) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
      
      // for(obj = 0; obj < NOobj; obj++)
    }
    
    std::cout << "* done." << std::endl; 
    
    //if(NOobj > 0)
  } 
  
  return check_obj_ids.size();
  
}

int AddObjsToChunk(long int * flag_vals, vector<object_props_dbl *> & detections, long int NOobj, int obj_limit, int chunk_x_start, int chunk_y_start, int chunk_z_start, int chunk_x_size, int chunk_y_size, int chunk_z_size, vector<long int> & check_obj_ids, int * data_metric, int * xyz_order){

  float progress;
  long int obj, obj_batch;
  int x, y, g, f;
  int temp_x[2], temp_y[2], temp_z[2];

  // reorder the datacube and subcube limits to be in x,y,z order
  temp_x[0] = chunk_x_start;
  temp_y[0] = chunk_y_start;
  temp_z[0] = chunk_z_start;
  temp_x[1] = chunk_x_size;
  temp_y[1] = chunk_y_size;
  temp_z[1] = chunk_z_size;
  switch(xyz_order[0]){
  case 1:
    chunk_x_start = temp_x[0];
    chunk_x_size = temp_x[1];
    break;
  case 2:
    chunk_x_start = temp_y[0];
    chunk_x_size = temp_y[1];
    break;
  case 3:
    chunk_x_start = temp_z[0];
    chunk_x_size = temp_z[1];
    break;
  default:
    chunk_x_start = temp_x[0];
    chunk_x_size = temp_x[1];
    break;
  }
  switch(xyz_order[1]){
  case 1:
    chunk_y_start = temp_x[0];
    chunk_y_size = temp_x[1];
    break;
  case 2:
    chunk_y_start = temp_y[0];
    chunk_y_size = temp_y[1];
    break;
  case 3:
    chunk_y_start = temp_z[0];
    chunk_y_size = temp_z[1];
    break;
  default:
    chunk_y_start = temp_y[0];
    chunk_y_size = temp_y[1];
    break;
  }
  switch(xyz_order[2]){
  case 1:
    chunk_z_start = temp_x[0];
    chunk_z_size = temp_x[1];
    break;
  case 2:
    chunk_z_start = temp_y[0];
    chunk_z_size = temp_y[1];
    break;
  case 3:
    chunk_z_start = temp_z[0];
    chunk_z_size = temp_z[1];
    break;
  default:
    chunk_z_start = temp_z[0];
    chunk_z_size = temp_z[1];
    break;
  }

  // add sources into the flag_vals array
  check_obj_ids.resize(0);
  if(NOobj > 0){
    
    std::cout << "Updating flag_vals array for existing objects . . . " << std::endl;
    progress = 0.0;
    std::cout << "0 | |:| | : | |:| | 100% complete" << std::endl;
    
    for(obj = 0; obj < NOobj; ++obj){
      
      // calculate obj_batch number for this object
      obj_batch = (long int) floor(((double) obj / (double) obj_limit));
      
      // move to the next object if it's been re-initialised
      if(detections[obj_batch][(obj - (obj_batch * obj_limit))].ShowVoxels() < 1){ 
	
	while(progress <= (((double) (obj + 1)) / ((double) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
	continue; 
	
      }
      
      // move to the next object if the bounding box is outside of the current chunk
      if((detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmax() < chunk_x_start) || (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmax() < chunk_y_start) || (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmin() >= (chunk_x_start + chunk_x_size)) || (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmin() >= (chunk_y_start + chunk_y_size))){ 
	
	while(progress <= (((double) (obj + 1)) / ((double) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
	continue; 
	
      }
      
      // add the object to the initial list of object ids to be checked
      check_obj_ids.push_back(obj);

      // reconstruct the object from the sparse representation and update the flag_vals array
      for(x = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0); x <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(1); ++x){
	
	for(y = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(2); y <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(3); ++y){
	  
	  // x,y position of object lies inside boundaries of chunk, so flag the flag_vals array with the obj number for every object string
	  if((x >= chunk_x_start) && (x < (chunk_x_start + chunk_x_size)) && (y >= chunk_y_start) && (y < (chunk_y_start + chunk_y_size))){
	    
	    for(g = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((((y - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + x - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0)))); g < detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((((y - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + x - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0) + 1))); ++g){
	      
	      for(f = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings((2 * g)); ((f <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings(((2 * g) + 1))) && (f < (chunk_z_size + chunk_z_start))); ++f){
		
		//flag_vals[(((f - chunk_z_start) * (data_x_size - chunk_x_start) * (data_y_size - chunk_y_start)) + ((y - chunk_y_start) * (data_x_size - chunk_x_start)) + x - chunk_x_start)] = obj;
		flag_vals[(((f - chunk_z_start) * data_metric[2]) + ((y - chunk_y_start) * data_metric[1]) + ((x - chunk_x_start) * data_metric[0]))] = obj;

		// for(f . . . )
	      }
	      
	      // for(g . . . )
	    }
	    
	    // if((x >= chunk_x_start) && (x < (chunk_x_start + chunk_x_size)) && (y >= chunk_y_start) && (y < (chunk_y_start + chunk_y_size)))
	  }
	  
	  // for(y = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmin(); y <= detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmax(); y++)
	}
	
	// for(x = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmin(); x <= detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmax(); x++)
      }
      
      while(progress <= (((double) (obj + 1)) / ((double) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
      
      // for(obj = 0; obj < NOobj; obj++)
    }
    
    std::cout << "* done." << std::endl; 
    
    //if(NOobj > 0)
  } 
  
  return check_obj_ids.size();
  
}


