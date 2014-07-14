#include<iostream>
#include<RJJ_ObjGen.h>

using namespace std;

// functions using floats

int CreateObjects(float * data_vals, int * flag_vals, int size_x, int size_y, int size_z, int chunk_x_start, int chunk_y_start, int chunk_z_start, int merge_x, int merge_y, int merge_z, int min_x_size, int min_y_size, int min_z_size, int min_v_size, float intens_thresh_min, float intens_thresh_max, int flag_value, int start_obj, vector<object_props *> & detections, vector<int> & obj_ids, vector<int> & check_obj_ids, int obj_limit, int max_x_val, int max_y_val, int max_z_val, int ss_mode, int * data_metric, int * xyz_order){
 
  int x,y,z,obj,existing,sx,sy,sz,sx_start,sy_start,sz_start,sx_finish,sy_finish,sz_finish,init_limit,g;
  int i,j,k,NOi, obj_batch, obj_batch_2,NO_obj_ids,NO_check_obj_ids;
  int prev,x_start,y_start,temp_x[3],temp_y[3],temp_z[3];
  float progress;
  vector<float> temp_mom0, temp_RAPV, temp_DECPV, temp_obj_spec, temp_ref_spec, temp_vfield;
  vector<int> match_init, temp_sparse_reps_grid, temp_sparse_reps_strings;

  // convert the merging lengths from distances (in pixels/voxels) to empty space between linked detections (in pixels/voxels) --- subtract 1 --- this is a more naturual
  // definition when processing a large file in chunks, because it makes it easier to build up sources that need to be retained for the rest of the input file
  --merge_x;
  --merge_y;
  --merge_z;

  // reorder the datacube and subcube limits to be in x,y,z order
  temp_x[0] = chunk_x_start;
  temp_y[0] = chunk_y_start;
  temp_z[0] = chunk_z_start;
  temp_x[1] = size_x;
  temp_y[1] = size_y;
  temp_z[1] = size_z;
  temp_x[2] = max_x_val;
  temp_y[2] = max_y_val;
  temp_z[2] = max_z_val;
  switch(xyz_order[0]){
  case 1:
    chunk_x_start = temp_x[0];
    size_x = temp_x[1];
    max_x_val = temp_x[2];
    break;
  case 2:
    chunk_x_start = temp_y[0];
    size_x = temp_y[1];
    max_x_val = temp_y[2];
    break;
  case 3:
    chunk_x_start = temp_z[0];
    size_x = temp_z[1];
    max_x_val = temp_z[2];
    break;
  default:
    chunk_x_start = temp_x[0];
    size_x = temp_x[1];
    max_x_val = temp_x[2];
    break;
  }
  switch(xyz_order[1]){
  case 1:
    chunk_y_start = temp_x[0];
    size_y = temp_x[1];
    max_y_val = temp_x[2];
    break;
  case 2:
    chunk_y_start = temp_y[0];
    size_y = temp_y[1];
    max_y_val = temp_y[2];
    break;
  case 3:
    chunk_y_start = temp_z[0];
    size_y = temp_z[1];
    max_y_val = temp_z[2];
    break;
  default:
    chunk_y_start = temp_y[0];
    size_y = temp_y[1];
    max_y_val = temp_y[2];
    break;
  }
  switch(xyz_order[2]){
  case 1:
    chunk_z_start = temp_x[0];
    size_z = temp_x[1];
    max_z_val = temp_x[2];
    break;
  case 2:
    chunk_z_start = temp_y[0];
    size_z = temp_y[1];
    max_z_val = temp_y[2];
    break;
  case 3:
    chunk_z_start = temp_z[0];
    size_z = temp_z[1];
    max_z_val = temp_z[2];
    break;
  default:
    chunk_z_start = temp_z[0];
    size_z = temp_z[1];
    max_z_val = temp_z[2];
    break;
  }

  // create temporary arrays
  if(10000 > (((2 * merge_x) + 3) * ((2 * merge_y) + 3) * ((2 * merge_z) + 3))){
    
    init_limit = 10000;
    
  } else {
    
    init_limit = ((2 * merge_x) + 3) * ((2 * merge_y) + 3) * ((2 * merge_z) + 3) + 10;
    
  }
  match_init.reserve(init_limit);
  match_init.resize(0);

  temp_sparse_reps_grid.reserve(1000000);
  temp_sparse_reps_grid.resize(0);
  temp_sparse_reps_strings.reserve(1000000);
  temp_sparse_reps_strings.resize(0);
  temp_mom0.reserve(1000000);
  temp_mom0.resize(0);
  temp_RAPV.reserve(1000000);
  temp_RAPV.resize(0);
  temp_DECPV.reserve(1000000);
  temp_DECPV.resize(0);
  temp_obj_spec.reserve(1000000);
  temp_obj_spec.resize(0);
  temp_ref_spec.reserve(1000000);
  temp_ref_spec.resize(0);
  temp_vfield.reserve(1000000);
  temp_vfield.resize(0);
 
  // 0. initialise variables and arrays
  obj = start_obj;
  x_start = 0;
  y_start = 0;
  if(chunk_x_start > 0){ x_start = merge_x + 1; }
  if(chunk_y_start > 0){ y_start = merge_y + 1; }
  NO_obj_ids = obj_ids.size();
  NO_check_obj_ids = check_obj_ids.size();

  // 1. Create list of `coherent' objects from neighbouring voxels
  
  // for each grid point, check if it has been flagged as a source voxel, and assign
  // it an object ID if it is. If a neighbouring voxel has already been flagged, assign 
  // that ID to the grid point, otherwise, assign current value of obj and increment obj.
  std::cout << "0 | |:| | : | |:| | 100% complete" << std::endl;
  progress = 0.0;
  for(z = 0; ((z < size_z) && ((z + chunk_z_start) < max_z_val)); ++z){

    // if the z value is sufficiently large that objects have started to pop out of the merging box,
    // then check if the objects outside of the merging box are sufficiently large
    if(z > (merge_z + 1)){
      
      // initialise NOi to store the number of current check_obj_ids values that haven't moved outside of 
      // the merging volume, and that form the basis for the new check_obj_ids list
      NOi = 0;
      
      // apply size threshold to objects added since last size thresholding
      for(i = 0; i < NO_check_obj_ids; ++i){
	
	// calculate the obj_batch value for the existing object
	obj_batch = (int) floorf(((float) check_obj_ids[i] / (float) obj_limit));
	
	// move to the next object if this one has been re-initialised
	if(detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	
	// check if this object has moved outside of the merging box available to the next plane, if it 
	// has apply a size threshold, otherwise use it to seed the new list of objects to be size thresholded
	// the next time
	if(detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmax() < (chunk_z_start + z - merge_z - 1)){
	  
	  // apply the size threshold, and if it fails re-initialise the object and pop its id to the list
	  // of available obj_ids
	  if((((detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetRAmin() + 1) < min_x_size) || ((detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetDECmin() + 1) < min_y_size) || ((detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmax() - detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmin() + 1) < min_z_size) || (detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].ShowVoxels() < min_v_size)) && (((detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetRAmin() - chunk_x_start) > merge_x) || (chunk_x_start == 0)) && ((((chunk_x_start + size_x - 1 - detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetRAmax()) > merge_x) || ((chunk_x_start + size_x) >= max_x_val)) && (((detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetDECmin() - chunk_y_start) > merge_y) || (chunk_y_start == 0)) && (((chunk_y_start + size_y - 1 - detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetDECmax()) > merge_y) || ((chunk_y_start + size_y) >= max_y_val)) && (((detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmin() - chunk_z_start) > merge_z) || (chunk_z_start == 0)) && (((chunk_z_start + size_z - 1 - detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmax()) > merge_z) || ((chunk_z_start + size_z) >= max_z_val)))){
	    
	    // re-initialise object
	    detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].ReInit();
	    if((detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].Get_srep_update() != 0) && (detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].Get_srep_size(0) >= 0)){
	      
	      detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].ReInit_srep();
	      detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].ReInit_mini();

	    }
	    detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].ReInit_size();
	    detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].Set_srep_update(0);

	    // flag the `object' values within the bounding box and remove it from the flag_vals array
	    sx_start = detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetRAmin() - chunk_x_start;
	    if(sx_start < 0){ sx_start = 0; }
	    sy_start = detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetDECmin() - chunk_y_start;
	    if(sy_start < 0){ sy_start = 0; }
	    sz_start = detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmin() - chunk_z_start;
	    if(sz_start < 0){ sz_start = 0; }
	    for(sz = sz_start; ((sz <= (detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmax() - chunk_z_start)) && (sz < size_z)); ++sz){
	      
	      for(sy = sy_start; ((sy <= (detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetDECmax() - chunk_y_start)) && (sy < size_y)); ++sy){
		
		for(sx = sx_start; ((sx <= (detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetRAmax() - chunk_x_start)) && (sx < size_x)); ++sx){
		  
		  if(flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx * data_metric[0]))] == check_obj_ids[i]){ flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx * data_metric[0]))] = -99; }
		  
		  // for(sx = x_start[(match_init[i] - (obj_batch_2 * obj_limit))]; sx <= x_finish[(match_init[i] - (obj_batch_2 * obj_limit))]; ++sx)
		}
		
		// for(sy = y_start[(match_init[i] - (obj_batch_2 * obj_limit))]; sy <= y_finish[(match_init[i] - (obj_batch_2 * obj_limit))]; ++sy)
	      }
	      
	      // for(sz = z_start[(match_init[i] - (obj_batch_2 * obj_limit))]; sz <= z_finish[(match_init[i] - (obj_batch_2 * obj_limit))]; ++sz)
	    }
	    
	    // add object id to list of available ids
	    j = -1;
	    for(g = 0; g < NO_obj_ids; ++g){ if(obj_ids[g] == check_obj_ids[i]){ j = 1; break; } }
	    if(j == -1){
	      ++NO_obj_ids;
	      obj_ids.push_back(check_obj_ids[i]);
	    }

	  }
	  
	  // if(z_finish[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))] < (z - merge_z - 1))
	} else {
	  
	  check_obj_ids[NOi] = check_obj_ids[i];
	  ++NOi;
	  
	  // else . . . if(z_finish[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))] < (z - merge_z - 1))
	}
	
	// for(i = 0; i < NO_check_obj_ids; ++i)
      }

      // update NO_check_obj_ids
      NO_check_obj_ids = NOi;
      check_obj_ids.resize(NOi);
      
      // if(z > (merge_z + 2))
    }
    
    for(y = y_start; (y < size_y) && ((chunk_y_start + y) < max_y_val); ++y){
      
      for(x = x_start; (x < size_x) && ((chunk_x_start + x) < max_x_val); ++x){
	
	// if this is a source voxel, check if it is associated with a previously identified
	// source voxel
	if(flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (x*data_metric[0]))] == flag_value){ 
	  
	  // 1.  initialise variables
	  existing = flag_value;
	  NOi = 0;
	  match_init.resize(0);
	  
	  // 3. initial pass, check the voxels within the merging distance to see if any have been
	  // identified as previous objects
	  
	  // 3a. set limits of region to be searched
	  sx_start = x - merge_x - 1;
	  if(sx_start < 0){ sx_start = 0; }
	  sy_start = y - merge_y - 1;
	  if(sy_start < 0){ sy_start = 0; }
	  sz_start = z - merge_z - 1;
	  if(sz_start < 0){ sz_start = 0; }
	  
	  // 3b. search through previous region and modify existing flag if required

	  if(ss_mode == 1){
	    
	    // use a rectangular paralleliped to do the linking

	    for(sz = sz_start; (sz < z); ++sz){
	    
	      for(sy = sy_start; (sy <= (y + merge_y + 1)) && (sy < size_y); ++sy){
	      
		for(sx = sx_start; (sx <= (x + merge_x + 1)) && (sx < size_x); ++sx){
		
		  // check if this test voxel is source
		  if((flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] == flag_value) || (flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] >= 0)){
		    
		    // check if this is a new objid
		    prev = -1;
		    for(i = 0; (i < NOi) && (i < init_limit); ++i){
		      
		      if(match_init[i] == flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]){
			
			prev = 1;
			break;
			
		      }
		      
		    }
		    if(prev == -1){ 
		      
		      // add current object to list of previous matches
		      match_init.push_back(flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]);
		      ++NOi;

		    }
		    
		    // update existing flag
		    if(((flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] < existing) && (flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] >= 0)) || (existing == flag_value)){ existing = flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]; }
		    
		  }
		  
		  // for(sx = x_start; (sx < (x + merge_x + 1)) && (sx < size_x); ++sx)
		}
		
		// for(sy = y_start; (sy < (y + merge_y + 1)) && (sy < size_y); ++sy)
	      }
	      
	      // for(sz = z_start; (sz < z); ++sz)
	    }
	    for(sy = sy_start; (sy < y); ++sy){
	      
	      for(sx = sx_start; (sx <= (x + merge_x + 1)) && (sx < size_x); ++sx){
		
		// check if this test voxel is source
		if((flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] == flag_value) || (flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] >= 0)){
		  
		  // check if this is a new objid, move on if not
		  prev = -1;
		  for(i = 0; (i < NOi) && (i < init_limit); ++i){
		    
		    if(match_init[i] == flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]){
		      
		      prev = 1;
		      break;
		      
		    }
		    
		  }
		  if(prev == -1){ 
		    
		    // add current object to list of previous matches
		    ++NOi;
		    match_init.push_back(flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]);

		  }
		  
		  // update existing flag
		  if(((flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] < existing) && (flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] >= 0)) || (existing == flag_value)){ existing = flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]; }
		  
		}		
		
		// for(sx = x_start; (sx < = (x + merge_x + 1)) && (sx < size_x); ++sx)
	      }
	      
	      // for(sy = y_start; (sy < y) && (sy < size_y); ++sy)
	    }
	    for(sx = sx_start; sx < x; ++sx){
	      
	      // check if this test voxel is source
	      if((flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))] == flag_value) || (flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))] >= 0)){
		
		// check if this is a new objid, move on if not
		prev = -1;
		for(i = 0; (i < NOi) && (i < init_limit); ++i){
		  
		  if(match_init[i] == flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))]){
		    
		    prev = 1;
		    break;
		    
		  }
		  
		}
		if(prev == -1){ 
		  
		  // add current object to list of previous matches
		  ++NOi;
		  match_init.push_back(flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))]);

		}
		
		// update existing flag
		if(((flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))] < existing) && (flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))] >= 0)) || (existing == flag_value)){ existing = flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))]; }
		
	      }	  
	      
	      // for(sx = x_start; sx <= x; ++sx)
	    }
	  
	    // if(ss_mode == 1)
	  } else {

	    // use a cylinder for the linking lengths

	    for(sz = sz_start; (sz < z); ++sz){
	    
	      for(sy = sy_start; (sy <= (y + merge_y + 1)) && (sy < size_y); ++sy){
	      
		for(sx = sx_start; (sx <= (x + merge_x + 1)) && (sx < size_x); ++sx){
		
		  // check if the test voxel is within the boundaries of the spatial ellipse
		  if(((((float) (x - sx)) * ((float) (x - sx)) / (((float) (merge_x + 1)) * ((float) (merge_x + 1)))) + (((float) (y - sy)) * ((float) (y - sy)) / (((float) (merge_y + 1)) * ((float) (merge_y + 1))))) > 1.0){ continue; }

		  // check if this test voxel is source
		  if((flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] == flag_value) || (flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] >= 0)){

		    // check if this is a new objid
		    prev = -1;
		    for(i = 0; (i < NOi) && (i < init_limit); ++i){
		      
		      if(match_init[i] == flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]){
			
			prev = 1;
			break;
			
		      }
		      
		    }
		    if(prev == -1){ 
		      
		      // add current object to list of previous matches
		      ++NOi;
		      match_init.push_back(flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]);

		    }
		    
		    // update existing flag
		    if(((flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] < existing) && (flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] >= 0)) || (existing == flag_value)){ existing = flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]; }
		    
		  }
		  
		  // for(sx = x_start; (sx < (x + merge_x + 1)) && (sx < size_x); ++sx)
		}
		
		// for(sy = y_start; (sy < (y + merge_y + 1)) && (sy < size_y); ++sy)
	      }
	      
	      // for(sz = z_start; (sz < z); ++sz)
	    }
	    for(sy = sy_start; (sy < y); ++sy){
	      
	      for(sx = sx_start; (sx <= (x + merge_x + 1)) && (sx < size_x); ++sx){
		
		// check if the test voxel is within the boundaries of the spatial ellipse
		if(((((float) (x - sx)) * ((float) (x - sx)) / (((float) (merge_x + 1)) * ((float) (merge_x + 1)))) + (((float) (y - sy)) * ((float) (y - sy)) / (((float) (merge_y + 1)) * ((float) (merge_y + 1))))) > 1.0){ continue; }

		// check if this test voxel is source
		if((flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] == flag_value) || (flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] >= 0)){
		  
		  // check if this is a new objid, move on if not
		  prev = -1;
		  for(i = 0; (i < NOi) && (i < init_limit); ++i){
		    
		    if(match_init[i] == flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]){
		      
		      prev = 1;
		      break;
		      
		    }
		    
		  }
		  if(prev == -1){ 
		    
		    // add current object to list of previous matches
		    ++NOi;
		    match_init.push_back(flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]);

		  }
		  
		  // update existing flag
		  if(((flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] < existing) && (flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] >= 0)) || (existing == flag_value)){ existing = flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]; }
		  
		}		
		
		// for(sx = x_start; (sx < = (x + merge_x + 1)) && (sx < size_x); ++sx)
	      }
	      
	      // for(sy = y_start; (sy < y) && (sy < size_y); ++sy)
	    }
	    for(sx = sx_start; sx < x; ++sx){
	      
	      // check if this test voxel is source
	      if((flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))] == flag_value) || (flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))] >= 0)){
		
		// check if this is a new objid, move on if not
		prev = -1;
		for(i = 0; (i < NOi) && (i < init_limit); ++i){
		  
		  if(match_init[i] == flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))]){
		    
		    prev = 1;
		    break;
		    
		  }
		  
		}
		if(prev == -1){ 
		  
		  // add current object to list of previous matches
		  ++NOi;
		  match_init.push_back(flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))]);
		  
		}
		
		// update existing flag
		if(((flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))] < existing) && (flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))] >= 0)) || (existing == flag_value)){ existing = flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))]; }
		
	      }	  
	      
	      // for(sx = x_start; sx <= x; ++sx)
	    }

	    // else . . . if(ss_mode == 1)
	  }

	  // 3c. assign an object number to this voxel, depending upon the value of the existing flag,
	  // and if this is part of an existing object or objects, then daisy chain from this voxel to all the others
	  if(existing == flag_value){
	    
	    // get obj value from array of obj_ids, then replace it with an incremented obj value
	    
	    // ensure that the obj_ids array is sorted
	    sort(obj_ids.begin(),obj_ids.end());
	    
	    // assign value to array
	    flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (x*data_metric[0]))] = obj_ids[0]; 
	    
	    // calculate the obj_batch value for this obj_id
	    obj_batch = (int) floorf(((float) obj_ids[0] / (float) obj_limit));
	    
	    // add this initial point to the object
	    detections[obj_batch][(obj_ids[0] - (obj_batch * obj_limit))].AddPoint(((float) (chunk_x_start + x)),((float) (chunk_y_start + y)),((float) (chunk_z_start + z)),data_vals[((z*data_metric[2]) + (y*data_metric[1]) + (x*data_metric[0]))]);

	    // push id to list of objects that need to be size thresholded
	    i = -1;
	    for(g = 0; g < NO_check_obj_ids; ++g){ if(check_obj_ids[g] == obj_ids[0]){ i = 1; break; } }
	    if(i == -1){

	      ++NO_check_obj_ids;
	      check_obj_ids.push_back(obj_ids[0]);
	      
	    }

	    // adjust obj and obj_ids
	    if(NO_obj_ids > 1){
	      
	      for(i = 0; i < (NO_obj_ids - 1); ++i){ obj_ids[i] = obj_ids[(i + 1)]; }
	      --NO_obj_ids;
	      obj_ids.resize(NO_obj_ids);

	    } else {
	      
	      ++obj;

	      obj_ids[0] = obj;
	      
	      // if this obj_id modulo obj_limit == 0, then it's the first obj of a new batch
	      // so create a new batch of objects
	      if((obj % obj_limit) == 0){ 
				
		detections.push_back(new object_props[obj_limit]);
				
		// if((obj % obj_limit) == 0)
	      }
			      
	      // if(NO_obj_ids > 1) { } else 
	    }
	    
	    // if(existing == -1)
	  } else if(existing >= 0){
	    
	    // assign the current voxel being processed the correct flag value, and update
	    // the bounding box
	    
	    // add object to list of objects within the merging box, provided that it hasn't been added already
	    i = -1;
	    for(g = 0; g < NO_check_obj_ids; ++g){ if(check_obj_ids[g] == existing){ i = 1; break; } }
	    if(i == -1){

	      ++NO_check_obj_ids;
	      check_obj_ids.push_back(existing);

	    }

	    // assign value to array
	    flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (x*data_metric[0]))] = existing;
	    
	    // calculate the obj_batch value for the existing object
	    obj_batch = (int) floorf(((float) existing / (float) obj_limit));
	    
	    // add this point to the existing object
	    detections[obj_batch][(existing - (obj_batch * obj_limit))].AddPoint((chunk_x_start + x),(chunk_y_start + y),(chunk_z_start + z),data_vals[((z*data_metric[2]) + (y*data_metric[1]) + (x*data_metric[0]))]);
	    
	    // for each object detected within the merging volume, update all of the values within the bounding box of this object
	    // to be the same as the other part of this object, and update the bounding box at the same time
	    if(NOi > 1){
	      
	      for(i = 0; i < NOi; ++i){
		
		// move on if this is the existing object
		if(match_init[i] == existing){ 

		  continue; 

		}
		
		// calculate a temporary obj_batch value
		obj_batch_2 = (int) floorf(((float) match_init[i] / (float) obj_limit));		
		
		// add this object to an existing object
		//detections[obj_batch][(existing - (obj_batch * obj_limit))].AddObject(detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))]);
		detections[obj_batch][(existing - (obj_batch * obj_limit))].AddObject(detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))],temp_sparse_reps_grid,temp_sparse_reps_strings,temp_mom0,temp_RAPV,temp_DECPV,temp_ref_spec,temp_obj_spec,temp_vfield);
		
		// flag the `object' values within the bounding box to merge it with the base object
		sx_start = detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetRAmin() - chunk_x_start;
		if(sx_start < 0){ sx_start = 0; }
		sy_start = detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetDECmin() - chunk_y_start;
		if(sy_start < 0){ sy_start = 0; }
		sz_start = detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetFREQmin() - chunk_z_start;
		if(sz_start < 0){ sz_start = 0; }
		for(sz = sz_start; (sz <= detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetFREQmax() - chunk_z_start) && (sz < size_z); ++sz){
		  
		  for(sy = sy_start; (sy <= detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetDECmax() - chunk_y_start) && (sy < size_y); ++sy){
		    
		    for(sx = sx_start; (sx <= detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetRAmax() - chunk_x_start) && (sx < size_x); ++sx){
		      
		      if(flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] == match_init[i]){ flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] = existing; }
		      
		      // for(sx = x_start[(match_init[i] - (obj_batch_2 * obj_limit))]; sx <= x_finish[(match_init[i] - (obj_batch_2 * obj_limit))]; ++sx)
		    }
		    
		    // for(sy = y_start[(match_init[i] - (obj_batch_2 * obj_limit))]; sy <= y_finish[(match_init[i] - (obj_batch_2 * obj_limit))]; ++sy)
		  }
		  
		  // for(sz = z_start[(match_init[i] - (obj_batch_2 * obj_limit))]; sz <= z_finish[(match_init[i] - (obj_batch_2 * obj_limit))]; ++sz)
		}
		
		// re-initialise the `object' values for the match_init[i] object that has just been merged into the existing object, 
		// and push the obj value (match_init[i]) to the top of the array of values
		detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].ReInit();
		if((detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].Get_srep_update() != 0) && (detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].Get_srep_size(0) >= 0)){
		  
		  detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].ReInit_srep();
		  detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].ReInit_mini();

		}
		detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].ReInit_size();
		detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].Set_srep_update(0);

		j = -1;
		for(g = 0; g < NO_obj_ids; ++g){ if(obj_ids[g] == match_init[i]){ j = 1; break; } }
		if(j == -1){
		  
		  ++NO_obj_ids;
		  obj_ids.push_back(match_init[i]);
		  		  
		}
		
		// for(i = 0; i < NOi; ++i)
	      }
	      
	      // if(NOi > 1)
	    }
	    
	    // else . . . if(existing == -1)
	  }
	  
	  // if(flag_vals[((z*size_x*size_y) + (y*size_x) + x)] == -1)
	}
	
	// update progress on display
	while(progress <= (((float) ((z * size_x * size_y) + (y * size_x) + x + 1)) / ((float) (size_x * size_y * size_z)))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
	
	// for(x = 0; x < size_x; ++x)
      }
      
      // for(y = 0; y < size_y; ++y)
    }
    
    // for(z = 0; z < size_z; ++z)
  }
  std::cout << "* done." << std::endl;
  
  // 4. generate sparse representations of coherent objects that are retained after
  // this `chunk' of the datacube is processed
  std::cout << "Generating/updating sparse representations of sources . . . " << std::endl;
  std::cout << "0 | |:| | : | |:| | 100% complete" << std::endl;
  progress = 0.0;
  for(i = 0; i < obj; ++i){
    
    // calculate the obj_batch value for this object
    obj_batch = floorf(((float) i / (float) obj_limit));
    
    // move on if this object has been re-initialised
    if(detections[obj_batch][(i - (obj_batch * obj_limit))].ShowVoxels() < 1){ 
      
      // update progress on display
      while(progress <= (((float) (i + 1) / ((float) obj)))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }  
      continue; 
      
    }
    
    // move on if this object is not within the boundaries of the current chunk, in other words,
    // move on if this object doesn't need to have it's sparse representation and postage stamp images updated
    if((detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmax() < chunk_x_start) || (detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmax() < chunk_y_start) || (detections[obj_batch][(i - (obj_batch * obj_limit))].GetFREQmax() < chunk_z_start) || (detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmin() >= (chunk_x_start + size_x)) || (detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmin() >= (chunk_y_start + size_y)) || (detections[obj_batch][(i - (obj_batch * obj_limit))].GetFREQmin() >= (chunk_z_start + size_z))){ 

      // update progress on display
      while(progress <= (((float) (i + 1) / ((float) obj)))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }  
      continue; 
      
    }
    
    // move on if this object's sparse representation and postage stamp images don't need updating
    if((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) >= 0) && (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_update() == -1)){
      
      // update progress on display
      while(progress <= (((float) (i + 1) / ((float) obj)))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }  
      continue; 
	
    }
          
    // update the sparse representation for this object
    
    // 0. if there is no existing sparse representation and postage stamp images, create them from scratch, 
    // otherwise concatenate the existing sparse representation and postage stamp images with the updates
    if((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) < 0) && (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_update() >= 0)){
      
      // branch for creating a brand new sparse representation
      
      // update the sparse_reps_update array
      detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_update(-1);

      // initialise the number of object sections that make up this object
      NOi = 0;   
      
      // calculate the edge of the bounding box in `chunk' co-ordinates
      sx_start = detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmin() - chunk_x_start;
      if(sx_start < 0){ sx_start = 0; }
      sx_finish = detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmax() - chunk_x_start;
      sy_start = detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmin() - chunk_y_start;
      if(sy_start < 0){ sy_start = 0; }
      sy_finish = detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmax() - chunk_y_start;
      sz_start = detections[obj_batch][(i - (obj_batch * obj_limit))].GetFREQmin() - chunk_z_start;
      if(sz_start < 0){ sz_start = 0; }
      sz_finish = detections[obj_batch][(i - (obj_batch * obj_limit))].GetFREQmax() - chunk_z_start;
      
      // count the number of object sections that make up this object
      for(sy = sy_start; (sy <= sy_finish) && (sy < size_y); ++sy){
	
	for(sx = sx_start; (sx <= sx_finish) && (sx < size_x); ++sx){
	  
	  // initialise the dummy integer j
	  j = -1;
	  
	  for(sz = sz_start; (sz <= sz_finish) && (sz < size_z); ++sz){
	    
	    // change j to reflect if this voxel in the flag_vals array belongs to the source
	    if(flag_vals[((sz * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))] == i){ 
	      
	      // if previously this was not part of an object string, j == -1, then increment the number of object strings and
	      // adjust the value of j
	      if(j == -1){ ++NOi; }
	      
	      // change j to reflect that this voxel belongs to the object
	      j = 1; 
	      
	    } else { j = -1; }
	    
	    // for(sz = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetFREQmin() - chunk_z_size; sz <= (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetFREQmax() - chunk_z_size); ++sz)
	  }
	  
	  // for(sx = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmin() - chunk_x_size; sx <= (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmax() - chunk_x_size); ++sx)
	}
	
	// for(sy = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmin() - chunk_y_size; sy <= (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmax() - chunk_y_size); ++sy)
      }

      // create an array element to store this sparse representation
      detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_size(0,detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmin());
      detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_size(1,detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmax());
      detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_size(2,detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmin());
      detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_size(3,detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmax());
      detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_size(4,detections[obj_batch][(i - (obj_batch * obj_limit))].GetFREQmin());
      detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_size(5,detections[obj_batch][(i - (obj_batch * obj_limit))].GetFREQmax());
      detections[obj_batch][(i - (obj_batch * obj_limit))].Create_srep_grid((((detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmin() + 1) * (detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmin() + 1)) + 1));
      detections[obj_batch][(i - (obj_batch * obj_limit))].Create_srep_strings((2 * NOi));
      
      // calculate the extent of the sparse_rep volume
      sx_finish = detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmin() + 1;
      sy_finish = detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmin() + 1;
      sz_finish = detections[obj_batch][(i - (obj_batch * obj_limit))].GetFREQmax() - detections[obj_batch][(i - (obj_batch * obj_limit))].GetFREQmin() + 1;
      
      // create arrays to store the postage stamp images and initialise them
      detections[obj_batch][(i - (obj_batch * obj_limit))].Create_mom0((sx_finish * sy_finish));
      detections[obj_batch][(i - (obj_batch * obj_limit))].Create_RAPV((sx_finish * sz_finish));
      detections[obj_batch][(i - (obj_batch * obj_limit))].Create_DECPV((sy_finish * sz_finish));
      detections[obj_batch][(i - (obj_batch * obj_limit))].Create_ospec(sz_finish);
      if(sz_finish >= 10){

	detections[obj_batch][(i - (obj_batch * obj_limit))].Create_rspec((2 * sz_finish));

      } else {

	detections[obj_batch][(i - (obj_batch * obj_limit))].Create_rspec((10 + sz_finish));

      }
      detections[obj_batch][(i - (obj_batch * obj_limit))].Create_vfield((sx_finish * sy_finish));
      
      // initialise the mini_mom0 values
      for(sy = 0; sy < sy_finish; ++sy){
	for(sx = 0; sx < sx_finish; ++sx){
	  detections[obj_batch][(i - (obj_batch * obj_limit))].Set_mom0(((sy * sx_finish) + sx),0.0);
	}
      }

      // initialise the mini_RAPV values
      for(sz = 0; sz < sz_finish; ++sz){
	for(sx = 0; sx < sx_finish; ++sx){
	  detections[obj_batch][(i - (obj_batch * obj_limit))].Set_RAPV(((sz * sx_finish) + sx),0.0);
	}
      }
      
      // initialise the mini_DECPV values
      for(sz = 0; sz < sz_finish; ++sz){
	for(sy = 0; sy < sy_finish; ++sy){
	  detections[obj_batch][(i - (obj_batch * obj_limit))].Set_DECPV(((sz * sy_finish) + sy),0.0);
	}
      }
      
      // initialise the mini_spec values
      for(sz = 0; sz < sz_finish; ++sz){ detections[obj_batch][(i - (obj_batch * obj_limit))].Set_ospec(sz,0.0); }
      if(sz_finish >= 10){

	for(sz = 0; sz < (2 * sz_finish); ++sz){ detections[obj_batch][(i - (obj_batch * obj_limit))].Set_rspec(sz,0.0); }

      } else {

	for(sz = 0; sz < (10 + sz_finish); ++sz){ detections[obj_batch][(i - (obj_batch * obj_limit))].Set_rspec(sz,0.0); }

      }
            
      // initialise the mini_vfield values
      for(sy = 0; sy < sy_finish; ++sy){
	for(sx = 0; sx < sx_finish; ++sx){
	  detections[obj_batch][(i - (obj_batch * obj_limit))].Set_vfield(((sy * sx_finish) + sx),0.0);
	}
      }

      // populate the sparse array representation and postage stamp images
      NOi = 0;
      for(sy = 0; sy < sy_finish; ++sy){
	
	for(sx = 0; sx < sx_finish; ++sx){
	  
	  // write the current number of object strings to the grid
	  detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_grid(((sy * sx_finish) + sx),NOi);
	  
	  // initialise the dummy integer j 
	  j = -1;
	  
	  for(sz = 0; sz < sz_finish; ++sz){
	    
	    // change j to reflect if this voxel in the flag_vals array belongs to the source
	    if(flag_vals[(((sz + sz_start) * data_metric[2]) + ((sy + sy_start) * data_metric[1]) + ((sx + sx_start) * data_metric[0]))] == i){ 
	      
	      // update mini_mom0
	      detections[obj_batch][(i - (obj_batch * obj_limit))].Add_mom0(((sy * sx_finish) + sx),data_vals[(((sz + sz_start) * data_metric[2]) + ((sy + sy_start) * data_metric[1]) + ((sx + sx_start) * data_metric[0]))]);
	      
	      // update mini_RAPV
	      detections[obj_batch][(i - (obj_batch * obj_limit))].Add_RAPV(((sz * sx_finish) + sx),data_vals[(((sz + sz_start) * data_metric[2]) + ((sy + sy_start) * data_metric[1]) + ((sx + sx_start) * data_metric[0]))]);
	      
	      // update mini_DECPV
	      detections[obj_batch][(i - (obj_batch * obj_limit))].Add_DECPV(((sz * sy_finish) + sy),data_vals[(((sz + sz_start) * data_metric[2]) + ((sy + sy_start) * data_metric[1]) + ((sx + sx_start) * data_metric[0]))]);
	      
	      // update mini_obj_spec
	      detections[obj_batch][(i - (obj_batch * obj_limit))].Add_ospec(sz,data_vals[(((sz + sz_start) * data_metric[2]) + ((sy + sy_start) * data_metric[1]) + ((sx + sx_start) * data_metric[0]))]);
  
	      // update mini_vfield
	      detections[obj_batch][(i - (obj_batch * obj_limit))].Add_vfield(((sy * sx_finish) + sx),(((float) (sz + sz_start + chunk_z_start)) * data_vals[(((sz + sz_start) * data_metric[2]) + ((sy + sy_start) * data_metric[1]) + ((sx + sx_start) * data_metric[0]))]));
	      
	      // if previously this was not part of an object string, j == -1, then increment the number of object strings and
	      // write the position to the sparse_reps array
	      if(j == -1){ 
		
		// write the beginning of the object string to sparse_reps
		detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_strings((2 * NOi),(sz_start + chunk_z_start + sz));
		
		// increment NOi
		++NOi; 
		
		// initialise the end of the object string in sparse_reps
		detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_strings(((2 * NOi) - 1),(sz_start + chunk_z_start + sz));

	      } else {
		
		// update the end of the object string in sparse_reps
		detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_strings(((2 * NOi) - 1),(sz_start + chunk_z_start + sz));
		
	      }
	      
	      // change j to reflect that this voxel belongs to the object
	      j = 1; 
	      
	    } else { 
	      
	      // change j to reflect that this voxel doesn't belong to the object
	      j = -1; 
	      
	    }
	    
	    // for(sz = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetFREmin() - chunk_z_start + 1; sz <= (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetFREQmax() - chunk_z_start + 1); ++sz)
	  }
	  
	  // for(sx = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmin() - chunk_x_start + 1; sx <= (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmax() - chunk_x_start + 1); ++sx)
	}
	
	//for(sy = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmin() - chunk_y_start + 1; sy <= (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmax() - chunk_y_start + 1); ++sy)
      }
            
      // update sparse_reps with the final value of NOi
      detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_grid((sx_finish * sy_finish),NOi);
            
      // create the reference spectrum
      if(sz_finish >= 10){
	
	for(sz = -1 * (int) floorf((0.5 * (float) sz_finish)); sz < (sz_finish + ((int) floorf((0.5 * (float) sz_finish)))); ++sz){
	  
	  for(sy = 0; sy < sy_finish; ++sy){
	  
	    for(sx = 0; sx < sx_finish; ++sx){
	    
	      if(((sz + sz_start) >= 0) && ((sz + sz_start) < size_z)){
		
		// update mini_ref_spec	
		detections[obj_batch][(i - (obj_batch * obj_limit))].Add_rspec((sz + ((int) floorf((0.5 * (float) sz_finish)))),data_vals[(((sz + sz_start) * data_metric[2]) + ((sy + sy_start) * data_metric[1]) + ((sx + sx_start) * data_metric[0]))]);	
		
	      }
	      
	    }
	    
	  }
      
	}

	// if(sz_finish >= 10)
      } else {

	for(sz = -5; sz < (sz_finish + 5); ++sz){
	  
	  for(sy = 0; sy < sy_finish; ++sy){
	  
	    for(sx = 0; sx < sx_finish; ++sx){
	    
	      if(((sz + sz_start) >= 0) && ((sz + sz_start) < size_z)){
		
		// update mini_ref_spec
		detections[obj_batch][(i - (obj_batch * obj_limit))].Add_rspec((sz + 5),data_vals[(((sz + sz_start) * data_metric[2]) + ((sy + sy_start) * data_metric[1]) + ((sx + sx_start) * data_metric[0]))]);	
		
	      }
	      
	    }
	    
	  }
      
	}

	// else . . . if(sz_finish >= 10)
      }
      
     // if(sparse_reps[obj_batch][(obj - (obj_batch * obj_limit))][0] < 0) . . . else
    } else {
            
      // branch for updating an existing sparse representation
      
      // count the number of object strings comprising the object
      
      // 1. initialise the number of object strings
      NOi = 0;
      
      // 2. convert the new bounding box of the object into `chunk' co-ordinates
      sx_start = detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmin() - chunk_x_start;
      sx_finish = detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmax() - chunk_x_start;
      sy_start = detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmin() - chunk_y_start;
      sy_finish = detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmax() - chunk_y_start;
      sz_start = detections[obj_batch][(i - (obj_batch * obj_limit))].GetFREQmin() - chunk_z_start;
      if(sz_start < 0){ sz_start = 0; }
      sz_finish = detections[obj_batch][(i - (obj_batch * obj_limit))].GetFREQmax() - chunk_z_start;
      
      // 3. search through the bounding box, and do one of the following
      // 3.a if inside the chunk and outside the existing sparse rep grid, run the routine to count the number
      //     of object strings along the current line of sight
      // 3.b otherwise use the grid value to increment NOi
      for(sy = sy_start; sy <= sy_finish; ++sy){
	
	for(sx = sx_start; sx <= sx_finish; ++sx){
	  	  
	  // case 3.a triggered
	  if(((sx > merge_x) || (chunk_x_start == 0)) && ((sy > merge_y) || (chunk_y_start == 0))){
	    	    
	    // initialise the dummy integer j 
	    j = -1;
	    
	    for(sz = sz_start; sz <= sz_finish; ++sz){
	      
	      // change i to reflect if this voxel in the flag_vals array belongs to the source
	      if(flag_vals[((sz * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))] == i){ 
		
		// if previously this was not part of an object string, i == -1, then increment the number of object strings and
		// write the position to the sparse_reps array
		if(j == -1){ 
		  
		  // increment NOi
		  ++NOi; 
		  
		} 
		
		// change j to reflect that this voxel belongs to the object
		j = 1; 
		
	      } else { j = -1; }
	      
	      // for(sz = sz_start; sz <= sz_finish; ++sz)
	    }
	    
	  } else {
	    
	    // case 3.b triggered  
	    NOi+=(detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_grid((((sy + chunk_y_start - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx + chunk_x_start - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_grid((((sy + chunk_y_start - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx + chunk_x_start - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0))));
	    
	  }     
	  
	  // for(sx = sx_start; sx <= sx_finish; ++sx)
	}
	
	// for(sy = sy_start; sy <= sy_finish; ++sy)
      }
      
      // 4. allocate memory to store temporary sparse representation and postage stamp images, then
      //    initialise the postage stamp image arrays
      temp_mom0.resize(0);
      for(sy = 0; sy < (sy_finish - sy_start + 1); ++sy){
	for(sx = 0; sx < (sx_finish - sx_start + 1); ++sx){
	  temp_mom0.push_back(0.0);
	}
      }
      temp_RAPV.resize(0);
      for(sz = 0; sz < (sz_finish - sz_start + 1); ++sz){
	for(sx = 0; sx < (sx_finish - sx_start + 1); ++sx){
	  temp_RAPV.push_back(0.0);
	}
      }
      temp_DECPV.resize(0);
      for(sz = 0; sz < (sz_finish - sz_start + 1); ++sz){
	for(sy = 0; sy < (sy_finish - sy_start + 1); ++sy){
	  temp_DECPV.push_back(0.0);
	}
      }
      temp_obj_spec.resize(0);
      for(sz = 0; sz < (sz_finish - sz_start + 1); ++sz){ temp_obj_spec.push_back(0.0); }
      if((sz_finish - sz_start + 1) >= 10){

	temp_ref_spec.resize(0);
	for(sz = 0; sz < (2 * (sz_finish - sz_start + 1)); ++sz){ temp_ref_spec.push_back(0.0); }

      } else {

	temp_ref_spec.resize(0);
	for(sz = 0; sz < (sz_finish - sz_start + 11); ++sz){ temp_ref_spec.push_back(0.0); }

      }
      temp_vfield.resize(0);
      for(sy = 0; sy < (sy_finish - sy_start + 1); ++sy){
	for(sx = 0; sx < (sx_finish - sx_start + 1); ++sx){
	  temp_vfield.push_back(0.0);
	}
      }
     
      // 5. populate temp_sparse_rep arrays
      // 5.a if inside the chunk and outside the existing sparse rep grid, run the routine to count the number
      //     of object strings along the current line of sight, populate new strings array and add values to new 
      //     postage stamp images
      // 5.b if outside of the chunk, use the grid value to increment NOi, populate strings array from previous version and 
      //     fill in valuese for postage stamp images

      // initialise number of object strings
      NOi = 0;
      
      // initialise RA PV image
      for(sz = 0; sz < (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1); ++sz){

	for(sx = 0; sx < (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1); ++sx){
	
	  temp_RAPV[(((sz + detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) - sz_start - chunk_z_start) * (sx_finish - sx_start + 1)) + sx + detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) - sx_start - chunk_x_start)]+=detections[obj_batch][(i - (obj_batch * obj_limit))].Get_RAPV(((sz * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx));
	  
	}

      }

      // initialise Dec PV image
      for(sz = 0; sz < (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1); ++sz){

	for(sy = 0; sy < (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2) + 1); ++sy){
	
	  temp_DECPV[(((sz + detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) - sz_start - chunk_z_start) * (sy_finish - sy_start + 1)) + sy + detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2) - sy_start - chunk_y_start)]+=detections[obj_batch][(i - (obj_batch * obj_limit))].Get_DECPV(((sz * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2) + 1)) + sy));
	  
	}

      }

      // initialise mini spec's
      for(sz = 0; sz < (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1); ++sz){ temp_obj_spec[(sz + detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) - sz_start - chunk_z_start)]+=detections[obj_batch][(i - (obj_batch * obj_limit))].Get_ospec(sz); }

      if((sz_finish - sz_start + 1) >= 10){

	j = sz_start + chunk_z_start - ((int) floorf((0.5 * (float) (sz_finish - sz_start + 1))));

	if((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1) >= 10){

	  k = detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) - ((int) floorf((0.5 * (float) (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1))));

	  for(sz = 0; sz < (2 * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); ++sz){

	    temp_ref_spec[(sz + k - j)]+=detections[obj_batch][(i - (obj_batch * obj_limit))].Get_rspec(sz);

	  }

	} else {
	  
	  k = detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) - 5;

	  for(sz = 0; sz < (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 11); ++sz){
	    
	    temp_ref_spec[(sz + k - j)]+=detections[obj_batch][(i - (obj_batch * obj_limit))].Get_rspec(sz);

	  }
	  
	}

      } else {

	j = sz_start + chunk_z_start - 5;

	if((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1) >= 10){

	  k = detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) - ((int) floorf((0.5 * (float) (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1))));

	  for(sz = 0; sz < (2 * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); ++sz){

	    temp_ref_spec[(sz + k - j)]+=detections[obj_batch][(i - (obj_batch * obj_limit))].Get_rspec(sz);

	  }

	} else {
	  
	  k = detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) - 5;

	  for(sz = 0; sz < (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 11); ++sz){
	    
	    temp_ref_spec[(sz + k - j)]+=detections[obj_batch][(i - (obj_batch * obj_limit))].Get_rspec(sz);

	  }
	  
	}

      }
      
      // search through LoS's of updated bounding box and take action as dictated by case 5.a or 5.b
      for(sy = sy_start; sy <= sy_finish; ++sy){
	
	for(sx = sx_start; sx <= sx_finish; ++sx){
	  	  
	  // add value to temp_sparse_rep_grid
	  temp_sparse_reps_grid[(((sy - sy_start) * (sx_finish - sx_start + 1)) + sx - sx_start)] = NOi;
	  
	  // case 5.a triggered for moment-0 image and sparse representation
	  if(((sx > merge_x) || (chunk_x_start == 0)) && ((sy > merge_y) || (chunk_y_start == 0))){
	    
	    // initialise the dummy integer j 
	    j = -1;
	    
	    for(sz = sz_start; sz <= sz_finish; ++sz){
	      
	      // change j to reflect if this voxel in the flag_vals array belongs to the source
	      if(flag_vals[((sz * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))] == i){ 
		
		// update mini_mom0
		temp_mom0[(((sy - sy_start) * (sx_finish - sx_start + 1)) + sx - sx_start)]+=data_vals[((sz * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))];
		
		// update mini_RAPV
		temp_RAPV[(((sz - sz_start) * (sx_finish - sx_start + 1)) + sx - sx_start)]+=data_vals[((sz * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))];
		
		// update mini_DECPV
		temp_DECPV[(((sz - sz_start) * (sy_finish - sy_start + 1)) + sy - sy_start)]+=data_vals[((sz * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))];
		
		// update mini_obj_spec
		temp_obj_spec[(sz - sz_start)]+=data_vals[((sz * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))];
		
		// update mini_vfield
		temp_vfield[(((sy - sy_start) * (sx_finish - sx_start + 1)) + sx - sx_start)]+=(((float) (sz + chunk_z_start)) * data_vals[((sz * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))]);
		
		// if previously this was not part of an object string, j == -1, then increment the number of object strings and
		// write the position to the sparse_reps array
		if(j == -1){ 
		  
		  // write the beginning of the object string to sparse_reps
		  temp_sparse_reps_strings[(2 * NOi)] = chunk_z_start + sz;
		  
		  // increment NOi
		  ++NOi; 
		  		  
		  // initialise the end of the object string in sparse_reps
		  temp_sparse_reps_strings[((2 * NOi) - 1)] = chunk_z_start + sz;
		  
		} else {
		  
		  // update the end of the object string in sparse_reps
		  temp_sparse_reps_strings[((2 * NOi) - 1)] = chunk_z_start + sz;
		  
		}
		
		// change j to reflect that this voxel belongs to the object
		j = 1; 
		
	      } else { 
		
		// change j to reflect that this voxel doesn't belong to the object
		j = -1; 
		
	      }
	      
	      // for(sz = sz_start; sz <= sz_finish; ++sz)
	    }
	    
	    // add values to temp_ref_spec
	    if((sz_finish - sz_start + 1) >= 10){
	      
	      k = (int) floorf((0.5 * (float) (sz_finish - sz_start + 1)));
	      
	      for(sz = 0; sz < (2 * (sz_finish - sz_start + 1)); ++sz){
	  	
		if(((sz + sz_start - k) >= 0) && ((sz + sz_start - k) < size_z)){
		  
		  // update mini_ref_spec
		  temp_ref_spec[(sz + sz_start - k)]+=data_vals[(((sz + sz_start - k) * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))];	
		  
		}
		
	      }
	      
	      // if(sz_finish >= 10)
	    } else {
	      
	      for(sz = 0; sz < (sz_finish - sz_start + 11); ++sz){
		
		if(((sz + sz_start - 5) >= 0) && ((sz + sz_start - 5) < size_z)){
		  
		  // update mini_ref_spec
		  temp_ref_spec[(sz + sz_start - 5)]+=data_vals[(((sz + sz_start - 5) * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))];	
		  
		}
		
	      }
	      
	      // else . . . if(sz_finish >= 10)
	    }
	    
	  } else {
	    // case 5.b triggered  
	    	    
	    // update temp_sparse_rep_strings
	    j = 0;
	    if(((sx + chunk_x_start) >= detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0)) && ((sx + chunk_x_start) <= (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) + 1)) && ((sy + chunk_y_start) >= detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2)) && ((sy + chunk_y_start) <= (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(3) + 1))){
	      
	      for(g = detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_grid((((sy + chunk_y_start - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx + chunk_x_start - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0))); g < detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_grid((((sy + chunk_y_start - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx + chunk_x_start - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1)); ++g){
		
		temp_sparse_reps_strings[(NOi + j)] = detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_strings((2 * g));
		temp_sparse_reps_strings[(NOi + j + 1)] = detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_strings(((2 * g) + 1));
		++j;
		
	      }
	      
	      // update moment 0 postage stamp image
	      temp_mom0[(((sy - sy_start) * (detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmin() + 1)) + sx - sx_start)] = detections[obj_batch][(i - (obj_batch * obj_limit))].Get_mom0((((sy + chunk_y_start - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx + chunk_x_start - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0)));
	      
	      // update vfield postage stamp image
	      temp_vfield[(((sy - sy_start) * (detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmin() + 1)) + sx - sx_start)] = detections[obj_batch][(i - (obj_batch * obj_limit))].Get_vfield((((sy + chunk_y_start - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx + chunk_x_start - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0)));

	      // update NOi
	      NOi+=(detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_grid((((sy + chunk_y_start - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx + chunk_x_start - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_grid((((sy + chunk_y_start - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx + chunk_x_start - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0))));
	      	      
	    }
	    
	    // else . . . case 5.b
	  }
	    
	  // for(sx = sx_start; sx <= sx_finish; ++sx)
	}
		
	// for(sy = sy_start; sy <= sy_finish; ++sy)
      }
            
      temp_sparse_reps_grid[((sx_finish - sx_start + 1) * (sy_finish - sy_start + 1))] = NOi;

      // 6. replace the current sparse representation with the new sparse representation
      detections[obj_batch][(i - (obj_batch * obj_limit))].ReInit_srep();
      detections[obj_batch][(i - (obj_batch * obj_limit))].ReInit_mini();
 
      // update the sparse_reps_size array with the new bounding box
      detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_size(0,detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmin());
      detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_size(1,detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmax());
      detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_size(2,detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmin());
      detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_size(3,detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmax());
      detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_size(4,detections[obj_batch][(i - (obj_batch * obj_limit))].GetFREQmin());
      detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_size(5,detections[obj_batch][(i - (obj_batch * obj_limit))].GetFREQmax());

      // update the sparse_reps_update array
      detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_update(-1);
          
      // populate the new sparse_reps and postage stamps
      detections[obj_batch][(i - (obj_batch * obj_limit))].Create_srep_grid((1 + ((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2) + 1))));
      detections[obj_batch][(i - (obj_batch * obj_limit))].Create_srep_strings((2 * NOi));
      detections[obj_batch][(i - (obj_batch * obj_limit))].Create_mom0(((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2) + 1)));
      detections[obj_batch][(i - (obj_batch * obj_limit))].Create_RAPV(((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1)));
      detections[obj_batch][(i - (obj_batch * obj_limit))].Create_DECPV(((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2) + 1) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1)));
      detections[obj_batch][(i - (obj_batch * obj_limit))].Create_vfield(((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2) + 1)));
  
      detections[obj_batch][(i - (obj_batch * obj_limit))].Create_ospec((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1));
      if((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1) >= 10){
       
	detections[obj_batch][(i - (obj_batch * obj_limit))].Create_rspec((2 * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1)));
	
      } else {
	
	detections[obj_batch][(i - (obj_batch * obj_limit))].Create_rspec((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 11));

      }
   
      for(j = 0; j < (1 + ((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2) + 1))); ++j){ detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_grid(j,temp_sparse_reps_grid[j]); }
      for(j = 0; j < (2 * NOi); ++j){ detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_strings(j,temp_sparse_reps_strings[j]); }
      for(j = 0; j < ((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); ++j){ detections[obj_batch][(i - (obj_batch * obj_limit))].Set_mom0(j,temp_mom0[j]); }
      for(j = 0; j < ((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); ++j){ detections[obj_batch][(i - (obj_batch * obj_limit))].Set_RAPV(j,temp_RAPV[j]); }
      for(j = 0; j < ((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2) + 1) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); ++j){ detections[obj_batch][(i - (obj_batch * obj_limit))].Set_DECPV(j,temp_DECPV[j]); }
      for(j = 0; j < ((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); ++j){ detections[obj_batch][(i - (obj_batch * obj_limit))].Set_vfield(j,temp_vfield[j]); }

      for(j = 0; j < (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1); ++j){ detections[obj_batch][(i - (obj_batch * obj_limit))].Set_ospec(j,temp_obj_spec[j]); }
      if((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1) >= 10){

	for(j = 0; j < (2 * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); ++j){ detections[obj_batch][(i - (obj_batch * obj_limit))].Set_rspec(j,temp_ref_spec[j]); }
	
      } else {

	for(j = 0; j < (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 11); ++j){ detections[obj_batch][(i - (obj_batch * obj_limit))].Set_rspec(j,temp_ref_spec[j]); }

      }
      
      // else . . . if(sparse_reps[obj_batch][(obj - (obj_batch * obj_limit))][0] < 0)
    }
    
    // update progress on display
    while(progress <= (((float) (i + 1) / ((float) obj)))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }  
   
    // for(i = 0; i < obj; ++i)
  }
  std::cout << "* done." << std::endl;
    
  // 5. return the number of coherent objects 
  return obj;
  
}

long int CreateObjects(float * data_vals, long int * flag_vals, int size_x, int size_y, int size_z, int chunk_x_start, int chunk_y_start, int chunk_z_start, int merge_x, int merge_y, int merge_z, int min_x_size, int min_y_size, int min_z_size, int min_v_size, float intens_thresh_min, float intens_thresh_max, long int flag_value, long int start_obj, vector<object_props *> & detections, vector<long int> & obj_ids, vector<long int> & check_obj_ids, int obj_limit, int max_x_val, int max_y_val, int max_z_val, int ss_mode, int * data_metric, int * xyz_order){
 
  long int o,obj,existing,obj_batch,obj_batch_2;
  int x,y,z,sx,sy,sz,sx_start,sy_start,sz_start,sx_finish,sy_finish,sz_finish;
  int g,i,j,k,NO_obj_ids,NO_check_obj_ids,NOi;
  int init_limit,prev,x_start,y_start,temp_x[3],temp_y[3],temp_z[3];
  float progress;
  vector<float> temp_mom0, temp_RAPV, temp_DECPV, temp_obj_spec, temp_ref_spec, temp_vfield;
  vector<long int> match_init;
  vector<int> temp_sparse_reps_grid, temp_sparse_reps_strings;

  // convert the merging lengths from distances (in pixels/voxels) to empty space between linked detections (in pixels/voxels) --- subtract 1 --- this is a more naturual
  // definition when processing a large file in chunks, because it makes it easier to build up sources that need to be retained for the rest of the input file
  --merge_x;
  --merge_y;
  --merge_z;

  // reorder the datacube and subcube limits to be in x,y,z order
  temp_x[0] = chunk_x_start;
  temp_y[0] = chunk_y_start;
  temp_z[0] = chunk_z_start;
  temp_x[1] = size_x;
  temp_y[1] = size_y;
  temp_z[1] = size_z;
  temp_x[2] = max_x_val;
  temp_y[2] = max_y_val;
  temp_z[2] = max_z_val;
  switch(xyz_order[0]){
  case 1:
    chunk_x_start = temp_x[0];
    size_x = temp_x[1];
    max_x_val = temp_x[2];
    break;
  case 2:
    chunk_x_start = temp_y[0];
    size_x = temp_y[1];
    max_x_val = temp_y[2];
    break;
  case 3:
    chunk_x_start = temp_z[0];
    size_x = temp_z[1];
    max_x_val = temp_z[2];
    break;
  default:
    chunk_x_start = temp_x[0];
    size_x = temp_x[1];
    max_x_val = temp_x[2];
    break;
  }
  switch(xyz_order[1]){
  case 1:
    chunk_y_start = temp_x[0];
    size_y = temp_x[1];
    max_y_val = temp_x[2];
    break;
  case 2:
    chunk_y_start = temp_y[0];
    size_y = temp_y[1];
    max_y_val = temp_y[2];
    break;
  case 3:
    chunk_y_start = temp_z[0];
    size_y = temp_z[1];
    max_y_val = temp_z[2];
    break;
  default:
    chunk_y_start = temp_y[0];
    size_y = temp_y[1];
    max_y_val = temp_y[2];
    break;
  }
  switch(xyz_order[2]){
  case 1:
    chunk_z_start = temp_x[0];
    size_z = temp_x[1];
    max_z_val = temp_x[2];
    break;
  case 2:
    chunk_z_start = temp_y[0];
    size_z = temp_y[1];
    max_z_val = temp_y[2];
    break;
  case 3:
    chunk_z_start = temp_z[0];
    size_z = temp_z[1];
    max_z_val = temp_z[2];
    break;
  default:
    chunk_z_start = temp_z[0];
    size_z = temp_z[1];
    max_z_val = temp_z[2];
    break;
  }

  // create temporary arrays
  if(10000 > (((2 * merge_x) + 3) * ((2 * merge_y) + 3) * ((2 * merge_z) + 3))){
    
    init_limit = 10000;
    
  } else {
    
    init_limit = ((2 * merge_x) + 3) * ((2 * merge_y) + 3) * ((2 * merge_z) + 3) + 10;
    
  }
  match_init.reserve(init_limit);
  match_init.resize(0);

  temp_sparse_reps_grid.reserve(1000000);
  temp_sparse_reps_grid.resize(0);
  temp_sparse_reps_strings.reserve(1000000);
  temp_sparse_reps_strings.resize(0);
  temp_mom0.reserve(1000000);
  temp_mom0.resize(0);
  temp_RAPV.reserve(1000000);
  temp_RAPV.resize(0);
  temp_DECPV.reserve(1000000);
  temp_DECPV.resize(0);
  temp_obj_spec.reserve(1000000);
  temp_obj_spec.resize(0);
  temp_ref_spec.reserve(1000000);
  temp_ref_spec.resize(0);
  temp_vfield.reserve(1000000);
  temp_vfield.resize(0);
 
  // 0. initialise variables and arrays
  obj = start_obj;
  x_start = 0;
  y_start = 0;
  if(chunk_x_start > 0){ x_start = merge_x + 1; }
  if(chunk_y_start > 0){ y_start = merge_y + 1; }
  NO_obj_ids = obj_ids.size();
  NO_check_obj_ids = check_obj_ids.size();

  // 1. Create list of `coherent' objects from neighbouring voxels
  
  // for each grid point, check if it has been flagged as a source voxel, and assign
  // it an object ID if it is. If a neighbouring voxel has already been flagged, assign 
  // that ID to the grid point, otherwise, assign current value of obj and increment obj.
  std::cout << "0 | |:| | : | |:| | 100% complete" << std::endl;
  progress = 0.0;
  for(z = 0; ((z < size_z) && ((z + chunk_z_start) < max_z_val)); ++z){

    // if the z value is sufficiently large that objects have started to pop out of the merging box,
    // then check if the objects outside of the merging box are sufficiently large
    if(z > (merge_z + 1)){
      
      // initialise NOi to store the number of current check_obj_ids values that haven't moved outside of 
      // the merging volume, and that form the basis for the new check_obj_ids list
      NOi = 0;
      
      // apply size threshold to objects added since last size thresholding
      for(i = 0; i < NO_check_obj_ids; ++i){
	
	// calculate the obj_batch value for the existing object
	obj_batch = (long int) floor(((double) check_obj_ids[i] / (double) obj_limit));
	
	// move to the next object if this one has been re-initialised
	if(detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	
	// check if this object has moved outside of the merging box available to the next plane, if it 
	// has apply a size threshold, otherwise use it to seed the new list of objects to be size thresholded
	// the next time
	if(detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmax() < (chunk_z_start + z - merge_z - 1)){
	  
	  // apply the size threshold, and if it fails re-initialise the object and pop its id to the list
	  // of available obj_ids
	  if((((detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetRAmin() + 1) < min_x_size) || ((detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetDECmin() + 1) < min_y_size) || ((detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmax() - detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmin() + 1) < min_z_size) || (detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].ShowVoxels() < min_v_size)) && (((detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetRAmin() - chunk_x_start) > merge_x) || (chunk_x_start == 0)) && ((((chunk_x_start + size_x - 1 - detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetRAmax()) > merge_x) || ((chunk_x_start + size_x) >= max_x_val)) && (((detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetDECmin() - chunk_y_start) > merge_y) || (chunk_y_start == 0)) && (((chunk_y_start + size_y - 1 - detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetDECmax()) > merge_y) || ((chunk_y_start + size_y) >= max_y_val)) && (((detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmin() - chunk_z_start) > merge_z) || (chunk_z_start == 0)) && (((chunk_z_start + size_z - 1 - detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmax()) > merge_z) || ((chunk_z_start + size_z) >= max_z_val)))){
	    
	    // re-initialise object
	    detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].ReInit();
	    if((detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].Get_srep_update() != 0) && (detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].Get_srep_size(0) >= 0)){
	      
	      detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].ReInit_srep();
	      detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].ReInit_mini();

	    }
	    detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].ReInit_size();
	    detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].Set_srep_update(0);

	    // flag the `object' values within the bounding box and remove it from the flag_vals array
	    sx_start = detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetRAmin() - chunk_x_start;
	    if(sx_start < 0){ sx_start = 0; }
	    sy_start = detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetDECmin() - chunk_y_start;
	    if(sy_start < 0){ sy_start = 0; }
	    sz_start = detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmin() - chunk_z_start;
	    if(sz_start < 0){ sz_start = 0; }
	    for(sz = sz_start; ((sz <= (detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmax() - chunk_z_start)) && (sz < size_z)); ++sz){
	      
	      for(sy = sy_start; ((sy <= (detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetDECmax() - chunk_y_start)) && (sy < size_y)); ++sy){
		
		for(sx = sx_start; ((sx <= (detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetRAmax() - chunk_x_start)) && (sx < size_x)); ++sx){
		  
		  if(flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx * data_metric[0]))] == check_obj_ids[i]){ flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx * data_metric[0]))] = -99; }
		  
		  // for(sx = x_start[(match_init[i] - (obj_batch_2 * obj_limit))]; sx <= x_finish[(match_init[i] - (obj_batch_2 * obj_limit))]; ++sx)
		}
		
		// for(sy = y_start[(match_init[i] - (obj_batch_2 * obj_limit))]; sy <= y_finish[(match_init[i] - (obj_batch_2 * obj_limit))]; ++sy)
	      }
	      
	      // for(sz = z_start[(match_init[i] - (obj_batch_2 * obj_limit))]; sz <= z_finish[(match_init[i] - (obj_batch_2 * obj_limit))]; ++sz)
	    }
	    
	    // add object id to list of available ids
	    j = -1;
	    for(g = 0; g < NO_obj_ids; ++g){ if(obj_ids[g] == check_obj_ids[i]){ j = 1; break; } }
	    if(j == -1){
	      ++NO_obj_ids;
	      obj_ids.push_back(check_obj_ids[i]);
	    }

	  }
	  
	  // if(z_finish[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))] < (z - merge_z - 1))
	} else {
	  
	  check_obj_ids[NOi] = check_obj_ids[i];
	  ++NOi;
	  
	  // else . . . if(z_finish[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))] < (z - merge_z - 1))
	}
	
	// for(i = 0; i < NO_check_obj_ids; ++i)
      }

      // update NO_check_obj_ids
      NO_check_obj_ids = NOi;
      check_obj_ids.resize(NOi);
      
      // if(z > (merge_z + 2))
    }
    
    for(y = y_start; (y < size_y) && ((chunk_y_start + y) < max_y_val); ++y){
      
      for(x = x_start; (x < size_x) && ((chunk_x_start + x) < max_x_val); ++x){
	
	// if this is a source voxel, check if it is associated with a previously identified
	// source voxel
	if(flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (x*data_metric[0]))] == flag_value){ 
	  
	  // 1.  initialise variables
	  existing = flag_value;
	  NOi = 0;
	  match_init.resize(0);
	  
	  // 3. initial pass, check the voxels within the merging distance to see if any have been
	  // identified as previous objects
	  
	  // 3a. set limits of region to be searched
	  sx_start = x - merge_x - 1;
	  if(sx_start < 0){ sx_start = 0; }
	  sy_start = y - merge_y - 1;
	  if(sy_start < 0){ sy_start = 0; }
	  sz_start = z - merge_z - 1;
	  if(sz_start < 0){ sz_start = 0; }
	  
	  // 3b. search through previous region and modify existing flag if required

	  if(ss_mode == 1){
	    
	    // use a rectangular paralleliped to do the linking

	    for(sz = sz_start; (sz < z); ++sz){
	    
	      for(sy = sy_start; (sy <= (y + merge_y + 1)) && (sy < size_y); ++sy){
	      
		for(sx = sx_start; (sx <= (x + merge_x + 1)) && (sx < size_x); ++sx){
		
		  // check if this test voxel is source
		  if((flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] == flag_value) || (flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] >= 0)){
		    
		    // check if this is a new objid
		    prev = -1;
		    for(i = 0; (i < NOi) && (i < init_limit); ++i){
		      
		      if(match_init[i] == flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]){
			
			prev = 1;
			break;
			
		      }
		      
		    }
		    if(prev == -1){ 
		      
		      // add current object to list of previous matches
		      match_init.push_back(flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]);
		      ++NOi;

		    }
		    
		    // update existing flag
		    if(((flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] < existing) && (flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] >= 0)) || (existing == flag_value)){ existing = flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]; }
		    
		  }
		  
		  // for(sx = x_start; (sx < (x + merge_x + 1)) && (sx < size_x); ++sx)
		}
		
		// for(sy = y_start; (sy < (y + merge_y + 1)) && (sy < size_y); ++sy)
	      }
	      
	      // for(sz = z_start; (sz < z); ++sz)
	    }
	    for(sy = sy_start; (sy < y); ++sy){
	      
	      for(sx = sx_start; (sx <= (x + merge_x + 1)) && (sx < size_x); ++sx){
		
		// check if this test voxel is source
		if((flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] == flag_value) || (flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] >= 0)){
		  
		  // check if this is a new objid, move on if not
		  prev = -1;
		  for(i = 0; (i < NOi) && (i < init_limit); ++i){
		    
		    if(match_init[i] == flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]){
		      
		      prev = 1;
		      break;
		      
		    }
		    
		  }
		  if(prev == -1){ 
		    
		    // add current object to list of previous matches
		    ++NOi;
		    match_init.push_back(flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]);

		  }
		  
		  // update existing flag
		  if(((flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] < existing) && (flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] >= 0)) || (existing == flag_value)){ existing = flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]; }
		  
		}		
		
		// for(sx = x_start; (sx < = (x + merge_x + 1)) && (sx < size_x); ++sx)
	      }
	      
	      // for(sy = y_start; (sy < y) && (sy < size_y); ++sy)
	    }
	    for(sx = sx_start; sx < x; ++sx){
	      
	      // check if this test voxel is source
	      if((flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))] == flag_value) || (flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))] >= 0)){
		
		// check if this is a new objid, move on if not
		prev = -1;
		for(i = 0; (i < NOi) && (i < init_limit); ++i){
		  
		  if(match_init[i] == flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))]){
		    
		    prev = 1;
		    break;
		    
		  }
		  
		}
		if(prev == -1){ 
		  
		  // add current object to list of previous matches
		  ++NOi;
		  match_init.push_back(flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))]);

		}
		
		// update existing flag
		if(((flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))] < existing) && (flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))] >= 0)) || (existing == flag_value)){ existing = flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))]; }
		
	      }	  
	      
	      // for(sx = x_start; sx <= x; ++sx)
	    }
	  
	    // if(ss_mode == 1)
	  } else {

	    // use a cylinder for the linking lengths

	    for(sz = sz_start; (sz < z); ++sz){
	    
	      for(sy = sy_start; (sy <= (y + merge_y + 1)) && (sy < size_y); ++sy){
	      
		for(sx = sx_start; (sx <= (x + merge_x + 1)) && (sx < size_x); ++sx){
		
		  // check if the test voxel is within the boundaries of the spatial ellipse
		  if(((((float) (x - sx)) * ((float) (x - sx)) / (((float) (merge_x + 1)) * ((float) (merge_x + 1)))) + (((float) (y - sy)) * ((float) (y - sy)) / (((float) (merge_y + 1)) * ((float) (merge_y + 1))))) > 1.0){ continue; }

		  // check if this test voxel is source
		  if((flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] == flag_value) || (flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] >= 0)){

		    // check if this is a new objid
		    prev = -1;
		    for(i = 0; (i < NOi) && (i < init_limit); ++i){
		      
		      if(match_init[i] == flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]){
			
			prev = 1;
			break;
			
		      }
		      
		    }
		    if(prev == -1){ 
		      
		      // add current object to list of previous matches
		      ++NOi;
		      match_init.push_back(flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]);

		    }
		    
		    // update existing flag
		    if(((flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] < existing) && (flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] >= 0)) || (existing == flag_value)){ existing = flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]; }
		    
		  }
		  
		  // for(sx = x_start; (sx < (x + merge_x + 1)) && (sx < size_x); ++sx)
		}
		
		// for(sy = y_start; (sy < (y + merge_y + 1)) && (sy < size_y); ++sy)
	      }
	      
	      // for(sz = z_start; (sz < z); ++sz)
	    }
	    for(sy = sy_start; (sy < y); ++sy){
	      
	      for(sx = sx_start; (sx <= (x + merge_x + 1)) && (sx < size_x); ++sx){
		
		// check if the test voxel is within the boundaries of the spatial ellipse
		if(((((float) (x - sx)) * ((float) (x - sx)) / (((float) (merge_x + 1)) * ((float) (merge_x + 1)))) + (((float) (y - sy)) * ((float) (y - sy)) / (((float) (merge_y + 1)) * ((float) (merge_y + 1))))) > 1.0){ continue; }

		// check if this test voxel is source
		if((flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] == flag_value) || (flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] >= 0)){
		  
		  // check if this is a new objid, move on if not
		  prev = -1;
		  for(i = 0; (i < NOi) && (i < init_limit); ++i){
		    
		    if(match_init[i] == flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]){
		      
		      prev = 1;
		      break;
		      
		    }
		    
		  }
		  if(prev == -1){ 
		    
		    // add current object to list of previous matches
		    ++NOi;
		    match_init.push_back(flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]);

		  }
		  
		  // update existing flag
		  if(((flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] < existing) && (flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] >= 0)) || (existing == flag_value)){ existing = flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]; }
		  
		}		
		
		// for(sx = x_start; (sx < = (x + merge_x + 1)) && (sx < size_x); ++sx)
	      }
	      
	      // for(sy = y_start; (sy < y) && (sy < size_y); ++sy)
	    }
	    for(sx = sx_start; sx < x; ++sx){
	      
	      // check if this test voxel is source
	      if((flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))] == flag_value) || (flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))] >= 0)){
		
		// check if this is a new objid, move on if not
		prev = -1;
		for(i = 0; (i < NOi) && (i < init_limit); ++i){
		  
		  if(match_init[i] == flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))]){
		    
		    prev = 1;
		    break;
		    
		  }
		  
		}
		if(prev == -1){ 
		  
		  // add current object to list of previous matches
		  ++NOi;
		  match_init.push_back(flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))]);
		  
		}
		
		// update existing flag
		if(((flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))] < existing) && (flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))] >= 0)) || (existing == flag_value)){ existing = flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))]; }
		
	      }	  
	      
	      // for(sx = x_start; sx <= x; ++sx)
	    }

	    // else . . . if(ss_mode == 1)
	  }

	  // 3c. assign an object number to this voxel, depending upon the value of the existing flag,
	  // and if this is part of an existing object or objects, then daisy chain from this voxel to all the others
	  if(existing == flag_value){
	    
	    // get obj value from array of obj_ids, then replace it with an incremented obj value
	    
	    // ensure that the obj_ids array is sorted
	    sort(obj_ids.begin(),obj_ids.end());
	    
	    // assign value to array
	    flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (x*data_metric[0]))] = obj_ids[0]; 
	    
	    // calculate the obj_batch value for this obj_id
	    obj_batch = (long int) floor(((double) obj_ids[0] / (double) obj_limit));
	    
	    // add this initial point to the object
	    detections[obj_batch][(obj_ids[0] - (obj_batch * obj_limit))].AddPoint(((float) (chunk_x_start + x)),((float) (chunk_y_start + y)),((float) (chunk_z_start + z)),data_vals[((z*data_metric[2]) + (y*data_metric[1]) + (x*data_metric[0]))]);

	    // push id to list of objects that need to be size thresholded
	    i = -1;
	    for(g = 0; g < NO_check_obj_ids; ++g){ if(check_obj_ids[g] == obj_ids[0]){ i = 1; break; } }
	    if(i == -1){

	      ++NO_check_obj_ids;
	      check_obj_ids.push_back(obj_ids[0]);
	      
	    }

	    // adjust obj and obj_ids
	    if(NO_obj_ids > 1){
	      
	      for(i = 0; i < (NO_obj_ids - 1); ++i){ obj_ids[i] = obj_ids[(i + 1)]; }
	      --NO_obj_ids;
	      obj_ids.resize(NO_obj_ids);

	    } else {
	      
	      ++obj;

	      obj_ids[0] = obj;
	      
	      // if this obj_id modulo obj_limit == 0, then it's the first obj of a new batch
	      // so create a new batch of objects
	      if((obj % obj_limit) == 0){ 
				
		detections.push_back(new object_props[obj_limit]);
				
		// if((obj % obj_limit) == 0)
	      }
			      
	      // if(NO_obj_ids > 1) { } else 
	    }
	    
	    // if(existing == -1)
	  } else if(existing >= 0){
	    
	    // assign the current voxel being processed the correct flag value, and update
	    // the bounding box
	    
	    // add object to list of objects within the merging box, provided that it hasn't been added already
	    i = -1;
	    for(g = 0; g < NO_check_obj_ids; ++g){ if(check_obj_ids[g] == existing){ i = 1; break; } }
	    if(i == -1){

	      ++NO_check_obj_ids;
	      check_obj_ids.push_back(existing);

	    }

	    // assign value to array
	    flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (x*data_metric[0]))] = existing;
	    
	    // calculate the obj_batch value for the existing object
	    obj_batch = (long int) floor(((double) existing / (double) obj_limit));
	    
	    // add this point to the existing object
	    detections[obj_batch][(existing - (obj_batch * obj_limit))].AddPoint((chunk_x_start + x),(chunk_y_start + y),(chunk_z_start + z),data_vals[((z*data_metric[2]) + (y*data_metric[1]) + (x*data_metric[0]))]);
	    
	    // for each object detected within the merging volume, update all of the values within the bounding box of this object
	    // to be the same as the other part of this object, and update the bounding box at the same time
	    if(NOi > 1){
	      
	      for(i = 0; i < NOi; ++i){
		
		// move on if this is the existing object
		if(match_init[i] == existing){ 

		  continue; 

		}
		
		// calculate a temporary obj_batch value
		obj_batch_2 = (long int) floor(((double) match_init[i] / (double) obj_limit));		
		
		// add this object to an existing object
		//detections[obj_batch][(existing - (obj_batch * obj_limit))].AddObject(detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))]);
		detections[obj_batch][(existing - (obj_batch * obj_limit))].AddObject(detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))],temp_sparse_reps_grid,temp_sparse_reps_strings,temp_mom0,temp_RAPV,temp_DECPV,temp_ref_spec,temp_obj_spec,temp_vfield);
		
		// flag the `object' values within the bounding box to merge it with the base object
		sx_start = detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetRAmin() - chunk_x_start;
		if(sx_start < 0){ sx_start = 0; }
		sy_start = detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetDECmin() - chunk_y_start;
		if(sy_start < 0){ sy_start = 0; }
		sz_start = detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetFREQmin() - chunk_z_start;
		if(sz_start < 0){ sz_start = 0; }
		for(sz = sz_start; (sz <= detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetFREQmax() - chunk_z_start) && (sz < size_z); ++sz){
		  
		  for(sy = sy_start; (sy <= detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetDECmax() - chunk_y_start) && (sy < size_y); ++sy){
		    
		    for(sx = sx_start; (sx <= detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetRAmax() - chunk_x_start) && (sx < size_x); ++sx){
		      
		      if(flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] == match_init[i]){ flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] = existing; }
		      
		      // for(sx = x_start[(match_init[i] - (obj_batch_2 * obj_limit))]; sx <= x_finish[(match_init[i] - (obj_batch_2 * obj_limit))]; ++sx)
		    }
		    
		    // for(sy = y_start[(match_init[i] - (obj_batch_2 * obj_limit))]; sy <= y_finish[(match_init[i] - (obj_batch_2 * obj_limit))]; ++sy)
		  }
		  
		  // for(sz = z_start[(match_init[i] - (obj_batch_2 * obj_limit))]; sz <= z_finish[(match_init[i] - (obj_batch_2 * obj_limit))]; ++sz)
		}
		
		// re-initialise the `object' values for the match_init[i] object that has just been merged into the existing object, 
		// and push the obj value (match_init[i]) to the top of the array of values
		detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].ReInit();
		if((detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].Get_srep_update() != 0) && (detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].Get_srep_size(0) >= 0)){
		  
		  detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].ReInit_srep();
		  detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].ReInit_mini();

		}
		detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].ReInit_size();
		detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].Set_srep_update(0);

		j = -1;
		for(g = 0; g < NO_obj_ids; ++g){ if(obj_ids[g] == match_init[i]){ j = 1; break; } }
		if(j == -1){
		  
		  ++NO_obj_ids;
		  obj_ids.push_back(match_init[i]);
		  		  
		}
		
		// for(i = 0; i < NOi; ++i)
	      }
	      
	      // if(NOi > 1)
	    }
	    
	    // else . . . if(existing == -1)
	  }
	  
	  // if(flag_vals[((z*size_x*size_y) + (y*size_x) + x)] == -1)
	}
	
	// update progress on display
	while(progress <= (((float) ((z * size_x * size_y) + (y * size_x) + x + 1)) / ((float) (size_x * size_y * size_z)))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
	
	// for(x = 0; x < size_x; ++x)
      }
      
      // for(y = 0; y < size_y; ++y)
    }
    
    // for(z = 0; z < size_z; ++z)
  }
  std::cout << "* done." << std::endl;
  
  // 4. generate sparse representations of coherent objects that are retained after
  // this `chunk' of the datacube is processed
  std::cout << "Generating/updating sparse representations of sources . . . " << std::endl;
  std::cout << "0 | |:| | : | |:| | 100% complete" << std::endl;
  progress = 0.0;
  for(o = 0; o < obj; ++o){
    
    // calculate the obj_batch value for this object
    obj_batch = (long int) floor(((double) o / (double) obj_limit));
    
    // move on if this object has been re-initialised
    if(detections[obj_batch][(o - (obj_batch * obj_limit))].ShowVoxels() < 1){ 
      
      // update progress on display
      while(progress <= (((double) (o + 1) / ((double) obj)))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }  
      continue; 
      
    }
    
    // move on if this object is not within the boundaries of the current chunk, in other words,
    // move on if this object doesn't need to have it's sparse representation and postage stamp images updated
    if((detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmax() < chunk_x_start) || (detections[obj_batch][(o - (obj_batch * obj_limit))].GetDECmax() < chunk_y_start) || (detections[obj_batch][(o - (obj_batch * obj_limit))].GetFREQmax() < chunk_z_start) || (detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmin() >= (chunk_x_start + size_x)) || (detections[obj_batch][(o - (obj_batch * obj_limit))].GetDECmin() >= (chunk_y_start + size_y)) || (detections[obj_batch][(o - (obj_batch * obj_limit))].GetFREQmin() >= (chunk_z_start + size_z))){ 

      // update progress on display
      while(progress <= (((double) (o + 1) / ((double) obj)))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }  
      continue; 
      
    }
    
    // move on if this object's sparse representation and postage stamp images don't need updating
    if((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) >= 0) && (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_update() == -1)){
      
      // update progress on display
      while(progress <= (((double) (o + 1) / ((double) obj)))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }  
      continue; 
	
    }
          
    // update the sparse representation for this object
    
    // 0. if there is no existing sparse representation and postage stamp images, create them from scratch, 
    // otherwise concatenate the existing sparse representation and postage stamp images with the updates
    if((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) < 0) && (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_update() >= 0)){
      
      // branch for creating a brand new sparse representation
      
      // update the sparse_reps_update array
      detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_update(-1);

      // initialise the number of object sections that make up this object
      NOi = 0;   
      
      // calculate the edge of the bounding box in `chunk' co-ordinates
      sx_start = detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmin() - chunk_x_start;
      if(sx_start < 0){ sx_start = 0; }
      sx_finish = detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmax() - chunk_x_start;
      sy_start = detections[obj_batch][(o - (obj_batch * obj_limit))].GetDECmin() - chunk_y_start;
      if(sy_start < 0){ sy_start = 0; }
      sy_finish = detections[obj_batch][(o - (obj_batch * obj_limit))].GetDECmax() - chunk_y_start;
      sz_start = detections[obj_batch][(o - (obj_batch * obj_limit))].GetFREQmin() - chunk_z_start;
      if(sz_start < 0){ sz_start = 0; }
      sz_finish = detections[obj_batch][(o - (obj_batch * obj_limit))].GetFREQmax() - chunk_z_start;
      
      // count the number of object sections that make up this object
      for(sy = sy_start; (sy <= sy_finish) && (sy < size_y); ++sy){
	
	for(sx = sx_start; (sx <= sx_finish) && (sx < size_x); ++sx){
	  
	  // initialise the dummy integer j
	  j = -1;
	  
	  for(sz = sz_start; (sz <= sz_finish) && (sz < size_z); ++sz){
	    
	    // change j to reflect if this voxel in the flag_vals array belongs to the source
	    if(flag_vals[((sz * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))] == o){ 
	      
	      // if previously this was not part of an object string, j == -1, then increment the number of object strings and
	      // adjust the value of j
	      if(j == -1){ ++NOi; }
	      
	      // change j to reflect that this voxel belongs to the object
	      j = 1; 
	      
	    } else { j = -1; }
	    
	    // for(sz = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetFREQmin() - chunk_z_size; sz <= (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetFREQmax() - chunk_z_size); ++sz)
	  }
	  
	  // for(sx = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmin() - chunk_x_size; sx <= (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmax() - chunk_x_size); ++sx)
	}
	
	// for(sy = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmin() - chunk_y_size; sy <= (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmax() - chunk_y_size); ++sy)
      }

      // create an array element to store this sparse representation
      detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_size(0,detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmin());
      detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_size(1,detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmax());
      detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_size(2,detections[obj_batch][(o - (obj_batch * obj_limit))].GetDECmin());
      detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_size(3,detections[obj_batch][(o - (obj_batch * obj_limit))].GetDECmax());
      detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_size(4,detections[obj_batch][(o - (obj_batch * obj_limit))].GetFREQmin());
      detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_size(5,detections[obj_batch][(o - (obj_batch * obj_limit))].GetFREQmax());
      detections[obj_batch][(o - (obj_batch * obj_limit))].Create_srep_grid((((detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmin() + 1) * (detections[obj_batch][(o - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(o - (obj_batch * obj_limit))].GetDECmin() + 1)) + 1));
      detections[obj_batch][(o - (obj_batch * obj_limit))].Create_srep_strings((2 * NOi));
      
      // calculate the extent of the sparse_rep volume
      sx_finish = detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmin() + 1;
      sy_finish = detections[obj_batch][(o - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(o - (obj_batch * obj_limit))].GetDECmin() + 1;
      sz_finish = detections[obj_batch][(o - (obj_batch * obj_limit))].GetFREQmax() - detections[obj_batch][(o - (obj_batch * obj_limit))].GetFREQmin() + 1;
      
      // create arrays to store the postage stamp images and initialise them
      detections[obj_batch][(o - (obj_batch * obj_limit))].Create_mom0((sx_finish * sy_finish));
      detections[obj_batch][(o - (obj_batch * obj_limit))].Create_RAPV((sx_finish * sz_finish));
      detections[obj_batch][(o - (obj_batch * obj_limit))].Create_DECPV((sy_finish * sz_finish));
      detections[obj_batch][(o - (obj_batch * obj_limit))].Create_ospec(sz_finish);
      if(sz_finish >= 10){

	detections[obj_batch][(o - (obj_batch * obj_limit))].Create_rspec((2 * sz_finish));

      } else {

	detections[obj_batch][(o - (obj_batch * obj_limit))].Create_rspec((10 + sz_finish));

      }
      detections[obj_batch][(o - (obj_batch * obj_limit))].Create_vfield((sx_finish * sy_finish));
      
      // initialise the mini_mom0 values
      for(sy = 0; sy < sy_finish; ++sy){
	for(sx = 0; sx < sx_finish; ++sx){
	  detections[obj_batch][(o - (obj_batch * obj_limit))].Set_mom0(((sy * sx_finish) + sx),0.0);
	}
      }

      // initialise the mini_RAPV values
      for(sz = 0; sz < sz_finish; ++sz){
	for(sx = 0; sx < sx_finish; ++sx){
	  detections[obj_batch][(o - (obj_batch * obj_limit))].Set_RAPV(((sz * sx_finish) + sx),0.0);
	}
      }
      
      // initialise the mini_DECPV values
      for(sz = 0; sz < sz_finish; ++sz){
	for(sy = 0; sy < sy_finish; ++sy){
	  detections[obj_batch][(o - (obj_batch * obj_limit))].Set_DECPV(((sz * sy_finish) + sy),0.0);
	}
      }
      
      // initialise the mini_spec values
      for(sz = 0; sz < sz_finish; ++sz){ detections[obj_batch][(o - (obj_batch * obj_limit))].Set_ospec(sz,0.0); }
      if(sz_finish >= 10){

	for(sz = 0; sz < (2 * sz_finish); ++sz){ detections[obj_batch][(o - (obj_batch * obj_limit))].Set_rspec(sz,0.0); }

      } else {

	for(sz = 0; sz < (10 + sz_finish); ++sz){ detections[obj_batch][(o - (obj_batch * obj_limit))].Set_rspec(sz,0.0); }

      }
            
      // initialise the mini_vfield values
      for(sy = 0; sy < sy_finish; ++sy){
	for(sx = 0; sx < sx_finish; ++sx){
	  detections[obj_batch][(o - (obj_batch * obj_limit))].Set_vfield(((sy * sx_finish) + sx),0.0);
	}
      }

      // populate the sparse array representation and postage stamp images
      NOi = 0;
      for(sy = 0; sy < sy_finish; ++sy){
	
	for(sx = 0; sx < sx_finish; ++sx){
	  
	  // write the current number of object strings to the grid
	  detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_grid(((sy * sx_finish) + sx),NOi);
	  
	  // initialise the dummy integer j 
	  j = -1;
	  
	  for(sz = 0; sz < sz_finish; ++sz){
	    
	    // change j to reflect if this voxel in the flag_vals array belongs to the source
	    if(flag_vals[(((sz + sz_start) * data_metric[2]) + ((sy + sy_start) * data_metric[1]) + ((sx + sx_start) * data_metric[0]))] == o){ 
	      
	      // update mini_mom0
	      detections[obj_batch][(o - (obj_batch * obj_limit))].Add_mom0(((sy * sx_finish) + sx),data_vals[(((sz + sz_start) * data_metric[2]) + ((sy + sy_start) * data_metric[1]) + ((sx + sx_start) * data_metric[0]))]);
	      
	      // update mini_RAPV
	      detections[obj_batch][(o - (obj_batch * obj_limit))].Add_RAPV(((sz * sx_finish) + sx),data_vals[(((sz + sz_start) * data_metric[2]) + ((sy + sy_start) * data_metric[1]) + ((sx + sx_start) * data_metric[0]))]);
	      
	      // update mini_DECPV
	      detections[obj_batch][(o - (obj_batch * obj_limit))].Add_DECPV(((sz * sy_finish) + sy),data_vals[(((sz + sz_start) * data_metric[2]) + ((sy + sy_start) * data_metric[1]) + ((sx + sx_start) * data_metric[0]))]);
	      
	      // update mini_obj_spec
	      detections[obj_batch][(o - (obj_batch * obj_limit))].Add_ospec(sz,data_vals[(((sz + sz_start) * data_metric[2]) + ((sy + sy_start) * data_metric[1]) + ((sx + sx_start) * data_metric[0]))]);
  
	      // update mini_vfield
	      detections[obj_batch][(o - (obj_batch * obj_limit))].Add_vfield(((sy * sx_finish) + sx),(((float) (sz + sz_start + chunk_z_start)) * data_vals[(((sz + sz_start) * data_metric[2]) + ((sy + sy_start) * data_metric[1]) + ((sx + sx_start) * data_metric[0]))]));
	      
	      // if previously this was not part of an object string, j == -1, then increment the number of object strings and
	      // write the position to the sparse_reps array
	      if(j == -1){ 
		
		// write the beginning of the object string to sparse_reps
		detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_strings((2 * NOi),(sz_start + chunk_z_start + sz));
		
		// increment NOi
		++NOi; 
		
		// initialise the end of the object string in sparse_reps
		detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_strings(((2 * NOi) - 1),(sz_start + chunk_z_start + sz));

	      } else {
		
		// update the end of the object string in sparse_reps
		detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_strings(((2 * NOi) - 1),(sz_start + chunk_z_start + sz));
		
	      }
	      
	      // change j to reflect that this voxel belongs to the object
	      j = 1; 
	      
	    } else { 
	      
	      // change j to reflect that this voxel doesn't belong to the object
	      j = -1; 
	      
	    }
	    
	    // for(sz = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetFREmin() - chunk_z_start + 1; sz <= (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetFREQmax() - chunk_z_start + 1); ++sz)
	  }
	  
	  // for(sx = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmin() - chunk_x_start + 1; sx <= (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmax() - chunk_x_start + 1); ++sx)
	}
	
	//for(sy = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmin() - chunk_y_start + 1; sy <= (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmax() - chunk_y_start + 1); ++sy)
      }
            
      // update sparse_reps with the final value of NOi
      detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_grid((sx_finish * sy_finish),NOi);
            
      // create the reference spectrum
      if(sz_finish >= 10){
	
	for(sz = -1 * (int) floorf((0.5 * (float) sz_finish)); sz < (sz_finish + ((int) floorf((0.5 * (float) sz_finish)))); ++sz){
	  
	  for(sy = 0; sy < sy_finish; ++sy){
	  
	    for(sx = 0; sx < sx_finish; ++sx){
	    
	      if(((sz + sz_start) >= 0) && ((sz + sz_start) < size_z)){
		
		// update mini_ref_spec	
		detections[obj_batch][(o - (obj_batch * obj_limit))].Add_rspec((sz + ((int) floorf((0.5 * (float) sz_finish)))),data_vals[(((sz + sz_start) * data_metric[2]) + ((sy + sy_start) * data_metric[1]) + ((sx + sx_start) * data_metric[0]))]);	
		
	      }
	      
	    }
	    
	  }
      
	}

	// if(sz_finish >= 10)
      } else {

	for(sz = -5; sz < (sz_finish + 5); ++sz){
	  
	  for(sy = 0; sy < sy_finish; ++sy){
	  
	    for(sx = 0; sx < sx_finish; ++sx){
	    
	      if(((sz + sz_start) >= 0) && ((sz + sz_start) < size_z)){
		
		// update mini_ref_spec
		detections[obj_batch][(o - (obj_batch * obj_limit))].Add_rspec((sz + 5),data_vals[(((sz + sz_start) * data_metric[2]) + ((sy + sy_start) * data_metric[1]) + ((sx + sx_start) * data_metric[0]))]);	
		
	      }
	      
	    }
	    
	  }
      
	}

	// else . . . if(sz_finish >= 10)
      }
      
     // if(sparse_reps[obj_batch][(obj - (obj_batch * obj_limit))][0] < 0) . . . else
    } else {
            
      // branch for updating an existing sparse representation
      
      // count the number of object strings comprising the object
      
      // 1. initialise the number of object strings
      NOi = 0;
      
      // 2. convert the new bounding box of the object into `chunk' co-ordinates
      sx_start = detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmin() - chunk_x_start;
      sx_finish = detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmax() - chunk_x_start;
      sy_start = detections[obj_batch][(o - (obj_batch * obj_limit))].GetDECmin() - chunk_y_start;
      sy_finish = detections[obj_batch][(o - (obj_batch * obj_limit))].GetDECmax() - chunk_y_start;
      sz_start = detections[obj_batch][(o - (obj_batch * obj_limit))].GetFREQmin() - chunk_z_start;
      if(sz_start < 0){ sz_start = 0; }
      sz_finish = detections[obj_batch][(o - (obj_batch * obj_limit))].GetFREQmax() - chunk_z_start;
      
      // 3. search through the bounding box, and do one of the following
      // 3.a if inside the chunk and outside the existing sparse rep grid, run the routine to count the number
      //     of object strings along the current line of sight
      // 3.b otherwise use the grid value to increment NOi
      for(sy = sy_start; sy <= sy_finish; ++sy){
	
	for(sx = sx_start; sx <= sx_finish; ++sx){
	  	  
	  // case 3.a triggered
	  if(((sx > merge_x) || (chunk_x_start == 0)) && ((sy > merge_y) || (chunk_y_start == 0))){
	    	    
	    // initialise the dummy integer j 
	    j = -1;
	    
	    for(sz = sz_start; sz <= sz_finish; ++sz){
	      
	      // change i to reflect if this voxel in the flag_vals array belongs to the source
	      if(flag_vals[((sz * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))] == o){ 
		
		// if previously this was not part of an object string, i == -1, then increment the number of object strings and
		// write the position to the sparse_reps array
		if(j == -1){ 
		  
		  // increment NOi
		  ++NOi; 
		  
		} 
		
		// change j to reflect that this voxel belongs to the object
		j = 1; 
		
	      } else { j = -1; }
	      
	      // for(sz = sz_start; sz <= sz_finish; ++sz)
	    }
	    
	  } else {
	    
	    // case 3.b triggered  
	    NOi+=(detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_grid((((sy + chunk_y_start - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx + chunk_x_start - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_grid((((sy + chunk_y_start - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx + chunk_x_start - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0))));
	    
	  }     
	  
	  // for(sx = sx_start; sx <= sx_finish; ++sx)
	}
	
	// for(sy = sy_start; sy <= sy_finish; ++sy)
      }
      
      // 4. allocate memory to store temporary sparse representation and postage stamp images, then
      //    initialise the postage stamp image arrays
      temp_mom0.resize(0);
      for(sy = 0; sy < (sy_finish - sy_start + 1); ++sy){
	for(sx = 0; sx < (sx_finish - sx_start + 1); ++sx){
	  temp_mom0.push_back(0.0);
	}
      }
      temp_RAPV.resize(0);
      for(sz = 0; sz < (sz_finish - sz_start + 1); ++sz){
	for(sx = 0; sx < (sx_finish - sx_start + 1); ++sx){
	  temp_RAPV.push_back(0.0);
	}
      }
      temp_DECPV.resize(0);
      for(sz = 0; sz < (sz_finish - sz_start + 1); ++sz){
	for(sy = 0; sy < (sy_finish - sy_start + 1); ++sy){
	  temp_DECPV.push_back(0.0);
	}
      }
      temp_obj_spec.resize(0);
      for(sz = 0; sz < (sz_finish - sz_start + 1); ++sz){ temp_obj_spec.push_back(0.0); }
      if((sz_finish - sz_start + 1) >= 10){

	temp_ref_spec.resize(0);
	for(sz = 0; sz < (2 * (sz_finish - sz_start + 1)); ++sz){ temp_ref_spec.push_back(0.0); }

      } else {

	temp_ref_spec.resize(0);
	for(sz = 0; sz < (sz_finish - sz_start + 11); ++sz){ temp_ref_spec.push_back(0.0); }

      }
      temp_vfield.resize(0);
      for(sy = 0; sy < (sy_finish - sy_start + 1); ++sy){
	for(sx = 0; sx < (sx_finish - sx_start + 1); ++sx){
	  temp_vfield.push_back(0.0);
	}
      }
     
      // 5. populate temp_sparse_rep arrays
      // 5.a if inside the chunk and outside the existing sparse rep grid, run the routine to count the number
      //     of object strings along the current line of sight, populate new strings array and add values to new 
      //     postage stamp images
      // 5.b if outside of the chunk, use the grid value to increment NOi, populate strings array from previous version and 
      //     fill in valuese for postage stamp images

      // initialise number of object strings
      NOi = 0;
      
      // initialise RA PV image
      for(sz = 0; sz < (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1); ++sz){

	for(sx = 0; sx < (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1); ++sx){
	
	  temp_RAPV[(((sz + detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) - sz_start - chunk_z_start) * (sx_finish - sx_start + 1)) + sx + detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) - sx_start - chunk_x_start)]+=detections[obj_batch][(o - (obj_batch * obj_limit))].Get_RAPV(((sz * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx));
	  
	}

      }

      // initialise Dec PV image
      for(sz = 0; sz < (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1); ++sz){

	for(sy = 0; sy < (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2) + 1); ++sy){
	
	  temp_DECPV[(((sz + detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) - sz_start - chunk_z_start) * (sy_finish - sy_start + 1)) + sy + detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2) - sy_start - chunk_y_start)]+=detections[obj_batch][(o - (obj_batch * obj_limit))].Get_DECPV(((sz * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2) + 1)) + sy));
	  
	}

      }

      // initialise mini spec's
      for(sz = 0; sz < (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1); ++sz){ temp_obj_spec[(sz + detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) - sz_start - chunk_z_start)]+=detections[obj_batch][(o - (obj_batch * obj_limit))].Get_ospec(sz); }

      if((sz_finish - sz_start + 1) >= 10){

	j = sz_start + chunk_z_start - ((int) floorf((0.5 * (float) (sz_finish - sz_start + 1))));

	if((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1) >= 10){

	  k = detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) - ((int) floorf((0.5 * (float) (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1))));

	  for(sz = 0; sz < (2 * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); ++sz){

	    temp_ref_spec[(sz + k - j)]+=detections[obj_batch][(o - (obj_batch * obj_limit))].Get_rspec(sz);

	  }

	} else {
	  
	  k = detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) - 5;

	  for(sz = 0; sz < (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 11); ++sz){
	    
	    temp_ref_spec[(sz + k - j)]+=detections[obj_batch][(o - (obj_batch * obj_limit))].Get_rspec(sz);

	  }
	  
	}

      } else {

	j = sz_start + chunk_z_start - 5;

	if((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1) >= 10){

	  k = detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) - ((int) floorf((0.5 * (float) (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1))));

	  for(sz = 0; sz < (2 * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); ++sz){

	    temp_ref_spec[(sz + k - j)]+=detections[obj_batch][(o - (obj_batch * obj_limit))].Get_rspec(sz);

	  }

	} else {
	  
	  k = detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) - 5;

	  for(sz = 0; sz < (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 11); ++sz){
	    
	    temp_ref_spec[(sz + k - j)]+=detections[obj_batch][(o - (obj_batch * obj_limit))].Get_rspec(sz);

	  }
	  
	}

      }
      
      // search through LoS's of updated bounding box and take action as dictated by case 5.a or 5.b
      for(sy = sy_start; sy <= sy_finish; ++sy){
	
	for(sx = sx_start; sx <= sx_finish; ++sx){
	  	  
	  // add value to temp_sparse_rep_grid
	  temp_sparse_reps_grid[(((sy - sy_start) * (sx_finish - sx_start + 1)) + sx - sx_start)] = NOi;
	  
	  // case 5.a triggered for moment-0 image and sparse representation
	  if(((sx > merge_x) || (chunk_x_start == 0)) && ((sy > merge_y) || (chunk_y_start == 0))){
	    
	    // initialise the dummy integer j 
	    j = -1;
	    
	    for(sz = sz_start; sz <= sz_finish; ++sz){
	      
	      // change j to reflect if this voxel in the flag_vals array belongs to the source
	      if(flag_vals[((sz * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))] == o){ 
		
		// update mini_mom0
		temp_mom0[(((sy - sy_start) * (sx_finish - sx_start + 1)) + sx - sx_start)]+=data_vals[((sz * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))];
		
		// update mini_RAPV
		temp_RAPV[(((sz - sz_start) * (sx_finish - sx_start + 1)) + sx - sx_start)]+=data_vals[((sz * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))];
		
		// update mini_DECPV
		temp_DECPV[(((sz - sz_start) * (sy_finish - sy_start + 1)) + sy - sy_start)]+=data_vals[((sz * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))];
		
		// update mini_obj_spec
		temp_obj_spec[(sz - sz_start)]+=data_vals[((sz * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))];
		
		// update mini_vfield
		temp_vfield[(((sy - sy_start) * (sx_finish - sx_start + 1)) + sx - sx_start)]+=(((float) (sz + chunk_z_start)) * data_vals[((sz * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))]);
		
		// if previously this was not part of an object string, j == -1, then increment the number of object strings and
		// write the position to the sparse_reps array
		if(j == -1){ 
		  
		  // write the beginning of the object string to sparse_reps
		  temp_sparse_reps_strings[(2 * NOi)] = chunk_z_start + sz;
		  
		  // increment NOi
		  ++NOi; 
		  		  
		  // initialise the end of the object string in sparse_reps
		  temp_sparse_reps_strings[((2 * NOi) - 1)] = chunk_z_start + sz;
		  
		} else {
		  
		  // update the end of the object string in sparse_reps
		  temp_sparse_reps_strings[((2 * NOi) - 1)] = chunk_z_start + sz;
		  
		}
		
		// change j to reflect that this voxel belongs to the object
		j = 1; 
		
	      } else { 
		
		// change j to reflect that this voxel doesn't belong to the object
		j = -1; 
		
	      }
	      
	      // for(sz = sz_start; sz <= sz_finish; ++sz)
	    }
	    
	    // add values to temp_ref_spec
	    if((sz_finish - sz_start + 1) >= 10){
	      
	      k = (int) floorf((0.5 * (float) (sz_finish - sz_start + 1)));
	      
	      for(sz = 0; sz < (2 * (sz_finish - sz_start + 1)); ++sz){
	  	
		if(((sz + sz_start - k) >= 0) && ((sz + sz_start - k) < size_z)){
		  
		  // update mini_ref_spec
		  temp_ref_spec[(sz + sz_start - k)]+=data_vals[(((sz + sz_start - k) * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))];	
		  
		}
		
	      }
	      
	      // if(sz_finish >= 10)
	    } else {
	      
	      for(sz = 0; sz < (sz_finish - sz_start + 11); ++sz){
		
		if(((sz + sz_start - 5) >= 0) && ((sz + sz_start - 5) < size_z)){
		  
		  // update mini_ref_spec
		  temp_ref_spec[(sz + sz_start - 5)]+=data_vals[(((sz + sz_start - 5) * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))];	
		  
		}
		
	      }
	      
	      // else . . . if(sz_finish >= 10)
	    }
	    
	  } else {
	    // case 5.b triggered  
	    	    
	    // update temp_sparse_rep_strings
	    j = 0;
	    if(((sx + chunk_x_start) >= detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0)) && ((sx + chunk_x_start) <= (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) + 1)) && ((sy + chunk_y_start) >= detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2)) && ((sy + chunk_y_start) <= (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(3) + 1))){
	      
	      for(g = detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_grid((((sy + chunk_y_start - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx + chunk_x_start - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0))); g < detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_grid((((sy + chunk_y_start - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx + chunk_x_start - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1)); ++g){
		
		temp_sparse_reps_strings[(NOi + j)] = detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_strings((2 * g));
		temp_sparse_reps_strings[(NOi + j + 1)] = detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_strings(((2 * g) + 1));
		++j;
		
	      }
	      
	      // update moment 0 postage stamp image
	      temp_mom0[(((sy - sy_start) * (detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmin() + 1)) + sx - sx_start)] = detections[obj_batch][(o - (obj_batch * obj_limit))].Get_mom0((((sy + chunk_y_start - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx + chunk_x_start - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0)));
	      
	      // update vfield postage stamp image
	      temp_vfield[(((sy - sy_start) * (detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmin() + 1)) + sx - sx_start)] = detections[obj_batch][(o - (obj_batch * obj_limit))].Get_vfield((((sy + chunk_y_start - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx + chunk_x_start - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0)));

	      // update NOi
	      NOi+=(detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_grid((((sy + chunk_y_start - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx + chunk_x_start - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_grid((((sy + chunk_y_start - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx + chunk_x_start - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0))));
	      	      
	    }
	    
	    // else . . . case 5.b
	  }
	    
	  // for(sx = sx_start; sx <= sx_finish; ++sx)
	}
		
	// for(sy = sy_start; sy <= sy_finish; ++sy)
      }
            
      temp_sparse_reps_grid[((sx_finish - sx_start + 1) * (sy_finish - sy_start + 1))] = NOi;

      // 6. replace the current sparse representation with the new sparse representation
      detections[obj_batch][(o - (obj_batch * obj_limit))].ReInit_srep();
      detections[obj_batch][(o - (obj_batch * obj_limit))].ReInit_mini();
 
      // update the sparse_reps_size array with the new bounding box
      detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_size(0,detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmin());
      detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_size(1,detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmax());
      detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_size(2,detections[obj_batch][(o - (obj_batch * obj_limit))].GetDECmin());
      detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_size(3,detections[obj_batch][(o - (obj_batch * obj_limit))].GetDECmax());
      detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_size(4,detections[obj_batch][(o - (obj_batch * obj_limit))].GetFREQmin());
      detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_size(5,detections[obj_batch][(o - (obj_batch * obj_limit))].GetFREQmax());

      // update the sparse_reps_update array
      detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_update(-1);
          
      // populate the new sparse_reps and postage stamps
      detections[obj_batch][(o - (obj_batch * obj_limit))].Create_srep_grid((1 + ((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2) + 1))));
      detections[obj_batch][(o - (obj_batch * obj_limit))].Create_srep_strings((2 * NOi));
      detections[obj_batch][(o - (obj_batch * obj_limit))].Create_mom0(((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2) + 1)));
      detections[obj_batch][(o - (obj_batch * obj_limit))].Create_RAPV(((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1)));
      detections[obj_batch][(o - (obj_batch * obj_limit))].Create_DECPV(((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2) + 1) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1)));
      detections[obj_batch][(o - (obj_batch * obj_limit))].Create_vfield(((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2) + 1)));
  
      detections[obj_batch][(o - (obj_batch * obj_limit))].Create_ospec((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1));
      if((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1) >= 10){
       
	detections[obj_batch][(o - (obj_batch * obj_limit))].Create_rspec((2 * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1)));
	
      } else {
	
	detections[obj_batch][(o - (obj_batch * obj_limit))].Create_rspec((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 11));

      }
   
      for(j = 0; j < (1 + ((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2) + 1))); ++j){ detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_grid(j,temp_sparse_reps_grid[j]); }
      for(j = 0; j < (2 * NOi); ++j){ detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_strings(j,temp_sparse_reps_strings[j]); }
      for(j = 0; j < ((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); ++j){ detections[obj_batch][(o - (obj_batch * obj_limit))].Set_mom0(j,temp_mom0[j]); }
      for(j = 0; j < ((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); ++j){ detections[obj_batch][(o - (obj_batch * obj_limit))].Set_RAPV(j,temp_RAPV[j]); }
      for(j = 0; j < ((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2) + 1) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); ++j){ detections[obj_batch][(o - (obj_batch * obj_limit))].Set_DECPV(j,temp_DECPV[j]); }
      for(j = 0; j < ((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); ++j){ detections[obj_batch][(o - (obj_batch * obj_limit))].Set_vfield(j,temp_vfield[j]); }

      for(j = 0; j < (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1); ++j){ detections[obj_batch][(o - (obj_batch * obj_limit))].Set_ospec(j,temp_obj_spec[j]); }
      if((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1) >= 10){

	for(j = 0; j < (2 * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); ++j){ detections[obj_batch][(o - (obj_batch * obj_limit))].Set_rspec(j,temp_ref_spec[j]); }
	
      } else {

	for(j = 0; j < (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 11); ++j){ detections[obj_batch][(o - (obj_batch * obj_limit))].Set_rspec(j,temp_ref_spec[j]); }

      }
      
      // else . . . if(sparse_reps[obj_batch][(obj - (obj_batch * obj_limit))][0] < 0)
    }
    
    // update progress on display
    while(progress <= (((double) (o + 1) / ((double) obj)))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }  
   
    // for(o = 0; o < obj; ++o)
  }
  std::cout << "* done." << std::endl;
    
  // 5. return the number of coherent objects 
  return obj;
  
}

// functions using doubles

int CreateObjects(double * data_vals, int * flag_vals, int size_x, int size_y, int size_z, int chunk_x_start, int chunk_y_start, int chunk_z_start, int merge_x, int merge_y, int merge_z, int min_x_size, int min_y_size, int min_z_size, int min_v_size, double intens_thresh_min, double intens_thresh_max, int flag_value, int start_obj, vector<object_props_dbl *> & detections, vector<int> & obj_ids, vector<int> & check_obj_ids, int obj_limit, int max_x_val, int max_y_val, int max_z_val, int ss_mode, int * data_metric, int * xyz_order){
 
  int x,y,z,obj,existing,sx,sy,sz,sx_start,sy_start,sz_start,sx_finish,sy_finish,sz_finish,init_limit,g;
  int i,j,k,NOi, obj_batch, obj_batch_2,NO_obj_ids,NO_check_obj_ids;
  int prev,x_start,y_start,temp_x[3],temp_y[3],temp_z[3];
  float progress;
  vector<double> temp_mom0, temp_RAPV, temp_DECPV, temp_obj_spec, temp_ref_spec, temp_vfield;
  vector<int> match_init, temp_sparse_reps_grid, temp_sparse_reps_strings;

  // convert the merging lengths from distances (in pixels/voxels) to empty space between linked detections (in pixels/voxels) --- subtract 1 --- this is a more naturual
  // definition when processing a large file in chunks, because it makes it easier to build up sources that need to be retained for the rest of the input file
  --merge_x;
  --merge_y;
  --merge_z;

  // reorder the datacube and subcube limits to be in x,y,z order
  temp_x[0] = chunk_x_start;
  temp_y[0] = chunk_y_start;
  temp_z[0] = chunk_z_start;
  temp_x[1] = size_x;
  temp_y[1] = size_y;
  temp_z[1] = size_z;
  temp_x[2] = max_x_val;
  temp_y[2] = max_y_val;
  temp_z[2] = max_z_val;
  switch(xyz_order[0]){
  case 1:
    chunk_x_start = temp_x[0];
    size_x = temp_x[1];
    max_x_val = temp_x[2];
    break;
  case 2:
    chunk_x_start = temp_y[0];
    size_x = temp_y[1];
    max_x_val = temp_y[2];
    break;
  case 3:
    chunk_x_start = temp_z[0];
    size_x = temp_z[1];
    max_x_val = temp_z[2];
    break;
  default:
    chunk_x_start = temp_x[0];
    size_x = temp_x[1];
    max_x_val = temp_x[2];
    break;
  }
  switch(xyz_order[1]){
  case 1:
    chunk_y_start = temp_x[0];
    size_y = temp_x[1];
    max_y_val = temp_x[2];
    break;
  case 2:
    chunk_y_start = temp_y[0];
    size_y = temp_y[1];
    max_y_val = temp_y[2];
    break;
  case 3:
    chunk_y_start = temp_z[0];
    size_y = temp_z[1];
    max_y_val = temp_z[2];
    break;
  default:
    chunk_y_start = temp_y[0];
    size_y = temp_y[1];
    max_y_val = temp_y[2];
    break;
  }
  switch(xyz_order[2]){
  case 1:
    chunk_z_start = temp_x[0];
    size_z = temp_x[1];
    max_z_val = temp_x[2];
    break;
  case 2:
    chunk_z_start = temp_y[0];
    size_z = temp_y[1];
    max_z_val = temp_y[2];
    break;
  case 3:
    chunk_z_start = temp_z[0];
    size_z = temp_z[1];
    max_z_val = temp_z[2];
    break;
  default:
    chunk_z_start = temp_z[0];
    size_z = temp_z[1];
    max_z_val = temp_z[2];
    break;
  }

  // create temporary arrays
  if(10000 > (((2 * merge_x) + 3) * ((2 * merge_y) + 3) * ((2 * merge_z) + 3))){
    
    init_limit = 10000;
    
  } else {
    
    init_limit = ((2 * merge_x) + 3) * ((2 * merge_y) + 3) * ((2 * merge_z) + 3) + 10;
    
  }
  match_init.reserve(init_limit);
  match_init.resize(0);

  temp_sparse_reps_grid.reserve(1000000);
  temp_sparse_reps_grid.resize(0);
  temp_sparse_reps_strings.reserve(1000000);
  temp_sparse_reps_strings.resize(0);
  temp_mom0.reserve(1000000);
  temp_mom0.resize(0);
  temp_RAPV.reserve(1000000);
  temp_RAPV.resize(0);
  temp_DECPV.reserve(1000000);
  temp_DECPV.resize(0);
  temp_obj_spec.reserve(1000000);
  temp_obj_spec.resize(0);
  temp_ref_spec.reserve(1000000);
  temp_ref_spec.resize(0);
  temp_vfield.reserve(1000000);
  temp_vfield.resize(0);
 
  // 0. initialise variables and arrays
  obj = start_obj;
  x_start = 0;
  y_start = 0;
  if(chunk_x_start > 0){ x_start = merge_x + 1; }
  if(chunk_y_start > 0){ y_start = merge_y + 1; }
  NO_obj_ids = obj_ids.size();
  NO_check_obj_ids = check_obj_ids.size();

  // 1. Create list of `coherent' objects from neighbouring voxels
  
  // for each grid point, check if it has been flagged as a source voxel, and assign
  // it an object ID if it is. If a neighbouring voxel has already been flagged, assign 
  // that ID to the grid point, otherwise, assign current value of obj and increment obj.
  std::cout << "0 | |:| | : | |:| | 100% complete" << std::endl;
  progress = 0.0;
  for(z = 0; ((z < size_z) && ((z + chunk_z_start) < max_z_val)); ++z){

    // if the z value is sufficiently large that objects have started to pop out of the merging box,
    // then check if the objects outside of the merging box are sufficiently large
    if(z > (merge_z + 1)){
      
      // initialise NOi to store the number of current check_obj_ids values that haven't moved outside of 
      // the merging volume, and that form the basis for the new check_obj_ids list
      NOi = 0;
      
      // apply size threshold to objects added since last size thresholding
      for(i = 0; i < NO_check_obj_ids; ++i){
	
	// calculate the obj_batch value for the existing object
	obj_batch = (int) floorf(((float) check_obj_ids[i] / (float) obj_limit));
	
	// move to the next object if this one has been re-initialised
	if(detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	
	// check if this object has moved outside of the merging box available to the next plane, if it 
	// has apply a size threshold, otherwise use it to seed the new list of objects to be size thresholded
	// the next time
	if(detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmax() < (chunk_z_start + z - merge_z - 1)){
	  
	  // apply the size threshold, and if it fails re-initialise the object and pop its id to the list
	  // of available obj_ids
	  if((((detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetRAmin() + 1) < min_x_size) || ((detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetDECmin() + 1) < min_y_size) || ((detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmax() - detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmin() + 1) < min_z_size) || (detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].ShowVoxels() < min_v_size)) && (((detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetRAmin() - chunk_x_start) > merge_x) || (chunk_x_start == 0)) && ((((chunk_x_start + size_x - 1 - detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetRAmax()) > merge_x) || ((chunk_x_start + size_x) >= max_x_val)) && (((detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetDECmin() - chunk_y_start) > merge_y) || (chunk_y_start == 0)) && (((chunk_y_start + size_y - 1 - detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetDECmax()) > merge_y) || ((chunk_y_start + size_y) >= max_y_val)) && (((detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmin() - chunk_z_start) > merge_z) || (chunk_z_start == 0)) && (((chunk_z_start + size_z - 1 - detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmax()) > merge_z) || ((chunk_z_start + size_z) >= max_z_val)))){
	    
	    // re-initialise object
	    detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].ReInit();
	    if((detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].Get_srep_update() != 0) && (detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].Get_srep_size(0) >= 0)){
	      
	      detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].ReInit_srep();
	      detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].ReInit_mini();

	    }
	    detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].ReInit_size();
	    detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].Set_srep_update(0);

	    // flag the `object' values within the bounding box and remove it from the flag_vals array
	    sx_start = detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetRAmin() - chunk_x_start;
	    if(sx_start < 0){ sx_start = 0; }
	    sy_start = detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetDECmin() - chunk_y_start;
	    if(sy_start < 0){ sy_start = 0; }
	    sz_start = detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmin() - chunk_z_start;
	    if(sz_start < 0){ sz_start = 0; }
	    for(sz = sz_start; ((sz <= (detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmax() - chunk_z_start)) && (sz < size_z)); ++sz){
	      
	      for(sy = sy_start; ((sy <= (detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetDECmax() - chunk_y_start)) && (sy < size_y)); ++sy){
		
		for(sx = sx_start; ((sx <= (detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetRAmax() - chunk_x_start)) && (sx < size_x)); ++sx){
		  
		  if(flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx * data_metric[0]))] == check_obj_ids[i]){ flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx * data_metric[0]))] = -99; }
		  
		  // for(sx = x_start[(match_init[i] - (obj_batch_2 * obj_limit))]; sx <= x_finish[(match_init[i] - (obj_batch_2 * obj_limit))]; ++sx)
		}
		
		// for(sy = y_start[(match_init[i] - (obj_batch_2 * obj_limit))]; sy <= y_finish[(match_init[i] - (obj_batch_2 * obj_limit))]; ++sy)
	      }
	      
	      // for(sz = z_start[(match_init[i] - (obj_batch_2 * obj_limit))]; sz <= z_finish[(match_init[i] - (obj_batch_2 * obj_limit))]; ++sz)
	    }
	    
	    // add object id to list of available ids
	    j = -1;
	    for(g = 0; g < NO_obj_ids; ++g){ if(obj_ids[g] == check_obj_ids[i]){ j = 1; break; } }
	    if(j == -1){
	      ++NO_obj_ids;
	      obj_ids.push_back(check_obj_ids[i]);
	    }

	  }
	  
	  // if(z_finish[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))] < (z - merge_z - 1))
	} else {
	  
	  check_obj_ids[NOi] = check_obj_ids[i];
	  ++NOi;
	  
	  // else . . . if(z_finish[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))] < (z - merge_z - 1))
	}
	
	// for(i = 0; i < NO_check_obj_ids; ++i)
      }

      // update NO_check_obj_ids
      NO_check_obj_ids = NOi;
      check_obj_ids.resize(NOi);
      
      // if(z > (merge_z + 2))
    }
    
    for(y = y_start; (y < size_y) && ((chunk_y_start + y) < max_y_val); ++y){
      
      for(x = x_start; (x < size_x) && ((chunk_x_start + x) < max_x_val); ++x){
	
	// if this is a source voxel, check if it is associated with a previously identified
	// source voxel
	if(flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (x*data_metric[0]))] == flag_value){ 
	  
	  // 1.  initialise variables
	  existing = flag_value;
	  NOi = 0;
	  match_init.resize(0);
	  
	  // 3. initial pass, check the voxels within the merging distance to see if any have been
	  // identified as previous objects
	  
	  // 3a. set limits of region to be searched
	  sx_start = x - merge_x - 1;
	  if(sx_start < 0){ sx_start = 0; }
	  sy_start = y - merge_y - 1;
	  if(sy_start < 0){ sy_start = 0; }
	  sz_start = z - merge_z - 1;
	  if(sz_start < 0){ sz_start = 0; }
	  
	  // 3b. search through previous region and modify existing flag if required

	  if(ss_mode == 1){
	    
	    // use a rectangular paralleliped to do the linking

	    for(sz = sz_start; (sz < z); ++sz){
	    
	      for(sy = sy_start; (sy <= (y + merge_y + 1)) && (sy < size_y); ++sy){
	      
		for(sx = sx_start; (sx <= (x + merge_x + 1)) && (sx < size_x); ++sx){
		
		  // check if this test voxel is source
		  if((flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] == flag_value) || (flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] >= 0)){
		    
		    // check if this is a new objid
		    prev = -1;
		    for(i = 0; (i < NOi) && (i < init_limit); ++i){
		      
		      if(match_init[i] == flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]){
			
			prev = 1;
			break;
			
		      }
		      
		    }
		    if(prev == -1){ 
		      
		      // add current object to list of previous matches
		      match_init.push_back(flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]);
		      ++NOi;

		    }
		    
		    // update existing flag
		    if(((flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] < existing) && (flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] >= 0)) || (existing == flag_value)){ existing = flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]; }
		    
		  }
		  
		  // for(sx = x_start; (sx < (x + merge_x + 1)) && (sx < size_x); ++sx)
		}
		
		// for(sy = y_start; (sy < (y + merge_y + 1)) && (sy < size_y); ++sy)
	      }
	      
	      // for(sz = z_start; (sz < z); ++sz)
	    }
	    for(sy = sy_start; (sy < y); ++sy){
	      
	      for(sx = sx_start; (sx <= (x + merge_x + 1)) && (sx < size_x); ++sx){
		
		// check if this test voxel is source
		if((flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] == flag_value) || (flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] >= 0)){
		  
		  // check if this is a new objid, move on if not
		  prev = -1;
		  for(i = 0; (i < NOi) && (i < init_limit); ++i){
		    
		    if(match_init[i] == flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]){
		      
		      prev = 1;
		      break;
		      
		    }
		    
		  }
		  if(prev == -1){ 
		    
		    // add current object to list of previous matches
		    ++NOi;
		    match_init.push_back(flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]);

		  }
		  
		  // update existing flag
		  if(((flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] < existing) && (flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] >= 0)) || (existing == flag_value)){ existing = flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]; }
		  
		}		
		
		// for(sx = x_start; (sx < = (x + merge_x + 1)) && (sx < size_x); ++sx)
	      }
	      
	      // for(sy = y_start; (sy < y) && (sy < size_y); ++sy)
	    }
	    for(sx = sx_start; sx < x; ++sx){
	      
	      // check if this test voxel is source
	      if((flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))] == flag_value) || (flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))] >= 0)){
		
		// check if this is a new objid, move on if not
		prev = -1;
		for(i = 0; (i < NOi) && (i < init_limit); ++i){
		  
		  if(match_init[i] == flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))]){
		    
		    prev = 1;
		    break;
		    
		  }
		  
		}
		if(prev == -1){ 
		  
		  // add current object to list of previous matches
		  ++NOi;
		  match_init.push_back(flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))]);

		}
		
		// update existing flag
		if(((flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))] < existing) && (flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))] >= 0)) || (existing == flag_value)){ existing = flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))]; }
		
	      }	  
	      
	      // for(sx = x_start; sx <= x; ++sx)
	    }
	  
	    // if(ss_mode == 1)
	  } else {

	    // use a cylinder for the linking lengths

	    for(sz = sz_start; (sz < z); ++sz){
	    
	      for(sy = sy_start; (sy <= (y + merge_y + 1)) && (sy < size_y); ++sy){
	      
		for(sx = sx_start; (sx <= (x + merge_x + 1)) && (sx < size_x); ++sx){
		
		  // check if the test voxel is within the boundaries of the spatial ellipse
		  if(((((float) (x - sx)) * ((float) (x - sx)) / (((float) (merge_x + 1)) * ((float) (merge_x + 1)))) + (((float) (y - sy)) * ((float) (y - sy)) / (((float) (merge_y + 1)) * ((float) (merge_y + 1))))) > 1.0){ continue; }

		  // check if this test voxel is source
		  if((flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] == flag_value) || (flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] >= 0)){

		    // check if this is a new objid
		    prev = -1;
		    for(i = 0; (i < NOi) && (i < init_limit); ++i){
		      
		      if(match_init[i] == flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]){
			
			prev = 1;
			break;
			
		      }
		      
		    }
		    if(prev == -1){ 
		      
		      // add current object to list of previous matches
		      ++NOi;
		      match_init.push_back(flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]);

		    }
		    
		    // update existing flag
		    if(((flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] < existing) && (flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] >= 0)) || (existing == flag_value)){ existing = flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]; }
		    
		  }
		  
		  // for(sx = x_start; (sx < (x + merge_x + 1)) && (sx < size_x); ++sx)
		}
		
		// for(sy = y_start; (sy < (y + merge_y + 1)) && (sy < size_y); ++sy)
	      }
	      
	      // for(sz = z_start; (sz < z); ++sz)
	    }
	    for(sy = sy_start; (sy < y); ++sy){
	      
	      for(sx = sx_start; (sx <= (x + merge_x + 1)) && (sx < size_x); ++sx){
		
		// check if the test voxel is within the boundaries of the spatial ellipse
		if(((((float) (x - sx)) * ((float) (x - sx)) / (((float) (merge_x + 1)) * ((float) (merge_x + 1)))) + (((float) (y - sy)) * ((float) (y - sy)) / (((float) (merge_y + 1)) * ((float) (merge_y + 1))))) > 1.0){ continue; }

		// check if this test voxel is source
		if((flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] == flag_value) || (flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] >= 0)){
		  
		  // check if this is a new objid, move on if not
		  prev = -1;
		  for(i = 0; (i < NOi) && (i < init_limit); ++i){
		    
		    if(match_init[i] == flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]){
		      
		      prev = 1;
		      break;
		      
		    }
		    
		  }
		  if(prev == -1){ 
		    
		    // add current object to list of previous matches
		    ++NOi;
		    match_init.push_back(flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]);

		  }
		  
		  // update existing flag
		  if(((flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] < existing) && (flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] >= 0)) || (existing == flag_value)){ existing = flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]; }
		  
		}		
		
		// for(sx = x_start; (sx < = (x + merge_x + 1)) && (sx < size_x); ++sx)
	      }
	      
	      // for(sy = y_start; (sy < y) && (sy < size_y); ++sy)
	    }
	    for(sx = sx_start; sx < x; ++sx){
	      
	      // check if this test voxel is source
	      if((flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))] == flag_value) || (flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))] >= 0)){
		
		// check if this is a new objid, move on if not
		prev = -1;
		for(i = 0; (i < NOi) && (i < init_limit); ++i){
		  
		  if(match_init[i] == flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))]){
		    
		    prev = 1;
		    break;
		    
		  }
		  
		}
		if(prev == -1){ 
		  
		  // add current object to list of previous matches
		  ++NOi;
		  match_init.push_back(flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))]);
		  
		}
		
		// update existing flag
		if(((flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))] < existing) && (flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))] >= 0)) || (existing == flag_value)){ existing = flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))]; }
		
	      }	  
	      
	      // for(sx = x_start; sx <= x; ++sx)
	    }

	    // else . . . if(ss_mode == 1)
	  }

	  // 3c. assign an object number to this voxel, depending upon the value of the existing flag,
	  // and if this is part of an existing object or objects, then daisy chain from this voxel to all the others
	  if(existing == flag_value){
	    
	    // get obj value from array of obj_ids, then replace it with an incremented obj value
	    
	    // ensure that the obj_ids array is sorted
	    sort(obj_ids.begin(),obj_ids.end());
	    
	    // assign value to array
	    flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (x*data_metric[0]))] = obj_ids[0]; 
	    
	    // calculate the obj_batch value for this obj_id
	    obj_batch = (int) floorf(((float) obj_ids[0] / (float) obj_limit));
	    
	    // add this initial point to the object
	    detections[obj_batch][(obj_ids[0] - (obj_batch * obj_limit))].AddPoint(((double) (chunk_x_start + x)),((double) (chunk_y_start + y)),((double) (chunk_z_start + z)),data_vals[((z*data_metric[2]) + (y*data_metric[1]) + (x*data_metric[0]))]);

	    // push id to list of objects that need to be size thresholded
	    i = -1;
	    for(g = 0; g < NO_check_obj_ids; ++g){ if(check_obj_ids[g] == obj_ids[0]){ i = 1; break; } }
	    if(i == -1){

	      ++NO_check_obj_ids;
	      check_obj_ids.push_back(obj_ids[0]);
	      
	    }

	    // adjust obj and obj_ids
	    if(NO_obj_ids > 1){
	      
	      for(i = 0; i < (NO_obj_ids - 1); ++i){ obj_ids[i] = obj_ids[(i + 1)]; }
	      --NO_obj_ids;
	      obj_ids.resize(NO_obj_ids);

	    } else {
	      
	      ++obj;

	      obj_ids[0] = obj;
	      
	      // if this obj_id modulo obj_limit == 0, then it's the first obj of a new batch
	      // so create a new batch of objects
	      if((obj % obj_limit) == 0){ 
				
		detections.push_back(new object_props_dbl[obj_limit]);
				
		// if((obj % obj_limit) == 0)
	      }
			      
	      // if(NO_obj_ids > 1) { } else 
	    }
	    
	    // if(existing == -1)
	  } else if(existing >= 0){
	    
	    // assign the current voxel being processed the correct flag value, and update
	    // the bounding box
	    
	    // add object to list of objects within the merging box, provided that it hasn't been added already
	    i = -1;
	    for(g = 0; g < NO_check_obj_ids; ++g){ if(check_obj_ids[g] == existing){ i = 1; break; } }
	    if(i == -1){

	      ++NO_check_obj_ids;
	      check_obj_ids.push_back(existing);

	    }

	    // assign value to array
	    flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (x*data_metric[0]))] = existing;
	    
	    // calculate the obj_batch value for the existing object
	    obj_batch = (int) floorf(((float) existing / (float) obj_limit));
	    
	    // add this point to the existing object
	    detections[obj_batch][(existing - (obj_batch * obj_limit))].AddPoint((chunk_x_start + x),(chunk_y_start + y),(chunk_z_start + z),data_vals[((z*data_metric[2]) + (y*data_metric[1]) + (x*data_metric[0]))]);
	    
	    // for each object detected within the merging volume, update all of the values within the bounding box of this object
	    // to be the same as the other part of this object, and update the bounding box at the same time
	    if(NOi > 1){
	      
	      for(i = 0; i < NOi; ++i){
		
		// move on if this is the existing object
		if(match_init[i] == existing){ 

		  continue; 

		}
		
		// calculate a temporary obj_batch value
		obj_batch_2 = (int) floorf(((float) match_init[i] / (float) obj_limit));		
		
		// add this object to an existing object
		//detections[obj_batch][(existing - (obj_batch * obj_limit))].AddObject(detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))]);
		detections[obj_batch][(existing - (obj_batch * obj_limit))].AddObject(detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))],temp_sparse_reps_grid,temp_sparse_reps_strings,temp_mom0,temp_RAPV,temp_DECPV,temp_ref_spec,temp_obj_spec,temp_vfield);
		
		// flag the `object' values within the bounding box to merge it with the base object
		sx_start = detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetRAmin() - chunk_x_start;
		if(sx_start < 0){ sx_start = 0; }
		sy_start = detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetDECmin() - chunk_y_start;
		if(sy_start < 0){ sy_start = 0; }
		sz_start = detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetFREQmin() - chunk_z_start;
		if(sz_start < 0){ sz_start = 0; }
		for(sz = sz_start; (sz <= detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetFREQmax() - chunk_z_start) && (sz < size_z); ++sz){
		  
		  for(sy = sy_start; (sy <= detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetDECmax() - chunk_y_start) && (sy < size_y); ++sy){
		    
		    for(sx = sx_start; (sx <= detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetRAmax() - chunk_x_start) && (sx < size_x); ++sx){
		      
		      if(flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] == match_init[i]){ flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] = existing; }
		      
		      // for(sx = x_start[(match_init[i] - (obj_batch_2 * obj_limit))]; sx <= x_finish[(match_init[i] - (obj_batch_2 * obj_limit))]; ++sx)
		    }
		    
		    // for(sy = y_start[(match_init[i] - (obj_batch_2 * obj_limit))]; sy <= y_finish[(match_init[i] - (obj_batch_2 * obj_limit))]; ++sy)
		  }
		  
		  // for(sz = z_start[(match_init[i] - (obj_batch_2 * obj_limit))]; sz <= z_finish[(match_init[i] - (obj_batch_2 * obj_limit))]; ++sz)
		}
			
		// re-initialise the `object' values for the match_init[i] object that has just been merged into the existing object, 
		// and push the obj value (match_init[i]) to the top of the array of values
		detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].ReInit();
		if((detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].Get_srep_update() != 0) && (detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].Get_srep_size(0) >= 0)){
		  
		  detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].ReInit_srep();
		  detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].ReInit_mini();

		}
		detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].ReInit_size();
		detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].Set_srep_update(0);

		j = -1;
		for(g = 0; g < NO_obj_ids; ++g){ if(obj_ids[g] == match_init[i]){ j = 1; break; } }
		if(j == -1){
		  
		  ++NO_obj_ids;
		  obj_ids.push_back(match_init[i]);
		  		  
		}
		
		// for(i = 0; i < NOi; ++i)
	      }
	      
	      // if(NOi > 1)
	    }
	    
	    // else . . . if(existing == -1)
	  }
	  
	  // if(flag_vals[((z*size_x*size_y) + (y*size_x) + x)] == -1)
	}
	
	// update progress on display
	while(progress <= (((float) ((z * size_x * size_y) + (y * size_x) + x + 1)) / ((float) (size_x * size_y * size_z)))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
	
	// for(x = 0; x < size_x; ++x)
      }
      
      // for(y = 0; y < size_y; ++y)
    }
    
    // for(z = 0; z < size_z; ++z)
  }
  std::cout << "* done." << std::endl;
  
  // 4. generate sparse representations of coherent objects that are retained after
  // this `chunk' of the datacube is processed
  std::cout << "Generating/updating sparse representations of sources . . . " << std::endl;
  std::cout << "0 | |:| | : | |:| | 100% complete" << std::endl;
  progress = 0.0;
  for(i = 0; i < obj; ++i){
    
    // calculate the obj_batch value for this object
    obj_batch = floorf(((float) i / (float) obj_limit));
    
    // move on if this object has been re-initialised
    if(detections[obj_batch][(i - (obj_batch * obj_limit))].ShowVoxels() < 1){ 
      
      // update progress on display
      while(progress <= (((float) (i + 1) / ((float) obj)))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }  
      continue; 
      
    }
    
    // move on if this object is not within the boundaries of the current chunk, in other words,
    // move on if this object doesn't need to have it's sparse representation and postage stamp images updated
    if((detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmax() < chunk_x_start) || (detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmax() < chunk_y_start) || (detections[obj_batch][(i - (obj_batch * obj_limit))].GetFREQmax() < chunk_z_start) || (detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmin() >= (chunk_x_start + size_x)) || (detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmin() >= (chunk_y_start + size_y)) || (detections[obj_batch][(i - (obj_batch * obj_limit))].GetFREQmin() >= (chunk_z_start + size_z))){ 

      // update progress on display
      while(progress <= (((float) (i + 1) / ((float) obj)))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }  
      continue; 
      
    }
    
    // move on if this object's sparse representation and postage stamp images don't need updating
    if((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) >= 0) && (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_update() == -1)){
      
      // update progress on display
      while(progress <= (((float) (i + 1) / ((float) obj)))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }  
      continue; 
	
    }
          
    // update the sparse representation for this object
    
    // 0. if there is no existing sparse representation and postage stamp images, create them from scratch, 
    // otherwise concatenate the existing sparse representation and postage stamp images with the updates
    if((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) < 0) && (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_update() >= 0)){
      
      // branch for creating a brand new sparse representation
      
      // update the sparse_reps_update array
      detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_update(-1);

      // initialise the number of object sections that make up this object
      NOi = 0;   
      
      // calculate the edge of the bounding box in `chunk' co-ordinates
      sx_start = detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmin() - chunk_x_start;
      if(sx_start < 0){ sx_start = 0; }
      sx_finish = detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmax() - chunk_x_start;
      sy_start = detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmin() - chunk_y_start;
      if(sy_start < 0){ sy_start = 0; }
      sy_finish = detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmax() - chunk_y_start;
      sz_start = detections[obj_batch][(i - (obj_batch * obj_limit))].GetFREQmin() - chunk_z_start;
      if(sz_start < 0){ sz_start = 0; }
      sz_finish = detections[obj_batch][(i - (obj_batch * obj_limit))].GetFREQmax() - chunk_z_start;
      
      // count the number of object sections that make up this object
      for(sy = sy_start; (sy <= sy_finish) && (sy < size_y); ++sy){
	
	for(sx = sx_start; (sx <= sx_finish) && (sx < size_x); ++sx){
	  
	  // initialise the dummy integer j
	  j = -1;
	  
	  for(sz = sz_start; (sz <= sz_finish) && (sz < size_z); ++sz){
	    
	    // change j to reflect if this voxel in the flag_vals array belongs to the source
	    if(flag_vals[((sz * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))] == i){ 
	      
	      // if previously this was not part of an object string, j == -1, then increment the number of object strings and
	      // adjust the value of j
	      if(j == -1){ ++NOi; }
	      
	      // change j to reflect that this voxel belongs to the object
	      j = 1; 
	      
	    } else { j = -1; }
	    
	    // for(sz = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetFREQmin() - chunk_z_size; sz <= (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetFREQmax() - chunk_z_size); ++sz)
	  }
	  
	  // for(sx = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmin() - chunk_x_size; sx <= (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmax() - chunk_x_size); ++sx)
	}
	
	// for(sy = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmin() - chunk_y_size; sy <= (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmax() - chunk_y_size); ++sy)
      }

      // create an array element to store this sparse representation
      detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_size(0,detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmin());
      detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_size(1,detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmax());
      detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_size(2,detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmin());
      detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_size(3,detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmax());
      detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_size(4,detections[obj_batch][(i - (obj_batch * obj_limit))].GetFREQmin());
      detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_size(5,detections[obj_batch][(i - (obj_batch * obj_limit))].GetFREQmax());
      detections[obj_batch][(i - (obj_batch * obj_limit))].Create_srep_grid((((detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmin() + 1) * (detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmin() + 1)) + 1));
      detections[obj_batch][(i - (obj_batch * obj_limit))].Create_srep_strings((2 * NOi));
      
      // calculate the extent of the sparse_rep volume
      sx_finish = detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmin() + 1;
      sy_finish = detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmin() + 1;
      sz_finish = detections[obj_batch][(i - (obj_batch * obj_limit))].GetFREQmax() - detections[obj_batch][(i - (obj_batch * obj_limit))].GetFREQmin() + 1;
      
      // create arrays to store the postage stamp images and initialise them
      detections[obj_batch][(i - (obj_batch * obj_limit))].Create_mom0((sx_finish * sy_finish));
      detections[obj_batch][(i - (obj_batch * obj_limit))].Create_RAPV((sx_finish * sz_finish));
      detections[obj_batch][(i - (obj_batch * obj_limit))].Create_DECPV((sy_finish * sz_finish));
      detections[obj_batch][(i - (obj_batch * obj_limit))].Create_ospec(sz_finish);
      if(sz_finish >= 10){

	detections[obj_batch][(i - (obj_batch * obj_limit))].Create_rspec((2 * sz_finish));

      } else {

	detections[obj_batch][(i - (obj_batch * obj_limit))].Create_rspec((10 + sz_finish));

      }
      detections[obj_batch][(i - (obj_batch * obj_limit))].Create_vfield((sx_finish * sy_finish));
      
      // initialise the mini_mom0 values
      for(sy = 0; sy < sy_finish; ++sy){
	for(sx = 0; sx < sx_finish; ++sx){
	  detections[obj_batch][(i - (obj_batch * obj_limit))].Set_mom0(((sy * sx_finish) + sx),0.0);
	}
      }

      // initialise the mini_RAPV values
      for(sz = 0; sz < sz_finish; ++sz){
	for(sx = 0; sx < sx_finish; ++sx){
	  detections[obj_batch][(i - (obj_batch * obj_limit))].Set_RAPV(((sz * sx_finish) + sx),0.0);
	}
      }
      
      // initialise the mini_DECPV values
      for(sz = 0; sz < sz_finish; ++sz){
	for(sy = 0; sy < sy_finish; ++sy){
	  detections[obj_batch][(i - (obj_batch * obj_limit))].Set_DECPV(((sz * sy_finish) + sy),0.0);
	}
      }
      
      // initialise the mini_spec values
      for(sz = 0; sz < sz_finish; ++sz){ detections[obj_batch][(i - (obj_batch * obj_limit))].Set_ospec(sz,0.0); }
      if(sz_finish >= 10){

	for(sz = 0; sz < (2 * sz_finish); ++sz){ detections[obj_batch][(i - (obj_batch * obj_limit))].Set_rspec(sz,0.0); }

      } else {

	for(sz = 0; sz < (10 + sz_finish); ++sz){ detections[obj_batch][(i - (obj_batch * obj_limit))].Set_rspec(sz,0.0); }

      }
            
      // initialise the mini_vfield values
      for(sy = 0; sy < sy_finish; ++sy){
	for(sx = 0; sx < sx_finish; ++sx){
	  detections[obj_batch][(i - (obj_batch * obj_limit))].Set_vfield(((sy * sx_finish) + sx),0.0);
	}
      }

      // populate the sparse array representation and postage stamp images
      NOi = 0;
      for(sy = 0; sy < sy_finish; ++sy){
	
	for(sx = 0; sx < sx_finish; ++sx){
	  
	  // write the current number of object strings to the grid
	  detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_grid(((sy * sx_finish) + sx),NOi);
	  
	  // initialise the dummy integer j 
	  j = -1;
	  
	  for(sz = 0; sz < sz_finish; ++sz){
	    
	    // change j to reflect if this voxel in the flag_vals array belongs to the source
	    if(flag_vals[(((sz + sz_start) * data_metric[2]) + ((sy + sy_start) * data_metric[1]) + ((sx + sx_start) * data_metric[0]))] == i){ 
	      
	      // update mini_mom0
	      detections[obj_batch][(i - (obj_batch * obj_limit))].Add_mom0(((sy * sx_finish) + sx),data_vals[(((sz + sz_start) * data_metric[2]) + ((sy + sy_start) * data_metric[1]) + ((sx + sx_start) * data_metric[0]))]);
	      
	      // update mini_RAPV
	      detections[obj_batch][(i - (obj_batch * obj_limit))].Add_RAPV(((sz * sx_finish) + sx),data_vals[(((sz + sz_start) * data_metric[2]) + ((sy + sy_start) * data_metric[1]) + ((sx + sx_start) * data_metric[0]))]);
	      
	      // update mini_DECPV
	      detections[obj_batch][(i - (obj_batch * obj_limit))].Add_DECPV(((sz * sy_finish) + sy),data_vals[(((sz + sz_start) * data_metric[2]) + ((sy + sy_start) * data_metric[1]) + ((sx + sx_start) * data_metric[0]))]);
	      
	      // update mini_obj_spec
	      detections[obj_batch][(i - (obj_batch * obj_limit))].Add_ospec(sz,data_vals[(((sz + sz_start) * data_metric[2]) + ((sy + sy_start) * data_metric[1]) + ((sx + sx_start) * data_metric[0]))]);
  
	      // update mini_vfield
	      detections[obj_batch][(i - (obj_batch * obj_limit))].Add_vfield(((sy * sx_finish) + sx),(((double) (sz + sz_start + chunk_z_start)) * data_vals[(((sz + sz_start) * data_metric[2]) + ((sy + sy_start) * data_metric[1]) + ((sx + sx_start) * data_metric[0]))]));
	      
	      // if previously this was not part of an object string, j == -1, then increment the number of object strings and
	      // write the position to the sparse_reps array
	      if(j == -1){ 
		
		// write the beginning of the object string to sparse_reps
		detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_strings((2 * NOi),(sz_start + chunk_z_start + sz));
		
		// increment NOi
		++NOi; 
		
		// initialise the end of the object string in sparse_reps
		detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_strings(((2 * NOi) - 1),(sz_start + chunk_z_start + sz));

	      } else {
		
		// update the end of the object string in sparse_reps
		detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_strings(((2 * NOi) - 1),(sz_start + chunk_z_start + sz));
		
	      }
	      
	      // change j to reflect that this voxel belongs to the object
	      j = 1; 
	      
	    } else { 
	      
	      // change j to reflect that this voxel doesn't belong to the object
	      j = -1; 
	      
	    }
	    
	    // for(sz = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetFREmin() - chunk_z_start + 1; sz <= (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetFREQmax() - chunk_z_start + 1); ++sz)
	  }
	  
	  // for(sx = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmin() - chunk_x_start + 1; sx <= (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmax() - chunk_x_start + 1); ++sx)
	}
	
	//for(sy = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmin() - chunk_y_start + 1; sy <= (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmax() - chunk_y_start + 1); ++sy)
      }
            
      // update sparse_reps with the final value of NOi
      detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_grid((sx_finish * sy_finish),NOi);
            
      // create the reference spectrum
      if(sz_finish >= 10){
	
	for(sz = -1 * (int) floorf((0.5 * (float) sz_finish)); sz < (sz_finish + ((int) floorf((0.5 * (float) sz_finish)))); ++sz){
	  
	  for(sy = 0; sy < sy_finish; ++sy){
	  
	    for(sx = 0; sx < sx_finish; ++sx){
	    
	      if(((sz + sz_start) >= 0) && ((sz + sz_start) < size_z)){
		
		// update mini_ref_spec	
		detections[obj_batch][(i - (obj_batch * obj_limit))].Add_rspec((sz + ((int) floorf((0.5 * (double) sz_finish)))),data_vals[(((sz + sz_start) * data_metric[2]) + ((sy + sy_start) * data_metric[1]) + ((sx + sx_start) * data_metric[0]))]);	
		
	      }
	      
	    }
	    
	  }
      
	}

	// if(sz_finish >= 10)
      } else {

	for(sz = -5; sz < (sz_finish + 5); ++sz){
	  
	  for(sy = 0; sy < sy_finish; ++sy){
	  
	    for(sx = 0; sx < sx_finish; ++sx){
	    
	      if(((sz + sz_start) >= 0) && ((sz + sz_start) < size_z)){
		
		// update mini_ref_spec
		detections[obj_batch][(i - (obj_batch * obj_limit))].Add_rspec((sz + 5),data_vals[(((sz + sz_start) * data_metric[2]) + ((sy + sy_start) * data_metric[1]) + ((sx + sx_start) * data_metric[0]))]);	
		
	      }
	      
	    }
	    
	  }
      
	}

	// else . . . if(sz_finish >= 10)
      }
      
     // if(sparse_reps[obj_batch][(obj - (obj_batch * obj_limit))][0] < 0) . . . else
    } else {
            
      // branch for updating an existing sparse representation
      
      // count the number of object strings comprising the object
      
      // 1. initialise the number of object strings
      NOi = 0;
      
      // 2. convert the new bounding box of the object into `chunk' co-ordinates
      sx_start = detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmin() - chunk_x_start;
      sx_finish = detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmax() - chunk_x_start;
      sy_start = detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmin() - chunk_y_start;
      sy_finish = detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmax() - chunk_y_start;
      sz_start = detections[obj_batch][(i - (obj_batch * obj_limit))].GetFREQmin() - chunk_z_start;
      if(sz_start < 0){ sz_start = 0; }
      sz_finish = detections[obj_batch][(i - (obj_batch * obj_limit))].GetFREQmax() - chunk_z_start;
      
      // 3. search through the bounding box, and do one of the following
      // 3.a if inside the chunk and outside the existing sparse rep grid, run the routine to count the number
      //     of object strings along the current line of sight
      // 3.b otherwise use the grid value to increment NOi
      for(sy = sy_start; sy <= sy_finish; ++sy){
	
	for(sx = sx_start; sx <= sx_finish; ++sx){
	  	  
	  // case 3.a triggered
	  if(((sx > merge_x) || (chunk_x_start == 0)) && ((sy > merge_y) || (chunk_y_start == 0))){
	    	    
	    // initialise the dummy integer j 
	    j = -1;
	    
	    for(sz = sz_start; sz <= sz_finish; ++sz){
	      
	      // change i to reflect if this voxel in the flag_vals array belongs to the source
	      if(flag_vals[((sz * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))] == i){ 
		
		// if previously this was not part of an object string, i == -1, then increment the number of object strings and
		// write the position to the sparse_reps array
		if(j == -1){ 
		  
		  // increment NOi
		  ++NOi; 
		  
		} 
		
		// change j to reflect that this voxel belongs to the object
		j = 1; 
		
	      } else { j = -1; }
	      
	      // for(sz = sz_start; sz <= sz_finish; ++sz)
	    }
	    
	  } else {
	    
	    // case 3.b triggered  
	    NOi+=(detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_grid((((sy + chunk_y_start - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx + chunk_x_start - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_grid((((sy + chunk_y_start - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx + chunk_x_start - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0))));
	    
	  }     
	  
	  // for(sx = sx_start; sx <= sx_finish; ++sx)
	}
	
	// for(sy = sy_start; sy <= sy_finish; ++sy)
      }
      
      // 4. allocate memory to store temporary sparse representation and postage stamp images, then
      //    initialise the postage stamp image arrays
      temp_mom0.resize(0);
      for(sy = 0; sy < (sy_finish - sy_start + 1); ++sy){
	for(sx = 0; sx < (sx_finish - sx_start + 1); ++sx){
	  temp_mom0.push_back(0.0);
	}
      }
      temp_RAPV.resize(0);
      for(sz = 0; sz < (sz_finish - sz_start + 1); ++sz){
	for(sx = 0; sx < (sx_finish - sx_start + 1); ++sx){
	  temp_RAPV.push_back(0.0);
	}
      }
      temp_DECPV.resize(0);
      for(sz = 0; sz < (sz_finish - sz_start + 1); ++sz){
	for(sy = 0; sy < (sy_finish - sy_start + 1); ++sy){
	  temp_DECPV.push_back(0.0);
	}
      }
      temp_obj_spec.resize(0);
      for(sz = 0; sz < (sz_finish - sz_start + 1); ++sz){ temp_obj_spec.push_back(0.0); }
      if((sz_finish - sz_start + 1) >= 10){

	temp_ref_spec.resize(0);
	for(sz = 0; sz < (2 * (sz_finish - sz_start + 1)); ++sz){ temp_ref_spec.push_back(0.0); }

      } else {

	temp_ref_spec.resize(0);
	for(sz = 0; sz < (sz_finish - sz_start + 11); ++sz){ temp_ref_spec.push_back(0.0); }

      }
      temp_vfield.resize(0);
      for(sy = 0; sy < (sy_finish - sy_start + 1); ++sy){
	for(sx = 0; sx < (sx_finish - sx_start + 1); ++sx){
	  temp_vfield.push_back(0.0);
	}
      }
     
      // 5. populate temp_sparse_rep arrays
      // 5.a if inside the chunk and outside the existing sparse rep grid, run the routine to count the number
      //     of object strings along the current line of sight, populate new strings array and add values to new 
      //     postage stamp images
      // 5.b if outside of the chunk, use the grid value to increment NOi, populate strings array from previous version and 
      //     fill in valuese for postage stamp images

      // initialise number of object strings
      NOi = 0;
      
      // initialise RA PV image
      for(sz = 0; sz < (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1); ++sz){

	for(sx = 0; sx < (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1); ++sx){
	
	  temp_RAPV[(((sz + detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) - sz_start - chunk_z_start) * (sx_finish - sx_start + 1)) + sx + detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) - sx_start - chunk_x_start)]+=detections[obj_batch][(i - (obj_batch * obj_limit))].Get_RAPV(((sz * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx));
	  
	}

      }

      // initialise Dec PV image
      for(sz = 0; sz < (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1); ++sz){

	for(sy = 0; sy < (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2) + 1); ++sy){
	
	  temp_DECPV[(((sz + detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) - sz_start - chunk_z_start) * (sy_finish - sy_start + 1)) + sy + detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2) - sy_start - chunk_y_start)]+=detections[obj_batch][(i - (obj_batch * obj_limit))].Get_DECPV(((sz * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2) + 1)) + sy));
	  
	}

      }

      // initialise mini spec's
      for(sz = 0; sz < (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1); ++sz){ temp_obj_spec[(sz + detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) - sz_start - chunk_z_start)]+=detections[obj_batch][(i - (obj_batch * obj_limit))].Get_ospec(sz); }

      if((sz_finish - sz_start + 1) >= 10){

	j = sz_start + chunk_z_start - ((int) floorf((0.5 * (float) (sz_finish - sz_start + 1))));

	if((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1) >= 10){

	  k = detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) - ((int) floorf((0.5 * (double) (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1))));

	  for(sz = 0; sz < (2 * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); ++sz){

	    temp_ref_spec[(sz + k - j)]+=detections[obj_batch][(i - (obj_batch * obj_limit))].Get_rspec(sz);

	  }

	} else {
	  
	  k = detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) - 5;

	  for(sz = 0; sz < (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 11); ++sz){
	    
	    temp_ref_spec[(sz + k - j)]+=detections[obj_batch][(i - (obj_batch * obj_limit))].Get_rspec(sz);

	  }
	  
	}

      } else {

	j = sz_start + chunk_z_start - 5;

	if((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1) >= 10){

	  k = detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) - ((int) floorf((0.5 * (double) (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1))));

	  for(sz = 0; sz < (2 * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); ++sz){

	    temp_ref_spec[(sz + k - j)]+=detections[obj_batch][(i - (obj_batch * obj_limit))].Get_rspec(sz);

	  }

	} else {
	  
	  k = detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) - 5;

	  for(sz = 0; sz < (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 11); ++sz){
	    
	    temp_ref_spec[(sz + k - j)]+=detections[obj_batch][(i - (obj_batch * obj_limit))].Get_rspec(sz);

	  }
	  
	}

      }
      
      // search through LoS's of updated bounding box and take action as dictated by case 5.a or 5.b
      for(sy = sy_start; sy <= sy_finish; ++sy){
	
	for(sx = sx_start; sx <= sx_finish; ++sx){
	  	  
	  // add value to temp_sparse_rep_grid
	  temp_sparse_reps_grid[(((sy - sy_start) * (sx_finish - sx_start + 1)) + sx - sx_start)] = NOi;
	  
	  // case 5.a triggered for moment-0 image and sparse representation
	  if(((sx > merge_x) || (chunk_x_start == 0)) && ((sy > merge_y) || (chunk_y_start == 0))){
	    
	    // initialise the dummy integer j 
	    j = -1;
	    
	    for(sz = sz_start; sz <= sz_finish; ++sz){
	      
	      // change j to reflect if this voxel in the flag_vals array belongs to the source
	      if(flag_vals[((sz * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))] == i){ 
		
		// update mini_mom0
		temp_mom0[(((sy - sy_start) * (sx_finish - sx_start + 1)) + sx - sx_start)]+=data_vals[((sz * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))];
		
		// update mini_RAPV
		temp_RAPV[(((sz - sz_start) * (sx_finish - sx_start + 1)) + sx - sx_start)]+=data_vals[((sz * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))];
		
		// update mini_DECPV
		temp_DECPV[(((sz - sz_start) * (sy_finish - sy_start + 1)) + sy - sy_start)]+=data_vals[((sz * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))];
		
		// update mini_obj_spec
		temp_obj_spec[(sz - sz_start)]+=data_vals[((sz * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))];
		
		// update mini_vfield
		temp_vfield[(((sy - sy_start) * (sx_finish - sx_start + 1)) + sx - sx_start)]+=(((double) (sz + chunk_z_start)) * data_vals[((sz * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))]);
		
		// if previously this was not part of an object string, j == -1, then increment the number of object strings and
		// write the position to the sparse_reps array
		if(j == -1){ 
		  
		  // write the beginning of the object string to sparse_reps
		  temp_sparse_reps_strings[(2 * NOi)] = chunk_z_start + sz;
		  
		  // increment NOi
		  ++NOi; 
		  		  
		  // initialise the end of the object string in sparse_reps
		  temp_sparse_reps_strings[((2 * NOi) - 1)] = chunk_z_start + sz;
		  
		} else {
		  
		  // update the end of the object string in sparse_reps
		  temp_sparse_reps_strings[((2 * NOi) - 1)] = chunk_z_start + sz;
		  
		}
		
		// change j to reflect that this voxel belongs to the object
		j = 1; 
		
	      } else { 
		
		// change j to reflect that this voxel doesn't belong to the object
		j = -1; 
		
	      }
	      
	      // for(sz = sz_start; sz <= sz_finish; ++sz)
	    }
	    
	    // add values to temp_ref_spec
	    if((sz_finish - sz_start + 1) >= 10){
	      
	      k = (int) floorf((0.5 * (float) (sz_finish - sz_start + 1)));
	      
	      for(sz = 0; sz < (2 * (sz_finish - sz_start + 1)); ++sz){
	  	
		if(((sz + sz_start - k) >= 0) && ((sz + sz_start - k) < size_z)){
		  
		  // update mini_ref_spec
		  temp_ref_spec[(sz + sz_start - k)]+=data_vals[(((sz + sz_start - k) * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))];	
		  
		}
		
	      }
	      
	      // if(sz_finish >= 10)
	    } else {
	      
	      for(sz = 0; sz < (sz_finish - sz_start + 11); ++sz){
		
		if(((sz + sz_start - 5) >= 0) && ((sz + sz_start - 5) < size_z)){
		  
		  // update mini_ref_spec
		  temp_ref_spec[(sz + sz_start - 5)]+=data_vals[(((sz + sz_start - 5) * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))];	
		  
		}
		
	      }
	      
	      // else . . . if(sz_finish >= 10)
	    }
	    
	  } else {
	    // case 5.b triggered  
	    	    
	    // update temp_sparse_rep_strings
	    j = 0;
	    if(((sx + chunk_x_start) >= detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0)) && ((sx + chunk_x_start) <= (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) + 1)) && ((sy + chunk_y_start) >= detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2)) && ((sy + chunk_y_start) <= (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(3) + 1))){
	      
	      for(g = detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_grid((((sy + chunk_y_start - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx + chunk_x_start - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0))); g < detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_grid((((sy + chunk_y_start - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx + chunk_x_start - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1)); ++g){
		
		temp_sparse_reps_strings[(NOi + j)] = detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_strings((2 * g));
		temp_sparse_reps_strings[(NOi + j + 1)] = detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_strings(((2 * g) + 1));
		++j;
		
	      }
	      
	      // update moment 0 postage stamp image
	      temp_mom0[(((sy - sy_start) * (detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmin() + 1)) + sx - sx_start)] = detections[obj_batch][(i - (obj_batch * obj_limit))].Get_mom0((((sy + chunk_y_start - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx + chunk_x_start - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0)));
	      
	      // update vfield postage stamp image
	      temp_vfield[(((sy - sy_start) * (detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmin() + 1)) + sx - sx_start)] = detections[obj_batch][(i - (obj_batch * obj_limit))].Get_vfield((((sy + chunk_y_start - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx + chunk_x_start - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0)));

	      // update NOi
	      NOi+=(detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_grid((((sy + chunk_y_start - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx + chunk_x_start - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_grid((((sy + chunk_y_start - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx + chunk_x_start - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0))));
	      	      
	    }
	    
	    // else . . . case 5.b
	  }
	    
	  // for(sx = sx_start; sx <= sx_finish; ++sx)
	}
		
	// for(sy = sy_start; sy <= sy_finish; ++sy)
      }
            
      temp_sparse_reps_grid[((sx_finish - sx_start + 1) * (sy_finish - sy_start + 1))] = NOi;

      // 6. replace the current sparse representation with the new sparse representation
      detections[obj_batch][(i - (obj_batch * obj_limit))].ReInit_srep();
      detections[obj_batch][(i - (obj_batch * obj_limit))].ReInit_mini();
 
      // update the sparse_reps_size array with the new bounding box
      detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_size(0,detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmin());
      detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_size(1,detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmax());
      detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_size(2,detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmin());
      detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_size(3,detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmax());
      detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_size(4,detections[obj_batch][(i - (obj_batch * obj_limit))].GetFREQmin());
      detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_size(5,detections[obj_batch][(i - (obj_batch * obj_limit))].GetFREQmax());

      // update the sparse_reps_update array
      detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_update(-1);
          
      // populate the new sparse_reps and postage stamps
      detections[obj_batch][(i - (obj_batch * obj_limit))].Create_srep_grid((1 + ((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2) + 1))));
      detections[obj_batch][(i - (obj_batch * obj_limit))].Create_srep_strings((2 * NOi));
      detections[obj_batch][(i - (obj_batch * obj_limit))].Create_mom0(((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2) + 1)));
      detections[obj_batch][(i - (obj_batch * obj_limit))].Create_RAPV(((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1)));
      detections[obj_batch][(i - (obj_batch * obj_limit))].Create_DECPV(((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2) + 1) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1)));
      detections[obj_batch][(i - (obj_batch * obj_limit))].Create_vfield(((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2) + 1)));
  
      detections[obj_batch][(i - (obj_batch * obj_limit))].Create_ospec((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1));
      if((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1) >= 10){
       
	detections[obj_batch][(i - (obj_batch * obj_limit))].Create_rspec((2 * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1)));
	
      } else {
	
	detections[obj_batch][(i - (obj_batch * obj_limit))].Create_rspec((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 11));

      }
   
      for(j = 0; j < (1 + ((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2) + 1))); ++j){ detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_grid(j,temp_sparse_reps_grid[j]); }
      for(j = 0; j < (2 * NOi); ++j){ detections[obj_batch][(i - (obj_batch * obj_limit))].Set_srep_strings(j,temp_sparse_reps_strings[j]); }
      for(j = 0; j < ((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); ++j){ detections[obj_batch][(i - (obj_batch * obj_limit))].Set_mom0(j,temp_mom0[j]); }
      for(j = 0; j < ((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); ++j){ detections[obj_batch][(i - (obj_batch * obj_limit))].Set_RAPV(j,temp_RAPV[j]); }
      for(j = 0; j < ((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2) + 1) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); ++j){ detections[obj_batch][(i - (obj_batch * obj_limit))].Set_DECPV(j,temp_DECPV[j]); }
      for(j = 0; j < ((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); ++j){ detections[obj_batch][(i - (obj_batch * obj_limit))].Set_vfield(j,temp_vfield[j]); }

      for(j = 0; j < (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1); ++j){ detections[obj_batch][(i - (obj_batch * obj_limit))].Set_ospec(j,temp_obj_spec[j]); }
      if((detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1) >= 10){

	for(j = 0; j < (2 * (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); ++j){ detections[obj_batch][(i - (obj_batch * obj_limit))].Set_rspec(j,temp_ref_spec[j]); }
	
      } else {

	for(j = 0; j < (detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(i - (obj_batch * obj_limit))].Get_srep_size(4) + 11); ++j){ detections[obj_batch][(i - (obj_batch * obj_limit))].Set_rspec(j,temp_ref_spec[j]); }

      }
      
      // else . . . if(sparse_reps[obj_batch][(obj - (obj_batch * obj_limit))][0] < 0)
    }
    
    // update progress on display
    while(progress <= (((float) (i + 1) / ((float) obj)))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }  
   
    // for(i = 0; i < obj; ++i)
  }
  std::cout << "* done." << std::endl;
    
  // 5. return the number of coherent objects 
  return obj;
  
}

long int CreateObjects(double * data_vals, long int * flag_vals, int size_x, int size_y, int size_z, int chunk_x_start, int chunk_y_start, int chunk_z_start, int merge_x, int merge_y, int merge_z, int min_x_size, int min_y_size, int min_z_size, int min_v_size, double intens_thresh_min, double intens_thresh_max, long int flag_value, long int start_obj, vector<object_props_dbl *> & detections, vector<long int> & obj_ids, vector<long int> & check_obj_ids, int obj_limit, int max_x_val, int max_y_val, int max_z_val, int ss_mode, int * data_metric, int * xyz_order){
 
  long int o,obj,existing,obj_batch,obj_batch_2;
  int x,y,z,sx,sy,sz,sx_start,sy_start,sz_start,sx_finish,sy_finish,sz_finish;
  int g,i,j,k,NO_obj_ids,NO_check_obj_ids,NOi;
  int init_limit,prev,x_start,y_start,temp_x[3],temp_y[3],temp_z[3];
  float progress;
  vector<double> temp_mom0, temp_RAPV, temp_DECPV, temp_obj_spec, temp_ref_spec, temp_vfield;
  vector<long int> match_init;
  vector<int> temp_sparse_reps_grid, temp_sparse_reps_strings;

  // convert the merging lengths from distances (in pixels/voxels) to empty space between linked detections (in pixels/voxels) --- subtract 1 --- this is a more naturual
  // definition when processing a large file in chunks, because it makes it easier to build up sources that need to be retained for the rest of the input file
  --merge_x;
  --merge_y;
  --merge_z;

  // reorder the datacube and subcube limits to be in x,y,z order
  temp_x[0] = chunk_x_start;
  temp_y[0] = chunk_y_start;
  temp_z[0] = chunk_z_start;
  temp_x[1] = size_x;
  temp_y[1] = size_y;
  temp_z[1] = size_z;
  temp_x[2] = max_x_val;
  temp_y[2] = max_y_val;
  temp_z[2] = max_z_val;
  switch(xyz_order[0]){
  case 1:
    chunk_x_start = temp_x[0];
    size_x = temp_x[1];
    max_x_val = temp_x[2];
    break;
  case 2:
    chunk_x_start = temp_y[0];
    size_x = temp_y[1];
    max_x_val = temp_y[2];
    break;
  case 3:
    chunk_x_start = temp_z[0];
    size_x = temp_z[1];
    max_x_val = temp_z[2];
    break;
  default:
    chunk_x_start = temp_x[0];
    size_x = temp_x[1];
    max_x_val = temp_x[2];
    break;
  }
  switch(xyz_order[1]){
  case 1:
    chunk_y_start = temp_x[0];
    size_y = temp_x[1];
    max_y_val = temp_x[2];
    break;
  case 2:
    chunk_y_start = temp_y[0];
    size_y = temp_y[1];
    max_y_val = temp_y[2];
    break;
  case 3:
    chunk_y_start = temp_z[0];
    size_y = temp_z[1];
    max_y_val = temp_z[2];
    break;
  default:
    chunk_y_start = temp_y[0];
    size_y = temp_y[1];
    max_y_val = temp_y[2];
    break;
  }
  switch(xyz_order[2]){
  case 1:
    chunk_z_start = temp_x[0];
    size_z = temp_x[1];
    max_z_val = temp_x[2];
    break;
  case 2:
    chunk_z_start = temp_y[0];
    size_z = temp_y[1];
    max_z_val = temp_y[2];
    break;
  case 3:
    chunk_z_start = temp_z[0];
    size_z = temp_z[1];
    max_z_val = temp_z[2];
    break;
  default:
    chunk_z_start = temp_z[0];
    size_z = temp_z[1];
    max_z_val = temp_z[2];
    break;
  }

  // create temporary arrays
  if(10000 > (((2 * merge_x) + 3) * ((2 * merge_y) + 3) * ((2 * merge_z) + 3))){
    
    init_limit = 10000;
    
  } else {
    
    init_limit = ((2 * merge_x) + 3) * ((2 * merge_y) + 3) * ((2 * merge_z) + 3) + 10;
    
  }
  match_init.reserve(init_limit);
  match_init.resize(0);

  temp_sparse_reps_grid.reserve(1000000);
  temp_sparse_reps_grid.resize(0);
  temp_sparse_reps_strings.reserve(1000000);
  temp_sparse_reps_strings.resize(0);
  temp_mom0.reserve(1000000);
  temp_mom0.resize(0);
  temp_RAPV.reserve(1000000);
  temp_RAPV.resize(0);
  temp_DECPV.reserve(1000000);
  temp_DECPV.resize(0);
  temp_obj_spec.reserve(1000000);
  temp_obj_spec.resize(0);
  temp_ref_spec.reserve(1000000);
  temp_ref_spec.resize(0);
  temp_vfield.reserve(1000000);
  temp_vfield.resize(0);
 
  // 0. initialise variables and arrays
  obj = start_obj;
  x_start = 0;
  y_start = 0;
  if(chunk_x_start > 0){ x_start = merge_x + 1; }
  if(chunk_y_start > 0){ y_start = merge_y + 1; }
  NO_obj_ids = obj_ids.size();
  NO_check_obj_ids = check_obj_ids.size();

  // 1. Create list of `coherent' objects from neighbouring voxels
  
  // for each grid point, check if it has been flagged as a source voxel, and assign
  // it an object ID if it is. If a neighbouring voxel has already been flagged, assign 
  // that ID to the grid point, otherwise, assign current value of obj and increment obj.
  std::cout << "0 | |:| | : | |:| | 100% complete" << std::endl;
  progress = 0.0;
  for(z = 0; ((z < size_z) && ((z + chunk_z_start) < max_z_val)); ++z){

    // if the z value is sufficiently large that objects have started to pop out of the merging box,
    // then check if the objects outside of the merging box are sufficiently large
    if(z > (merge_z + 1)){
      
      // initialise NOi to store the number of current check_obj_ids values that haven't moved outside of 
      // the merging volume, and that form the basis for the new check_obj_ids list
      NOi = 0;
      
      // apply size threshold to objects added since last size thresholding
      for(i = 0; i < NO_check_obj_ids; ++i){
	
	// calculate the obj_batch value for the existing object
	obj_batch = (long int) floor(((double) check_obj_ids[i] / (double) obj_limit));
	
	// move to the next object if this one has been re-initialised
	if(detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	
	// check if this object has moved outside of the merging box available to the next plane, if it 
	// has apply a size threshold, otherwise use it to seed the new list of objects to be size thresholded
	// the next time
	if(detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmax() < (chunk_z_start + z - merge_z - 1)){
	  
	  // apply the size threshold, and if it fails re-initialise the object and pop its id to the list
	  // of available obj_ids
	  if((((detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetRAmin() + 1) < min_x_size) || ((detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetDECmin() + 1) < min_y_size) || ((detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmax() - detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmin() + 1) < min_z_size) || (detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].ShowVoxels() < min_v_size)) && (((detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetRAmin() - chunk_x_start) > merge_x) || (chunk_x_start == 0)) && ((((chunk_x_start + size_x - 1 - detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetRAmax()) > merge_x) || ((chunk_x_start + size_x) >= max_x_val)) && (((detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetDECmin() - chunk_y_start) > merge_y) || (chunk_y_start == 0)) && (((chunk_y_start + size_y - 1 - detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetDECmax()) > merge_y) || ((chunk_y_start + size_y) >= max_y_val)) && (((detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmin() - chunk_z_start) > merge_z) || (chunk_z_start == 0)) && (((chunk_z_start + size_z - 1 - detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmax()) > merge_z) || ((chunk_z_start + size_z) >= max_z_val)))){
	    
	    // re-initialise object
	    detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].ReInit();
	    if((detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].Get_srep_update() != 0) && (detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].Get_srep_size(0) >= 0)){
	      
	      detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].ReInit_srep();
	      detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].ReInit_mini();

	    }
	    detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].ReInit_size();
	    detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].Set_srep_update(0);

	    // flag the `object' values within the bounding box and remove it from the flag_vals array
	    sx_start = detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetRAmin() - chunk_x_start;
	    if(sx_start < 0){ sx_start = 0; }
	    sy_start = detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetDECmin() - chunk_y_start;
	    if(sy_start < 0){ sy_start = 0; }
	    sz_start = detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmin() - chunk_z_start;
	    if(sz_start < 0){ sz_start = 0; }
	    for(sz = sz_start; ((sz <= (detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmax() - chunk_z_start)) && (sz < size_z)); ++sz){
	      
	      for(sy = sy_start; ((sy <= (detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetDECmax() - chunk_y_start)) && (sy < size_y)); ++sy){
		
		for(sx = sx_start; ((sx <= (detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetRAmax() - chunk_x_start)) && (sx < size_x)); ++sx){
		  
		  if(flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx * data_metric[0]))] == check_obj_ids[i]){ flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx * data_metric[0]))] = -99; }
		  
		  // for(sx = x_start[(match_init[i] - (obj_batch_2 * obj_limit))]; sx <= x_finish[(match_init[i] - (obj_batch_2 * obj_limit))]; ++sx)
		}
		
		// for(sy = y_start[(match_init[i] - (obj_batch_2 * obj_limit))]; sy <= y_finish[(match_init[i] - (obj_batch_2 * obj_limit))]; ++sy)
	      }
	      
	      // for(sz = z_start[(match_init[i] - (obj_batch_2 * obj_limit))]; sz <= z_finish[(match_init[i] - (obj_batch_2 * obj_limit))]; ++sz)
	    }
	    
	    // add object id to list of available ids
	    j = -1;
	    for(g = 0; g < NO_obj_ids; ++g){ if(obj_ids[g] == check_obj_ids[i]){ j = 1; break; } }
	    if(j == -1){
	      ++NO_obj_ids;
	      obj_ids.push_back(check_obj_ids[i]);
	    }

	  }
	  
	  // if(z_finish[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))] < (z - merge_z - 1))
	} else {
	  
	  check_obj_ids[NOi] = check_obj_ids[i];
	  ++NOi;
	  
	  // else . . . if(z_finish[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))] < (z - merge_z - 1))
	}
	
	// for(i = 0; i < NO_check_obj_ids; ++i)
      }

      // update NO_check_obj_ids
      NO_check_obj_ids = NOi;
      check_obj_ids.resize(NOi);
      
      // if(z > (merge_z + 2))
    }
    
    for(y = y_start; (y < size_y) && ((chunk_y_start + y) < max_y_val); ++y){
      
      for(x = x_start; (x < size_x) && ((chunk_x_start + x) < max_x_val); ++x){
	
	// if this is a source voxel, check if it is associated with a previously identified
	// source voxel
	if(flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (x*data_metric[0]))] == flag_value){ 
	  
	  // 1.  initialise variables
	  existing = flag_value;
	  NOi = 0;
	  match_init.resize(0);
	  
	  // 3. initial pass, check the voxels within the merging distance to see if any have been
	  // identified as previous objects
	  
	  // 3a. set limits of region to be searched
	  sx_start = x - merge_x - 1;
	  if(sx_start < 0){ sx_start = 0; }
	  sy_start = y - merge_y - 1;
	  if(sy_start < 0){ sy_start = 0; }
	  sz_start = z - merge_z - 1;
	  if(sz_start < 0){ sz_start = 0; }
	  
	  // 3b. search through previous region and modify existing flag if required

	  if(ss_mode == 1){
	    
	    // use a rectangular paralleliped to do the linking

	    for(sz = sz_start; (sz < z); ++sz){
	    
	      for(sy = sy_start; (sy <= (y + merge_y + 1)) && (sy < size_y); ++sy){
	      
		for(sx = sx_start; (sx <= (x + merge_x + 1)) && (sx < size_x); ++sx){
		
		  // check if this test voxel is source
		  if((flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] == flag_value) || (flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] >= 0)){
		    
		    // check if this is a new objid
		    prev = -1;
		    for(i = 0; (i < NOi) && (i < init_limit); ++i){
		      
		      if(match_init[i] == flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]){
			
			prev = 1;
			break;
			
		      }
		      
		    }
		    if(prev == -1){ 
		      
		      // add current object to list of previous matches
		      match_init.push_back(flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]);
		      ++NOi;

		    }
		    
		    // update existing flag
		    if(((flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] < existing) && (flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] >= 0)) || (existing == flag_value)){ existing = flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]; }
		    
		  }
		  
		  // for(sx = x_start; (sx < (x + merge_x + 1)) && (sx < size_x); ++sx)
		}
		
		// for(sy = y_start; (sy < (y + merge_y + 1)) && (sy < size_y); ++sy)
	      }
	      
	      // for(sz = z_start; (sz < z); ++sz)
	    }
	    for(sy = sy_start; (sy < y); ++sy){
	      
	      for(sx = sx_start; (sx <= (x + merge_x + 1)) && (sx < size_x); ++sx){
		
		// check if this test voxel is source
		if((flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] == flag_value) || (flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] >= 0)){
		  
		  // check if this is a new objid, move on if not
		  prev = -1;
		  for(i = 0; (i < NOi) && (i < init_limit); ++i){
		    
		    if(match_init[i] == flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]){
		      
		      prev = 1;
		      break;
		      
		    }
		    
		  }
		  if(prev == -1){ 
		    
		    // add current object to list of previous matches
		    ++NOi;
		    match_init.push_back(flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]);

		  }
		  
		  // update existing flag
		  if(((flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] < existing) && (flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] >= 0)) || (existing == flag_value)){ existing = flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]; }
		  
		}		
		
		// for(sx = x_start; (sx < = (x + merge_x + 1)) && (sx < size_x); ++sx)
	      }
	      
	      // for(sy = y_start; (sy < y) && (sy < size_y); ++sy)
	    }
	    for(sx = sx_start; sx < x; ++sx){
	      
	      // check if this test voxel is source
	      if((flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))] == flag_value) || (flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))] >= 0)){
		
		// check if this is a new objid, move on if not
		prev = -1;
		for(i = 0; (i < NOi) && (i < init_limit); ++i){
		  
		  if(match_init[i] == flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))]){
		    
		    prev = 1;
		    break;
		    
		  }
		  
		}
		if(prev == -1){ 
		  
		  // add current object to list of previous matches
		  ++NOi;
		  match_init.push_back(flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))]);

		}
		
		// update existing flag
		if(((flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))] < existing) && (flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))] >= 0)) || (existing == flag_value)){ existing = flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))]; }
		
	      }	  
	      
	      // for(sx = x_start; sx <= x; ++sx)
	    }
	  
	    // if(ss_mode == 1)
	  } else {

	    // use a cylinder for the linking lengths

	    for(sz = sz_start; (sz < z); ++sz){
	    
	      for(sy = sy_start; (sy <= (y + merge_y + 1)) && (sy < size_y); ++sy){
	      
		for(sx = sx_start; (sx <= (x + merge_x + 1)) && (sx < size_x); ++sx){
		
		  // check if the test voxel is within the boundaries of the spatial ellipse
		  if(((((float) (x - sx)) * ((float) (x - sx)) / (((float) (merge_x + 1)) * ((float) (merge_x + 1)))) + (((float) (y - sy)) * ((float) (y - sy)) / (((float) (merge_y + 1)) * ((float) (merge_y + 1))))) > 1.0){ continue; }

		  // check if this test voxel is source
		  if((flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] == flag_value) || (flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] >= 0)){

		    // check if this is a new objid
		    prev = -1;
		    for(i = 0; (i < NOi) && (i < init_limit); ++i){
		      
		      if(match_init[i] == flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]){
			
			prev = 1;
			break;
			
		      }
		      
		    }
		    if(prev == -1){ 
		      
		      // add current object to list of previous matches
		      ++NOi;
		      match_init.push_back(flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]);

		    }
		    
		    // update existing flag
		    if(((flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] < existing) && (flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] >= 0)) || (existing == flag_value)){ existing = flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]; }
		    
		  }
		  
		  // for(sx = x_start; (sx < (x + merge_x + 1)) && (sx < size_x); ++sx)
		}
		
		// for(sy = y_start; (sy < (y + merge_y + 1)) && (sy < size_y); ++sy)
	      }
	      
	      // for(sz = z_start; (sz < z); ++sz)
	    }
	    for(sy = sy_start; (sy < y); ++sy){
	      
	      for(sx = sx_start; (sx <= (x + merge_x + 1)) && (sx < size_x); ++sx){
		
		// check if the test voxel is within the boundaries of the spatial ellipse
		if(((((float) (x - sx)) * ((float) (x - sx)) / (((float) (merge_x + 1)) * ((float) (merge_x + 1)))) + (((float) (y - sy)) * ((float) (y - sy)) / (((float) (merge_y + 1)) * ((float) (merge_y + 1))))) > 1.0){ continue; }

		// check if this test voxel is source
		if((flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] == flag_value) || (flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] >= 0)){
		  
		  // check if this is a new objid, move on if not
		  prev = -1;
		  for(i = 0; (i < NOi) && (i < init_limit); ++i){
		    
		    if(match_init[i] == flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]){
		      
		      prev = 1;
		      break;
		      
		    }
		    
		  }
		  if(prev == -1){ 
		    
		    // add current object to list of previous matches
		    ++NOi;
		    match_init.push_back(flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]);

		  }
		  
		  // update existing flag
		  if(((flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] < existing) && (flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] >= 0)) || (existing == flag_value)){ existing = flag_vals[((z*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))]; }
		  
		}		
		
		// for(sx = x_start; (sx < = (x + merge_x + 1)) && (sx < size_x); ++sx)
	      }
	      
	      // for(sy = y_start; (sy < y) && (sy < size_y); ++sy)
	    }
	    for(sx = sx_start; sx < x; ++sx){
	      
	      // check if this test voxel is source
	      if((flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))] == flag_value) || (flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))] >= 0)){
		
		// check if this is a new objid, move on if not
		prev = -1;
		for(i = 0; (i < NOi) && (i < init_limit); ++i){
		  
		  if(match_init[i] == flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))]){
		    
		    prev = 1;
		    break;
		    
		  }
		  
		}
		if(prev == -1){ 
		  
		  // add current object to list of previous matches
		  ++NOi;
		  match_init.push_back(flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))]);
		  
		}
		
		// update existing flag
		if(((flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))] < existing) && (flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))] >= 0)) || (existing == flag_value)){ existing = flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (sx*data_metric[0]))]; }
		
	      }	  
	      
	      // for(sx = x_start; sx <= x; ++sx)
	    }

	    // else . . . if(ss_mode == 1)
	  }

	  // 3c. assign an object number to this voxel, depending upon the value of the existing flag,
	  // and if this is part of an existing object or objects, then daisy chain from this voxel to all the others
	  if(existing == flag_value){
	    
	    // get obj value from array of obj_ids, then replace it with an incremented obj value
	    
	    // ensure that the obj_ids array is sorted
	    sort(obj_ids.begin(),obj_ids.end());
	    
	    // assign value to array
	    flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (x*data_metric[0]))] = obj_ids[0]; 
	    
	    // calculate the obj_batch value for this obj_id
	    obj_batch = (long int) floor(((double) obj_ids[0] / (double) obj_limit));
	    
	    // add this initial point to the object
	    detections[obj_batch][(obj_ids[0] - (obj_batch * obj_limit))].AddPoint(((double) (chunk_x_start + x)),((double) (chunk_y_start + y)),((double) (chunk_z_start + z)),data_vals[((z*data_metric[2]) + (y*data_metric[1]) + (x*data_metric[0]))]);

	    // push id to list of objects that need to be size thresholded
	    i = -1;
	    for(g = 0; g < NO_check_obj_ids; ++g){ if(check_obj_ids[g] == obj_ids[0]){ i = 1; break; } }
	    if(i == -1){

	      ++NO_check_obj_ids;
	      check_obj_ids.push_back(obj_ids[0]);
	      
	    }

	    // adjust obj and obj_ids
	    if(NO_obj_ids > 1){
	      
	      for(i = 0; i < (NO_obj_ids - 1); ++i){ obj_ids[i] = obj_ids[(i + 1)]; }
	      --NO_obj_ids;
	      obj_ids.resize(NO_obj_ids);

	    } else {
	      
	      ++obj;

	      obj_ids[0] = obj;
	      
	      // if this obj_id modulo obj_limit == 0, then it's the first obj of a new batch
	      // so create a new batch of objects
	      if((obj % obj_limit) == 0){ 
				
		detections.push_back(new object_props_dbl[obj_limit]);
				
		// if((obj % obj_limit) == 0)
	      }
			      
	      // if(NO_obj_ids > 1) { } else 
	    }
	    
	    // if(existing == -1)
	  } else if(existing >= 0){
	    
	    // assign the current voxel being processed the correct flag value, and update
	    // the bounding box
	    
	    // add object to list of objects within the merging box, provided that it hasn't been added already
	    i = -1;
	    for(g = 0; g < NO_check_obj_ids; ++g){ if(check_obj_ids[g] == existing){ i = 1; break; } }
	    if(i == -1){

	      ++NO_check_obj_ids;
	      check_obj_ids.push_back(existing);

	    }

	    // assign value to array
	    flag_vals[((z*data_metric[2]) + (y*data_metric[1]) + (x*data_metric[0]))] = existing;
	    
	    // calculate the obj_batch value for the existing object
	    obj_batch = (long int) floor(((double) existing / (double) obj_limit));
	    
	    // add this point to the existing object
	    detections[obj_batch][(existing - (obj_batch * obj_limit))].AddPoint((chunk_x_start + x),(chunk_y_start + y),(chunk_z_start + z),data_vals[((z*data_metric[2]) + (y*data_metric[1]) + (x*data_metric[0]))]);
	    
	    // for each object detected within the merging volume, update all of the values within the bounding box of this object
	    // to be the same as the other part of this object, and update the bounding box at the same time
	    if(NOi > 1){
	      
	      for(i = 0; i < NOi; ++i){
		
		// move on if this is the existing object
		if(match_init[i] == existing){ 

		  continue; 

		}
		
		// calculate a temporary obj_batch value
		obj_batch_2 = (long int) floor(((double) match_init[i] / (double) obj_limit));		
		
		// add this object to an existing object
		//detections[obj_batch][(existing - (obj_batch * obj_limit))].AddObject(detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))]);
		detections[obj_batch][(existing - (obj_batch * obj_limit))].AddObject(detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))],temp_sparse_reps_grid,temp_sparse_reps_strings,temp_mom0,temp_RAPV,temp_DECPV,temp_ref_spec,temp_obj_spec,temp_vfield);
		
		// flag the `object' values within the bounding box to merge it with the base object
		sx_start = detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetRAmin() - chunk_x_start;
		if(sx_start < 0){ sx_start = 0; }
		sy_start = detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetDECmin() - chunk_y_start;
		if(sy_start < 0){ sy_start = 0; }
		sz_start = detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetFREQmin() - chunk_z_start;
		if(sz_start < 0){ sz_start = 0; }
		for(sz = sz_start; (sz <= detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetFREQmax() - chunk_z_start) && (sz < size_z); ++sz){
		  
		  for(sy = sy_start; (sy <= detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetDECmax() - chunk_y_start) && (sy < size_y); ++sy){
		    
		    for(sx = sx_start; (sx <= detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetRAmax() - chunk_x_start) && (sx < size_x); ++sx){
		      
		      if(flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] == match_init[i]){ flag_vals[((sz*data_metric[2]) + (sy*data_metric[1]) + (sx*data_metric[0]))] = existing; }
		      
		      // for(sx = x_start[(match_init[i] - (obj_batch_2 * obj_limit))]; sx <= x_finish[(match_init[i] - (obj_batch_2 * obj_limit))]; ++sx)
		    }
		    
		    // for(sy = y_start[(match_init[i] - (obj_batch_2 * obj_limit))]; sy <= y_finish[(match_init[i] - (obj_batch_2 * obj_limit))]; ++sy)
		  }
		  
		  // for(sz = z_start[(match_init[i] - (obj_batch_2 * obj_limit))]; sz <= z_finish[(match_init[i] - (obj_batch_2 * obj_limit))]; ++sz)
		}
			
		// re-initialise the `object' values for the match_init[i] object that has just been merged into the existing object, 
		// and push the obj value (match_init[i]) to the top of the array of values
		detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].ReInit();
		if((detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].Get_srep_update() != 0) && (detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].Get_srep_size(0) >= 0)){
		  
		  detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].ReInit_srep();
		  detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].ReInit_mini();

		}
		detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].ReInit_size();
		detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].Set_srep_update(0);

		j = -1;
		for(g = 0; g < NO_obj_ids; ++g){ if(obj_ids[g] == match_init[i]){ j = 1; break; } }
		if(j == -1){
		  
		  ++NO_obj_ids;
		  obj_ids.push_back(match_init[i]);
		  		  
		}
		
		// for(i = 0; i < NOi; ++i)
	      }
	      
	      // if(NOi > 1)
	    }
	    
	    // else . . . if(existing == -1)
	  }
	  
	  // if(flag_vals[((z*size_x*size_y) + (y*size_x) + x)] == -1)
	}
	
	// update progress on display
	while(progress <= (((float) ((z * size_x * size_y) + (y * size_x) + x + 1)) / ((float) (size_x * size_y * size_z)))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
	
	// for(x = 0; x < size_x; ++x)
      }
      
      // for(y = 0; y < size_y; ++y)
    }
    
    // for(z = 0; z < size_z; ++z)
  }
  std::cout << "* done." << std::endl;
  
  // 4. generate sparse representations of coherent objects that are retained after
  // this `chunk' of the datacube is processed
  std::cout << "Generating/updating sparse representations of sources . . . " << std::endl;
  std::cout << "0 | |:| | : | |:| | 100% complete" << std::endl;
  progress = 0.0;
  for(o = 0; o < obj; ++o){
    
    // calculate the obj_batch value for this object
    obj_batch = (long int) floor(((double) o / (double) obj_limit));
    
    // move on if this object has been re-initialised
    if(detections[obj_batch][(o - (obj_batch * obj_limit))].ShowVoxels() < 1){ 
      
      // update progress on display
      while(progress <= (((double) (o + 1) / ((double) obj)))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }  
      continue; 
      
    }
    
    // move on if this object is not within the boundaries of the current chunk, in other words,
    // move on if this object doesn't need to have it's sparse representation and postage stamp images updated
    if((detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmax() < chunk_x_start) || (detections[obj_batch][(o - (obj_batch * obj_limit))].GetDECmax() < chunk_y_start) || (detections[obj_batch][(o - (obj_batch * obj_limit))].GetFREQmax() < chunk_z_start) || (detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmin() >= (chunk_x_start + size_x)) || (detections[obj_batch][(o - (obj_batch * obj_limit))].GetDECmin() >= (chunk_y_start + size_y)) || (detections[obj_batch][(o - (obj_batch * obj_limit))].GetFREQmin() >= (chunk_z_start + size_z))){ 

      // update progress on display
      while(progress <= (((double) (o + 1) / ((double) obj)))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }  
      continue; 
      
    }
    
    // move on if this object's sparse representation and postage stamp images don't need updating
    if((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) >= 0) && (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_update() == -1)){
      
      // update progress on display
      while(progress <= (((double) (i + 1) / ((double) obj)))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }  
      continue; 
	
    }
          
    // update the sparse representation for this object
    
    // 0. if there is no existing sparse representation and postage stamp images, create them from scratch, 
    // otherwise concatenate the existing sparse representation and postage stamp images with the updates
    if((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) < 0) && (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_update() >= 0)){
      
      // branch for creating a brand new sparse representation
      
      // update the sparse_reps_update array
      detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_update(-1);

      // initialise the number of object sections that make up this object
      NOi = 0;   
      
      // calculate the edge of the bounding box in `chunk' co-ordinates
      sx_start = detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmin() - chunk_x_start;
      if(sx_start < 0){ sx_start = 0; }
      sx_finish = detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmax() - chunk_x_start;
      sy_start = detections[obj_batch][(o - (obj_batch * obj_limit))].GetDECmin() - chunk_y_start;
      if(sy_start < 0){ sy_start = 0; }
      sy_finish = detections[obj_batch][(o - (obj_batch * obj_limit))].GetDECmax() - chunk_y_start;
      sz_start = detections[obj_batch][(o - (obj_batch * obj_limit))].GetFREQmin() - chunk_z_start;
      if(sz_start < 0){ sz_start = 0; }
      sz_finish = detections[obj_batch][(o - (obj_batch * obj_limit))].GetFREQmax() - chunk_z_start;
      
      // count the number of object sections that make up this object
      for(sy = sy_start; (sy <= sy_finish) && (sy < size_y); ++sy){
	
	for(sx = sx_start; (sx <= sx_finish) && (sx < size_x); ++sx){
	  
	  // initialise the dummy integer j
	  j = -1;
	  
	  for(sz = sz_start; (sz <= sz_finish) && (sz < size_z); ++sz){
	    
	    // change j to reflect if this voxel in the flag_vals array belongs to the source
	    if(flag_vals[((sz * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))] == o){ 
	      
	      // if previously this was not part of an object string, j == -1, then increment the number of object strings and
	      // adjust the value of j
	      if(j == -1){ ++NOi; }
	      
	      // change j to reflect that this voxel belongs to the object
	      j = 1; 
	      
	    } else { j = -1; }
	    
	    // for(sz = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetFREQmin() - chunk_z_size; sz <= (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetFREQmax() - chunk_z_size); ++sz)
	  }
	  
	  // for(sx = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmin() - chunk_x_size; sx <= (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmax() - chunk_x_size); ++sx)
	}
	
	// for(sy = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmin() - chunk_y_size; sy <= (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmax() - chunk_y_size); ++sy)
      }

      // create an array element to store this sparse representation
      detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_size(0,detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmin());
      detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_size(1,detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmax());
      detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_size(2,detections[obj_batch][(o - (obj_batch * obj_limit))].GetDECmin());
      detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_size(3,detections[obj_batch][(o - (obj_batch * obj_limit))].GetDECmax());
      detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_size(4,detections[obj_batch][(o - (obj_batch * obj_limit))].GetFREQmin());
      detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_size(5,detections[obj_batch][(o - (obj_batch * obj_limit))].GetFREQmax());
      detections[obj_batch][(o - (obj_batch * obj_limit))].Create_srep_grid((((detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmin() + 1) * (detections[obj_batch][(o - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(o - (obj_batch * obj_limit))].GetDECmin() + 1)) + 1));
      detections[obj_batch][(o - (obj_batch * obj_limit))].Create_srep_strings((2 * NOi));
      
      // calculate the extent of the sparse_rep volume
      sx_finish = detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmin() + 1;
      sy_finish = detections[obj_batch][(o - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(o - (obj_batch * obj_limit))].GetDECmin() + 1;
      sz_finish = detections[obj_batch][(o - (obj_batch * obj_limit))].GetFREQmax() - detections[obj_batch][(o - (obj_batch * obj_limit))].GetFREQmin() + 1;
      
      // create arrays to store the postage stamp images and initialise them
      detections[obj_batch][(o - (obj_batch * obj_limit))].Create_mom0((sx_finish * sy_finish));
      detections[obj_batch][(o - (obj_batch * obj_limit))].Create_RAPV((sx_finish * sz_finish));
      detections[obj_batch][(o - (obj_batch * obj_limit))].Create_DECPV((sy_finish * sz_finish));
      detections[obj_batch][(o - (obj_batch * obj_limit))].Create_ospec(sz_finish);
      if(sz_finish >= 10){

	detections[obj_batch][(o - (obj_batch * obj_limit))].Create_rspec((2 * sz_finish));

      } else {

	detections[obj_batch][(o - (obj_batch * obj_limit))].Create_rspec((10 + sz_finish));

      }
      detections[obj_batch][(o - (obj_batch * obj_limit))].Create_vfield((sx_finish * sy_finish));
      
      // initialise the mini_mom0 values
      for(sy = 0; sy < sy_finish; ++sy){
	for(sx = 0; sx < sx_finish; ++sx){
	  detections[obj_batch][(o - (obj_batch * obj_limit))].Set_mom0(((sy * sx_finish) + sx),0.0);
	}
      }

      // initialise the mini_RAPV values
      for(sz = 0; sz < sz_finish; ++sz){
	for(sx = 0; sx < sx_finish; ++sx){
	  detections[obj_batch][(o - (obj_batch * obj_limit))].Set_RAPV(((sz * sx_finish) + sx),0.0);
	}
      }
      
      // initialise the mini_DECPV values
      for(sz = 0; sz < sz_finish; ++sz){
	for(sy = 0; sy < sy_finish; ++sy){
	  detections[obj_batch][(o - (obj_batch * obj_limit))].Set_DECPV(((sz * sy_finish) + sy),0.0);
	}
      }
      
      // initialise the mini_spec values
      for(sz = 0; sz < sz_finish; ++sz){ detections[obj_batch][(o - (obj_batch * obj_limit))].Set_ospec(sz,0.0); }
      if(sz_finish >= 10){

	for(sz = 0; sz < (2 * sz_finish); ++sz){ detections[obj_batch][(o - (obj_batch * obj_limit))].Set_rspec(sz,0.0); }

      } else {

	for(sz = 0; sz < (10 + sz_finish); ++sz){ detections[obj_batch][(o - (obj_batch * obj_limit))].Set_rspec(sz,0.0); }

      }
            
      // initialise the mini_vfield values
      for(sy = 0; sy < sy_finish; ++sy){
	for(sx = 0; sx < sx_finish; ++sx){
	  detections[obj_batch][(o - (obj_batch * obj_limit))].Set_vfield(((sy * sx_finish) + sx),0.0);
	}
      }

      // populate the sparse array representation and postage stamp images
      NOi = 0;
      for(sy = 0; sy < sy_finish; ++sy){
	
	for(sx = 0; sx < sx_finish; ++sx){
	  
	  // write the current number of object strings to the grid
	  detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_grid(((sy * sx_finish) + sx),NOi);
	  
	  // initialise the dummy integer j 
	  j = -1;
	  
	  for(sz = 0; sz < sz_finish; ++sz){
	    
	    // change j to reflect if this voxel in the flag_vals array belongs to the source
	    if(flag_vals[(((sz + sz_start) * data_metric[2]) + ((sy + sy_start) * data_metric[1]) + ((sx + sx_start) * data_metric[0]))] == o){ 
	      
	      // update mini_mom0
	      detections[obj_batch][(o - (obj_batch * obj_limit))].Add_mom0(((sy * sx_finish) + sx),data_vals[(((sz + sz_start) * data_metric[2]) + ((sy + sy_start) * data_metric[1]) + ((sx + sx_start) * data_metric[0]))]);
	      
	      // update mini_RAPV
	      detections[obj_batch][(o - (obj_batch * obj_limit))].Add_RAPV(((sz * sx_finish) + sx),data_vals[(((sz + sz_start) * data_metric[2]) + ((sy + sy_start) * data_metric[1]) + ((sx + sx_start) * data_metric[0]))]);
	      
	      // update mini_DECPV
	      detections[obj_batch][(o - (obj_batch * obj_limit))].Add_DECPV(((sz * sy_finish) + sy),data_vals[(((sz + sz_start) * data_metric[2]) + ((sy + sy_start) * data_metric[1]) + ((sx + sx_start) * data_metric[0]))]);
	      
	      // update mini_obj_spec
	      detections[obj_batch][(o - (obj_batch * obj_limit))].Add_ospec(sz,data_vals[(((sz + sz_start) * data_metric[2]) + ((sy + sy_start) * data_metric[1]) + ((sx + sx_start) * data_metric[0]))]);
  
	      // update mini_vfield
	      detections[obj_batch][(o - (obj_batch * obj_limit))].Add_vfield(((sy * sx_finish) + sx),(((double) (sz + sz_start + chunk_z_start)) * data_vals[(((sz + sz_start) * data_metric[2]) + ((sy + sy_start) * data_metric[1]) + ((sx + sx_start) * data_metric[0]))]));
	      
	      // if previously this was not part of an object string, j == -1, then increment the number of object strings and
	      // write the position to the sparse_reps array
	      if(j == -1){ 
		
		// write the beginning of the object string to sparse_reps
		detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_strings((2 * NOi),(sz_start + chunk_z_start + sz));
		
		// increment NOi
		++NOi; 
		
		// initialise the end of the object string in sparse_reps
		detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_strings(((2 * NOi) - 1),(sz_start + chunk_z_start + sz));

	      } else {
		
		// update the end of the object string in sparse_reps
		detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_strings(((2 * NOi) - 1),(sz_start + chunk_z_start + sz));
		
	      }
	      
	      // change j to reflect that this voxel belongs to the object
	      j = 1; 
	      
	    } else { 
	      
	      // change j to reflect that this voxel doesn't belong to the object
	      j = -1; 
	      
	    }
	    
	    // for(sz = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetFREmin() - chunk_z_start + 1; sz <= (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetFREQmax() - chunk_z_start + 1); ++sz)
	  }
	  
	  // for(sx = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmin() - chunk_x_start + 1; sx <= (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmax() - chunk_x_start + 1); ++sx)
	}
	
	//for(sy = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmin() - chunk_y_start + 1; sy <= (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmax() - chunk_y_start + 1); ++sy)
      }
            
      // update sparse_reps with the final value of NOi
      detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_grid((sx_finish * sy_finish),NOi);
            
      // create the reference spectrum
      if(sz_finish >= 10){
	
	for(sz = -1 * (int) floorf((0.5 * (float) sz_finish)); sz < (sz_finish + ((int) floorf((0.5 * (float) sz_finish)))); ++sz){
	  
	  for(sy = 0; sy < sy_finish; ++sy){
	  
	    for(sx = 0; sx < sx_finish; ++sx){
	    
	      if(((sz + sz_start) >= 0) && ((sz + sz_start) < size_z)){
		
		// update mini_ref_spec	
		detections[obj_batch][(o - (obj_batch * obj_limit))].Add_rspec((sz + ((int) floorf((0.5 * (double) sz_finish)))),data_vals[(((sz + sz_start) * data_metric[2]) + ((sy + sy_start) * data_metric[1]) + ((sx + sx_start) * data_metric[0]))]);	
		
	      }
	      
	    }
	    
	  }
      
	}

	// if(sz_finish >= 10)
      } else {

	for(sz = -5; sz < (sz_finish + 5); ++sz){
	  
	  for(sy = 0; sy < sy_finish; ++sy){
	  
	    for(sx = 0; sx < sx_finish; ++sx){
	    
	      if(((sz + sz_start) >= 0) && ((sz + sz_start) < size_z)){
		
		// update mini_ref_spec
		detections[obj_batch][(o - (obj_batch * obj_limit))].Add_rspec((sz + 5),data_vals[(((sz + sz_start) * data_metric[2]) + ((sy + sy_start) * data_metric[1]) + ((sx + sx_start) * data_metric[0]))]);	
		
	      }
	      
	    }
	    
	  }
      
	}

	// else . . . if(sz_finish >= 10)
      }
      
     // if(sparse_reps[obj_batch][(obj - (obj_batch * obj_limit))][0] < 0) . . . else
    } else {
            
      // branch for updating an existing sparse representation
      
      // count the number of object strings comprising the object
      
      // 1. initialise the number of object strings
      NOi = 0;
      
      // 2. convert the new bounding box of the object into `chunk' co-ordinates
      sx_start = detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmin() - chunk_x_start;
      sx_finish = detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmax() - chunk_x_start;
      sy_start = detections[obj_batch][(o - (obj_batch * obj_limit))].GetDECmin() - chunk_y_start;
      sy_finish = detections[obj_batch][(o - (obj_batch * obj_limit))].GetDECmax() - chunk_y_start;
      sz_start = detections[obj_batch][(o - (obj_batch * obj_limit))].GetFREQmin() - chunk_z_start;
      if(sz_start < 0){ sz_start = 0; }
      sz_finish = detections[obj_batch][(o - (obj_batch * obj_limit))].GetFREQmax() - chunk_z_start;
      
      // 3. search through the bounding box, and do one of the following
      // 3.a if inside the chunk and outside the existing sparse rep grid, run the routine to count the number
      //     of object strings along the current line of sight
      // 3.b otherwise use the grid value to increment NOi
      for(sy = sy_start; sy <= sy_finish; ++sy){
	
	for(sx = sx_start; sx <= sx_finish; ++sx){
	  	  
	  // case 3.a triggered
	  if(((sx > merge_x) || (chunk_x_start == 0)) && ((sy > merge_y) || (chunk_y_start == 0))){
	    	    
	    // initialise the dummy integer j 
	    j = -1;
	    
	    for(sz = sz_start; sz <= sz_finish; ++sz){
	      
	      // change i to reflect if this voxel in the flag_vals array belongs to the source
	      if(flag_vals[((sz * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))] == o){ 
		
		// if previously this was not part of an object string, i == -1, then increment the number of object strings and
		// write the position to the sparse_reps array
		if(j == -1){ 
		  
		  // increment NOi
		  ++NOi; 
		  
		} 
		
		// change j to reflect that this voxel belongs to the object
		j = 1; 
		
	      } else { j = -1; }
	      
	      // for(sz = sz_start; sz <= sz_finish; ++sz)
	    }
	    
	  } else {
	    
	    // case 3.b triggered  
	    NOi+=(detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_grid((((sy + chunk_y_start - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx + chunk_x_start - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_grid((((sy + chunk_y_start - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx + chunk_x_start - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0))));
	    
	  }     
	  
	  // for(sx = sx_start; sx <= sx_finish; ++sx)
	}
	
	// for(sy = sy_start; sy <= sy_finish; ++sy)
      }
      
      // 4. allocate memory to store temporary sparse representation and postage stamp images, then
      //    initialise the postage stamp image arrays
      temp_mom0.resize(0);
      for(sy = 0; sy < (sy_finish - sy_start + 1); ++sy){
	for(sx = 0; sx < (sx_finish - sx_start + 1); ++sx){
	  temp_mom0.push_back(0.0);
	}
      }
      temp_RAPV.resize(0);
      for(sz = 0; sz < (sz_finish - sz_start + 1); ++sz){
	for(sx = 0; sx < (sx_finish - sx_start + 1); ++sx){
	  temp_RAPV.push_back(0.0);
	}
      }
      temp_DECPV.resize(0);
      for(sz = 0; sz < (sz_finish - sz_start + 1); ++sz){
	for(sy = 0; sy < (sy_finish - sy_start + 1); ++sy){
	  temp_DECPV.push_back(0.0);
	}
      }
      temp_obj_spec.resize(0);
      for(sz = 0; sz < (sz_finish - sz_start + 1); ++sz){ temp_obj_spec.push_back(0.0); }
      if((sz_finish - sz_start + 1) >= 10){

	temp_ref_spec.resize(0);
	for(sz = 0; sz < (2 * (sz_finish - sz_start + 1)); ++sz){ temp_ref_spec.push_back(0.0); }

      } else {

	temp_ref_spec.resize(0);
	for(sz = 0; sz < (sz_finish - sz_start + 11); ++sz){ temp_ref_spec.push_back(0.0); }

      }
      temp_vfield.resize(0);
      for(sy = 0; sy < (sy_finish - sy_start + 1); ++sy){
	for(sx = 0; sx < (sx_finish - sx_start + 1); ++sx){
	  temp_vfield.push_back(0.0);
	}
      }
     
      // 5. populate temp_sparse_rep arrays
      // 5.a if inside the chunk and outside the existing sparse rep grid, run the routine to count the number
      //     of object strings along the current line of sight, populate new strings array and add values to new 
      //     postage stamp images
      // 5.b if outside of the chunk, use the grid value to increment NOi, populate strings array from previous version and 
      //     fill in valuese for postage stamp images

      // initialise number of object strings
      NOi = 0;
      
      // initialise RA PV image
      for(sz = 0; sz < (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1); ++sz){

	for(sx = 0; sx < (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1); ++sx){
	
	  temp_RAPV[(((sz + detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) - sz_start - chunk_z_start) * (sx_finish - sx_start + 1)) + sx + detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) - sx_start - chunk_x_start)]+=detections[obj_batch][(o - (obj_batch * obj_limit))].Get_RAPV(((sz * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx));
	  
	}

      }

      // initialise Dec PV image
      for(sz = 0; sz < (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1); ++sz){

	for(sy = 0; sy < (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2) + 1); ++sy){
	
	  temp_DECPV[(((sz + detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) - sz_start - chunk_z_start) * (sy_finish - sy_start + 1)) + sy + detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2) - sy_start - chunk_y_start)]+=detections[obj_batch][(o - (obj_batch * obj_limit))].Get_DECPV(((sz * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2) + 1)) + sy));
	  
	}

      }

      // initialise mini spec's
      for(sz = 0; sz < (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1); ++sz){ temp_obj_spec[(sz + detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) - sz_start - chunk_z_start)]+=detections[obj_batch][(o - (obj_batch * obj_limit))].Get_ospec(sz); }

      if((sz_finish - sz_start + 1) >= 10){

	j = sz_start + chunk_z_start - ((int) floorf((0.5 * (float) (sz_finish - sz_start + 1))));

	if((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1) >= 10){

	  k = detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) - ((int) floorf((0.5 * (double) (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1))));

	  for(sz = 0; sz < (2 * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); ++sz){

	    temp_ref_spec[(sz + k - j)]+=detections[obj_batch][(o - (obj_batch * obj_limit))].Get_rspec(sz);

	  }

	} else {
	  
	  k = detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) - 5;

	  for(sz = 0; sz < (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 11); ++sz){
	    
	    temp_ref_spec[(sz + k - j)]+=detections[obj_batch][(o - (obj_batch * obj_limit))].Get_rspec(sz);

	  }
	  
	}

      } else {

	j = sz_start + chunk_z_start - 5;

	if((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1) >= 10){

	  k = detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) - ((int) floorf((0.5 * (double) (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1))));

	  for(sz = 0; sz < (2 * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); ++sz){

	    temp_ref_spec[(sz + k - j)]+=detections[obj_batch][(o - (obj_batch * obj_limit))].Get_rspec(sz);

	  }

	} else {
	  
	  k = detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) - 5;

	  for(sz = 0; sz < (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 11); ++sz){
	    
	    temp_ref_spec[(sz + k - j)]+=detections[obj_batch][(o - (obj_batch * obj_limit))].Get_rspec(sz);

	  }
	  
	}

      }
      
      // search through LoS's of updated bounding box and take action as dictated by case 5.a or 5.b
      for(sy = sy_start; sy <= sy_finish; ++sy){
	
	for(sx = sx_start; sx <= sx_finish; ++sx){
	  	  
	  // add value to temp_sparse_rep_grid
	  temp_sparse_reps_grid[(((sy - sy_start) * (sx_finish - sx_start + 1)) + sx - sx_start)] = NOi;
	  
	  // case 5.a triggered for moment-0 image and sparse representation
	  if(((sx > merge_x) || (chunk_x_start == 0)) && ((sy > merge_y) || (chunk_y_start == 0))){
	    
	    // initialise the dummy integer j 
	    j = -1;
	    
	    for(sz = sz_start; sz <= sz_finish; ++sz){
	      
	      // change j to reflect if this voxel in the flag_vals array belongs to the source
	      if(flag_vals[((sz * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))] == o){ 
		
		// update mini_mom0
		temp_mom0[(((sy - sy_start) * (sx_finish - sx_start + 1)) + sx - sx_start)]+=data_vals[((sz * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))];
		
		// update mini_RAPV
		temp_RAPV[(((sz - sz_start) * (sx_finish - sx_start + 1)) + sx - sx_start)]+=data_vals[((sz * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))];
		
		// update mini_DECPV
		temp_DECPV[(((sz - sz_start) * (sy_finish - sy_start + 1)) + sy - sy_start)]+=data_vals[((sz * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))];
		
		// update mini_obj_spec
		temp_obj_spec[(sz - sz_start)]+=data_vals[((sz * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))];
		
		// update mini_vfield
		temp_vfield[(((sy - sy_start) * (sx_finish - sx_start + 1)) + sx - sx_start)]+=(((double) (sz + chunk_z_start)) * data_vals[((sz * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))]);
		
		// if previously this was not part of an object string, j == -1, then increment the number of object strings and
		// write the position to the sparse_reps array
		if(j == -1){ 
		  
		  // write the beginning of the object string to sparse_reps
		  temp_sparse_reps_strings[(2 * NOi)] = chunk_z_start + sz;
		  
		  // increment NOi
		  ++NOi; 
		  		  
		  // initialise the end of the object string in sparse_reps
		  temp_sparse_reps_strings[((2 * NOi) - 1)] = chunk_z_start + sz;
		  
		} else {
		  
		  // update the end of the object string in sparse_reps
		  temp_sparse_reps_strings[((2 * NOi) - 1)] = chunk_z_start + sz;
		  
		}
		
		// change j to reflect that this voxel belongs to the object
		j = 1; 
		
	      } else { 
		
		// change j to reflect that this voxel doesn't belong to the object
		j = -1; 
		
	      }
	      
	      // for(sz = sz_start; sz <= sz_finish; ++sz)
	    }
	    
	    // add values to temp_ref_spec
	    if((sz_finish - sz_start + 1) >= 10){
	      
	      k = (int) floorf((0.5 * (float) (sz_finish - sz_start + 1)));
	      
	      for(sz = 0; sz < (2 * (sz_finish - sz_start + 1)); ++sz){
	  	
		if(((sz + sz_start - k) >= 0) && ((sz + sz_start - k) < size_z)){
		  
		  // update mini_ref_spec
		  temp_ref_spec[(sz + sz_start - k)]+=data_vals[(((sz + sz_start - k) * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))];	
		  
		}
		
	      }
	      
	      // if(sz_finish >= 10)
	    } else {
	      
	      for(sz = 0; sz < (sz_finish - sz_start + 11); ++sz){
		
		if(((sz + sz_start - 5) >= 0) && ((sz + sz_start - 5) < size_z)){
		  
		  // update mini_ref_spec
		  temp_ref_spec[(sz + sz_start - 5)]+=data_vals[(((sz + sz_start - 5) * data_metric[2]) + (sy * data_metric[1]) + (sx * data_metric[0]))];	
		  
		}
		
	      }
	      
	      // else . . . if(sz_finish >= 10)
	    }
	    
	  } else {
	    // case 5.b triggered  
	    	    
	    // update temp_sparse_rep_strings
	    j = 0;
	    if(((sx + chunk_x_start) >= detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0)) && ((sx + chunk_x_start) <= (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) + 1)) && ((sy + chunk_y_start) >= detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2)) && ((sy + chunk_y_start) <= (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(3) + 1))){
	      
	      for(g = detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_grid((((sy + chunk_y_start - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx + chunk_x_start - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0))); g < detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_grid((((sy + chunk_y_start - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx + chunk_x_start - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1)); ++g){
		
		temp_sparse_reps_strings[(NOi + j)] = detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_strings((2 * g));
		temp_sparse_reps_strings[(NOi + j + 1)] = detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_strings(((2 * g) + 1));
		++j;
		
	      }
	      
	      // update moment 0 postage stamp image
	      temp_mom0[(((sy - sy_start) * (detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmin() + 1)) + sx - sx_start)] = detections[obj_batch][(o - (obj_batch * obj_limit))].Get_mom0((((sy + chunk_y_start - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx + chunk_x_start - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0)));
	      
	      // update vfield postage stamp image
	      temp_vfield[(((sy - sy_start) * (detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmin() + 1)) + sx - sx_start)] = detections[obj_batch][(o - (obj_batch * obj_limit))].Get_vfield((((sy + chunk_y_start - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx + chunk_x_start - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0)));

	      // update NOi
	      NOi+=(detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_grid((((sy + chunk_y_start - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx + chunk_x_start - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_grid((((sy + chunk_y_start - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + sx + chunk_x_start - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0))));
	      	      
	    }
	    
	    // else . . . case 5.b
	  }
	    
	  // for(sx = sx_start; sx <= sx_finish; ++sx)
	}
		
	// for(sy = sy_start; sy <= sy_finish; ++sy)
      }
            
      temp_sparse_reps_grid[((sx_finish - sx_start + 1) * (sy_finish - sy_start + 1))] = NOi;

      // 6. replace the current sparse representation with the new sparse representation
      detections[obj_batch][(o - (obj_batch * obj_limit))].ReInit_srep();
      detections[obj_batch][(o - (obj_batch * obj_limit))].ReInit_mini();
 
      // update the sparse_reps_size array with the new bounding box
      detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_size(0,detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmin());
      detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_size(1,detections[obj_batch][(o - (obj_batch * obj_limit))].GetRAmax());
      detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_size(2,detections[obj_batch][(o - (obj_batch * obj_limit))].GetDECmin());
      detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_size(3,detections[obj_batch][(o - (obj_batch * obj_limit))].GetDECmax());
      detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_size(4,detections[obj_batch][(o - (obj_batch * obj_limit))].GetFREQmin());
      detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_size(5,detections[obj_batch][(o - (obj_batch * obj_limit))].GetFREQmax());

      // update the sparse_reps_update array
      detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_update(-1);
          
      // populate the new sparse_reps and postage stamps
      detections[obj_batch][(o - (obj_batch * obj_limit))].Create_srep_grid((1 + ((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2) + 1))));
      detections[obj_batch][(o - (obj_batch * obj_limit))].Create_srep_strings((2 * NOi));
      detections[obj_batch][(o - (obj_batch * obj_limit))].Create_mom0(((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2) + 1)));
      detections[obj_batch][(o - (obj_batch * obj_limit))].Create_RAPV(((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1)));
      detections[obj_batch][(o - (obj_batch * obj_limit))].Create_DECPV(((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2) + 1) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1)));
      detections[obj_batch][(o - (obj_batch * obj_limit))].Create_vfield(((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2) + 1)));
  
      detections[obj_batch][(o - (obj_batch * obj_limit))].Create_ospec((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1));
      if((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1) >= 10){
       
	detections[obj_batch][(o - (obj_batch * obj_limit))].Create_rspec((2 * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1)));
	
      } else {
	
	detections[obj_batch][(o - (obj_batch * obj_limit))].Create_rspec((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 11));

      }
   
      for(j = 0; j < (1 + ((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2) + 1))); ++j){ detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_grid(j,temp_sparse_reps_grid[j]); }
      for(j = 0; j < (2 * NOi); ++j){ detections[obj_batch][(o - (obj_batch * obj_limit))].Set_srep_strings(j,temp_sparse_reps_strings[j]); }
      for(j = 0; j < ((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); ++j){ detections[obj_batch][(o - (obj_batch * obj_limit))].Set_mom0(j,temp_mom0[j]); }
      for(j = 0; j < ((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); ++j){ detections[obj_batch][(o - (obj_batch * obj_limit))].Set_RAPV(j,temp_RAPV[j]); }
      for(j = 0; j < ((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2) + 1) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); ++j){ detections[obj_batch][(o - (obj_batch * obj_limit))].Set_DECPV(j,temp_DECPV[j]); }
      for(j = 0; j < ((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); ++j){ detections[obj_batch][(o - (obj_batch * obj_limit))].Set_vfield(j,temp_vfield[j]); }

      for(j = 0; j < (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1); ++j){ detections[obj_batch][(o - (obj_batch * obj_limit))].Set_ospec(j,temp_obj_spec[j]); }
      if((detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1) >= 10){

	for(j = 0; j < (2 * (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); ++j){ detections[obj_batch][(o - (obj_batch * obj_limit))].Set_rspec(j,temp_ref_spec[j]); }
	
      } else {

	for(j = 0; j < (detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(o - (obj_batch * obj_limit))].Get_srep_size(4) + 11); ++j){ detections[obj_batch][(o - (obj_batch * obj_limit))].Set_rspec(j,temp_ref_spec[j]); }

      }
      
      // else . . . if(sparse_reps[obj_batch][(obj - (obj_batch * obj_limit))][0] < 0)
    }
    
    // update progress on display
    while(progress <= (((double) (o + 1) / ((double) obj)))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }  
   
    // for(o = 0; o < obj; ++o)
  }
  std::cout << "* done." << std::endl;
    
  // 5. return the number of coherent objects 
  return obj;
  
}

