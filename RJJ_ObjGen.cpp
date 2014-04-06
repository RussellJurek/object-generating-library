#include<iostream>
#include<RJJ_ObjGen.h>

object_props::object_props(){ 
    
  NOvox = 0; 
  ra = dec = freq = ra_i = dec_i = freq_i = tot_intens = avg_intens = sigma_intens = rms = 0.0;
  ra_min = dec_min = freq_min = min_intens = 1E10;
  ra_max = dec_max = freq_max = max_intens = -1E10;
  
}

object_props::~object_props(){ ; }

void object_props::ReInit(){

  NOvox = 0;
  ra = dec = freq = ra_i = dec_i = freq_i = tot_intens = avg_intens = sigma_intens = rms = 0.0;
  ra_min = dec_min = freq_min = min_intens = 1E10;
  ra_max = dec_max = freq_max = max_intens = -1E10;

}

void object_props::AddVoxel(int value){ NOvox+=value; }

void object_props::AddRa(float value){ ra+=value; }

void object_props::AddDec(float value){ dec+=value; }

void object_props::AddFreq(float value){ freq+=value; }

void object_props::AddRA_i(float value){ ra_i+=value; }

void object_props::AddDec_i(float value){ dec_i+=value; }

void object_props::AddFreq_i(float value){ freq_i+=value; }

void object_props::AddTotIntens(float value){ tot_intens+=value; }

void object_props::AddAvgIntens(float value){ avg_intens+=value; }

void object_props::AddSigmaItens(float value){ sigma_intens+=(value * value); }

void object_props::AdjustRange(float value){
    
  if(value <= min_intens){ min_intens = value; }
  if(value >= max_intens){ max_intens = value; }

}

void object_props::CalcProps(){
  
  ra = ra / NOvox; 
  dec = dec / NOvox; 
  freq = freq / NOvox;
  ra_i = ra_i / tot_intens;
  dec_i = dec_i / tot_intens;
  freq_i = freq_i / tot_intens;
  avg_intens = avg_intens / NOvox;
  rms = sqrtf((sigma_intens / NOvox));
  sigma_intens = sqrtf(((sigma_intens / NOvox) - (avg_intens * avg_intens)));
  
}

void object_props::ShowProps(int id){
    
  std::cout << "ID: " << id << " #voxels: " << NOvox << " ra: " << ra << " dec: " << dec << " freq: " << freq << " ra_i: " << ra_i << " dec_i: " << dec_i << " freq_i: " << freq_i<< " Tot_i: " << tot_intens << " Mean_i: " << avg_intens << " Sig_i: " << sigma_intens << " RMS_i: " << rms << " Min_i: " << min_intens << " Max_i: " << max_intens << " RA_rng: " << (abs(ra_max - ra_min) + 1) << " Dec_rng: " << (abs(dec_max - dec_min) + 1) << " F_rng: " << (abs(freq_max - freq_min) + 1) << std::endl;
  
}

void object_props::ShowProps_file(int id, std::fstream& output_file){
    
  output_file << id << " " << NOvox << " " << ra << " " << dec << " " << freq << " " << ra_i << " " << dec_i << " " << freq_i << " " << tot_intens << " " << avg_intens << " " << sigma_intens << " " << rms << " " << min_intens << " " << max_intens << " " << ra_min << " " << ra_max << " " << dec_min << " " << dec_max << " " << freq_min << " " << freq_max << std::flush;
  
}

void object_props::AdjustRArange(float value){

  if(value <= ra_min){ ra_min = value; }
  if(value >= ra_max){ ra_max = value; }

}

void object_props::AdjustDECrange(float value){

  if(value <= dec_min){ dec_min = value; }
  if(value >= dec_max){ dec_max = value; }

}

void object_props::AdjustFREQrange(float value){

  if(value <= freq_min){ freq_min = value; }
  if(value >= freq_max){ freq_max = value; }

}

int object_props::ShowArea(){ return (int) ((abs(ra_max - ra_min) + 1) * (abs(dec_max - dec_min) + 1)); }

int object_props::ShowVoxels(){ return NOvox; }

int object_props::ShowRArange(){ return (int) (abs(ra_max - ra_min) + 1); }

int object_props::ShowDECrange(){ return (int) (abs(dec_max - dec_min) + 1); }

int object_props::ShowFREQrange(){ return (int) (abs(freq_max - freq_min) + 1); }

int object_props::GetRAmin(){ return ra_min; }

int object_props::GetRAmax(){ return ra_max; }

int object_props::GetDECmin(){ return dec_min; }

int object_props::GetDECmax(){ return dec_max; }

int object_props::GetFREQmin(){ return freq_min; }

int object_props::GetFREQmax(){ return freq_max; }

float object_props::GetRA(){ return ra; }

float object_props::GetDEC(){ return dec; }

float object_props::GetFREQ(){ return freq; }

float object_props::GetRAi(){ return ra_i; }

float object_props::GetDECi(){ return dec_i; }

float object_props::GetFREQi(){ return freq_i; }

float object_props::GetTI(){ return tot_intens; }

float object_props::GetSigmaI(){ return sigma_intens; }

float object_props::GetAvgI(){ return avg_intens; }

float object_props::GetMinI(){ return min_intens; }

float object_props::GetMaxI(){ return max_intens; }

void HeapSort_ids(int n, int ra[]){

  int i,ir,j,l;
  int rra;

  if (n<2) return;
  l = (n >> 1) + 1;
  ir = n;
  for(;;){

    if (l > 1) {
      rra = ra[--l];
    } else {
      rra=ra[ir];
      ra[ir]=ra[1];
      if (--ir == 1) {
	ra[1] = rra;
	break;
      }
    }
    i=l;
    j=l+l;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) j++;
      if (rra < ra[j]) {
	ra[i] = ra[j];
	i = j;
	j <<= 1;
      } else break;
    }
    ra[i] = rra;
  }
}
	    

int CreateObjects(float * data_vals, int * flag_vals, int size_x, int size_y, int size_z, int chunk_x_start, int chunk_y_start, int chunk_z_start, int merge_x, int merge_y, int merge_z, int min_x_size, int min_y_size, int min_z_size, int min_v_size, float intensity_threshold, int flag_value, int start_obj, object_props ** detections, int *** sparse_reps_size, int *** sparse_reps_grid, int *** sparse_reps_strings, int * obj_ids, int& NO_obj_ids, int * check_obj_ids, int& NO_check_obj_ids, int obj_limit, int obj_batch_limit, float *** mini_mom0, float *** mini_RAPV, float *** mini_DECPV, int max_x_val, int max_y_val, int max_z_val, float *** mini_obj_spec, float *** mini_ref_spec, int ** sparse_reps_update){
  
  void HeapSort_ids(int n, int ra[]);
  int x,y,z,obj,existing,sx,sy,sz,sx_start,sy_start,sz_start,sx_finish,sy_finish,sz_finish,init_limit,g;
  int i,j,k,NOi, * match_init, obj_batch, obj_batch_2;
  int * temp_sparse_reps_grid, * temp_sparse_reps_strings;
  int prev,x_start,y_start;
  float progress;
  float * temp_mom0, * temp_RAPV, * temp_DECPV, * temp_obj_spec, * temp_ref_spec;
  
  // create temporary arrays
  if(10000 > (((2 * merge_x) + 3) * ((2 * merge_y) + 3) * ((2 * merge_z) + 3))){
    
    init_limit = 10000;
    
  } else {
    
    init_limit = ((2 * merge_x) + 3) * ((2 * merge_y) + 3) * ((2 * merge_z) + 3) + 10;
    
  }
  match_init = new int [init_limit];
  
  // initialise the temporary arrays
  temp_sparse_reps_grid = new int [1000000];
  temp_sparse_reps_strings = new int [1000000];
  temp_mom0 = new float [1000000];
  temp_RAPV = new float [1000000];
  temp_DECPV = new float [1000000];
  temp_obj_spec = new float [1000000];
  temp_ref_spec = new float [1000000];
 
  // 0. initialise variables and arrays
  obj = start_obj;
  x_start = 0;
  y_start = 0;
  if(chunk_x_start > 0){ x_start = merge_x + 1; }
  if(chunk_y_start > 0){ y_start = merge_y + 1; }
  
  // 1. Create list of `coherent' objects from neighbouring voxels
  
  // for each grid point, check if it has been flagged as a source voxel, and assign
  // it an object ID if it is. If a neighbouring voxel has already been flagged, assign 
  // that ID to the grid point, otherwise, assign current value of obj and increment obj.
  std::cout << "0 | |:| | : | |:| | 100% complete" << std::endl;
  progress = 0.0;
  for(z = 0; ((z < size_z) && ((z + chunk_z_start) < max_z_val)); z++){

    // if the z value is sufficiently large that objects have started to pop out of the merging box,
    // then check if the objects outside of the merging box are sufficiently large
    if(z > (merge_z + 1)){
      
      // initialise NOi to store the number of current check_obj_ids values that haven't moved outside of 
      // the merging volume, and that form the basis for the new check_obj_ids list
      NOi = 0;
      
      // apply size threshold to objects added since last size thresholding
      for(i = 0; i < NO_check_obj_ids; i++){
	
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
	  if((((detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetRAmin() + 1) < min_x_size) || ((detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetDECmin() + 1) < min_y_size) || ((detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmax() - detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmin() + 1) < min_z_size) || (detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].ShowVoxels() < min_v_size)) && ((((chunk_x_start + size_x - 1 - detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetRAmax()) > merge_x) || ((chunk_x_start + size_x) >= max_x_val)) && (((chunk_y_start + size_y - 1 - detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetDECmax()) > merge_y) || ((chunk_y_start + size_y) >= max_y_val)) && (((chunk_z_start + size_z - 1 - detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmax()) > merge_z) || ((chunk_z_start + size_z) >= max_z_val)))){
	    
	    // re-initialise object
	    detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].ReInit();
	    if((sparse_reps_update[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))] != 0) && (sparse_reps_size[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))][0] >= 0)){
	      
	      delete [] sparse_reps_grid[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))];
	      sparse_reps_grid[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))] = NULL;
	      delete [] sparse_reps_strings[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))];
	      sparse_reps_strings[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))] = NULL;
	      delete [] mini_mom0[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))];
	      mini_mom0[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))] = NULL;
	      delete [] mini_RAPV[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))];
	      mini_RAPV[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))] = NULL;
	      delete [] mini_DECPV[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))];
	      mini_DECPV[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))] = NULL;
	      delete [] mini_obj_spec[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))];
	      mini_obj_spec[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))] = NULL;
	      delete [] mini_ref_spec[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))];
	      mini_ref_spec[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))] = NULL;

	    }
	    sparse_reps_size[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))][0] = -99;
	    sparse_reps_size[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))][1] = -99;
	    sparse_reps_size[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))][2] = -99;
	    sparse_reps_size[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))][3] = -99;
	    sparse_reps_size[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))][4] = -99;
	    sparse_reps_size[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))][5] = -99;
	    sparse_reps_update[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))] = 0;

	    // flag the `object' values within the bounding box and remove it from the flag_vals array
	    sx_start = detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetRAmin() - chunk_x_start;
	    if(sx_start < 0){ sx_start = 0; }
	    sy_start = detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetDECmin() - chunk_y_start;
	    if(sy_start < 0){ sy_start = 0; }
	    sz_start = detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmin() - chunk_z_start;
	    if(sz_start < 0){ sz_start = 0; }
	    for(sz = sz_start; ((sz <= (detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetFREQmax() - chunk_z_start)) && (sz < size_z)); sz++){
	      
	      for(sy = sy_start; ((sy <= (detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetDECmax() - chunk_y_start)) && (sy < size_y)); sy++){
		
		for(sx = sx_start; ((sx <= (detections[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))].GetRAmax() - chunk_x_start)) && (sx < size_x)); sx++){
		  
		  if(flag_vals[((sz*size_x*size_y) + (sy*size_x) + sx)] == check_obj_ids[i]){ flag_vals[((sz*size_x*size_y) + (sy*size_x) + sx)] = -99; }
		  
		  // for(sx = x_start[(match_init[i] - (obj_batch_2 * obj_limit))]; sx <= x_finish[(match_init[i] - (obj_batch_2 * obj_limit))]; sx++)
		}
		
		// for(sy = y_start[(match_init[i] - (obj_batch_2 * obj_limit))]; sy <= y_finish[(match_init[i] - (obj_batch_2 * obj_limit))]; sy++)
	      }
	      
	      // for(sz = z_start[(match_init[i] - (obj_batch_2 * obj_limit))]; sz <= z_finish[(match_init[i] - (obj_batch_2 * obj_limit))]; sz++)
	    }
	    
	    // add object id to list of available ids
	    j = -1;
	    for(g = 0; g < NO_obj_ids; g++){ if(obj_ids[g] == check_obj_ids[i]){ j = 1; break; } }
	    if(j == -1){
	      obj_ids[NO_obj_ids] = check_obj_ids[i];
	      NO_obj_ids++;
	    }

	  }
	  
	  // if(z_finish[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))] < (z - merge_z - 1))
	} else {
	  
	  check_obj_ids[NOi] = check_obj_ids[i];
	  NOi++;
	  
	  // else . . . if(z_finish[obj_batch][(check_obj_ids[i] - (obj_batch * obj_limit))] < (z - merge_z - 1))
	}
	
	// for(i = 0; i < NO_check_obj_ids; i++)
      }

      // update NO_check_obj_ids
      NO_check_obj_ids = NOi;
      
      // if(z > (merge_z + 2))
    }
    
    for(y = y_start; (y < size_y) && ((chunk_y_start + y) < max_y_val); y++){
      
      for(x = x_start; (x < size_x) && ((chunk_x_start + x) < max_x_val); x++){

	// if this is a source voxel, check if it is associated with a previously identified
	// source voxel
	if(flag_vals[((z*size_x*size_y) + (y*size_x) + x)] == flag_value){ 
	  
	  // 1.  initialise variables
	  existing = flag_value;
	  NOi = 0;
	  
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
	  for(sz = sz_start; (sz < z); sz++){
	    
	    for(sy = sy_start; (sy <= (y + merge_y + 1)) && (sy < size_y); sy++){
	      
	      for(sx = sx_start; (sx <= (x + merge_x + 1)) && (sx < size_x); sx++){
		
		// check if this test voxel is source
		if((flag_vals[((sz*size_x*size_y) + (sy*size_x) + sx)] == flag_value) || (flag_vals[((sz*size_x*size_y) + (sy*size_x) + sx)] >= 0)){
		  
		  // check if this is a new objid
		  prev = -1;
		  for(i = 0; (i < NOi) && (i < init_limit); i++){
		    
		    if(match_init[i] == flag_vals[((sz*size_x*size_y) + (sy*size_x) + sx)]){
		      
		      prev = 1;
		      break;
		      
		    }
		    
		  }
		  if(prev == -1){ 
		    
		    // add current object to list of previous matches
		    match_init[NOi] = flag_vals[((sz*size_x*size_y) + (sy*size_x) + sx)];
		    NOi++;
		    
		  }
  		  
		  // update existing flag
		  if(((flag_vals[((sz*size_x*size_y) + (sy*size_x) + sx)] < existing) && (flag_vals[((sz*size_x*size_y) + (sy*size_x) + sx)] >= 0)) || (existing == flag_value)){ existing = flag_vals[((sz*size_x*size_y) + (sy*size_x) + sx)]; }
		  
		}
		
		// for(sx = x_start; (sx < (x + merge_x + 1)) && (sx < size_x); sx++)
	      }
	      
	      // for(sy = y_start; (sy < (y + merge_y + 1)) && (sy < size_y); sy++)
	    }

	    // for(sz = z_start; (sz < z); sz++)
	  }
	  for(sy = sy_start; (sy < y); sy++){
	    
	    for(sx = sx_start; (sx <= (x + merge_x + 1)) && (sx < size_x); sx++){
	      
	      // check if this test voxel is source
	      if((flag_vals[((z*size_x*size_y) + (sy*size_x) + sx)] == flag_value) || (flag_vals[((z*size_x*size_y) + (sy*size_x) + sx)] >= 0)){
		
		// check if this is a new objid, move on if not
		prev = -1;
		for(i = 0; (i < NOi) && (i < init_limit); i++){
		  
		  if(match_init[i] == flag_vals[((z*size_x*size_y) + (sy*size_x) + sx)]){
		    
		    prev = 1;
		    break;
		    
		  }
		  
		}
		if(prev == -1){ 
		  
		  // add current object to list of previous matches
		  match_init[NOi] = flag_vals[((z*size_x*size_y) + (sy*size_x) + sx)];
		  NOi++;
		  
		}
		
		// update existing flag
		if(((flag_vals[((z*size_x*size_y) + (sy*size_x) + sx)] < existing) && (flag_vals[((z*size_x*size_y) + (sy*size_x) + sx)] >= 0)) || (existing == flag_value)){ existing = flag_vals[((z*size_x*size_y) + (sy*size_x) + sx)]; }
		
	      }		
	      
	      // for(sx = x_start; (sx < = (x + merge_x + 1)) && (sx < size_x); sx++)
	    }
	    
	    // for(sy = y_start; (sy < y) && (sy < size_y); sy++)
	  }
	  for(sx = sx_start; sx < x; sx++){
	    
	    // check if this test voxel is source
	    if((flag_vals[((z*size_x*size_y) + (y*size_x) + sx)] == flag_value) || (flag_vals[((z*size_x*size_y) + (y*size_x) + sx)] >= 0)){
	      
	      // check if this is a new objid, move on if not
	      prev = -1;
	      for(i = 0; (i < NOi) && (i < init_limit); i++){
		
		if(match_init[i] == flag_vals[((z*size_x*size_y) + (y*size_x) + sx)]){
		  
		  prev = 1;
		  break;
		  
		}
		
	      }
	      if(prev == -1){ 
		
		// add current object to list of previous matches
		match_init[NOi] = flag_vals[((z*size_x*size_y) + (y*size_x) + sx)];
		NOi++;
		
	      }
	      
	      // update existing flag
	      if(((flag_vals[((z*size_x*size_y) + (y*size_x) + sx)] < existing) && (flag_vals[((z*size_x*size_y) + (y*size_x) + sx)] >= 0)) || (existing == flag_value)){ existing = flag_vals[((z*size_x*size_y) + (y*size_x) + sx)]; }
	      
	    }	  
	    
	    // for(sx = x_start; sx <= x; sx++)
	  }
	  
	  // 3c. assign an object number to this voxel, depending upon the value of the existing flag,
	  // and if this is part of an existing object or objects, then daisy chain from this voxel to all the others
	  if(existing == flag_value){
	    
	    // get obj value from array of obj_ids, then replace it with an incremented obj value
	    
	    // ensure that the obj_ids array is sorted
	    HeapSort_ids(NO_obj_ids,obj_ids - 1);    
	    
	    // assign value to array
	    flag_vals[((z*size_x*size_y) + (y*size_x) + x)] = obj_ids[0]; 
	    
	    // calculate the obj_batch value for this obj_id
	    obj_batch = (int) floorf(((float) obj_ids[0] / (float) obj_limit));
	    
	    // create initial bounding box
	    detections[obj_batch][(obj_ids[0] - (obj_batch * obj_limit))].AdjustRArange((chunk_x_start + x));
	    detections[obj_batch][(obj_ids[0] - (obj_batch * obj_limit))].AdjustDECrange((chunk_y_start + y));
	    detections[obj_batch][(obj_ids[0] - (obj_batch * obj_limit))].AdjustFREQrange((chunk_z_start + z));
	    
	    // create initial voxel count
	    detections[obj_batch][(obj_ids[0] - (obj_batch * obj_limit))].AddVoxel(1);
	    
	    // call member functions for object_props to adjust appropriate values
	    detections[obj_batch][(obj_ids[0] - (obj_batch * obj_limit))].AddRa((float) (chunk_x_start + x));
	    detections[obj_batch][(obj_ids[0] - (obj_batch * obj_limit))].AddDec((float) (chunk_y_start + y));
	    detections[obj_batch][(obj_ids[0] - (obj_batch * obj_limit))].AddFreq((float) (chunk_z_start + z));
	    detections[obj_batch][(obj_ids[0] - (obj_batch * obj_limit))].AddRA_i(((float) (chunk_x_start + x) * data_vals[((z*size_x*size_y) + (y*size_x) + x)]));
	    detections[obj_batch][(obj_ids[0] - (obj_batch * obj_limit))].AddDec_i(((float) (chunk_y_start + y) * data_vals[((z*size_x*size_y) + (y*size_x) + x)]));
	    detections[obj_batch][(obj_ids[0] - (obj_batch * obj_limit))].AddFreq_i(((float) (chunk_z_start + z) * data_vals[((z*size_x*size_y) + (y*size_x) + x)]));
	    detections[obj_batch][(obj_ids[0] - (obj_batch * obj_limit))].AddTotIntens(data_vals[((z*size_x*size_y) + (y*size_x) + x)]);
	    detections[obj_batch][(obj_ids[0] - (obj_batch * obj_limit))].AddAvgIntens(data_vals[((z*size_x*size_y) + (y*size_x) + x)]);
	    detections[obj_batch][(obj_ids[0] - (obj_batch * obj_limit))].AddSigmaItens(data_vals[((z*size_x*size_y) + (y*size_x) + x)]);
	    detections[obj_batch][(obj_ids[0] - (obj_batch * obj_limit))].AdjustRange(data_vals[((z*size_x*size_y) + (y*size_x) + x)]);
	    
	    // change sparse_reps_update
	    sparse_reps_update[obj_batch][(obj_ids[0] - (obj_batch * obj_limit))] = 1;

	    // push id to list of objects that need to be size thresholded
	    i = -1;
	    for(g = 0; g < NO_check_obj_ids; g++){ if(check_obj_ids[g] == obj_ids[0]){ i = 1; break; } }
	    if(i == -1){

	      check_obj_ids[NO_check_obj_ids] = obj_ids[0];
	      NO_check_obj_ids++;
	    
	    }

	    // adjust obj and obj_ids
	    if(NO_obj_ids > 1){
	      
	      for(i = 0; i < (NO_obj_ids - 1); i++){ obj_ids[i] = obj_ids[(i + 1)]; }
	      NO_obj_ids--;
	      
	    } else {
	      
	      obj++;

	      if(obj < (obj_limit * obj_batch_limit)){ 

		obj_ids[0] = obj;
	      
		// if this obj_id modulo obj_limit == 0, then it's the first obj of a new batch
		// so create a new batch of objects
		if((obj % obj_limit) == 0){ 
		  
		  // calculate the obj_batch value for this obj_id
		  obj_batch = (int) floorf(((float) obj / (float) obj_limit));
		  
		  // create new batch of objects and associated arrays
		  detections[obj_batch] = new object_props [obj_limit];
		  sparse_reps_size[obj_batch] = new int * [obj_limit];
		  sparse_reps_grid[obj_batch] = new int * [obj_limit];
		  sparse_reps_strings[obj_batch] = new int * [obj_limit];
		  sparse_reps_update[obj_batch] = new int [obj_limit];
		  mini_mom0[obj_batch] = new float * [obj_limit];
		  mini_RAPV[obj_batch] = new float * [obj_limit];
		  mini_DECPV[obj_batch] = new float * [obj_limit];
		  mini_obj_spec[obj_batch] = new float * [obj_limit];
		  mini_ref_spec[obj_batch] = new float * [obj_limit];
		  for(i = 0; i < obj_limit; i ++){ 
		    
		    sparse_reps_size[obj_batch][i] = new int[6]; 
		    sparse_reps_size[obj_batch][i][0] = -99; 
		    sparse_reps_size[obj_batch][i][1] = -99; 
		    sparse_reps_size[obj_batch][i][2] = -99; 
		    sparse_reps_size[obj_batch][i][3] = -99; 
		    sparse_reps_size[obj_batch][i][4] = -99; 
		    sparse_reps_size[obj_batch][i][5] = -99; 
		    
		    sparse_reps_update[obj_batch][i] = 0;
		    sparse_reps_grid[obj_batch][i] = NULL;
		    sparse_reps_strings[obj_batch][i] = NULL;
		    mini_mom0[obj_batch][i] = NULL;
		    mini_RAPV[obj_batch][i] = NULL;
		    mini_DECPV[obj_batch][i] = NULL;
		    mini_obj_spec[obj_batch][i] = NULL;
		    mini_ref_spec[obj_batch][i] = NULL;
		   
		    // for(i = 0; i < obj_limit; i ++)
		  }
		  
		  // if((obj % obj_limit) == 0)
		}
		
		// if(obj < (obj_limit * obj_batch_limit))
	      }

	      // if(NO_obj_ids > 1) { } else 
	    }
	    
	    // if(existing == -1)
	  } else if(existing >= 0){
	    
	    // assign the current voxel being processed the correct flag value, and update
	    // the bounding box
	    
	    // add object to list of objects within the merging box, provided that it hasn't been added already
	    i = -1;
	    for(g = 0; g < NO_check_obj_ids; g++){ if(check_obj_ids[g] == existing){ i = 1; break; } }
	    if(i == -1){

	      check_obj_ids[NO_check_obj_ids] = existing;
	      NO_check_obj_ids++;
	    
	    }

	    // assign value to array
	    flag_vals[((z*size_x*size_y) + (y*size_x) + x)] = existing;
	    
	    // calculate the obj_batch value for the existing object
	    obj_batch = (int) floorf(((float) existing / (float) obj_limit));
	    
	    // modify sparse_reps_update
	    sparse_reps_update[obj_batch][(existing - (obj_batch * obj_limit))] = 1;

	    // update bounding box
	    detections[obj_batch][(existing - (obj_batch * obj_limit))].AdjustRArange((chunk_x_start + x));
	    detections[obj_batch][(existing - (obj_batch * obj_limit))].AdjustDECrange((chunk_y_start + y));
	    detections[obj_batch][(existing - (obj_batch * obj_limit))].AdjustFREQrange((chunk_z_start + z));
	    
	    // update properties of the existing object with the voxel being merged in
	    detections[obj_batch][(existing - (obj_batch * obj_limit))].AddVoxel(1);
	    detections[obj_batch][(existing - (obj_batch * obj_limit))].AddRa((float) (chunk_x_start + x));
	    detections[obj_batch][(existing - (obj_batch * obj_limit))].AddDec((float) (chunk_y_start + y));
	    detections[obj_batch][(existing - (obj_batch * obj_limit))].AddFreq((float) (chunk_z_start + z));
	    detections[obj_batch][(existing - (obj_batch * obj_limit))].AddRA_i(((float) (chunk_x_start + x) * data_vals[((z*size_x*size_y) + (y*size_x) + x)]));
	    detections[obj_batch][(existing - (obj_batch * obj_limit))].AddDec_i(((float) (chunk_y_start + y) * data_vals[((z*size_x*size_y) + (y*size_x) + x)]));
	    detections[obj_batch][(existing - (obj_batch * obj_limit))].AddFreq_i(((float) (chunk_z_start + z) * data_vals[((z*size_x*size_y) + (y*size_x) + x)]));
	    detections[obj_batch][(existing - (obj_batch * obj_limit))].AddTotIntens(data_vals[((z*size_x*size_y) + (y*size_x) + x)]);
	    detections[obj_batch][(existing - (obj_batch * obj_limit))].AddAvgIntens(data_vals[((z*size_x*size_y) + (y*size_x) + x)]);
	    detections[obj_batch][(existing - (obj_batch * obj_limit))].AddSigmaItens(data_vals[((z*size_x*size_y) + (y*size_x) + x)]);
	    detections[obj_batch][(existing - (obj_batch * obj_limit))].AdjustRange(data_vals[((z*size_x*size_y) + (y*size_x) + x)]);
	    
	    // for each object detected within the merging volume, update all of the values within the bounding box of this object
	    // to be the same as the other part of this object, and update the bounding box at the same time
	    if(NOi > 1){
	      
	      for(i = 0; i < NOi; i++){
		
		// move on if this is the existing object
		if(match_init[i] == existing){ 

		  continue; 

		}
		
		// calculate a temporary obj_batch value
		obj_batch_2 = (int) floorf(((float) match_init[i] / (float) obj_limit));		
		
		// update the properties of the existing object with the object being merged in
		detections[obj_batch][(existing - (obj_batch * obj_limit))].AddVoxel(detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].ShowVoxels());
		detections[obj_batch][(existing - (obj_batch * obj_limit))].AdjustRArange(detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetRAmin());
		detections[obj_batch][(existing - (obj_batch * obj_limit))].AdjustRArange(detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetRAmax());
		detections[obj_batch][(existing - (obj_batch * obj_limit))].AdjustDECrange(detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetDECmin());
		detections[obj_batch][(existing - (obj_batch * obj_limit))].AdjustDECrange(detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetDECmax());
		detections[obj_batch][(existing - (obj_batch * obj_limit))].AdjustFREQrange(detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetFREQmin());
		detections[obj_batch][(existing - (obj_batch * obj_limit))].AdjustFREQrange(detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetFREQmax());
		detections[obj_batch][(existing - (obj_batch * obj_limit))].AddRa(detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetRA());
		detections[obj_batch][(existing - (obj_batch * obj_limit))].AddDec(detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetDEC());
		detections[obj_batch][(existing - (obj_batch * obj_limit))].AddFreq(detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetFREQ());
		detections[obj_batch][(existing - (obj_batch * obj_limit))].AddRA_i(detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetRAi());
		detections[obj_batch][(existing - (obj_batch * obj_limit))].AddDec_i(detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetDECi());
		detections[obj_batch][(existing - (obj_batch * obj_limit))].AddFreq_i(detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetFREQi());
		detections[obj_batch][(existing - (obj_batch * obj_limit))].AddTotIntens(detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetTI());
		detections[obj_batch][(existing - (obj_batch * obj_limit))].AddAvgIntens(detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetAvgI());
		detections[obj_batch][(existing - (obj_batch * obj_limit))].AddSigmaItens(detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetSigmaI());
		detections[obj_batch][(existing - (obj_batch * obj_limit))].AdjustRange(detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetMinI());		
		detections[obj_batch][(existing - (obj_batch * obj_limit))].AdjustRange(detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetMaxI());		
		
		// flag the `object' values within the bounding box to merge it with the base object
		sx_start = detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetRAmin() - chunk_x_start;
		if(sx_start < 0){ sx_start = 0; }
		sy_start = detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetDECmin() - chunk_y_start;
		if(sy_start < 0){ sy_start = 0; }
		sz_start = detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetFREQmin() - chunk_z_start;
		if(sz_start < 0){ sz_start = 0; }
		for(sz = sz_start; (sz <= detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetFREQmax() - chunk_z_start) && (sz < size_z); sz++){
		  
		  for(sy = sy_start; (sy <= detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetDECmax() - chunk_y_start) && (sy < size_y); sy++){
		    
		    for(sx = sx_start; (sx <= detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetRAmax() - chunk_x_start) && (sx < size_x); sx++){
		      
		      if(flag_vals[((sz*size_x*size_y) + (sy*size_x) + sx)] == match_init[i]){ flag_vals[((sz*size_x*size_y) + (sy*size_x) + sx)] = existing; }
		      
		      // for(sx = x_start[(match_init[i] - (obj_batch_2 * obj_limit))]; sx <= x_finish[(match_init[i] - (obj_batch_2 * obj_limit))]; sx++)
		    }
		    
		    // for(sy = y_start[(match_init[i] - (obj_batch_2 * obj_limit))]; sy <= y_finish[(match_init[i] - (obj_batch_2 * obj_limit))]; sy++)
		  }
		  
		  // for(sz = z_start[(match_init[i] - (obj_batch_2 * obj_limit))]; sz <= z_finish[(match_init[i] - (obj_batch_2 * obj_limit))]; sz++)
		}
			
		// if a sparse representation for this object exists, merge it into the sparse representation of the base object
		if((sparse_reps_update[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))] != 0) && (sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][0] >= 0)){
		  
		  // write the existing object's sparse representations into temporary arrays and initialise temporary arrays at the same time, provided it exists
		  
		  // a. grid
		  for(g = 0; g < (1 + ((detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmin() + 1) * (detections[obj_batch][(existing - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetDECmin() + 1))); g++){ temp_sparse_reps_grid[g] = 0; }
		  
		  if((sparse_reps_update[obj_batch][(existing - (obj_batch * obj_limit))] != 0) && (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] >= 0)){
		    
		    for(sx = 0; sx < (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][1] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] + 1); sx++){
		      for(sy = 0; sy < (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][3] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][2] + 1); sy++){
			
			temp_sparse_reps_grid[(((sy + sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][2] - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetDECmin()) * (detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmin() + 1)) + sx + sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmin())] = sparse_reps_grid[obj_batch][(existing - (obj_batch * obj_limit))][((sy * (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][1] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] + 1)) + sx + 1)] - sparse_reps_grid[obj_batch][(existing - (obj_batch * obj_limit))][((sy * (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][1] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] + 1)) + sx)];
			
		      }
		      
		    }
		    
		  }
		  
		  // b. mom-0
		  for(g = 0; g < ((detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmin() + 1) * (detections[obj_batch][(existing - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetDECmin() + 1)); g++){ temp_mom0[g] = 0.0; }
		  
		  if((sparse_reps_update[obj_batch][(existing - (obj_batch * obj_limit))] != 0) && (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] >= 0)){
		    
		    for(sx = 0; sx < (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][1] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] + 1); sx++){
		   
		      for(sy = 0; sy < (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][3] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][2] + 1); sy++){
			
			temp_mom0[(((sy + sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][2] - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetDECmin()) * (detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmin() + 1)) + sx + sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmin())]+=mini_mom0[obj_batch][(existing - (obj_batch * obj_limit))][((sy * (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][1] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] + 1)) + sx)];
			
		      }
		    
		    }
		    
		  }
		  		  
		  // c. RAPV
		  for(g = 0; g < ((detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmin() + 1) * (detections[obj_batch][(existing - (obj_batch * obj_limit))].GetFREQmax() - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetFREQmin() + 1)); g++){ temp_RAPV[g] = 0.0; }
		  
		  if((sparse_reps_update[obj_batch][(existing - (obj_batch * obj_limit))] != 0) && (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] >= 0)){
		    
		    for(sx = 0; sx < (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][1] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] + 1); sx++){
		   
		      for(sz = 0; sz < (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][4] + 1); sz++){
			
			temp_RAPV[(((sz + sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][4] - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetFREQmin()) * (detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmin() + 1)) + sx + sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmin())]+=mini_RAPV[obj_batch][(existing - (obj_batch * obj_limit))][((sz * (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][1] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] + 1)) + sx)];
			
		      }
		      
		    }
		    
		  }
		  		  
		  // d. DECPV
		  for(g = 0; g < ((detections[obj_batch][(existing - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetDECmin() + 1) * (detections[obj_batch][(existing - (obj_batch * obj_limit))].GetFREQmax() - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetFREQmin() + 1)); g++){ temp_DECPV[g] = 0.0; }
		  
		  if((sparse_reps_update[obj_batch][(existing - (obj_batch * obj_limit))] != 0) && (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] >= 0)){
		    
		    for(sy = 0; sy < (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][3] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][2] + 1); sy++){
		 
		      for(sz = 0; sz < (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][4] + 1); sz++){
			
			temp_DECPV[(((sz + sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][4] - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetFREQmin()) * (detections[obj_batch][(existing - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetDECmin() + 1)) + sy + sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][2] - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetDECmin())]+=mini_DECPV[obj_batch][(existing - (obj_batch * obj_limit))][((sz * (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][3] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][2] + 1)) + sy)];
			
		      }
		      
		    }

		  }
		  
		  // e. ref_spec
		  if((sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][4] + 1) >= 10){
		    
		    for(g = 0; g < (2 * (detections[obj_batch][(existing - (obj_batch * obj_limit))].GetFREQmax() - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetFREQmin() + 1)); g++){ temp_ref_spec[g] = 0.0; }
		    
		    if((sparse_reps_update[obj_batch][(existing - (obj_batch * obj_limit))] != 0) && (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] >= 0)){
		      
		      for(sz = 0; sz < (2 * (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][4] + 1)); sz++){ temp_ref_spec[(sz + sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][4] - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetFREQmin())]+=mini_ref_spec[obj_batch][(existing - (obj_batch * obj_limit))][sz]; } 
		      
		    }
		    
		  } else {
		    
		    for(g = 0; g < (detections[obj_batch][(existing - (obj_batch * obj_limit))].GetFREQmax() - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetFREQmin() + 11); g++){ temp_ref_spec[g] = 0.0; }
		    
		    if((sparse_reps_update[obj_batch][(existing - (obj_batch * obj_limit))] != 0) && (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] >= 0)){
		      
		      for(sz = 0; sz < (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][4] + 11); sz++){ temp_ref_spec[(sz + sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][4] - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetFREQmin())]+=mini_ref_spec[obj_batch][(existing - (obj_batch * obj_limit))][sz]; } 		    
		      
		    }
		    
		  }
		  
		  // f. obj_spec
		  for(g = 0; g < (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][4] + 1); g++){ temp_obj_spec[g] = 0.0; }
		  
		  if((sparse_reps_update[obj_batch][(existing - (obj_batch * obj_limit))] != 0) && (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] >= 0)){
		    
		    for(sz = 0; sz < (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][4] + 1); sz++){
		      
		      temp_obj_spec[(sz + sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][4] - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetFREQmin())]+=mini_obj_spec[obj_batch][(existing - (obj_batch * obj_limit))][sz];
		      
		    }
		    
		  }
		  
		  // write the merged object's sparse representations into temporary arrays
		  		  
		  // a. grid
		  for(sx = 0; sx < (sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][1] - sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][0] + 1); sx++){
		    
		    for(sy = 0; sy < (sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][3] - sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][2] + 1); sy++){
		      
		      temp_sparse_reps_grid[(((sy + sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][2] - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetDECmin()) * (detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmin() + 1)) + sx + sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][0] - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmin())]+=(sparse_reps_grid[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][((sy * (sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][1] - sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][0] + 1)) + sx + 1)] - sparse_reps_grid[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][((sy * (sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][1] - sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][0] + 1)) + sx)]);
		      
		    }
		    
		  }
		  		  
		  // b. mom-0
		  for(sx = 0; sx < (sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][1] - sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][0] + 1); sx++){
		    
		    for(sy = 0; sy < (sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][3] - sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][2] + 1); sy++){
		      
		      temp_mom0[(((sy + sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][2] - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetDECmin()) * (detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmin() + 1)) + sx + sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][0] - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmin())]+=mini_mom0[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][((sy * (sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][1] - sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][0] + 1)) + sx)];
		      
		    }
		    
		  }
		  
		  // c. RAPV
		  for(sx = 0; sx < (sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][1] - sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][0] + 1); sx++){
		    
		    for(sz = 0; sz < (sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][5] - sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][4] + 1); sz++){
		      
		      temp_RAPV[(((sz + sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][4] - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetFREQmin()) * (detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmin() + 1)) + sx + sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][0] - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmin())]+=mini_RAPV[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][((sz * (sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][1] - sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][0] + 1)) + sx)];
		      
		    }
		    
		  }

		  // d. DECPV
		  for(sy = 0; sy < (sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][3] - sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][2] + 1); sy++){
		    
		    for(sz = 0; sz < (sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][5] - sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][4] + 1); sz++){
		      
		      temp_DECPV[(((sz + sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][4] - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetFREQmin()) * (detections[obj_batch][(existing - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetDECmin() + 1)) + sy + sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][2] - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetDECmin())]+=mini_DECPV[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][((sz * (sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][3] - sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][2] + 1)) + sy)];

		    }
		    
		  }
		  		  
		  // e. ref_spec
		  if((sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][5] - sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][4] + 1) >= 10){
		    
		    for(sz = 0; sz < (2 * (sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][5] - sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][4] + 1)); sz++){ 
		      
		      if((sz + sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][4] - detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].GetFREQmin()) >= 0){
			
			temp_ref_spec[(sz + sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][4] - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetFREQmin())]+=mini_ref_spec[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][sz]; 
			
		      }
		      
		    } 
		    
		  } else {
		    
		    for(sz = 0; sz < (sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][5] - sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][4] + 11); sz++){ 
		      
		      if((sz + sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][4] - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetFREQmin()) >= 0){
			
			temp_ref_spec[(sz + sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][4] - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetFREQmin())]+=mini_ref_spec[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][sz]; 
			
		      }
		      
		    } 		    
		    
		  }
		  		  
		  // f. obj_spec
		  for(sz = 0; sz < (sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][5] - sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][4] + 1); sz++){
		    
		    temp_obj_spec[(sz + sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][4] - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetFREQmin())]+=mini_obj_spec[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][sz];
		    
		  }
		  
		  // convert temp_sparse_reps_grid from differential to cumulative counts using temp_sparse_reps_string as an intermediary
		  temp_sparse_reps_strings[0] = 0;
		  for(g = 1; g < (((detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmin() + 1) * (detections[obj_batch][(existing - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetDECmin() + 1)) + 1); g++){ temp_sparse_reps_strings[g] = temp_sparse_reps_strings[(g - 1)] + temp_sparse_reps_grid[(g - 1)]; }
		  for(g = 0; g < (((detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmin() + 1) * (detections[obj_batch][(existing - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetDECmin() + 1)) + 1); g++){ temp_sparse_reps_grid[g] = temp_sparse_reps_strings[g]; }
		  
		  // write new sparse_reps_strings value to temp_sparse_reps_strings array, using various grids to achieve indexing
		  
		  // a. initialise temp_sparse_reps_strings
		  for(g = 0; g < (2 * temp_sparse_reps_grid[((detections[obj_batch][(existing - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetDECmin() + 1) * (detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmin() + 1))]); g++){ temp_sparse_reps_strings[g] = 0; }
		  

		  // b. for each line of sight through the new existing object bounding box, retrieve the channel range of each object string along this LoS
		  for(sy = detections[obj_batch][(existing - (obj_batch * obj_limit))].GetDECmin(); sy <= detections[obj_batch][(existing - (obj_batch * obj_limit))].GetDECmax(); sy++){
		    
		    for(sx = detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmin(); sx <= detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmax(); sx++){
		      
		      // initialise the number of object strings written to this LoS
		      k = 0;

		      // retrieve the starting index for object strings along this LoS
		      j = 2 * temp_sparse_reps_grid[(((sy - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetDECmin()) * (detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmin() + 1)) + sx - detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmin())];
		      		      
		      // write existing object's object strings to temp_strings_array
		      if((sparse_reps_update[obj_batch][(existing - (obj_batch * obj_limit))] != 0) && (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] >= 0)){
			
			if((sx >= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0]) && (sx <= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][1]) && (sy >= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][2]) && (sy <= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][3])){
			  
			  for(g = sparse_reps_grid[obj_batch][(existing - (obj_batch * obj_limit))][(((sy - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][2]) * (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][1] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] + 1)) + sx - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0])]; g < sparse_reps_grid[obj_batch][(existing - (obj_batch * obj_limit))][(((sy - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][2]) * (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][1] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] + 1)) + sx - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] + 1)]; g++){
			    
			    temp_sparse_reps_strings[(j + k)] = sparse_reps_strings[obj_batch][(existing - (obj_batch * obj_limit))][(2 * g)];
			    k++;
			    temp_sparse_reps_strings[(j + k)] = sparse_reps_strings[obj_batch][(existing - (obj_batch * obj_limit))][((2 * g) + 1)];
			    k++;
			    
			  }
			  
			  // if((sx >= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0]) && (sx <= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][1]) && (sy >= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][2]) && (sy <= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][3]))
			}
			
			// if((sparse_reps_update[obj_batch][(existing - (obj_batch * obj_limit))] != 0) && (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] >= 0))
		      }
		      
		      // write match_init[i] object strings to temp_strings_array
		      if((sx >= sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][0]) && (sx <= sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][1]) && (sy >= sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][2]) && (sy <= sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][3])){
			
			for(g = sparse_reps_grid[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][(((sy - sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][2]) * (sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][1] - sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][0] + 1)) + sx - sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][0])]; g < sparse_reps_grid[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][(((sy - sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][2]) * (sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][1] - sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][0] + 1)) + sx - sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][0] + 1)]; g++){
			  
			  temp_sparse_reps_strings[(j + k)] = sparse_reps_strings[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][(2 * g)];
			  k++;
			  temp_sparse_reps_strings[(j + k)] = sparse_reps_strings[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][((2 * g) + 1)];
			  k++;
			  
			}
			
			// if((sx >= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0]) && (sx <= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][1]) && (sy >= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][2]) && (sy <= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][3]))
		      }
		      
		      // for(sx = 0; sx < (); sx++)
		    }
		    
		    // for(sy = 0; sy < (); sy++)
		  }
		  		  
		  // over-write the existing object's sparse representations with the existing+merged sparse representations
		  sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] = detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmin();
		  sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][1] = detections[obj_batch][(existing - (obj_batch * obj_limit))].GetRAmax();
		  sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][2] = detections[obj_batch][(existing - (obj_batch * obj_limit))].GetDECmin();
		  sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][3] = detections[obj_batch][(existing - (obj_batch * obj_limit))].GetDECmax();
		  sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][4] = detections[obj_batch][(existing - (obj_batch * obj_limit))].GetFREQmin();
		  sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][5] = detections[obj_batch][(existing - (obj_batch * obj_limit))].GetFREQmax();
		  
		  // a. grid
		  delete [] sparse_reps_grid[obj_batch][(existing - (obj_batch * obj_limit))];
		  sparse_reps_grid[obj_batch][(existing - (obj_batch * obj_limit))] = new int[(((sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][1] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] + 1) * (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][3] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][2] + 1)) + 1)];
		  for(g = 0; g < (((sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][1] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] + 1) * (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][3] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][2] + 1)) + 1); g++){ sparse_reps_grid[obj_batch][(existing - (obj_batch * obj_limit))][g] = temp_sparse_reps_grid[g]; }
		  
		  // b. mini_mom0
		  delete [] mini_mom0[obj_batch][(existing - (obj_batch * obj_limit))];
		  mini_mom0[obj_batch][(existing - (obj_batch * obj_limit))] = new float[((sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][1] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] + 1) * (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][3] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][2] + 1))];
		  for(g = 0; g < ((sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][1] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] + 1) * (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][3] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][2] + 1)); g++){ mini_mom0[obj_batch][(existing - (obj_batch * obj_limit))][g] = temp_mom0[g]; }
		  
		  // c. mini_RAPV
		  delete [] mini_RAPV[obj_batch][(existing - (obj_batch * obj_limit))];
		  mini_RAPV[obj_batch][(existing - (obj_batch * obj_limit))] = new float[((sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][1] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] + 1) * (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][4] + 1))];
		  for(g = 0; g < ((sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][1] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] + 1) * (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][4] + 1)); g++){ mini_RAPV[obj_batch][(existing - (obj_batch * obj_limit))][g] = temp_RAPV[g]; }
		  
		  // d. mini_DECPV
		  delete [] mini_DECPV[obj_batch][(existing - (obj_batch * obj_limit))];
		  mini_DECPV[obj_batch][(existing - (obj_batch * obj_limit))] = new float[((sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][4] + 1) * (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][3] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][2] + 1))];
		  for(g = 0; g < ((sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][3] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][2] + 1) * (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][4] + 1)); g++){ mini_DECPV[obj_batch][(existing - (obj_batch * obj_limit))][g] = temp_DECPV[g]; }
		  
		  // e. mini_obj_spec
		  delete [] mini_obj_spec[obj_batch][(existing - (obj_batch * obj_limit))];
		  mini_obj_spec[obj_batch][(existing - (obj_batch * obj_limit))] = new float[(sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][4] + 1)];
		  for(g = 0; g < (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][4] + 1); g++){ mini_obj_spec[obj_batch][(existing - (obj_batch * obj_limit))][g] = temp_obj_spec[g]; }
		  
		  // f. mini_ref_spec
		  delete [] mini_ref_spec[obj_batch][(existing - (obj_batch * obj_limit))];
		  if((sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][4] + 1) >= 10){
		    
		    mini_ref_spec[obj_batch][(existing - (obj_batch * obj_limit))] = new float[(2 * (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][4] + 1))];
		    for(g = 0; g < (2 * (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][4] + 1)); g++){ mini_ref_spec[obj_batch][(existing - (obj_batch * obj_limit))][g] = temp_ref_spec[g]; }
		    
		  } else {
		    
		    mini_ref_spec[obj_batch][(existing - (obj_batch * obj_limit))] = new float[(sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][4] + 11)];
		    for(g = 0; g < (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][4] + 11); g++){ mini_ref_spec[obj_batch][(existing - (obj_batch * obj_limit))][g] = temp_ref_spec[g]; }
		    
		  }
		  
		  // g. sparse_reps_strings
		  delete [] sparse_reps_strings[obj_batch][(existing - (obj_batch * obj_limit))];
		  sparse_reps_strings[obj_batch][(existing - (obj_batch * obj_limit))] = new int[(2 * sparse_reps_grid[obj_batch][(existing - (obj_batch * obj_limit))][((sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][1] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] + 1) * (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][3] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][2] + 1))])];
		  for(g = 0; g < (2 * sparse_reps_grid[obj_batch][(existing - (obj_batch * obj_limit))][((sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][1] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] + 1) * (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][3] - sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][2] + 1))]); g++){ sparse_reps_strings[obj_batch][(existing - (obj_batch * obj_limit))][g] = temp_sparse_reps_strings[g]; }
		  
		}
		
		// re-initialise the `object' values for the match_init[i] object that has just been merged into the existing object, 
		// and push the obj value (match_init[i]) to the top of the array of values
		detections[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))].ReInit();
		if((sparse_reps_update[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))] != 0) && (sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][0] >= 0)){
		  
		  delete [] sparse_reps_grid[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))];
		  sparse_reps_grid[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))] = NULL;
		  delete [] sparse_reps_strings[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))];
		  sparse_reps_strings[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))] = NULL;
		  delete [] mini_mom0[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))];
		  mini_mom0[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))] = NULL;
		  delete [] mini_RAPV[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))];
		  mini_RAPV[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))] = NULL;
		  delete [] mini_DECPV[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))];
		  mini_DECPV[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))] = NULL;
		  delete [] mini_obj_spec[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))];
		  mini_obj_spec[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))] = NULL;
		  delete [] mini_ref_spec[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))];
		  mini_ref_spec[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))] = NULL;

		}
		sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][0] = -99;
		sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][1] = -99;
		sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][2] = -99;
		sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][3] = -99;
		sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][4] = -99;
		sparse_reps_size[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))][5] = -99;
		sparse_reps_update[obj_batch_2][(match_init[i] - (obj_batch_2 * obj_limit))] = 0;

		j = -1;
		for(g = 0; g < NO_obj_ids; g++){ if(obj_ids[g] == match_init[i]){ j = 1; break; } }
		if(j == -1){
		  
		  obj_ids[NO_obj_ids] = match_init[i];
		  NO_obj_ids++;
		  		  
		}
		
		// for(i = 0; i < NOi; i++)
	      }
	      
	      // if(NOi > 1)
	    }
	    
	    // else . . . if(existing == -1)
	  }
	  
	  // if(flag_vals[((z*size_x*size_y) + (y*size_x) + x)] == -1)
	}
	
	// update progress on display
	while(progress <= (((float) ((z * size_x * size_y) + (y * size_x) + x + 1)) / ((float) (size_x * size_y * size_z)))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
	
	// for(x = 0; x < size_x; x++)
      }
      
      // for(y = 0; y < size_y; y++)
    }
    
    // for(z = 0; z < size_z; z++)
  }
  std::cout << "* done." << std::endl;
  
  // 4. generate sparse representations of coherent objects that are retained after
  // this `chunk' of the datacube is processed
  std::cout << "Generating/updating sparse representations of sources . . . " << std::endl;
  std::cout << "0 | |:| | : | |:| | 100% complete" << std::endl;
  progress = 0.0;
  for(i = 0; i < obj; i++){
    
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
    if((sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][0] >= 0) && (sparse_reps_update[obj_batch][(i - (obj_batch * obj_limit))] == -1)){
      
      // update progress on display
      while(progress <= (((float) (i + 1) / ((float) obj)))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }  
      continue; 
	
    }
          
    // update the sparse representation for this object
    
    // 0. if there is no existing sparse representation and postage stamp images, create them from scratch, 
    // otherwise concatenate the existing sparse representation and postage stamp images with the updates
    if((sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][0] < 0) && (sparse_reps_update[obj_batch][(i - (obj_batch * obj_limit))] >= 0)){
      
      // branch for creating a brand new sparse representation
      
      // update the sparse_reps_update array
      sparse_reps_update[obj_batch][(i - (obj_batch * obj_limit))] = -1;

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
      for(sy = sy_start; (sy <= sy_finish) && (sy < size_y); sy++){
	
	for(sx = sx_start; (sx <= sx_finish) && (sx < size_x); sx++){
	  
	  // initialise the dummy integer j
	  j = -1;
	  
	  for(sz = sz_start; (sz <= sz_finish) && (sz < size_z); sz++){
	    
	    // change j to reflect if this voxel in the flag_vals array belongs to the source
	    if(flag_vals[((sz * size_x * size_y) + (sy * size_x) + sx)] == i){ 
	      
	      // if previously this was not part of an object string, j == -1, then increment the number of object strings and
	      // adjust the value of j
	      if(j == -1){ NOi++; }
	      
	      // change j to reflect that this voxel belongs to the object
	      j = 1; 
	      
	    } else { j = -1; }
	    
	    // for(sz = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetFREQmin() - chunk_z_size; sz <= (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetFREQmax() - chunk_z_size); sz++)
	  }
	  
	  // for(sx = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmin() - chunk_x_size; sx <= (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmax() - chunk_x_size); sx++)
	}
	
	// for(sy = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmin() - chunk_y_size; sy <= (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmax() - chunk_y_size); sy++)
      }

      // create an array element to store this sparse representation
      sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][0] = detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmin();
      sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][1] = detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmax();
      sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][2] = detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmin();
      sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][3] = detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmax();
      sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][4] = detections[obj_batch][(i - (obj_batch * obj_limit))].GetFREQmin();
      sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][5] = detections[obj_batch][(i - (obj_batch * obj_limit))].GetFREQmax();
      while(sparse_reps_grid[obj_batch][(i - (obj_batch * obj_limit))] == NULL){ sparse_reps_grid[obj_batch][(i - (obj_batch * obj_limit))] = new int [(((detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmin() + 1) * (detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmin() + 1)) + 1)]; }
      while(sparse_reps_strings[obj_batch][(i - (obj_batch * obj_limit))] == NULL){ sparse_reps_strings[obj_batch][(i - (obj_batch * obj_limit))] = new int [(2 * NOi)]; }
            
      // calculate the extent of the sparse_rep volume
      sx_finish = detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmin() + 1;
      sy_finish = detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmin() + 1;
      sz_finish = detections[obj_batch][(i - (obj_batch * obj_limit))].GetFREQmax() - detections[obj_batch][(i - (obj_batch * obj_limit))].GetFREQmin() + 1;

      // create arrays to store the postage stamp images and initialise them
      while(mini_mom0[obj_batch][(i - (obj_batch * obj_limit))] == NULL){ mini_mom0[obj_batch][(i - (obj_batch * obj_limit))] = new float [(sx_finish * sy_finish)]; }
      while(mini_RAPV[obj_batch][(i - (obj_batch * obj_limit))] == NULL){ mini_RAPV[obj_batch][(i - (obj_batch * obj_limit))] = new float [(sx_finish * sz_finish)]; }
      while(mini_DECPV[obj_batch][(i - (obj_batch * obj_limit))] == NULL){ mini_DECPV[obj_batch][(i - (obj_batch * obj_limit))] = new float [(sy_finish * sz_finish)]; }
      while(mini_obj_spec[obj_batch][(i - (obj_batch * obj_limit))] == NULL){ mini_obj_spec[obj_batch][(i - (obj_batch * obj_limit))] = new float [sz_finish]; }
      if(sz_finish >= 10){

	while(mini_ref_spec[obj_batch][(i - (obj_batch * obj_limit))] == NULL){ mini_ref_spec[obj_batch][(i - (obj_batch * obj_limit))] = new float [(2 * sz_finish)]; }  

      } else {

	while(mini_ref_spec[obj_batch][(i - (obj_batch * obj_limit))] == NULL){ mini_ref_spec[obj_batch][(i - (obj_batch * obj_limit))] = new float [(10 + sz_finish)]; } 

      }

      // initialise the mini_mom0 values
      for(sy = 0; sy < sy_finish; sy++){
	for(sx = 0; sx < sx_finish; sx++){
	  mini_mom0[obj_batch][(i - (obj_batch * obj_limit))][((sy * sx_finish) + sx)] = 0.0;
	}
      }

      // initialise the mini_RAPV values
      for(sz = 0; sz < sz_finish; sz++){
	for(sx = 0; sx < sx_finish; sx++){
	  mini_RAPV[obj_batch][(i - (obj_batch * obj_limit))][((sz * sx_finish) + sx)] = 0.0;
	}
      }
      
      // initialise the mini_DECPV values
      for(sz = 0; sz < sz_finish; sz++){
	for(sy = 0; sy < sy_finish; sy++){
	  mini_DECPV[obj_batch][(i - (obj_batch * obj_limit))][((sz * sy_finish) + sy)] = 0.0;
	}
      }
      
      // initialise the mini_spec values
      for(sz = 0; sz < sz_finish; sz++){ mini_obj_spec[obj_batch][(i - (obj_batch * obj_limit))][sz] = 0.0; }
      if(sz_finish >= 10){

	for(sz = 0; sz < (2 * sz_finish); sz++){ mini_ref_spec[obj_batch][(i - (obj_batch * obj_limit))][sz] = 0.0; }

      } else {

	for(sz = 0; sz < (10 + sz_finish); sz++){ mini_ref_spec[obj_batch][(i - (obj_batch * obj_limit))][sz] = 0.0; }

      }
            
      // populate the sparse array representation and postage stamp images
      NOi = 0;
      for(sy = 0; sy < sy_finish; sy++){
	
	for(sx = 0; sx < sx_finish; sx++){
	  
	  // write the current number of object strings to the grid
	  sparse_reps_grid[obj_batch][(i - (obj_batch * obj_limit))][((sy * sx_finish) + sx)] = NOi;
	  
	  // initialise the dummy integer j 
	  j = -1;
	  
	  for(sz = 0; sz < sz_finish; sz++){
	    
	    // change j to reflect if this voxel in the flag_vals array belongs to the source
	    if(flag_vals[(((sz + sz_start) * size_x * size_y) + ((sy + sy_start) * size_x) + sx + sx_start)] == i){ 
	      
	      // update mini_mom0
	      mini_mom0[obj_batch][(i - (obj_batch * obj_limit))][((sy * sx_finish) + sx)]+=data_vals[(((sz + sz_start) * size_x * size_y) + ((sy + sy_start) * size_x) + sx + sx_start)];
	      
	      // update mini_RAPV
	      mini_RAPV[obj_batch][(i - (obj_batch * obj_limit))][((sz * sx_finish) + sx)]+=data_vals[(((sz + sz_start) * size_x * size_y) + ((sy + sy_start) * size_x) + sx + sx_start)];
	      
	      // update mini_DECPV
	      mini_DECPV[obj_batch][(i - (obj_batch * obj_limit))][((sz * sy_finish) + sy)]+=data_vals[(((sz + sz_start) * size_x * size_y) + ((sy + sy_start) * size_x) + sx + sx_start)];
	      
	      // update mini_obj_spec
	      mini_obj_spec[obj_batch][(i - (obj_batch * obj_limit))][sz]+=data_vals[(((sz + sz_start) * size_x * size_y) + ((sy + sy_start) * size_x) + sx + sx_start)];
  
	      // if previously this was not part of an object string, j == -1, then increment the number of object strings and
	      // write the position to the sparse_reps array
	      if(j == -1){ 
		
		// write the beginning of the object string to sparse_reps
		sparse_reps_strings[obj_batch][(i - (obj_batch * obj_limit))][(2 * NOi)] = sz_start + chunk_z_start + sz;
		
		// increment NOi
		NOi++; 
		
		// initialise the end of the object string in sparse_reps
		sparse_reps_strings[obj_batch][(i - (obj_batch * obj_limit))][((2 * NOi) - 1)] = sz_start + chunk_z_start + sz;

	      } else {
		
		// update the end of the object string in sparse_reps
		sparse_reps_strings[obj_batch][(i - (obj_batch * obj_limit))][((2 * NOi) - 1)] = sz_start + chunk_z_start + sz;
		
	      }
	      
	      // change j to reflect that this voxel belongs to the object
	      j = 1; 
	      
	    } else { 
	      
	      // change j to reflect that this voxel doesn't belong to the object
	      j = -1; 
	      
	    }
	    
	    // for(sz = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetFREmin() - chunk_z_start + 1; sz <= (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetFREQmax() - chunk_z_start + 1); sz++)
	  }
	  
	  // for(sx = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmin() - chunk_x_start + 1; sx <= (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetRAmax() - chunk_x_start + 1); sx++)
	}
	
	//for(sy = detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmin() - chunk_y_start + 1; sy <= (detections[obj_batch][(obj - (obj_batch * obj_limit))].GetDECmax() - chunk_y_start + 1); sy++)
      }
            
      // update sparse_reps with the final value of NOi
      sparse_reps_grid[obj_batch][(i - (obj_batch * obj_limit))][(sx_finish * sy_finish)] = NOi;
            
      // create the reference spectrum
      if(sz_finish >= 10){

	for(sz = -1 * (int) floorf((0.5 * (float) sz_finish)); sz < (sz_finish + ((int) floorf((0.5 * (float) sz_finish)))); sz++){
	  
	  for(sy = 0; sy < sy_finish; sy++){
	  
	    for(sx = 0; sx < sx_finish; sx++){
	    
	      if(((sz + sz_start) >= 0) && ((sz + sz_start) < size_z)){
		
		// update mini_ref_spec
		mini_ref_spec[obj_batch][(i - (obj_batch * obj_limit))][(sz + ((int) floorf((0.5 * (float) sz_finish))))]+=data_vals[(((sz + sz_start) * size_x * size_y) + ((sy + sy_start) * size_x) + sx + sx_start)];	
		
	      }
	      
	    }
	    
	  }
      
	}

	// if(sz_finish >= 10)
      } else {

	for(sz = -5; sz < (sz_finish + 5); sz++){
	  
	  for(sy = 0; sy < sy_finish; sy++){
	  
	    for(sx = 0; sx < sx_finish; sx++){
	    
	      if(((sz + sz_start) >= 0) && ((sz + sz_start) < size_z)){
		
		// update mini_ref_spec
		mini_ref_spec[obj_batch][(i - (obj_batch * obj_limit))][(sz + 5)]+=data_vals[(((sz + sz_start) * size_x * size_y) + ((sy + sy_start) * size_x) + sx + sx_start)];	
		
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
      for(sy = sy_start; sy <= sy_finish; sy++){
	
	for(sx = sx_start; sx <= sx_finish; sx++){
	  	  
	  // case 3.a triggered
	  if(((sx > merge_x) || (chunk_x_start == 0)) && ((sy > merge_y) || (chunk_y_start == 0))){
	    	    
	    // initialise the dummy integer j 
	    j = -1;
	    
	    for(sz = sz_start; sz <= sz_finish; sz++){
	      
	      // change i to reflect if this voxel in the flag_vals array belongs to the source
	      if(flag_vals[((sz * size_x * size_y) + (sy * size_x) + sx)] == i){ 
		
		// if previously this was not part of an object string, i == -1, then increment the number of object strings and
		// write the position to the sparse_reps array
		if(j == -1){ 
		  
		  // increment NOi
		  NOi++; 
		  
		} 
		
		// change j to reflect that this voxel belongs to the object
		j = 1; 
		
	      } else { j = -1; }
	      
	      // for(sz = sz_start; sz <= sz_finish; sz++)
	    }
	    
	  } else {
	    
	    // case 3.b triggered  
	    NOi+=(sparse_reps_grid[obj_batch][(i - (obj_batch * obj_limit))][(((sy + chunk_y_start - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][2]) * (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][1] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][0] + 1)) + sx + chunk_x_start - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][0] + 1)] - sparse_reps_grid[obj_batch][(i - (obj_batch * obj_limit))][(((sy + chunk_y_start - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][2]) * (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][1] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][0] + 1)) + sx + chunk_x_start - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][0])]);
	    
	  }     
	  
	  // for(sx = sx_start; sx <= sx_finish; sx++)
	}
	
	// for(sy = sy_start; sy <= sy_finish; sy++)
      }
      
      // 4. allocate memory to store temporary sparse representation and postage stamp images, then
      //    initialise the postage stamp image arrays
      for(sy = 0; sy < (sy_finish - sy_start + 1); sy++){
	for(sx = 0; sx < (sx_finish - sx_start + 1); sx++){
	  temp_mom0[((sy * (sx_finish - sx_start + 1)) + sx)] = 0.0;
	}
      }
      for(sz = 0; sz < (sz_finish - sz_start + 1); sz++){
	for(sx = 0; sx < (sx_finish - sx_start + 1); sx++){
	  temp_RAPV[((sz * (sx_finish - sx_start + 1)) + sx)] = 0.0;
	}
      }
      for(sz = 0; sz < (sz_finish - sz_start + 1); sz++){
	for(sy = 0; sy < (sy_finish - sy_start + 1); sy++){
	  temp_DECPV[((sz * (sy_finish - sy_start + 1)) + sy)] = 0.0;
	}
      }
      for(sz = 0; sz < (sz_finish - sz_start + 1); sz++){ temp_obj_spec[sz] = 0.0; }
      if((sz_finish - sz_start + 1) >= 10){

	for(sz = 0; sz < (2 * (sz_finish - sz_start + 1)); sz++){ temp_ref_spec[sz] = 0.0; }

      } else {

	for(sz = 0; sz < (sz_finish - sz_start + 11); sz++){ temp_ref_spec[sz] = 0.0; }

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
      for(sz = 0; sz < (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][4] + 1); sz++){

	for(sx = 0; sx < (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][1] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][0] + 1); sx++){
	
	  temp_RAPV[(((sz + sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][4] - sz_start - chunk_z_start) * (sx_finish - sx_start + 1)) + sx + sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][0] - sx_start - chunk_x_start)]+=mini_RAPV[obj_batch][(i - (obj_batch * obj_limit))][((sz * (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][1] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][0] + 1)) + sx)];
	  
	}

      }

      // initialise Dec PV image
      for(sz = 0; sz < (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][4] + 1); sz++){

	for(sy = 0; sy < (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][3] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][2] + 1); sy++){
	
	  temp_DECPV[(((sz + sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][4] - sz_start - chunk_z_start) * (sy_finish - sy_start + 1)) + sy + sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][2] - sy_start - chunk_y_start)]+=mini_DECPV[obj_batch][(i - (obj_batch * obj_limit))][((sz * (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][3] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][2] + 1)) + sy)];
	  
	}

      }

      // initialise mini spec's
      for(sz = 0; sz < (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][4] + 1); sz++){ temp_obj_spec[(sz + sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][4] - sz_start - chunk_z_start)]+=mini_obj_spec[obj_batch][(i - (obj_batch * obj_limit))][sz]; }

      if((sz_finish - sz_start + 1) >= 10){

	j = sz_start + chunk_z_start - ((int) floorf((0.5 * (float) (sz_finish - sz_start + 1))));

	if((sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][4] + 1) >= 10){

	  k = sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][4] - ((int) floorf((0.5 * (float) (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][4] + 1))));

	  for(sz = 0; sz < (2 * (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][4] + 1)); sz++){

	    temp_ref_spec[(sz + k - j)]+=mini_ref_spec[obj_batch][(i - (obj_batch * obj_limit))][sz];

	  }

	} else {
	  
	  k = sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][4] - 5;

	  for(sz = 0; sz < (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][4] + 11); sz++){
	    
	    temp_ref_spec[(sz + k - j)]+=mini_ref_spec[obj_batch][(i - (obj_batch * obj_limit))][sz];

	  }
	  
	}

      } else {

	j = sz_start + chunk_z_start - 5;

	if((sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][4] + 1) >= 10){

	  k = sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][4] - ((int) floorf((0.5 * (float) (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][4] + 1))));

	  for(sz = 0; sz < (2 * (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][4] + 1)); sz++){

	    temp_ref_spec[(sz + k - j)]+=mini_ref_spec[obj_batch][(i - (obj_batch * obj_limit))][sz];

	  }

	} else {
	  
	  k = sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][4] - 5;

	  for(sz = 0; sz < (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][4] + 11); sz++){
	    
	    temp_ref_spec[(sz + k - j)]+=mini_ref_spec[obj_batch][(i - (obj_batch * obj_limit))][sz];

	  }
	  
	}

      }
      
      // search through LoS's of updated bounding box and take action as dictated by case 5.a or 5.b
      for(sy = sy_start; sy <= sy_finish; sy++){
	
	for(sx = sx_start; sx <= sx_finish; sx++){
	  	  
	  // add value to temp_sparse_rep_grid
	  temp_sparse_reps_grid[(((sy - sy_start) * (sx_finish - sx_start + 1)) + sx - sx_start)] = NOi;
	  
	  // case 5.a triggered for moment-0 image and sparse representation
	  if(((sx > merge_x) || (chunk_x_start == 0)) && ((sy > merge_y) || (chunk_y_start == 0))){
	    
	    // initialise the dummy integer j 
	    j = -1;
	    
	    for(sz = sz_start; sz <= sz_finish; sz++){
	      
	      // change j to reflect if this voxel in the flag_vals array belongs to the source
	      if(flag_vals[((sz * size_x * size_y) + (sy * size_x) + sx)] == i){ 
		
		// update mini_mom0
		temp_mom0[(((sy - sy_start) * (sx_finish - sx_start + 1)) + sx - sx_start)]+=data_vals[((sz * size_x * size_y) + (sy * size_x) + sx)];
		
		// update mini_RAPV
		temp_RAPV[(((sz - sz_start) * (sx_finish - sx_start + 1)) + sx - sx_start)]+=data_vals[((sz * size_x * size_y) + (sy * size_x) + sx)];
		
		// update mini_DECPV
		temp_DECPV[(((sz - sz_start) * (sy_finish - sy_start + 1)) + sy - sy_start)]+=data_vals[((sz * size_x * size_y) + (sy * size_x) + sx)];
		
		// update mini_obj_spect
		temp_obj_spec[(sz - sz_start)]+=data_vals[((sz * size_x * size_y) + (sy * size_x) + sx)];
		
		// if previously this was not part of an object string, j == -1, then increment the number of object strings and
		// write the position to the sparse_reps array
		if(j == -1){ 
		  
		  // write the beginning of the object string to sparse_reps
		  temp_sparse_reps_strings[(2 * NOi)] = chunk_z_start + sz;
		  
		  // increment NOi
		  NOi++; 
		  		  
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
	      
	      // for(sz = sz_start; sz <= sz_finish; sz++)
	    }
	    
	    // add values to temp_ref_spec
	    if((sz_finish - sz_start + 1) >= 10){
	      
	      k = (int) floorf((0.5 * (float) (sz_finish - sz_start + 1)));
	      
	      for(sz = 0; sz < (2 * (sz_finish - sz_start + 1)); sz++){
	  	
		if(((sz + sz_start - k) >= 0) && ((sz + sz_start - k) < size_z)){
		  
		  // update mini_ref_spec
		  temp_ref_spec[(sz + sz_start - k)]+=data_vals[(((sz + sz_start - k) * size_x * size_y) + (sy * size_x) + sx)];	
		  
		}
		
	      }
	      
	      // if(sz_finish >= 10)
	    } else {
	      
	      for(sz = 0; sz < (sz_finish - sz_start + 11); sz++){
		
		if(((sz + sz_start - 5) >= 0) && ((sz + sz_start - 5) < size_z)){
		  
		  // update mini_ref_spec
		  temp_ref_spec[(sz + sz_start - 5)]+=data_vals[(((sz + sz_start - 5) * size_x * size_y) + (sy * size_x) + sx)];	
		  
		}
		
	      }
	      
	      // else . . . if(sz_finish >= 10)
	    }
	    
	  } else {
	    // case 5.b triggered  
	    	    
	    // update temp_sparse_rep_strings
	    j = 0;
	    if(((sx + chunk_x_start) >= sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][0]) && ((sx + chunk_x_start) <= (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][1] + 1)) && ((sy + chunk_y_start) >= sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][2]) && ((sy + chunk_y_start) <= (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][3] + 1))){
	      
	      for(g = sparse_reps_grid[obj_batch][(i - (obj_batch * obj_limit))][(((sy + chunk_y_start - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][2]) * (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][1] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][0] + 1)) + sx + chunk_x_start - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][0])]; g < sparse_reps_grid[obj_batch][(i - (obj_batch * obj_limit))][(((sy + chunk_y_start - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][2]) * (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][1] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][0] + 1)) + sx + chunk_x_start - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][0] + 1)]; g++){
		
		temp_sparse_reps_strings[(NOi + j)] = sparse_reps_strings[obj_batch][(i - (obj_batch * obj_limit))][(2 * g)];
		temp_sparse_reps_strings[(NOi + j + 1)] = sparse_reps_strings[obj_batch][(i - (obj_batch * obj_limit))][((2 * g) + 1)];
		j++;
		
	      }
	      
	      // update moment 0 postage stamp image
	      temp_mom0[(((sy - sy_start) * (detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmin() + 1)) + sx - sx_start)] = mini_mom0[obj_batch][(i - (obj_batch * obj_limit))][(((sy + chunk_y_start - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][2]) * (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][1] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][0] + 1)) + sx + chunk_x_start - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][0])];
	      
	      // update NOi
	      NOi+=(sparse_reps_grid[obj_batch][(i - (obj_batch * obj_limit))][(((sy + chunk_y_start - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][2]) * (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][1] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][0] + 1)) + sx + chunk_x_start - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][0] + 1)] - sparse_reps_grid[obj_batch][(i - (obj_batch * obj_limit))][(((sy + chunk_y_start - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][2]) * (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][1] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][0] + 1)) + sx + chunk_x_start - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][0])]);
	      	      
	    }
	    
	    // else . . . case 5.b
	  }
	    
	  // for(sx = sx_start; sx <= sx_finish; sx++)
	}
		
	// for(sy = sy_start; sy <= sy_finish; sy++)
      }
            
      temp_sparse_reps_grid[((sx_finish - sx_start + 1) * (sy_finish - sy_start + 1))] = NOi;

      // 6. replace the current sparse representation with the new sparse representation
      delete [] sparse_reps_grid[obj_batch][(i - (obj_batch * obj_limit))];
      sparse_reps_grid[obj_batch][(i - (obj_batch * obj_limit))] = NULL;
      delete [] sparse_reps_strings[obj_batch][(i - (obj_batch * obj_limit))];
      sparse_reps_strings[obj_batch][(i - (obj_batch * obj_limit))] = NULL;
      delete [] mini_mom0[obj_batch][(i - (obj_batch * obj_limit))];
      mini_mom0[obj_batch][(i - (obj_batch * obj_limit))] = NULL;
      delete [] mini_RAPV[obj_batch][(i - (obj_batch * obj_limit))];
      mini_RAPV[obj_batch][(i - (obj_batch * obj_limit))] = NULL;
      delete [] mini_DECPV[obj_batch][(i - (obj_batch * obj_limit))];
      mini_DECPV[obj_batch][(i - (obj_batch * obj_limit))] = NULL;
      delete [] mini_obj_spec[obj_batch][(i - (obj_batch * obj_limit))];
      mini_obj_spec[obj_batch][(i - (obj_batch * obj_limit))] = NULL;
      delete [] mini_ref_spec[obj_batch][(i - (obj_batch * obj_limit))];
      mini_ref_spec[obj_batch][(i - (obj_batch * obj_limit))] = NULL;
 
      // update the sparse_reps_size array with the new bounding box
      sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][0] = detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmin();
      sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][1] = detections[obj_batch][(i - (obj_batch * obj_limit))].GetRAmax();
      sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][2] = detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmin();
      sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][3] = detections[obj_batch][(i - (obj_batch * obj_limit))].GetDECmax();
      sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][4] = detections[obj_batch][(i - (obj_batch * obj_limit))].GetFREQmin();
      sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][5] = detections[obj_batch][(i - (obj_batch * obj_limit))].GetFREQmax();

      // update the sparse_reps_update array
      sparse_reps_update[obj_batch][(i - (obj_batch * obj_limit))] = -1;
          
      // populate the new sparse_reps and postaage stamps
      sparse_reps_grid[obj_batch][(i - (obj_batch * obj_limit))] = new int [(1 + ((sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][1] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][0] + 1) * (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][3] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][2] + 1)))];
      sparse_reps_strings[obj_batch][(i - (obj_batch * obj_limit))] = new int [(2 * NOi)];
      mini_mom0[obj_batch][(i - (obj_batch * obj_limit))] = new float [((sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][1] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][0] + 1) * (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][3] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][2] + 1))];
      mini_RAPV[obj_batch][(i - (obj_batch * obj_limit))] = new float [((sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][1] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][0] + 1) * (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][4] + 1))];
      mini_DECPV[obj_batch][(i - (obj_batch * obj_limit))] = new float [((sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][3] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][2] + 1) * (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][4] + 1))];
  
      mini_obj_spec[obj_batch][(i - (obj_batch * obj_limit))] = new float [(sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][4] + 1)];
      if((sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][4] + 1) >= 10){
       
	mini_ref_spec[obj_batch][(i - (obj_batch * obj_limit))] = new float [(2 * (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][4] + 1))];
	
      } else {
	
	mini_ref_spec[obj_batch][(i - (obj_batch * obj_limit))] = new float [(sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][4] + 11)];

      }
   
      for(j = 0; j < (1 + ((sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][1] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][0] + 1) * (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][3] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][2] + 1))); j++){ sparse_reps_grid[obj_batch][(i - (obj_batch * obj_limit))][j] = temp_sparse_reps_grid[j]; }
      for(j = 0; j < (2 * NOi); j++){ sparse_reps_strings[obj_batch][(i - (obj_batch * obj_limit))][j] = temp_sparse_reps_strings[j]; }
      for(j = 0; j < ((sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][1] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][0] + 1) * (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][3] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][2] + 1)); j++){ mini_mom0[obj_batch][(i - (obj_batch * obj_limit))][j] = temp_mom0[j]; }
      for(j = 0; j < ((sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][1] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][0] + 1) * (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][4] + 1)); j++){ mini_RAPV[obj_batch][(i - (obj_batch * obj_limit))][j] = temp_RAPV[j]; }
      for(j = 0; j < ((sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][3] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][2] + 1) * (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][4] + 1)); j++){ mini_DECPV[obj_batch][(i - (obj_batch * obj_limit))][j] = temp_DECPV[j]; }

      for(j = 0; j < (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][4] + 1); j++){ mini_obj_spec[obj_batch][(i - (obj_batch * obj_limit))][j] = temp_obj_spec[j]; }
      if((sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][4] + 1) >= 10){

	for(j = 0; j < (2 * (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][4] + 1)); j++){ mini_ref_spec[obj_batch][(i - (obj_batch * obj_limit))][j] = temp_ref_spec[j]; }
	
      } else {

	for(j = 0; j < (sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][5] - sparse_reps_size[obj_batch][(i - (obj_batch * obj_limit))][4] + 11); j++){ mini_ref_spec[obj_batch][(i - (obj_batch * obj_limit))][j] = temp_ref_spec[j]; }

      }
      
      // else . . . if(sparse_reps[obj_batch][(obj - (obj_batch * obj_limit))][0] < 0)
    }
    
    // update progress on display
    while(progress <= (((float) (i + 1) / ((float) obj)))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }  
   
    // for(i = 0; i < obj; i++)
  }
  std::cout << "* done." << std::endl;
  
  // 5. free up memory
  delete [] match_init;
  match_init = NULL;
  delete [] temp_sparse_reps_grid;
  temp_sparse_reps_grid = NULL;
  delete [] temp_sparse_reps_strings;
  temp_sparse_reps_strings = NULL;
  delete [] temp_mom0;
  temp_mom0 = NULL;
  delete [] temp_RAPV;
  temp_RAPV = NULL;
  delete [] temp_DECPV;
  temp_DECPV = NULL;
  delete [] temp_obj_spec;
  temp_obj_spec = NULL;
  delete [] temp_ref_spec;
  temp_ref_spec = NULL;
  
  // 6. return the number of coherent objects 
  return obj;
  
}

float CreateMoment0Map(float * plot_array, int NOobj, object_props ** detections, float *** mini_mom0, int NOx, int NOy, int obj_limit, int obj_batch_limit){
  
  int i,j,k,sx_finish,sy_finish,obj_batch;
  float max = -99.0;
  
  // initialise plot_array
  for(i = 0; i < NOx; i++){

    for(j = 0; j < NOy; j++){
      
      plot_array[((j*NOx) + i)] = 0.0;

    }
    
  }
  
  // search through objects and construct moment 0 map
  for(k = 0; k < NOobj; k++){
    
    obj_batch = floorf((float) k / (float) obj_limit);

    if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() > 1){

      sx_finish = detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin() + 1;
      sy_finish = detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin() + 1;

      for(i = 0; i < sx_finish; i++){
	
	for(j = 0; j < sy_finish; j++){

	  plot_array[(((j + detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin()) * NOx) + i + detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin())]+=mini_mom0[obj_batch][(k - (obj_batch * obj_limit))][((j * sx_finish) + i)];
	 
	}

      }

    }

    // for(k = 0; k < NOobj; k++)
  }

  for(i = 0; i < NOx; i++){
    for(j = 0; j < NOy; j++){
      if(plot_array[((j * NOx) + i)] > max){ max = plot_array[((j * NOx) + i)]; }
    }
  }

  return max;

}

float CreateRAPVPlot(float * plot_array, int NOobj, object_props ** detections, float *** mini_RAPV, int NOx, int NOz, int obj_limit, int obj_batch_limit){
  
  int i,j,k,sx_finish,sz_finish,obj_batch;
  float max = -99.0;
  
  // initialise plot_array
  for(i = 0; i < NOx; i++){

    for(j = 0; j < NOz; j++){
      
      plot_array[((j*NOx) + i)] = 0.0;

    }
    
  }
  
  // search through objects and construct moment 0 map
  for(k = 0; k < NOobj; k++){
    
    obj_batch = floorf((float) k / (float) obj_limit);

    if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() > 1){

      sx_finish = detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin() + 1;
      sz_finish = detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmin() + 1;

      for(i = 0; i < sx_finish; i++){
	
	for(j = 0; j < sz_finish; j++){

	  plot_array[(((j + detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmin()) * NOx) + i + detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin())]+=mini_RAPV[obj_batch][(k - (obj_batch * obj_limit))][((j * sx_finish) + i)];
	 
	}

      }

    }

    // for(k = 0; k < NOobj; k++)
  }

  for(i = 0; i < NOx; i++){
    for(j = 0; j < NOz; j++){
      if(plot_array[((j * NOx) + i)] > max){ max = plot_array[((j * NOx) + i)]; }
    }
  }

  return max;

}

float CreateDecPVPlot(float * plot_array, int NOobj, object_props ** detections, float *** mini_DECPV, int NOy, int NOz, int obj_limit, int obj_batch_limit){
  
  int i,j,k,sy_finish,sz_finish,obj_batch;
  float max = -99.0;
  
  // initialise plot_array
  for(i = 0; i < NOy; i++){

    for(j = 0; j < NOz; j++){
      
      plot_array[((j*NOy) + i)] = 0.0;

    }
    
  }
  
  // search through objects and construct moment 0 map
  for(k = 0; k < NOobj; k++){
    
    obj_batch = floorf((float) k / (float) obj_limit);

    if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() > 1){

      sy_finish = detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin() + 1;
      sz_finish = detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmin() + 1;

      for(i = 0; i < sy_finish; i++){
	
	for(j = 0; j < sz_finish; j++){

	  plot_array[(((j + detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmin()) * NOy) + i + detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin())]+=mini_DECPV[obj_batch][(k - (obj_batch * obj_limit))][((j * sy_finish) + i)];
	 
	}

      }

    }

    // for(k = 0; k < NOobj; k++)
  }

  for(i = 0; i < NOy; i++){
    for(j = 0; j < NOz; j++){
      if(plot_array[((j * NOy) + i)] > max){ max = plot_array[((j * NOy) + i)]; }
    }
  }

  return max;

}

int CreateMoment0Bounds(int *** sparse_reps_grid, int size_x, int size_y, int min_x, int min_y, float * plot_x, float * plot_y, int obj, int obj_limit){

  int i, j, found, p = 0, obj_batch;

  // raster scan through the cube, and identify whether the source exists along a particular line of sight within the bounding volume
  // if it is, add entry to boundary list

  // 0. calculate obj_batch value
  obj_batch = floorf((float) obj / (float) obj_limit);
    
  // 1. identify left edge of source
  for(j = 0; j < size_y; j++){

    for(i = 0; i < size_x; i++){

      if(sparse_reps_grid[obj_batch][(obj - (obj_batch * obj_limit))][((j * size_x) + i)] != sparse_reps_grid[obj_batch][(obj - (obj_batch * obj_limit))][((j * size_x) + i + 1)]){

	plot_x[p] = (float) min_x + i - 1;
	plot_y[p] = (float) min_y + j;
	p++;
	break;

      }

    }

  }

  // 2. add top left point 
  plot_x[p] = plot_x[(p - 1)];
  plot_y[p] = plot_y[(p - 1)] + 1;
  p++;

  // 3. identify top right edge of source
  for(i = size_x - 1; i >= 0; i--){

    if(sparse_reps_grid[obj_batch][(obj - (obj_batch * obj_limit))][(((size_y - 1) * size_x) + i)] != sparse_reps_grid[obj_batch][(obj - (obj_batch * obj_limit))][(((size_y - 1) * size_x) + i + 1)]){ 

      found = min_x + i + 1;
      break; 

    }

  }

  // 4. add top right point, then just identified top right edge of source
  plot_x[p] = (float) found;
  plot_y[p] = (float) min_y + size_y;
  p++;
  plot_x[p] = (float) found;
  plot_y[p] = (float) min_y + size_y - 1;
  p++;

  // 5. identify remaining right edge
  for(j = size_y - 2; j >= 0; j--){

    for(i = size_x - 1; i >= 0; i--){

      if(sparse_reps_grid[obj_batch][(obj - (obj_batch * obj_limit))][((j * size_x) + i)] != sparse_reps_grid[obj_batch][(obj - (obj_batch * obj_limit))][((j * size_x) + i + 1)]){

	plot_x[p] = (float) min_x + i + 1;
	plot_y[p] = (float) min_y + j;
	p++;
	break;
	
      }
     
    }

  }

  // 6. add bottom two points
  plot_x[p] = plot_x[(p - 1)];
  plot_y[p] = plot_y[(p - 1)] - 1;
  p++;
  plot_x[p] = plot_x[0];
  plot_y[p] = plot_y[0] - 1;
  p++;

  // 7. add original boundary point to complete 
  plot_x[p] = plot_x[0];
  plot_y[p] = plot_y[0];
  p++;

  // 8. check that all of the boundary points are within the bounding box limits
  for(i = 0; i < p; i++){

    if(plot_x[i] < (min_x - 1)){ plot_x[i] = min_x - 1; }
    if(plot_y[i] < (min_y - 1)){ plot_y[i] = min_y - 1; }
    if(plot_x[i] > (min_x + size_x)){ plot_x[i] = min_x + size_x; }
    if(plot_y[i] > (min_y + size_y)){ plot_y[i] = min_y + size_y; }

  }

  // 9. return the number of points making up the boundary
  return p;

}

int CreateRAPVBounds(float *** mini_RAPV, int size_x, int size_z, int min_x, int min_z, float * plot_x, float * plot_y, int obj, int obj_limit){

  int i, j, found, p = 0, obj_batch;

  // raster scan through the cube, and identify whether the source exists along a particular line of sight within the bounding volume
  // if it is, add entry to boundary list

  // 0. calculate obj_batch value
  obj_batch = floorf((float) obj / (float) obj_limit);

  // 1. identify left edge of source
  for(j = 0; j < size_z; j++){
    
    for(i = 0; i < size_x; i++){

      if(mini_RAPV[obj_batch][(obj - (obj_batch * obj_limit))][((j * size_x) + i)] != 0.0){

	plot_x[p] = (float) min_x + i - 1;
	plot_y[p] = (float) min_z + j;
	p++;
	break;

      }

    }

  }

  // 2. add top left point 
  plot_x[p] = plot_x[(p - 1)];
  plot_y[p] = plot_y[(p - 1)] + 1;
  p++;

  // 3. identify top right edge of source
  for(i = size_x - 1; i >= 0; i--){

    if(mini_RAPV[obj_batch][(obj - (obj_batch * obj_limit))][(((size_z - 1) * size_x) + i)] != 0.0){ 

      found = min_x + i + 1;
      break; 

    }

  }

  // 4. add top right point, then just identified top right edge of source
  plot_x[p] = (float) found;
  plot_y[p] = (float) min_z + size_z;
  p++;
  plot_x[p] = (float) found;
  plot_y[p] = (float) min_z + size_z - 1;
  p++;

  // 5. identify remaining right edge
  for(j = size_z - 2; j >= 0; j--){

    for(i = size_x - 1; i >= 0; i--){

      if(mini_RAPV[obj_batch][(obj - (obj_batch * obj_limit))][((j * size_x) + i)] != 0.0){

	plot_x[p] = (float) min_x + i + 1;
	plot_y[p] = (float) min_z + j;
	p++;
	break;
	
      }
     
    }

  }

  // 6. add bottom two points
  plot_x[p] = plot_x[(p - 1)];
  plot_y[p] = plot_y[(p - 1)] - 1;
  p++;
  plot_x[p] = plot_x[0];
  plot_y[p] = plot_y[0] - 1;
  p++;

  // 7. add original boundary point to complete 
  plot_x[p] = plot_x[0];
  plot_y[p] = plot_y[0];
  p++;

  // 8. check that all of the boundary points are within the bounding box limits
  for(i = 0; i < p; i++){

    if(plot_x[i] < (min_x - 1)){ plot_x[i] = min_x - 1; }
    if(plot_y[i] < (min_z - 1)){ plot_y[i] = min_z - 1; }
    if(plot_x[i] > (min_x + size_x)){ plot_x[i] = min_x + size_x; }
    if(plot_y[i] > (min_z + size_z)){ plot_y[i] = min_z + size_z; }

  }

  // 9. return the number of points making up the boundary
  return p;

}

int CreateDecPVBounds(float *** mini_DECPV, int size_y, int size_z, int min_y, int min_z, float * plot_x, float * plot_y, int obj, int obj_limit){

  int i, j, found, p = 0, obj_batch;

  // raster scan through the cube, and identify whether the source exists along a particular line of sight within the bounding volume
  // if it is, add entry to boundary list

  // calculate obj_batch value
  obj_batch = floorf((float) obj / (float) obj_limit);

  // 1. identify left edge of source
  for(j = 0; j < size_z; j++){
    
    for(i = 0; i < size_y; i++){

      if(mini_DECPV[obj_batch][(obj - (obj_batch * obj_limit))][((j * size_y) + i)] != 0.0){

	plot_x[p] = (float) min_y + i - 1;
	plot_y[p] = (float) min_z + j;
	p++;
	break;

      }

    }

  }

  // 2. add top left point 
  plot_x[p] = plot_x[(p - 1)];
  plot_y[p] = plot_y[(p - 1)] + 1;
  p++;

  // 3. identify top right edge of source
  for(i = size_y - 1; i >= 0; i--){

    if(mini_DECPV[obj_batch][(obj - (obj_batch * obj_limit))][(((size_z - 1) * size_y) + i)] != 0.0){ 

      found = min_y + i + 1;
      break; 

    }

  }

  // 4. add top right point, then just identified top right edge of source
  plot_x[p] = (float) found;
  plot_y[p] = (float) min_z + size_z;
  p++;
  plot_x[p] = (float) found;
  plot_y[p] = (float) min_z + size_z - 1;
  p++;

  // 5. identify remaining right edge
  for(j = size_z - 2; j >= 0; j--){

    for(i = size_y - 1; i >= 0; i--){

      if(mini_DECPV[obj_batch][(obj - (obj_batch * obj_limit))][((j * size_y) + i)] != 0.0){

	plot_x[p] = (float) min_y + i + 1;
	plot_y[p] = (float) min_z + j;
	p++;
	break;
	
      }
     
    }

  }

  // 6. add bottom two points
  plot_x[p] = plot_x[(p - 1)];
  plot_y[p] = plot_y[(p - 1)] - 1;
  p++;
  plot_x[p] = plot_x[0];
  plot_y[p] = plot_y[0] - 1;
  p++;

  // 7. add original boundary point to complete 
  plot_x[p] = plot_x[0];
  plot_y[p] = plot_y[0];
  p++;

  // 8. check that all of the boundary points are within the bounding box limits
  for(i = 0; i < p; i++){

    if(plot_x[i] < (min_y - 1)){ plot_x[i] = min_y - 1; }
    if(plot_y[i] < (min_z - 1)){ plot_y[i] = min_z - 1; }
    if(plot_x[i] > (min_y + size_y)){ plot_x[i] = min_y + size_y; }
    if(plot_y[i] > (min_z + size_z)){ plot_y[i] = min_z + size_z; }

  }
  // 9. return the number of points making up the boundary
  return p;

}

