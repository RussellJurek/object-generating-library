#include<iostream>
#include<RJJ_ObjGen.h>

using namespace std;

// functions using floats

float CreateMoment0Map(float * plot_array, int NOobj, vector<object_props *> & detections, int NOx, int NOy, int obj_limit){
  
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

	  plot_array[(((j + detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin()) * NOx) + i + detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin())]+=detections[obj_batch][(k - (obj_batch * obj_limit))].Get_mom0(((j * sx_finish) + i));
	 
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

float CreateMoment0Map(float * plot_array, long int NOobj, vector<object_props *> & detections, int NOx, int NOy, int obj_limit){
  
  long int k,obj_batch;
  int i,j,sx_finish,sy_finish;
  float max = -99.0;
  
  // initialise plot_array
  for(i = 0; i < NOx; i++){

    for(j = 0; j < NOy; j++){
      
      plot_array[((j*NOx) + i)] = 0.0;

    }
    
  }
  
  // search through objects and construct moment 0 map
  for(k = 0; k < NOobj; k++){
    
    obj_batch = floorf((double) k / (double) obj_limit);

    if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() > 1){

      sx_finish = detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin() + 1;
      sy_finish = detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin() + 1;

      for(i = 0; i < sx_finish; i++){
	
	for(j = 0; j < sy_finish; j++){

	  plot_array[(((j + detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin()) * NOx) + i + detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin())]+=detections[obj_batch][(k - (obj_batch * obj_limit))].Get_mom0(((j * sx_finish) + i));
	 
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

float CreateRAPVPlot(float * plot_array, int NOobj, vector<object_props *> & detections, int NOx, int NOz, int obj_limit){
  
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

	  plot_array[(((j + detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmin()) * NOx) + i + detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin())]+=detections[obj_batch][(k - (obj_batch * obj_limit))].Get_RAPV(((j * sx_finish) + i));
	 
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

float CreateRAPVPlot(float * plot_array, long int NOobj, vector<object_props *> & detections, int NOx, int NOz, int obj_limit){
  
  long int k,obj_batch;
  int i,j,sx_finish,sz_finish;
  float max = -99.0;
  
  // initialise plot_array
  for(i = 0; i < NOx; i++){

    for(j = 0; j < NOz; j++){
      
      plot_array[((j*NOx) + i)] = 0.0;

    }
    
  }
  
  // search through objects and construct moment 0 map
  for(k = 0; k < NOobj; k++){
    
    obj_batch = floorf((double) k / (double) obj_limit);

    if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() > 1){

      sx_finish = detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin() + 1;
      sz_finish = detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmin() + 1;

      for(i = 0; i < sx_finish; i++){
	
	for(j = 0; j < sz_finish; j++){

	  plot_array[(((j + detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmin()) * NOx) + i + detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin())]+=detections[obj_batch][(k - (obj_batch * obj_limit))].Get_RAPV(((j * sx_finish) + i));
	 
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

float CreateDecPVPlot(float * plot_array, int NOobj, vector<object_props *> & detections, int NOy, int NOz, int obj_limit){
  
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

	  plot_array[(((j + detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmin()) * NOy) + i + detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin())]+=detections[obj_batch][(k - (obj_batch * obj_limit))].Get_DECPV(((j * sy_finish) + i));
	 
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

float CreateDecPVPlot(float * plot_array, long int NOobj, vector<object_props *> & detections, int NOy, int NOz, int obj_limit){
  
  long int k,obj_batch;
  int i,j,sy_finish,sz_finish;
  float max = -99.0;
  
  // initialise plot_array
  for(i = 0; i < NOy; i++){

    for(j = 0; j < NOz; j++){
      
      plot_array[((j*NOy) + i)] = 0.0;

    }
    
  }
  
  // search through objects and construct moment 0 map
  for(k = 0; k < NOobj; k++){
    
    obj_batch = floorf((double) k / (double) obj_limit);

    if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() > 1){

      sy_finish = detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin() + 1;
      sz_finish = detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmin() + 1;

      for(i = 0; i < sy_finish; i++){
	
	for(j = 0; j < sz_finish; j++){

	  plot_array[(((j + detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmin()) * NOy) + i + detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin())]+=detections[obj_batch][(k - (obj_batch * obj_limit))].Get_DECPV(((j * sy_finish) + i));
	 
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

int CreateMoment0Bounds(vector<object_props *> & detections, int size_x, int size_y, int min_x, int min_y, float * plot_x, float * plot_y, int obj, int obj_limit){

  int i, j, found, p = 0, obj_batch;

  // raster scan through the cube, and identify whether the source exists along a particular line of sight within the bounding volume
  // if it is, add entry to boundary list

  // 0. calculate obj_batch value
  obj_batch = floorf((float) obj / (float) obj_limit);
    
  // 1. identify left edge of source
  for(j = 0; j < size_y; j++){

    for(i = 0; i < size_x; i++){

      if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i + 1))){

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

    if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid((((size_y - 1) * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid((((size_y - 1) * size_x) + i + 1))){ 

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

      if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i + 1))){

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

int CreateMoment0Bounds(vector<object_props *> & detections, int size_x, int size_y, int min_x, int min_y, float * plot_x, float * plot_y, long int obj, int obj_limit){

  long int obj_batch;
  int i, j, found, p = 0;

  // raster scan through the cube, and identify whether the source exists along a particular line of sight within the bounding volume
  // if it is, add entry to boundary list

  // 0. calculate obj_batch value
  obj_batch = floorf((double) obj / (double) obj_limit);
    
  // 1. identify left edge of source
  for(j = 0; j < size_y; j++){

    for(i = 0; i < size_x; i++){

      if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i + 1))){

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

    if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid((((size_y - 1) * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid((((size_y - 1) * size_x) + i + 1))){ 

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

      if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i + 1))){

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

int CreateRAPVBounds(vector<object_props *> & detections, int size_x, int size_y, int size_z, int min_x, int min_z, float * plot_x, float * plot_y, int obj, int obj_limit){

  int i, j, k, f, found, p = 0, obj_batch;

  // raster scan through the cube, and identify whether the source exists along a particular line of sight and channel within the bounding volume
  // if it is, add entry to boundary list

  // 0. calculate obj_batch value
  obj_batch = floorf((float) obj / (float) obj_limit);
    
  // 1. identify left edge of source
  for(k = 0; k < size_z; k++){

    found = -1;
    
    for(i = 0; i < size_x; i++){

      for(j = 0; j < size_y; j++){
	
	if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i + 1))){
	  
	  for(f = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)); f < detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j  * size_x) + i + 1)); f++){
	    
	    if(((k + min_z) >= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings((2 * f))) && ((k + min_z) <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings(((2 * f) + 1)))){
	      
	      plot_x[p] = (float) min_x + i - 1;
	      plot_y[p] = (float) min_z + k;
	      p++;
	      found = 1;
	      break;
	      
	    }
	    
	    // for(f . . . )
	  }
	  
	  // if(sparse_reps_grid[obj_batch][(obj - (obj_batch * obj_limit))][((j * size_x) + i)] != sparse_reps_grid[obj_batch][(obj - (obj_batch * obj_limit))][((j * size_x) + i + 1)])
	}
	
	if(found == 1){ break; }
	
	// for(j = 0; j < size_y; j++)
      }
      
      if(found == 1){ break; }
      
      // for(i = 0; i < size_x; i++)
    }
    
    // for(k = 0; k < size_z; k++)
  }

  // 2. add top left point 
  plot_x[p] = plot_x[(p - 1)];
  plot_y[p] = plot_y[(p - 1)] + 1;
  p++;

  // 3. identify top right edge of source
  found = -1;
  for(i = size_x - 1; i >= 0; i--){
    
    for(j = 0; j < size_y; j++){
	
      if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i + 1))){
	  
	for(f = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)); f < detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j  * size_x) + i + 1)); f++){
	  
	  if(((min_z + size_z - 1) >= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings((2 * f))) && ((min_z + size_z - 1) <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings(((2 * f) + 1)))){
	    
	    found = min_x + i + 1;
	    break;
	    
	  }
	  
	  // for(f . . . )
	  }
	
	// if(sparse_reps_grid[obj_batch][(obj - (obj_batch * obj_limit))][((j * size_x) + i)] != sparse_reps_grid[obj_batch][(obj - (obj_batch * obj_limit))][((j * size_x) + i + 1)])
      }
      
      if(found >= 0){ break; }
      
      // for(j = 0; j < size_y; j++)
    }
    
    if(found >= 0){ break; }
      
  }

  // 4. add top right point, then just identified top right edge of source
  plot_x[p] = (float) found;
  plot_y[p] = (float) min_z + size_z;
  p++;
  plot_x[p] = (float) found;
  plot_y[p] = (float) min_z + size_z - 1;
  p++;

  // 5. identify remaining right edge
  for(k = size_z - 2; k >= 0; k--){

    found = -1;
    
    for(i = size_x - 1; i >= 0; i--){

      for(j = 0; j < size_y; j++){
	
	if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i + 1))){
	  
	  for(f = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)); f < detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j  * size_x) + i + 1)); f++){
	    
	    if(((min_z + k) >= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings((2 * f))) && ((min_z + k) <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings(((2 * f) + 1)))){
	      
	      plot_x[p] = (float) min_x + i + 1;
	      plot_y[p] = (float) min_z + k;
	      p++;
	      found = 1;
	      break;
	      
	    }
	    
	    // for(f . . . )
	  }
	  
	  // if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid[((j * size_x) + i)] != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid[((j * size_x) + i + 1)])
	}
	
	if(found == 1){ break; }
	
	// for(j = 0; j < size_y; j++)
      }
      
      if(found == 1){ break; }
      
      // for(i = 0; i < size_x; i++)
    }
    
    // for(k = 0; k < size_z; k++)
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

int CreateRAPVBounds(vector<object_props *> & detections, int size_x, int size_y, int size_z, int min_x, int min_z, float * plot_x, float * plot_y, long int obj, int obj_limit){

  long int obj_batch;
  int i, j, k, f, found, p = 0;

  // raster scan through the cube, and identify whether the source exists along a particular line of sight and channel within the bounding volume
  // if it is, add entry to boundary list

  // 0. calculate obj_batch value
  obj_batch = floorf((double) obj / (double) obj_limit);
    
  // 1. identify left edge of source
  for(k = 0; k < size_z; k++){

    found = -1;
    
    for(i = 0; i < size_x; i++){

      for(j = 0; j < size_y; j++){
	
	if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i + 1))){
	  
	  for(f = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)); f < detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j  * size_x) + i + 1)); f++){
	    
	    if(((k + min_z) >= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings((2 * f))) && ((k + min_z) <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings(((2 * f) + 1)))){
	      
	      plot_x[p] = (float) min_x + i - 1;
	      plot_y[p] = (float) min_z + k;
	      p++;
	      found = 1;
	      break;
	      
	    }
	    
	    // for(f . . . )
	  }
	  
	  // if(sparse_reps_grid[obj_batch][(obj - (obj_batch * obj_limit))][((j * size_x) + i)] != sparse_reps_grid[obj_batch][(obj - (obj_batch * obj_limit))][((j * size_x) + i + 1)])
	}
	
	if(found == 1){ break; }
	
	// for(j = 0; j < size_y; j++)
      }
      
      if(found == 1){ break; }
      
      // for(i = 0; i < size_x; i++)
    }
    
    // for(k = 0; k < size_z; k++)
  }

  // 2. add top left point 
  plot_x[p] = plot_x[(p - 1)];
  plot_y[p] = plot_y[(p - 1)] + 1;
  p++;

  // 3. identify top right edge of source
  found = -1;
  for(i = size_x - 1; i >= 0; i--){
    
    for(j = 0; j < size_y; j++){
	
      if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i + 1))){
	  
	for(f = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)); f < detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j  * size_x) + i + 1)); f++){
	  
	  if(((min_z + size_z - 1) >= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings((2 * f))) && ((min_z + size_z - 1) <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings(((2 * f) + 1)))){
	    
	    found = min_x + i + 1;
	    break;
	    
	  }
	  
	  // for(f . . . )
	  }
	
	// if(sparse_reps_grid[obj_batch][(obj - (obj_batch * obj_limit))][((j * size_x) + i)] != sparse_reps_grid[obj_batch][(obj - (obj_batch * obj_limit))][((j * size_x) + i + 1)])
      }
      
      if(found >= 0){ break; }
      
      // for(j = 0; j < size_y; j++)
    }
    
    if(found >= 0){ break; }
      
  }

  // 4. add top right point, then just identified top right edge of source
  plot_x[p] = (float) found;
  plot_y[p] = (float) min_z + size_z;
  p++;
  plot_x[p] = (float) found;
  plot_y[p] = (float) min_z + size_z - 1;
  p++;

  // 5. identify remaining right edge
  for(k = size_z - 2; k >= 0; k--){

    found = -1;
    
    for(i = size_x - 1; i >= 0; i--){

      for(j = 0; j < size_y; j++){
	
	if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i + 1))){
	  
	  for(f = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)); f < detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j  * size_x) + i + 1)); f++){
	    
	    if(((min_z + k) >= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings((2 * f))) && ((min_z + k) <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings(((2 * f) + 1)))){
	      
	      plot_x[p] = (float) min_x + i + 1;
	      plot_y[p] = (float) min_z + k;
	      p++;
	      found = 1;
	      break;
	      
	    }
	    
	    // for(f . . . )
	  }
	  
	  // if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid[((j * size_x) + i)] != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid[((j * size_x) + i + 1)])
	}
	
	if(found == 1){ break; }
	
	// for(j = 0; j < size_y; j++)
      }
      
      if(found == 1){ break; }
      
      // for(i = 0; i < size_x; i++)
    }
    
    // for(k = 0; k < size_z; k++)
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

int CreateDecPVBounds(vector<object_props *> & detections, int size_x, int size_y, int size_z, int min_y, int min_z, float * plot_x, float * plot_y, int obj, int obj_limit){
  
  int i, j, k, f, found, p = 0, obj_batch;
  
  // raster scan through the cube, and identify whether the source exists along a particular line of sight and channel within the bounding volume
  // if it is, add entry to boundary list
  
  // 0. calculate obj_batch value
  obj_batch = floorf((float) obj / (float) obj_limit);
  
  // 1. identify left edge of source
  for(k = 0; k < size_z; k++){

    found = -1;
    
    for(j = 0; j < size_y; j++){

      for(i = 0; i < size_x; i++){
	
	if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i + 1))){
	  
	  for(f = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)); f < detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j  * size_x) + i + 1)); f++){
	    
	    if(((min_z + k) >= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings((2 * f))) && ((min_z + k) <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings(((2 * f) + 1)))){
	      
	      plot_x[p] = (float) min_y + j - 1;
	      plot_y[p] = (float) min_z + k;
	      p++;
	      found = 1;
	      break;
	      
	    }
	    
	    // for(f . . . )
	  }
	  
	  // if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid[((j * size_x) + i)] != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid[((j * size_x) + i + 1)])
	}
	
	if(found == 1){ break; }
	
	// for(j = 0; j < size_y; j++)
      }
      
      if(found == 1){ break; }
      
      // for(i = 0; i < size_x; i++)
    }
    
    // for(k = 0; k < size_z; k++)
  }

  // 2. add top left point 
  plot_x[p] = plot_x[(p - 1)];
  plot_y[p] = plot_y[(p - 1)] + 1;
  p++;

  // 3. identify top right edge of source
  found = -1;
  for(j = size_y - 1; j >= 0; j--){
    
    for(i = 0; i < size_x; i++){
	
      if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i + 1))){
	  
	for(f = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)); f < detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j  * size_x) + i + 1)); f++){
	  
	  if(((min_z + size_z - 1) >= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings((2 * f))) && ((min_z + size_z - 1) <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings(((2 * f) + 1)))){
	    
	    found = min_y + j + 1;
	    break;
	    
	  }
	  
	  // for(f . . . )
	  }
	
	// if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid[((j * size_x) + i)] != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid[((j * size_x) + i + 1)])
      }
      
      if(found >= 0){ break; }
      
      // for(j = 0; j < size_y; j++)
    }
    
    if(found >= 0){ break; }
      
  }

  // 4. add top right point, then just identified top right edge of source
  plot_x[p] = (float) found;
  plot_y[p] = (float) min_z + size_z;
  p++;
  plot_x[p] = (float) found;
  plot_y[p] = (float) min_z + size_z - 1;
  p++;

  // 5. identify remaining right edge
  for(k = size_z - 2; k >= 0; k--){

    found = -1;
    
    for(j = size_y - 1; j >= 0; j--){

      for(i = 0; i < size_x; i++){
	
	if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i + 1))){
	  
	  for(f = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)); f < detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j  * size_x) + i + 1)); f++){
	    
	    if(((min_z + k) >= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings((2 * f))) && ((min_z + k) <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings(((2 * f) + 1)))){
	      
	      plot_x[p] = (float) min_y + j + 1;
	      plot_y[p] = (float) min_z + k;
	      p++;
	      found = 1;
	      break;
	      
	    }
	    
	    // for(f . . . )
	  }
	  
	  // if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid[((j * size_x) + i)] != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid[((j * size_x) + i + 1)])
	}
	
	if(found == 1){ break; }
	
	// for(j = 0; j < size_y; j++)
      }
      
      if(found == 1){ break; }
      
      // for(i = 0; i < size_x; i++)
    }
    
    // for(k = 0; k < size_z; k++)
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

int CreateDecPVBounds(vector<object_props *> & detections, int size_x, int size_y, int size_z, int min_y, int min_z, float * plot_x, float * plot_y, long int obj, int obj_limit){
  
  long int obj_batch;
  int i, j, k, f, found, p = 0;
  
  // raster scan through the cube, and identify whether the source exists along a particular line of sight and channel within the bounding volume
  // if it is, add entry to boundary list
  
  // 0. calculate obj_batch value
  obj_batch = floorf((double) obj / (double) obj_limit);
  
  // 1. identify left edge of source
  for(k = 0; k < size_z; k++){

    found = -1;
    
    for(j = 0; j < size_y; j++){

      for(i = 0; i < size_x; i++){
	
	if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i + 1))){
	  
	  for(f = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)); f < detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j  * size_x) + i + 1)); f++){
	    
	    if(((min_z + k) >= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings((2 * f))) && ((min_z + k) <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings(((2 * f) + 1)))){
	      
	      plot_x[p] = (float) min_y + j - 1;
	      plot_y[p] = (float) min_z + k;
	      p++;
	      found = 1;
	      break;
	      
	    }
	    
	    // for(f . . . )
	  }
	  
	  // if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid[((j * size_x) + i)] != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid[((j * size_x) + i + 1)])
	}
	
	if(found == 1){ break; }
	
	// for(j = 0; j < size_y; j++)
      }
      
      if(found == 1){ break; }
      
      // for(i = 0; i < size_x; i++)
    }
    
    // for(k = 0; k < size_z; k++)
  }

  // 2. add top left point 
  plot_x[p] = plot_x[(p - 1)];
  plot_y[p] = plot_y[(p - 1)] + 1;
  p++;

  // 3. identify top right edge of source
  found = -1;
  for(j = size_y - 1; j >= 0; j--){
    
    for(i = 0; i < size_x; i++){
	
      if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i + 1))){
	  
	for(f = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)); f < detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j  * size_x) + i + 1)); f++){
	  
	  if(((min_z + size_z - 1) >= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings((2 * f))) && ((min_z + size_z - 1) <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings(((2 * f) + 1)))){
	    
	    found = min_y + j + 1;
	    break;
	    
	  }
	  
	  // for(f . . . )
	  }
	
	// if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid[((j * size_x) + i)] != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid[((j * size_x) + i + 1)])
      }
      
      if(found >= 0){ break; }
      
      // for(j = 0; j < size_y; j++)
    }
    
    if(found >= 0){ break; }
      
  }

  // 4. add top right point, then just identified top right edge of source
  plot_x[p] = (float) found;
  plot_y[p] = (float) min_z + size_z;
  p++;
  plot_x[p] = (float) found;
  plot_y[p] = (float) min_z + size_z - 1;
  p++;

  // 5. identify remaining right edge
  for(k = size_z - 2; k >= 0; k--){

    found = -1;
    
    for(j = size_y - 1; j >= 0; j--){

      for(i = 0; i < size_x; i++){
	
	if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i + 1))){
	  
	  for(f = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)); f < detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j  * size_x) + i + 1)); f++){
	    
	    if(((min_z + k) >= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings((2 * f))) && ((min_z + k) <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings(((2 * f) + 1)))){
	      
	      plot_x[p] = (float) min_y + j + 1;
	      plot_y[p] = (float) min_z + k;
	      p++;
	      found = 1;
	      break;
	      
	    }
	    
	    // for(f . . . )
	  }
	  
	  // if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid[((j * size_x) + i)] != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid[((j * size_x) + i + 1)])
	}
	
	if(found == 1){ break; }
	
	// for(j = 0; j < size_y; j++)
      }
      
      if(found == 1){ break; }
      
      // for(i = 0; i < size_x; i++)
    }
    
    // for(k = 0; k < size_z; k++)
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

// functions using doubles

float CreateMoment0Map(float * plot_array, int NOobj, vector<object_props_dbl *> & detections, int NOx, int NOy, int obj_limit){
  
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

	  plot_array[(((j + detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin()) * NOx) + i + detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin())]+=detections[obj_batch][(k - (obj_batch * obj_limit))].Get_mom0(((j * sx_finish) + i));
	 
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

float CreateMoment0Map(float * plot_array, long int NOobj, vector<object_props_dbl *> & detections, int NOx, int NOy, int obj_limit){
  
  long int k,obj_batch;
  int i,j,sx_finish,sy_finish;
  float max = -99.0;
  
  // initialise plot_array
  for(i = 0; i < NOx; i++){

    for(j = 0; j < NOy; j++){
      
      plot_array[((j*NOx) + i)] = 0.0;

    }
    
  }
  
  // search through objects and construct moment 0 map
  for(k = 0; k < NOobj; k++){
    
    obj_batch = floorf((double) k / (double) obj_limit);

    if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() > 1){

      sx_finish = detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin() + 1;
      sy_finish = detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin() + 1;

      for(i = 0; i < sx_finish; i++){
	
	for(j = 0; j < sy_finish; j++){

	  plot_array[(((j + detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin()) * NOx) + i + detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin())]+=detections[obj_batch][(k - (obj_batch * obj_limit))].Get_mom0(((j * sx_finish) + i));
	 
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

float CreateRAPVPlot(float * plot_array, int NOobj, vector<object_props_dbl *> & detections, int NOx, int NOz, int obj_limit){
  
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

	  plot_array[(((j + detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmin()) * NOx) + i + detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin())]+=detections[obj_batch][(k - (obj_batch * obj_limit))].Get_RAPV(((j * sx_finish) + i));
	 
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

float CreateRAPVPlot(float * plot_array, long int NOobj, vector<object_props_dbl *> & detections, int NOx, int NOz, int obj_limit){
  
  long int k,obj_batch;
  int i,j,sx_finish,sz_finish;
  float max = -99.0;
  
  // initialise plot_array
  for(i = 0; i < NOx; i++){

    for(j = 0; j < NOz; j++){
      
      plot_array[((j*NOx) + i)] = 0.0;

    }
    
  }
  
  // search through objects and construct moment 0 map
  for(k = 0; k < NOobj; k++){
    
    obj_batch = floorf((double) k / (double) obj_limit);

    if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() > 1){

      sx_finish = detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin() + 1;
      sz_finish = detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmin() + 1;

      for(i = 0; i < sx_finish; i++){
	
	for(j = 0; j < sz_finish; j++){

	  plot_array[(((j + detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmin()) * NOx) + i + detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin())]+=detections[obj_batch][(k - (obj_batch * obj_limit))].Get_RAPV(((j * sx_finish) + i));
	 
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

float CreateDecPVPlot(float * plot_array, int NOobj, vector<object_props_dbl *> & detections, int NOy, int NOz, int obj_limit){
  
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

	  plot_array[(((j + detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmin()) * NOy) + i + detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin())]+=detections[obj_batch][(k - (obj_batch * obj_limit))].Get_DECPV(((j * sy_finish) + i));
	 
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

float CreateDecPVPlot(float * plot_array, long int NOobj, vector<object_props_dbl *> & detections, int NOy, int NOz, int obj_limit){
  
  long int k,obj_batch;
  int i,j,sy_finish,sz_finish;
  float max = -99.0;
  
  // initialise plot_array
  for(i = 0; i < NOy; i++){

    for(j = 0; j < NOz; j++){
      
      plot_array[((j*NOy) + i)] = 0.0;

    }
    
  }
  
  // search through objects and construct moment 0 map
  for(k = 0; k < NOobj; k++){
    
    obj_batch = floorf((double) k / (double) obj_limit);

    if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() > 1){

      sy_finish = detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin() + 1;
      sz_finish = detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmin() + 1;

      for(i = 0; i < sy_finish; i++){
	
	for(j = 0; j < sz_finish; j++){

	  plot_array[(((j + detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmin()) * NOy) + i + detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin())]+=detections[obj_batch][(k - (obj_batch * obj_limit))].Get_DECPV(((j * sy_finish) + i));
	 
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

int CreateMoment0Bounds(vector<object_props_dbl *> & detections, int size_x, int size_y, int min_x, int min_y, float * plot_x, float * plot_y, int obj, int obj_limit){

  int i, j, found, p = 0, obj_batch;

  // raster scan through the cube, and identify whether the source exists along a particular line of sight within the bounding volume
  // if it is, add entry to boundary list

  // 0. calculate obj_batch value
  obj_batch = floorf((float) obj / (float) obj_limit);
    
  // 1. identify left edge of source
  for(j = 0; j < size_y; j++){

    for(i = 0; i < size_x; i++){

      if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i + 1))){

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

    if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid((((size_y - 1) * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid((((size_y - 1) * size_x) + i + 1))){ 

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

      if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i + 1))){

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

int CreateMoment0Bounds(vector<object_props_dbl *> & detections, int size_x, int size_y, int min_x, int min_y, float * plot_x, float * plot_y, long int obj, int obj_limit){

  long int obj_batch;
  int i, j, found, p = 0;

  // raster scan through the cube, and identify whether the source exists along a particular line of sight within the bounding volume
  // if it is, add entry to boundary list

  // 0. calculate obj_batch value
  obj_batch = floorf((double) obj / (double) obj_limit);
    
  // 1. identify left edge of source
  for(j = 0; j < size_y; j++){

    for(i = 0; i < size_x; i++){

      if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i + 1))){

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

    if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid((((size_y - 1) * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid((((size_y - 1) * size_x) + i + 1))){ 

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

      if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i + 1))){

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

int CreateRAPVBounds(vector<object_props_dbl *> & detections, int size_x, int size_y, int size_z, int min_x, int min_z, float * plot_x, float * plot_y, int obj, int obj_limit){

  int i, j, k, f, found, p = 0, obj_batch;

  // raster scan through the cube, and identify whether the source exists along a particular line of sight and channel within the bounding volume
  // if it is, add entry to boundary list

  // 0. calculate obj_batch value
  obj_batch = floorf((float) obj / (float) obj_limit);
    
  // 1. identify left edge of source
  for(k = 0; k < size_z; k++){

    found = -1;
    
    for(i = 0; i < size_x; i++){

      for(j = 0; j < size_y; j++){
	
	if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i + 1))){
	  
	  for(f = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)); f < detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j  * size_x) + i + 1)); f++){
	    
	    if(((k + min_z) >= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings((2 * f))) && ((k + min_z) <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings(((2 * f) + 1)))){
	      
	      plot_x[p] = (float) min_x + i - 1;
	      plot_y[p] = (float) min_z + k;
	      p++;
	      found = 1;
	      break;
	      
	    }
	    
	    // for(f . . . )
	  }
	  
	  // if(sparse_reps_grid[obj_batch][(obj - (obj_batch * obj_limit))][((j * size_x) + i)] != sparse_reps_grid[obj_batch][(obj - (obj_batch * obj_limit))][((j * size_x) + i + 1)])
	}
	
	if(found == 1){ break; }
	
	// for(j = 0; j < size_y; j++)
      }
      
      if(found == 1){ break; }
      
      // for(i = 0; i < size_x; i++)
    }
    
    // for(k = 0; k < size_z; k++)
  }

  // 2. add top left point 
  plot_x[p] = plot_x[(p - 1)];
  plot_y[p] = plot_y[(p - 1)] + 1;
  p++;

  // 3. identify top right edge of source
  found = -1;
  for(i = size_x - 1; i >= 0; i--){
    
    for(j = 0; j < size_y; j++){
	
      if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i + 1))){
	  
	for(f = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)); f < detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j  * size_x) + i + 1)); f++){
	  
	  if(((min_z + size_z - 1) >= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings((2 * f))) && ((min_z + size_z - 1) <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings(((2 * f) + 1)))){
	    
	    found = min_x + i + 1;
	    break;
	    
	  }
	  
	  // for(f . . . )
	  }
	
	// if(sparse_reps_grid[obj_batch][(obj - (obj_batch * obj_limit))][((j * size_x) + i)] != sparse_reps_grid[obj_batch][(obj - (obj_batch * obj_limit))][((j * size_x) + i + 1)])
      }
      
      if(found >= 0){ break; }
      
      // for(j = 0; j < size_y; j++)
    }
    
    if(found >= 0){ break; }
      
  }

  // 4. add top right point, then just identified top right edge of source
  plot_x[p] = (float) found;
  plot_y[p] = (float) min_z + size_z;
  p++;
  plot_x[p] = (float) found;
  plot_y[p] = (float) min_z + size_z - 1;
  p++;

  // 5. identify remaining right edge
  for(k = size_z - 2; k >= 0; k--){

    found = -1;
    
    for(i = size_x - 1; i >= 0; i--){

      for(j = 0; j < size_y; j++){
	
	if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i + 1))){
	  
	  for(f = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)); f < detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j  * size_x) + i + 1)); f++){
	    
	    if(((min_z + k) >= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings((2 * f))) && ((min_z + k) <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings(((2 * f) + 1)))){
	      
	      plot_x[p] = (float) min_x + i + 1;
	      plot_y[p] = (float) min_z + k;
	      p++;
	      found = 1;
	      break;
	      
	    }
	    
	    // for(f . . . )
	  }
	  
	  // if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid[((j * size_x) + i)] != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid[((j * size_x) + i + 1)])
	}
	
	if(found == 1){ break; }
	
	// for(j = 0; j < size_y; j++)
      }
      
      if(found == 1){ break; }
      
      // for(i = 0; i < size_x; i++)
    }
    
    // for(k = 0; k < size_z; k++)
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

int CreateRAPVBounds(vector<object_props_dbl *> & detections, int size_x, int size_y, int size_z, int min_x, int min_z, float * plot_x, float * plot_y, long int obj, int obj_limit){

  long int obj_batch;
  int i, j, k, f, found, p = 0;

  // raster scan through the cube, and identify whether the source exists along a particular line of sight and channel within the bounding volume
  // if it is, add entry to boundary list

  // 0. calculate obj_batch value
  obj_batch = floorf((double) obj / (double) obj_limit);
    
  // 1. identify left edge of source
  for(k = 0; k < size_z; k++){

    found = -1;
    
    for(i = 0; i < size_x; i++){

      for(j = 0; j < size_y; j++){
	
	if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i + 1))){
	  
	  for(f = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)); f < detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j  * size_x) + i + 1)); f++){
	    
	    if(((k + min_z) >= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings((2 * f))) && ((k + min_z) <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings(((2 * f) + 1)))){
	      
	      plot_x[p] = (float) min_x + i - 1;
	      plot_y[p] = (float) min_z + k;
	      p++;
	      found = 1;
	      break;
	      
	    }
	    
	    // for(f . . . )
	  }
	  
	  // if(sparse_reps_grid[obj_batch][(obj - (obj_batch * obj_limit))][((j * size_x) + i)] != sparse_reps_grid[obj_batch][(obj - (obj_batch * obj_limit))][((j * size_x) + i + 1)])
	}
	
	if(found == 1){ break; }
	
	// for(j = 0; j < size_y; j++)
      }
      
      if(found == 1){ break; }
      
      // for(i = 0; i < size_x; i++)
    }
    
    // for(k = 0; k < size_z; k++)
  }

  // 2. add top left point 
  plot_x[p] = plot_x[(p - 1)];
  plot_y[p] = plot_y[(p - 1)] + 1;
  p++;

  // 3. identify top right edge of source
  found = -1;
  for(i = size_x - 1; i >= 0; i--){
    
    for(j = 0; j < size_y; j++){
	
      if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i + 1))){
	  
	for(f = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)); f < detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j  * size_x) + i + 1)); f++){
	  
	  if(((min_z + size_z - 1) >= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings((2 * f))) && ((min_z + size_z - 1) <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings(((2 * f) + 1)))){
	    
	    found = min_x + i + 1;
	    break;
	    
	  }
	  
	  // for(f . . . )
	  }
	
	// if(sparse_reps_grid[obj_batch][(obj - (obj_batch * obj_limit))][((j * size_x) + i)] != sparse_reps_grid[obj_batch][(obj - (obj_batch * obj_limit))][((j * size_x) + i + 1)])
      }
      
      if(found >= 0){ break; }
      
      // for(j = 0; j < size_y; j++)
    }
    
    if(found >= 0){ break; }
      
  }

  // 4. add top right point, then just identified top right edge of source
  plot_x[p] = (float) found;
  plot_y[p] = (float) min_z + size_z;
  p++;
  plot_x[p] = (float) found;
  plot_y[p] = (float) min_z + size_z - 1;
  p++;

  // 5. identify remaining right edge
  for(k = size_z - 2; k >= 0; k--){

    found = -1;
    
    for(i = size_x - 1; i >= 0; i--){

      for(j = 0; j < size_y; j++){
	
	if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i + 1))){
	  
	  for(f = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)); f < detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j  * size_x) + i + 1)); f++){
	    
	    if(((min_z + k) >= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings((2 * f))) && ((min_z + k) <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings(((2 * f) + 1)))){
	      
	      plot_x[p] = (float) min_x + i + 1;
	      plot_y[p] = (float) min_z + k;
	      p++;
	      found = 1;
	      break;
	      
	    }
	    
	    // for(f . . . )
	  }
	  
	  // if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid[((j * size_x) + i)] != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid[((j * size_x) + i + 1)])
	}
	
	if(found == 1){ break; }
	
	// for(j = 0; j < size_y; j++)
      }
      
      if(found == 1){ break; }
      
      // for(i = 0; i < size_x; i++)
    }
    
    // for(k = 0; k < size_z; k++)
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

int CreateDecPVBounds(vector<object_props_dbl *> & detections, int size_x, int size_y, int size_z, int min_y, int min_z, float * plot_x, float * plot_y, int obj, int obj_limit){
  
  int i, j, k, f, found, p = 0, obj_batch;
  
  // raster scan through the cube, and identify whether the source exists along a particular line of sight and channel within the bounding volume
  // if it is, add entry to boundary list
  
  // 0. calculate obj_batch value
  obj_batch = floorf((float) obj / (float) obj_limit);
  
  // 1. identify left edge of source
  for(k = 0; k < size_z; k++){

    found = -1;
    
    for(j = 0; j < size_y; j++){

      for(i = 0; i < size_x; i++){
	
	if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i + 1))){
	  
	  for(f = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)); f < detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j  * size_x) + i + 1)); f++){
	    
	    if(((min_z + k) >= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings((2 * f))) && ((min_z + k) <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings(((2 * f) + 1)))){
	      
	      plot_x[p] = (float) min_y + j - 1;
	      plot_y[p] = (float) min_z + k;
	      p++;
	      found = 1;
	      break;
	      
	    }
	    
	    // for(f . . . )
	  }
	  
	  // if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid[((j * size_x) + i)] != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid[((j * size_x) + i + 1)])
	}
	
	if(found == 1){ break; }
	
	// for(j = 0; j < size_y; j++)
      }
      
      if(found == 1){ break; }
      
      // for(i = 0; i < size_x; i++)
    }
    
    // for(k = 0; k < size_z; k++)
  }

  // 2. add top left point 
  plot_x[p] = plot_x[(p - 1)];
  plot_y[p] = plot_y[(p - 1)] + 1;
  p++;

  // 3. identify top right edge of source
  found = -1;
  for(j = size_y - 1; j >= 0; j--){
    
    for(i = 0; i < size_x; i++){
	
      if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i + 1))){
	  
	for(f = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)); f < detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j  * size_x) + i + 1)); f++){
	  
	  if(((min_z + size_z - 1) >= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings((2 * f))) && ((min_z + size_z - 1) <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings(((2 * f) + 1)))){
	    
	    found = min_y + j + 1;
	    break;
	    
	  }
	  
	  // for(f . . . )
	  }
	
	// if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid[((j * size_x) + i)] != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid[((j * size_x) + i + 1)])
      }
      
      if(found >= 0){ break; }
      
      // for(j = 0; j < size_y; j++)
    }
    
    if(found >= 0){ break; }
      
  }

  // 4. add top right point, then just identified top right edge of source
  plot_x[p] = (float) found;
  plot_y[p] = (float) min_z + size_z;
  p++;
  plot_x[p] = (float) found;
  plot_y[p] = (float) min_z + size_z - 1;
  p++;

  // 5. identify remaining right edge
  for(k = size_z - 2; k >= 0; k--){

    found = -1;
    
    for(j = size_y - 1; j >= 0; j--){

      for(i = 0; i < size_x; i++){
	
	if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i + 1))){
	  
	  for(f = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)); f < detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j  * size_x) + i + 1)); f++){
	    
	    if(((min_z + k) >= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings((2 * f))) && ((min_z + k) <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings(((2 * f) + 1)))){
	      
	      plot_x[p] = (float) min_y + j + 1;
	      plot_y[p] = (float) min_z + k;
	      p++;
	      found = 1;
	      break;
	      
	    }
	    
	    // for(f . . . )
	  }
	  
	  // if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid[((j * size_x) + i)] != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid[((j * size_x) + i + 1)])
	}
	
	if(found == 1){ break; }
	
	// for(j = 0; j < size_y; j++)
      }
      
      if(found == 1){ break; }
      
      // for(i = 0; i < size_x; i++)
    }
    
    // for(k = 0; k < size_z; k++)
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

int CreateDecPVBounds(vector<object_props_dbl *> & detections, int size_x, int size_y, int size_z, int min_y, int min_z, float * plot_x, float * plot_y, long int obj, int obj_limit){
  
  long int obj_batch;
  int i, j, k, f, found, p = 0;
  
  // raster scan through the cube, and identify whether the source exists along a particular line of sight and channel within the bounding volume
  // if it is, add entry to boundary list
  
  // 0. calculate obj_batch value
  obj_batch = floorf((double) obj / (double) obj_limit);
  
  // 1. identify left edge of source
  for(k = 0; k < size_z; k++){

    found = -1;
    
    for(j = 0; j < size_y; j++){

      for(i = 0; i < size_x; i++){
	
	if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i + 1))){
	  
	  for(f = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)); f < detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j  * size_x) + i + 1)); f++){
	    
	    if(((min_z + k) >= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings((2 * f))) && ((min_z + k) <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings(((2 * f) + 1)))){
	      
	      plot_x[p] = (float) min_y + j - 1;
	      plot_y[p] = (float) min_z + k;
	      p++;
	      found = 1;
	      break;
	      
	    }
	    
	    // for(f . . . )
	  }
	  
	  // if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid[((j * size_x) + i)] != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid[((j * size_x) + i + 1)])
	}
	
	if(found == 1){ break; }
	
	// for(j = 0; j < size_y; j++)
      }
      
      if(found == 1){ break; }
      
      // for(i = 0; i < size_x; i++)
    }
    
    // for(k = 0; k < size_z; k++)
  }

  // 2. add top left point 
  plot_x[p] = plot_x[(p - 1)];
  plot_y[p] = plot_y[(p - 1)] + 1;
  p++;

  // 3. identify top right edge of source
  found = -1;
  for(j = size_y - 1; j >= 0; j--){
    
    for(i = 0; i < size_x; i++){
	
      if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i + 1))){
	  
	for(f = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)); f < detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j  * size_x) + i + 1)); f++){
	  
	  if(((min_z + size_z - 1) >= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings((2 * f))) && ((min_z + size_z - 1) <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings(((2 * f) + 1)))){
	    
	    found = min_y + j + 1;
	    break;
	    
	  }
	  
	  // for(f . . . )
	  }
	
	// if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid[((j * size_x) + i)] != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid[((j * size_x) + i + 1)])
      }
      
      if(found >= 0){ break; }
      
      // for(j = 0; j < size_y; j++)
    }
    
    if(found >= 0){ break; }
      
  }

  // 4. add top right point, then just identified top right edge of source
  plot_x[p] = (float) found;
  plot_y[p] = (float) min_z + size_z;
  p++;
  plot_x[p] = (float) found;
  plot_y[p] = (float) min_z + size_z - 1;
  p++;

  // 5. identify remaining right edge
  for(k = size_z - 2; k >= 0; k--){

    found = -1;
    
    for(j = size_y - 1; j >= 0; j--){

      for(i = 0; i < size_x; i++){
	
	if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)) != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i + 1))){
	  
	  for(f = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j * size_x) + i)); f < detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((j  * size_x) + i + 1)); f++){
	    
	    if(((min_z + k) >= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings((2 * f))) && ((min_z + k) <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings(((2 * f) + 1)))){
	      
	      plot_x[p] = (float) min_y + j + 1;
	      plot_y[p] = (float) min_z + k;
	      p++;
	      found = 1;
	      break;
	      
	    }
	    
	    // for(f . . . )
	  }
	  
	  // if(detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid[((j * size_x) + i)] != detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid[((j * size_x) + i + 1)])
	}
	
	if(found == 1){ break; }
	
	// for(j = 0; j < size_y; j++)
      }
      
      if(found == 1){ break; }
      
      // for(i = 0; i < size_x; i++)
    }
    
    // for(k = 0; k < size_z; k++)
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

