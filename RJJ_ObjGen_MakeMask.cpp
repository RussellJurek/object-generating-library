#include<iostream>
#include<RJJ_ObjGen.h>
#include<string>
extern "C" {

#include<fitsio.h>
  
}

using namespace std;

// functions using floats

float AddIDsToChunk(float progress, int * temp_indices, int * flag_vals, vector<object_props *> & detections, int NOobj, int obj_limit, int chunk_x_start, int chunk_y_start, int chunk_z_start, int chunk_x_size, int chunk_y_size, int chunk_z_size, int * data_metric, int * xyz_order){

  int obj, obj_batch, x, y, g, f, temp_x[2], temp_y[2], temp_z[2];

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

  // add the objects within this datacube subregion to the output mask file
  for(obj = 0; obj < NOobj; obj++){
      
    // calculate obj_batch number for this object
    obj_batch = floorf(((float) obj / (float) obj_limit));
    
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
    
    // reconstruct the object from the sparse representation and update the flag_vals array
    for(x = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0); x <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(1); x++){
      
      for(y = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(2); y <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(3); y++){
	
	// x,y position of object lies inside boundaries of chunk, so flag the flag_vals array with the obj number for every object string
	if((x >= chunk_x_start) && (x < (chunk_x_start + chunk_x_size)) && (y >= chunk_y_start) && (y < (chunk_y_start + chunk_y_size))){
	  
	    for(g = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((((y - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + x - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0)))); g < detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((((y - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + x - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0) + 1))); g++){
	      
	      for(f = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings((2 * g)); ((f <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings(((2 * g) + 1))) && (f < (chunk_z_size + chunk_z_start))); f++){
		
		//flag_vals[(((f - chunk_z_start) * (data_x_size - chunk_x_start) * (data_y_size - chunk_y_start)) + ((y - chunk_y_start) * (data_x_size - chunk_x_start)) + x - chunk_x_start)] = temp_indices[obj];
		flag_vals[(((f - chunk_z_start) * data_metric[2]) + ((y - chunk_y_start) * data_metric[1]) + ((x - chunk_x_start) * data_metric[0]))] = temp_indices[obj];
		
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
  
  return progress;
  
}

float AddIDsToChunk(float progress, long int * temp_indices, long int * flag_vals, vector<object_props *> & detections, long int NOobj, int obj_limit, int chunk_x_start, int chunk_y_start, int chunk_z_start, int chunk_x_size, int chunk_y_size, int chunk_z_size, int * data_metric, int * xyz_order){

  long int obj, obj_batch;
  int x, y, g, f, temp_x[2], temp_y[2], temp_z[2];

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

  // add the objects within this datacube subregion to the output mask file
  for(obj = 0; obj < NOobj; obj++){
      
    // calculate obj_batch number for this object
    obj_batch = floorf(((double) obj / (double) obj_limit));
    
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
    
    // reconstruct the object from the sparse representation and update the flag_vals array
    for(x = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0); x <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(1); x++){
      
      for(y = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(2); y <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(3); y++){
	
	// x,y position of object lies inside boundaries of chunk, so flag the flag_vals array with the obj number for every object string
	if((x >= chunk_x_start) && (x < (chunk_x_start + chunk_x_size)) && (y >= chunk_y_start) && (y < (chunk_y_start + chunk_y_size))){
	  
	    for(g = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((((y - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + x - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0)))); g < detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((((y - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + x - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0) + 1))); g++){
	      
	      for(f = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings((2 * g)); ((f <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings(((2 * g) + 1))) && (f < (chunk_z_size + chunk_z_start))); f++){
		
		//flag_vals[(((f - chunk_z_start) * (data_x_size - chunk_x_start) * (data_y_size - chunk_y_start)) + ((y - chunk_y_start) * (data_x_size - chunk_x_start)) + x - chunk_x_start)] = temp_indices[obj];
		flag_vals[(((f - chunk_z_start) * data_metric[2]) + ((y - chunk_y_start) * data_metric[1]) + ((x - chunk_x_start) * data_metric[0]))] = temp_indices[obj];
		
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
  
  return progress;
  
}

void CreateFitsMask(std::string output_code, int NOx, int NOy, int NOf, vector<object_props *> & detections, int NOobj, int obj_limit, int * flag_vals, int chunk_x_size, int chunk_y_size, int temp_chunk_x_size, int temp_chunk_y_size, int * data_metric, int * xyz_order){

  fitsfile * mask_output;
  int * temp_indices; 
  long int fits_read_start[4], fits_read_finish[4], fits_read_inc[4] = {1,1,1,1};
  int obj, obj_batch, status, x, y, z, g, i, j;
  std::string mask_output_file;
  float progress;

  // create the mask_output file name
  mask_output_file = "!" + output_code + "_mask.fits";
  
  std::cout << "\nCreating output .fits file: " << mask_output_file << ", labelled with object indices." << std::endl;
  
  // create new file using the input mask file parameters as a template
  status = 0;
  fits_create_file(&mask_output,mask_output_file.c_str(),&status);
  fits_read_start[0] = NOx;
  fits_read_start[1] = NOy;
  fits_read_start[2] = NOf;
  fits_read_start[3] = 1;
  fits_create_img(mask_output,LONG_IMG,4,fits_read_start,&status);
  if(status != 0){
    
    std::cout << "WARNING!!! failed to create output mask .fits file. Not creating mask output file." << std::endl;
    if(status > 0){ fits_report_error(stderr,status); }

  } else {
    
    // create an array to store indices matching the catalogue IDs
    temp_indices = new int[NOobj];
    x = 1;
    for(obj = 0; obj < NOobj; obj++){
      
      // calculate obj_batch
      obj_batch = floorf(((float) obj / (float) obj_limit));
    
      // move on if this object has been re-initialised
      if(detections[obj_batch][(obj - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
           
      // write object properties to output file/terminal
      temp_indices[obj] = x;
      x++;
            
      //for(k = 0; k < NOobj; k++) 
    }
    
    progress = 0.0;
    std::cout << "0 | |:| | : | |:| | 100% complete" << std::endl;

    // set the constant fits_read_start and fits_read_finish array values
    fits_read_start[3] = fits_read_finish[3] = 1;
    fits_read_start[2] = 1;
    fits_read_finish[2] = NOf;
    
    // for each chunk, populate the flag_vals array with the catalogue IDs of the objects,
    // then write this array to the output file
    for(j = 0; j < temp_chunk_y_size; j++){

      for(i = 0; i < temp_chunk_x_size; i++){

	// set parameters for region to be processed	
	fits_read_start[0] = 1 + (i * chunk_x_size);
	fits_read_finish[0] = fits_read_start[0] + chunk_x_size - 1;
	if(fits_read_finish[0] > NOx){ fits_read_finish[0] = NOx; }
	
	fits_read_start[1] = 1 + (j * chunk_y_size);
	fits_read_finish[1] = fits_read_start[1] + chunk_y_size - 1;
	if(fits_read_finish[1] > NOy){ fits_read_finish[1] = NOy; }
	
	// initialise flag_vals array
	for(z = 0; z < NOf; z++){
	  for(y = 0; y < chunk_y_size; y++){
	    for(x = 0; x < chunk_x_size; x++){
	      flag_vals[((z * chunk_y_size * chunk_x_size) + (y * chunk_x_size) + x)] = 0;
	    }
	  }
	}
	
	// create metric for accessing this data chunk in x,y,z order
	CreateMetric(data_metric,xyz_order,(fits_read_finish[0] - fits_read_start[0] + 1),(fits_read_finish[1] - fits_read_start[1] + 1),(fits_read_finish[2] - fits_read_start[2] + 1));
	
	// add objects to flag_vals array for this chunk
	progress = AddIDsToChunk(progress,temp_indices,flag_vals,detections,NOobj,obj_limit,(fits_read_start[0] - 1),(fits_read_start[1] - 1),0,(fits_read_finish[0] - fits_read_start[0] + 1),(fits_read_finish[1] - fits_read_start[1] + 1),NOf,data_metric,xyz_order);

	// write updated flag_vals array to output file
	fits_write_subset(mask_output,TINT,fits_read_start,fits_read_finish,flag_vals,&status);	

	// for(i = 0; i < temp_chunk_x_size; i++)
      }

      // for(j = 0; j < temp_chunk_y_size; j++)
    }

    std::cout << "* done." << std::endl;
    
    // free up the temp_indices array
    delete [] temp_indices;

    // close the mask_output file
    fits_close_file(mask_output,&status);
    
    // else . . . if(status != 0)
  }
     
}

void CreateFitsMask(std::string output_code, int NOx, int NOy, int NOf, vector<object_props *> & detections, long int NOobj, int obj_limit, long int * flag_vals, int chunk_x_size, int chunk_y_size, int temp_chunk_x_size, int temp_chunk_y_size, int * data_metric, int * xyz_order){

  fitsfile * mask_output;
  long int * temp_indices; 
  long int fits_read_start[4], fits_read_finish[4], fits_read_inc[4] = {1,1,1,1};
  long int obj, obj_batch;
  int status, x, y, z, g, i, j;
  std::string mask_output_file;
  float progress;

  // create the mask_output file name
  mask_output_file = "!" + output_code + "_mask.fits";
  
  std::cout << "\nCreating output .fits file: " << mask_output_file << ", labelled with object indices." << std::endl;
  
  // create new file using the input mask file parameters as a template
  status = 0;
  fits_create_file(&mask_output,mask_output_file.c_str(),&status);
  fits_read_start[0] = NOx;
  fits_read_start[1] = NOy;
  fits_read_start[2] = NOf;
  fits_read_start[3] = 1;
  fits_create_img(mask_output,LONG_IMG,4,fits_read_start,&status);
  if(status != 0){
    
    std::cout << "WARNING!!! failed to create output mask .fits file. Not creating mask output file." << std::endl;
    if(status > 0){ fits_report_error(stderr,status); }

  } else {
    
    // create an array to store indices matching the catalogue IDs
    temp_indices = new long int[NOobj];
    x = 1;
    for(obj = 0; obj < NOobj; obj++){
      
      // calculate obj_batch
      obj_batch = floorf(((double) obj / (double) obj_limit));
    
      // move on if this object has been re-initialised
      if(detections[obj_batch][(obj - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
           
      // write object properties to output file/terminal
      temp_indices[obj] = x;
      x++;
            
      //for(k = 0; k < NOobj; k++) 
    }
    
    progress = 0.0;
    std::cout << "0 | |:| | : | |:| | 100% complete" << std::endl;

    // set the constant fits_read_start and fits_read_finish array values
    fits_read_start[3] = fits_read_finish[3] = 1;
    fits_read_start[2] = 1;
    fits_read_finish[2] = NOf;
    
    // for each chunk, populate the flag_vals array with the catalogue IDs of the objects,
    // then write this array to the output file
    for(j = 0; j < temp_chunk_y_size; j++){

      for(i = 0; i < temp_chunk_x_size; i++){

	// set parameters for region to be processed	
	fits_read_start[0] = 1 + (i * chunk_x_size);
	fits_read_finish[0] = fits_read_start[0] + chunk_x_size - 1;
	if(fits_read_finish[0] > NOx){ fits_read_finish[0] = NOx; }
	
	fits_read_start[1] = 1 + (j * chunk_y_size);
	fits_read_finish[1] = fits_read_start[1] + chunk_y_size - 1;
	if(fits_read_finish[1] > NOy){ fits_read_finish[1] = NOy; }
	
	// initialise flag_vals array
	for(z = 0; z < NOf; z++){
	  for(y = 0; y < chunk_y_size; y++){
	    for(x = 0; x < chunk_x_size; x++){
	      flag_vals[((z * chunk_y_size * chunk_x_size) + (y * chunk_x_size) + x)] = 0;
	    }
	  }
	}
	
	// create metric for accessing this data chunk in x,y,z order
	CreateMetric(data_metric,xyz_order,(fits_read_finish[0] - fits_read_start[0] + 1),(fits_read_finish[1] - fits_read_start[1] + 1),(fits_read_finish[2] - fits_read_start[2] + 1));
	
	// add objects to flag_vals array for this chunk
	progress = AddIDsToChunk(progress,temp_indices,flag_vals,detections,NOobj,obj_limit,(fits_read_start[0] - 1),(fits_read_start[1] - 1),0,(fits_read_finish[0] - fits_read_start[0] + 1),(fits_read_finish[1] - fits_read_start[1] + 1),NOf,data_metric,xyz_order);

	// write updated flag_vals array to output file
	fits_write_subset(mask_output,TINT,fits_read_start,fits_read_finish,flag_vals,&status);	

	// for(i = 0; i < temp_chunk_x_size; i++)
      }

      // for(j = 0; j < temp_chunk_y_size; j++)
    }

    std::cout << "* done." << std::endl;
    
    // free up the temp_indices array
    delete [] temp_indices;

    // close the mask_output file
    fits_close_file(mask_output,&status);
    
    // else . . . if(status != 0)
  }
     
}

// functions using doubles

float AddIDsToChunk(float progress, int * temp_indices, int * flag_vals, vector<object_props_dbl *> & detections, int NOobj, int obj_limit, int chunk_x_start, int chunk_y_start, int chunk_z_start, int chunk_x_size, int chunk_y_size, int chunk_z_size, int * data_metric, int * xyz_order){

  int obj, obj_batch, x, y, g, f, temp_x[2], temp_y[2], temp_z[2];

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

  // add the objects within this datacube subregion to the output mask file
  for(obj = 0; obj < NOobj; obj++){
      
    // calculate obj_batch number for this object
    obj_batch = floorf(((float) obj / (float) obj_limit));
    
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
    
    // reconstruct the object from the sparse representation and update the flag_vals array
    for(x = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0); x <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(1); x++){
      
      for(y = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(2); y <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(3); y++){
	
	// x,y position of object lies inside boundaries of chunk, so flag the flag_vals array with the obj number for every object string
	if((x >= chunk_x_start) && (x < (chunk_x_start + chunk_x_size)) && (y >= chunk_y_start) && (y < (chunk_y_start + chunk_y_size))){
	  
	    for(g = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((((y - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + x - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0)))); g < detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((((y - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + x - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0) + 1))); g++){
	      
	      for(f = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings((2 * g)); ((f <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings(((2 * g) + 1))) && (f < (chunk_z_size + chunk_z_start))); f++){
		
		//flag_vals[(((f - chunk_z_start) * (data_x_size - chunk_x_start) * (data_y_size - chunk_y_start)) + ((y - chunk_y_start) * (data_x_size - chunk_x_start)) + x - chunk_x_start)] = temp_indices[obj];
		flag_vals[(((f - chunk_z_start) * data_metric[2]) + ((y - chunk_y_start) * data_metric[1]) + ((x - chunk_x_start) * data_metric[0]))] = temp_indices[obj];
		
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
  
  return progress;
  
}

float AddIDsToChunk(float progress, long int * temp_indices, long int * flag_vals, vector<object_props_dbl *> & detections, long int NOobj, int obj_limit, int chunk_x_start, int chunk_y_start, int chunk_z_start, int chunk_x_size, int chunk_y_size, int chunk_z_size, int * data_metric, int * xyz_order){

  long int obj, obj_batch;
  int x, y, g, f, temp_x[2], temp_y[2], temp_z[2];

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

  // add the objects within this datacube subregion to the output mask file
  for(obj = 0; obj < NOobj; obj++){
      
    // calculate obj_batch number for this object
    obj_batch = floorf(((double) obj / (double) obj_limit));
    
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
    
    // reconstruct the object from the sparse representation and update the flag_vals array
    for(x = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0); x <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(1); x++){
      
      for(y = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(2); y <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(3); y++){
	
	// x,y position of object lies inside boundaries of chunk, so flag the flag_vals array with the obj number for every object string
	if((x >= chunk_x_start) && (x < (chunk_x_start + chunk_x_size)) && (y >= chunk_y_start) && (y < (chunk_y_start + chunk_y_size))){
	  
	    for(g = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((((y - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + x - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0)))); g < detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_grid(((((y - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(2)) * (detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0) + 1)) + x - detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_size(0) + 1))); g++){
	      
	      for(f = detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings((2 * g)); ((f <= detections[obj_batch][(obj - (obj_batch * obj_limit))].Get_srep_strings(((2 * g) + 1))) && (f < (chunk_z_size + chunk_z_start))); f++){
		
		//flag_vals[(((f - chunk_z_start) * (data_x_size - chunk_x_start) * (data_y_size - chunk_y_start)) + ((y - chunk_y_start) * (data_x_size - chunk_x_start)) + x - chunk_x_start)] = temp_indices[obj];
		flag_vals[(((f - chunk_z_start) * data_metric[2]) + ((y - chunk_y_start) * data_metric[1]) + ((x - chunk_x_start) * data_metric[0]))] = temp_indices[obj];
		
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
  
  return progress;
  
}

void CreateFitsMask(std::string output_code, int NOx, int NOy, int NOf, vector<object_props_dbl *> & detections, int NOobj, int obj_limit, int * flag_vals, int chunk_x_size, int chunk_y_size, int temp_chunk_x_size, int temp_chunk_y_size, int * data_metric, int * xyz_order){

  fitsfile * mask_output;
  int * temp_indices; 
  long int fits_read_start[4], fits_read_finish[4], fits_read_inc[4] = {1,1,1,1};
  int obj, obj_batch, status, x, y, z, g, i, j;
  std::string mask_output_file;
  float progress;

  // create the mask_output file name
  mask_output_file = "!" + output_code + "_mask.fits";
  
  std::cout << "\nCreating output .fits file: " << mask_output_file << ", labelled with object indices." << std::endl;
  
  // create new file using the input mask file parameters as a template
  status = 0;
  fits_create_file(&mask_output,mask_output_file.c_str(),&status);
  fits_read_start[0] = NOx;
  fits_read_start[1] = NOy;
  fits_read_start[2] = NOf;
  fits_read_start[3] = 1;
  fits_create_img(mask_output,LONG_IMG,4,fits_read_start,&status);
  if(status != 0){
    
    std::cout << "WARNING!!! failed to create output mask .fits file. Not creating mask output file." << std::endl;
    if(status > 0){ fits_report_error(stderr,status); }

  } else {
    
    // create an array to store indices matching the catalogue IDs
    temp_indices = new int[NOobj];
    x = 1;
    for(obj = 0; obj < NOobj; obj++){
      
      // calculate obj_batch
      obj_batch = floorf(((float) obj / (float) obj_limit));
    
      // move on if this object has been re-initialised
      if(detections[obj_batch][(obj - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
           
      // write object properties to output file/terminal
      temp_indices[obj] = x;
      x++;
            
      //for(k = 0; k < NOobj; k++) 
    }
    
    progress = 0.0;
    std::cout << "0 | |:| | : | |:| | 100% complete" << std::endl;

    // set the constant fits_read_start and fits_read_finish array values
    fits_read_start[3] = fits_read_finish[3] = 1;
    fits_read_start[2] = 1;
    fits_read_finish[2] = NOf;
    
    // for each chunk, populate the flag_vals array with the catalogue IDs of the objects,
    // then write this array to the output file
    for(j = 0; j < temp_chunk_y_size; j++){

      for(i = 0; i < temp_chunk_x_size; i++){

	// set parameters for region to be processed	
	fits_read_start[0] = 1 + (i * chunk_x_size);
	fits_read_finish[0] = fits_read_start[0] + chunk_x_size - 1;
	if(fits_read_finish[0] > NOx){ fits_read_finish[0] = NOx; }
	
	fits_read_start[1] = 1 + (j * chunk_y_size);
	fits_read_finish[1] = fits_read_start[1] + chunk_y_size - 1;
	if(fits_read_finish[1] > NOy){ fits_read_finish[1] = NOy; }
	
	// initialise flag_vals array
	for(z = 0; z < NOf; z++){
	  for(y = 0; y < chunk_y_size; y++){
	    for(x = 0; x < chunk_x_size; x++){
	      flag_vals[((z * chunk_y_size * chunk_x_size) + (y * chunk_x_size) + x)] = 0;
	    }
	  }
	}
	
	// create metric for accessing this data chunk in x,y,z order
	CreateMetric(data_metric,xyz_order,(fits_read_finish[0] - fits_read_start[0] + 1),(fits_read_finish[1] - fits_read_start[1] + 1),(fits_read_finish[2] - fits_read_start[2] + 1));
	
	// add objects to flag_vals array for this chunk
	progress = AddIDsToChunk(progress,temp_indices,flag_vals,detections,NOobj,obj_limit,(fits_read_start[0] - 1),(fits_read_start[1] - 1),0,(fits_read_finish[0] - fits_read_start[0] + 1),(fits_read_finish[1] - fits_read_start[1] + 1),NOf,data_metric,xyz_order);

	// write updated flag_vals array to output file
	fits_write_subset(mask_output,TINT,fits_read_start,fits_read_finish,flag_vals,&status);	

	// for(i = 0; i < temp_chunk_x_size; i++)
      }

      // for(j = 0; j < temp_chunk_y_size; j++)
    }

    std::cout << "* done." << std::endl;
    
    // free up the temp_indices array
    delete [] temp_indices;

    // close the mask_output file
    fits_close_file(mask_output,&status);
    
    // else . . . if(status != 0)
  }
     
}

void CreateFitsMask(std::string output_code, int NOx, int NOy, int NOf, vector<object_props_dbl *> & detections, long int NOobj, int obj_limit, long int * flag_vals, int chunk_x_size, int chunk_y_size, int temp_chunk_x_size, int temp_chunk_y_size, int * data_metric, int * xyz_order){

  fitsfile * mask_output;
  long int * temp_indices; 
  long int fits_read_start[4], fits_read_finish[4], fits_read_inc[4] = {1,1,1,1};
  long int obj, obj_batch;
  int status, x, y, z, g, i, j;
  std::string mask_output_file;
  float progress;

  // create the mask_output file name
  mask_output_file = "!" + output_code + "_mask.fits";
  
  std::cout << "\nCreating output .fits file: " << mask_output_file << ", labelled with object indices." << std::endl;
  
  // create new file using the input mask file parameters as a template
  status = 0;
  fits_create_file(&mask_output,mask_output_file.c_str(),&status);
  fits_read_start[0] = NOx;
  fits_read_start[1] = NOy;
  fits_read_start[2] = NOf;
  fits_read_start[3] = 1;
  fits_create_img(mask_output,LONG_IMG,4,fits_read_start,&status);
  if(status != 0){
    
    std::cout << "WARNING!!! failed to create output mask .fits file. Not creating mask output file." << std::endl;
    if(status > 0){ fits_report_error(stderr,status); }

  } else {
    
    // create an array to store indices matching the catalogue IDs
    temp_indices = new long int[NOobj];
    x = 1;
    for(obj = 0; obj < NOobj; obj++){
      
      // calculate obj_batch
      obj_batch = floorf(((double) obj / (double) obj_limit));
    
      // move on if this object has been re-initialised
      if(detections[obj_batch][(obj - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
           
      // write object properties to output file/terminal
      temp_indices[obj] = x;
      x++;
            
      //for(k = 0; k < NOobj; k++) 
    }
    
    progress = 0.0;
    std::cout << "0 | |:| | : | |:| | 100% complete" << std::endl;

    // set the constant fits_read_start and fits_read_finish array values
    fits_read_start[3] = fits_read_finish[3] = 1;
    fits_read_start[2] = 1;
    fits_read_finish[2] = NOf;
    
    // for each chunk, populate the flag_vals array with the catalogue IDs of the objects,
    // then write this array to the output file
    for(j = 0; j < temp_chunk_y_size; j++){

      for(i = 0; i < temp_chunk_x_size; i++){

	// set parameters for region to be processed	
	fits_read_start[0] = 1 + (i * chunk_x_size);
	fits_read_finish[0] = fits_read_start[0] + chunk_x_size - 1;
	if(fits_read_finish[0] > NOx){ fits_read_finish[0] = NOx; }
	
	fits_read_start[1] = 1 + (j * chunk_y_size);
	fits_read_finish[1] = fits_read_start[1] + chunk_y_size - 1;
	if(fits_read_finish[1] > NOy){ fits_read_finish[1] = NOy; }
	
	// initialise flag_vals array
	for(z = 0; z < NOf; z++){
	  for(y = 0; y < chunk_y_size; y++){
	    for(x = 0; x < chunk_x_size; x++){
	      flag_vals[((z * chunk_y_size * chunk_x_size) + (y * chunk_x_size) + x)] = 0;
	    }
	  }
	}
	
	// create metric for accessing this data chunk in x,y,z order
	CreateMetric(data_metric,xyz_order,(fits_read_finish[0] - fits_read_start[0] + 1),(fits_read_finish[1] - fits_read_start[1] + 1),(fits_read_finish[2] - fits_read_start[2] + 1));
	
	// add objects to flag_vals array for this chunk
	progress = AddIDsToChunk(progress,temp_indices,flag_vals,detections,NOobj,obj_limit,(fits_read_start[0] - 1),(fits_read_start[1] - 1),0,(fits_read_finish[0] - fits_read_start[0] + 1),(fits_read_finish[1] - fits_read_start[1] + 1),NOf,data_metric,xyz_order);

	// write updated flag_vals array to output file
	fits_write_subset(mask_output,TINT,fits_read_start,fits_read_finish,flag_vals,&status);	

	// for(i = 0; i < temp_chunk_x_size; i++)
      }

      // for(j = 0; j < temp_chunk_y_size; j++)
    }

    std::cout << "* done." << std::endl;
    
    // free up the temp_indices array
    delete [] temp_indices;

    // close the mask_output file
    fits_close_file(mask_output,&status);
    
    // else . . . if(status != 0)
  }
     
}

