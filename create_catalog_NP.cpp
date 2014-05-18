#include<fstream>
#include<iostream>
#include<string>
#include<sstream>
#include<iomanip>
#include<cmath>
#include<ctime>
#include<ctype.h>
#include<RJJ_ObjGen.h>
extern "C" {

#include<fitsio.h>

}

using namespace std;

// compile line
// c++ /Users/Jurek83/CatalogCode/latest_prototype/trunk/create_catalog_NP.cpp -o /Users/Jurek83/CatalogCode/latest_prototype/trunk/create_catalog_NP -O2 -lcfitsio -I/usr/local/pgplot -L/usr/local/pgplot -lcpgplot -lpgplot -I/Users/Jurek83/CatalogCode/latest_prototype/trunk -L/Users/Jurek83/CatalogCode/latest_prototype/trunk -lrjj_objgen

int main(int argc, char* argv[]){

  fitsfile * fits_input, * fits_input_source, * mask_input;
  stringstream dummy2;
  string dummy1, output_code, maskfile, sfile, snfile;
  char dummy3[100];
  int NOf,f,x,y,status,nkeys,i,g,NOg,j,k,m,length,NOx,NOy,xmid,ymid,xfinish,yfinish;
  float * data_vals, ratio, threshold, * signal_vals, fill_factor = 0.7;
  int * flag_vals, * data_metric, * xyz_order;
  int min_z_size = 1, max_z_size = 300, min_x_size = 3, min_y_size = 3, min_v_size, min_LoS_count;
  long fits_read_start[4], fits_read_finish[4], fits_read_inc[4] = {1,1,1,1};
  int chunk_size = 1024 * 1024 * 1024 / 4, NO_chunks, chunk_x_overlap = 10, chunk_y_overlap = 10, chunk_z_overlap = 10;
  int chunk_x_size, chunk_y_size, chunk_z_size, chunk_x_start, chunk_y_start, chunk_z_start;
  int best_NO_chunks, best_chunk_x_size, temp_chunk_x_size, best_chunk_y_size, temp_chunk_y_size;
  int SFmethod = 1, NOobj, obj, NOtrue, NOobj_positive;
  float label_y_shift = 0.0, progress, intens_thresh_min = -1E10, intens_thresh_max = 1E10;

  // variables unique to this code, used to specify 2 flag values
  int flag_val_em = -1;

  // variables for smoothing
  int used, scale2 = 7;

  // variables for merging
  int merge_x = 2, merge_y = 2, merge_z = 10, ss_mode = 1;

  // variables used for cataloguing
  vector<object_props *> detections;
  vector<int> obj_ids, check_obj_ids;
  int obj_limit = 1000, obj_batch_limit = 10000, obj_batch;
  int NO_check_obj_ids, NO_obj_ids, cat_mode;
  fstream outputfile;
  
  // variables used to scale datacube to Jy/beam
  float unit_scaling, beam_minor, beam_major, beam_SA, pixel_x_size, pixel_y_size;

  // variables used to scale velocity field postage stamps
  float cdelt3, crpix3, crval3, restfreq;
  int ctype3;

  // variables used to create output .fits file, which has flag_vals array replaced with object indices
  int mask_output_flag = -1;

  // de-bugging variables

  // check if the right number of command line parameters have been specified
  if((argc > 22) || ((argc < 15) && (argc != 2))){ cout << "WARNING: Incorrect number of arguments! Exiting.\nEnter Prototype_Catalog_v7 -h on the command line to see the command line input.\n"; return 1; }

  // get input parameters from command line
  dummy2.str("");
  dummy2.clear();
  dummy2 << argv[1];
  dummy2 >> output_code;

  // if the output_code is -h, display the command line parameters and exit
  if(output_code == "-h"){ cout << "\nUsage: ./Prototype_Catalog_v7_NP output_code mask_file data_file min_x_size min_y_size min_z_size min_LoS fill_factor flag_val Tflux_thresh_min Tflux_thresh_max merge_x merge_y merge_z [-SR] [-C] [-M] [-DO x-order y-order z-order] \n\nOutput: \n\"output_code\"_obj.cat ; file listing objects and calculated properties. \n\n" << endl; return 0; }

  // check if the right command line option was entered
  if((argc == 2) && (output_code != "-h")){ cout << "WARNING: Incorrect arguments! Exiting.\nEnter Prototype_Catalog_v6 -h on the command line to see the command line input.\n"; return 1; }

  // create initial arrays used for cataloguing and object construction
  InitObjGen(detections,NOobj,obj_limit,obj_ids,check_obj_ids,data_metric,xyz_order);

  dummy2.str("");
  dummy2.clear();
  dummy2 << argv[2];
  dummy2 >> maskfile;

  dummy2.str("");
  dummy2.clear();
  dummy2 << argv[3];
  dummy2 >> snfile;

  dummy2.str("");
  dummy2.clear();
  dummy2 << argv[4];
  dummy2 >> min_x_size;

  dummy2.str("");
  dummy2.clear();
  dummy2 << argv[5];
  dummy2 >> min_y_size;

  dummy2.str("");
  dummy2.clear();
  dummy2 << argv[6];
  dummy2 >> min_z_size;

  dummy2.str("");
  dummy2.clear();
  dummy2 << argv[7];
  dummy2 >> min_LoS_count;

  dummy2.str("");
  dummy2.clear();
  dummy2 << argv[8];
  dummy2 >> fill_factor;

  dummy2.str("");
  dummy2.clear();
  dummy2 << argv[9];
  dummy2 >> flag_val_em;

  dummy2.str("");
  dummy2.clear();
  dummy2 << argv[10];
  dummy2 >> intens_thresh_min;

  dummy2.str("");
  dummy2.clear();
  dummy2 << argv[11];
  dummy2 >> intens_thresh_max;

  dummy2.str("");
  dummy2.clear();
  dummy2 << argv[12];
  dummy2 >> merge_x;

  dummy2.str("");
  dummy2.clear();
  dummy2 << argv[13];
  dummy2 >> merge_y;

  dummy2.str("");
  dummy2.clear();
  dummy2 << argv[14];
  dummy2 >> merge_z;
 
  cat_mode = -1;
  ss_mode = 1;
  mask_output_flag = -1;
  xyz_order[0] = 1;
  xyz_order[1] = 2;
  xyz_order[2] = 3;
  for(i = 15; i < argc; i++){

    dummy2.str("");
    dummy2.clear();
    dummy2 << argv[i];
    dummy1.clear();
    dummy2 >> dummy1;
    if((dummy1 == "-sr") || (dummy1 == "-SR") || (dummy1 == "-sR") || (dummy1 == "-Sr")){ cat_mode = 1; }
    if((dummy1 == "-c") || (dummy1 == "-C")){ ss_mode = -1; }
    if((dummy1 == "-m") || (dummy1 == "-M")){ mask_output_flag = 1; }
    if((dummy1 == "-do") || (dummy1 == "-DO") || (dummy1 == "-dO") || (dummy1 == "-Do")){

      i++;
      dummy2.str("");
      dummy2.clear();
      dummy2 << argv[i];
      dummy1.clear();
      dummy2 >> xyz_order[0];
     
      i++;
      dummy2.str("");
      dummy2.clear();
      dummy2 << argv[i];
      dummy1.clear();
      dummy2 >> xyz_order[1];

      i++;
      dummy2.str("");
      dummy2.clear();
      dummy2 << argv[i];
      dummy1.clear();
      dummy2 >> xyz_order[2];

    }

  }
  
  // set overlap region of chunks
  chunk_x_overlap = merge_x + 1;
  chunk_y_overlap = merge_y + 1;
  chunk_z_overlap = merge_z + 1;

  // calculate min_v_size from filling factor and other minimum sizes of bounding volume
  if(fill_factor >= 0.0){ min_v_size = (int) ceilf(((float) min_x_size * (float) min_y_size * (float) min_z_size * fill_factor)); } else { min_v_size = (int) (-1.0 * fill_factor); }

  // display input parameters
  cout << "Input parameters are . . . " << endl;
  cout << "Output code: " << output_code << endl;
  cout << "Mask file: " << maskfile << endl;
  cout << "Source+noise cube: " << snfile << endl;
  cout << "Minimum extent in voxels along RA: " << min_x_size << endl;
  cout << "Minimum extent in voxels along Dec: " << min_y_size << endl;
  cout << "Minimum extent in voxels along frequency: " << min_z_size << endl;
  cout << "Sources must occupy at least " << min_LoS_count << " lines-of-sight through the datacube." << endl;
  cout << "For filling factor of " << fill_factor << ", minimum size in voxels: " << min_v_size << endl;
  cout << "Pseudo-total flux of sources must range from " << intens_thresh_min << " to " << intens_thresh_max << endl;
  cout << "Treating values of " << flag_val_em << " in mask file as emission source voxels." << endl;
  cout << "Using an overlap of: x = " << chunk_x_overlap << ", y = " << chunk_y_overlap << ", z = " << chunk_z_overlap << endl;
  if(ss_mode == 1){ 

    cout << "Using rectangle to link spatial dimensions." << endl;

  } else {

    cout << "Using ellipse to link spatial dimensions." << endl;

  }
  cout << "Merging objects separated by " << merge_x << " (x), " << merge_y << " (y), " << merge_z << " (z) voxels in respective dimensions." << endl;
  if(cat_mode > 0){ 

    cout << "Creating object catalogues containing sparse representations." << endl;

  } else {

    cout << "Creating object catalogues without sparse representations." << endl;

  }
  cout << "Mapping RA to data dimension: " << xyz_order[0] << ", Dec to data dimension: " << xyz_order[1] << " & Freq. to data dimension: " << xyz_order[2] << endl;

  // open input files
  status = 0;
  cout << "Attempting to open file: " << snfile << endl;
  fits_open_file(&fits_input,snfile.c_str(),READONLY,&status);
  cout << "File opened, status = " << status << endl;
  if(status > 0){ fits_report_error(stderr,status); return 1; }

  // get the size of the three dimensions of the input file
  fits_read_key(fits_input,TINT,"NAXIS1",&NOx,NULL,&status);
  fits_read_key(fits_input,TINT,"NAXIS2",&NOy,NULL,&status);
  fits_read_key(fits_input,TINT,"NAXIS3",&NOf,NULL,&status);
  if(status > 0){ 
    
    fits_report_error(stderr,status); 
    return 1; 
    
  }
  cout << NOx << " x " << NOy << " spatial elements, " << NOf << " frequency channels in data cube." << endl;
 
  // get the units of the datacube and calculate the scaling to Jy/beam
  unit_scaling = 1.0;
  fits_read_key(fits_input,TSTRING,"BUNIT",&dummy3,NULL,&status);
  if(status != 0){ 

    status = 0;
    cout << "WARNING!!! Could not find datacube units. Datacube and object properties will not be scaled. " << endl;

  } else {

    dummy1.clear();
    dummy1 = dummy3;
    for(i = 0; i < dummy1.length(); i++){ dummy1[i] = toupper(dummy1[i]); }
    cout << "Datacube units are " << dummy1 << ". Calculating conversion factor to mJy/beam." << endl;
    
    if(dummy1 == "JY/BEAM"){ 
      
      cout << "Setting unit_scaling = 1000.0 for mJy/beam." << endl;
      unit_scaling = 1000.0; 
      
    }
    if(dummy1 == "MJY/BEAM"){ 
      
      cout << "Setting unit_scaling = 1.0 for mJy/beam." << endl;
      unit_scaling = 1.0; 
      
    }
    if(dummy1 == "JY"){ 
      
      unit_scaling = 1000.0; 
      beam_SA = 1.0;
      fits_read_key(fits_input,TFLOAT,"BMAJ",&beam_major,NULL,&status);
      if(status != 0){ 
	
	status = 0;
	cout << "WARNING!!! Could not find beam semi-major axis in header. Using beam normalisation of 1." << endl;
	
      } else {
	
	fits_read_key(fits_input,TFLOAT,"BMIN",&beam_minor,NULL,&status);
	if(status != 0){
	  
	  status = 0;
	  cout << "WARNING!!! Could not find beam semi-minor axis in header. Using beam normalisation of 1." << endl;
	  
	} else {
	  
	  fits_read_key(fits_input,TFLOAT,"CDELT1",&pixel_x_size,NULL,&status);
	  pixel_x_size = fabs(pixel_x_size);
	  if(status != 0){
	    
	    status = 0;
	    cout << "WARNING!!! Could not find size of pixel in degrees along RA. Using beam normalisation of 1." << endl;
	    
	  } else {
	    
	    fits_read_key(fits_input,TFLOAT,"CDELT2",&pixel_y_size,NULL,&status);
	    pixel_y_size = fabs(pixel_y_size);
	    if(status != 0){
	      
	      status = 0;
	      cout << "WARNING!!! Could not find size of pixel in degrees along Dec. Using beam normalisation of 1." << endl;
	      
	    } else {
	      
	      beam_SA = 3.1415926535 * beam_minor * beam_major / (4.0 * logf(2.0) * pixel_x_size * pixel_y_size);
	      
	    }
	    
	  }
	  
	}
	
      }
      unit_scaling/=beam_SA;
      cout << "Setting unit_scaling = " << unit_scaling << " for mJy/beam. Unit scaling is using a beam size of " << beam_SA << " pixels = pi * (" << beam_minor <<  " degrees x " << beam_major << " degrees) / (4 * ln(2) * " << pixel_x_size << " deg/pixel x " << pixel_y_size << " deg/pixel)." << endl;
      
    }
    if(dummy1 == "JY"){ 
      
      unit_scaling = 1000.0; 
      beam_SA = 1.0;
      fits_read_key(fits_input,TFLOAT,"BMAJ",&beam_major,NULL,&status);
      if(status != 0){ 
	
	status = 0;
	cout << "WARNING!!! Could not find beam semi-major axis in header. Using beam normalisation of 1." << endl;
	
      } else {
	
	fits_read_key(fits_input,TFLOAT,"BMIN",&beam_minor,NULL,&status);
	if(status != 0){
	  
	  status = 0;
	  cout << "WARNING!!! Could not find beam semi-minor axis in header. Using beam normalisation of 1." << endl;
	  
	} else {
	  
	  fits_read_key(fits_input,TFLOAT,"CDELT1",&pixel_x_size,NULL,&status);
	  pixel_x_size = fabs(pixel_x_size);
	  if(status != 0){
	    
	    status = 0;
	    cout << "WARNING!!! Could not find size of pixel in degrees along RA. Using beam normalisation of 1." << endl;
	    
	  } else {
	    
	    fits_read_key(fits_input,TFLOAT,"CDELT2",&pixel_y_size,NULL,&status);
	    pixel_y_size = fabs(pixel_y_size);
	    if(status != 0){
	      
	      status = 0;
	      cout << "WARNING!!! Could not find size of pixel in degrees along Dec. Using beam normalisation of 1." << endl;
	      
	    } else {
	      
	      beam_SA = 3.1415926535 * beam_minor * beam_major / (4.0 * logf(2.0) * pixel_x_size * pixel_y_size);
	      
	    }
	    
	  }
	  
	}
	
      }
      unit_scaling/=beam_SA;
      cout << "Setting unit_scaling = " << unit_scaling << " for mJy/beam. Unit scaling uses a beam of size " << beam_SA << " pixels = pi * (" << beam_minor <<  " degrees x " << beam_major << " degrees) / (" << pixel_x_size << " deg/pixel x " << pixel_y_size << " deg/pixel)." << endl;
      
    }

  }

  // calculate the size of chunks, that cover the entire frequency range and are ideally square (shrink either axis if it's larger than the datacube)
  chunk_x_size = chunk_y_size = (int) floorf((sqrtf((floorf((float) (chunk_size / NOf))))));
  if(chunk_x_size > NOx){ chunk_x_size = NOx; }
  if(chunk_y_size > NOy){ chunk_y_size = NOy; }
  chunk_z_size = NOf;
  
  if(chunk_x_size < NOx){ temp_chunk_x_size = ((int) ceilf((((float) NOx / (float) (chunk_x_size - chunk_x_overlap))))); } else { temp_chunk_x_size = 1; }
  if(chunk_y_size < NOy){ temp_chunk_y_size = ((int) ceilf((((float) NOy / (float) (chunk_y_size - chunk_y_overlap))))); } else { temp_chunk_y_size = 1; }
  NO_chunks = temp_chunk_x_size * temp_chunk_y_size;

  // display the details of the 1GB chunks
  cout << "Processing file using overlapping " << (chunk_size * 4 / (1024 * 1024)) << "MB `chunks'. This is most efficiently achieved using " << NO_chunks << " chunk/s of size " << chunk_x_size << " (x) by " << chunk_y_size << " (y) by " << chunk_z_size << " (freq.) voxels. This amounts to " << temp_chunk_x_size << " (x) by " << temp_chunk_y_size << " (y) chunks." << endl;
  
  // create array to store signal
  data_vals = new float[(chunk_x_size * chunk_y_size * chunk_z_size)];
  
  // create array to store signal
  flag_vals = new int[(chunk_x_size * chunk_y_size * chunk_z_size)];
  
  // open file containing source mask
  status = 0;
  cout << "Attempting to open file: " << maskfile << endl;
  fits_open_file(&mask_input,maskfile.c_str(),READONLY,&status);
  cout << "File opened, status = " << status << endl;
  if(status > 0){ fits_report_error(stderr,status); return 1; }
  
  // for each chunk, load it into memory, then smooth it, updating progress in terms of 1GB chunk and overall
  for(j = 0; j < temp_chunk_y_size; j++){
    
    for(i = 0; i < temp_chunk_x_size; i++){
      
      cout << "\nProcessing chunk " << (1 + i + (j * temp_chunk_x_size)) << " of " << (temp_chunk_x_size * temp_chunk_y_size) << " . . . " << endl;
      
      // set parameters for region to be read from memory
      chunk_x_start = i * (chunk_x_size - chunk_x_overlap);
      if(chunk_x_start < 0){ chunk_x_start = 0; }
      chunk_y_start = j * (chunk_y_size - chunk_y_overlap);
      if(chunk_y_start < 0){ chunk_y_start = 0; }
      chunk_z_start = 0;

      fits_read_start[0] = 1 + chunk_x_start;
      fits_read_finish[0] = chunk_x_start + chunk_x_size;
      if(fits_read_finish[0] > NOx){ fits_read_finish[0] = NOx; }

      fits_read_start[1] = 1 + chunk_y_start;
      fits_read_finish[1] = chunk_y_start + chunk_y_size;
      if(fits_read_finish[1] > NOy){ fits_read_finish[1] = NOy; }
      
      fits_read_start[2] = 1;
      fits_read_finish[2] = chunk_z_size;
      
      fits_read_start[3] = 1;
      fits_read_finish[3] = 1;
      
      cout << "Covers datacube region: " << fits_read_start[0] << " <= x <= " << fits_read_finish[0] << ", " << fits_read_start[1] << " <= y <= " << fits_read_finish[1] << ", " << fits_read_start[2] << " <= z <= " << fits_read_finish[2] << endl; 

      // initialise mask array
      cout << "Initialising flag_vals array for this chunk . . . " << endl;
      for(f = 0; f < chunk_z_size; f++){
	for(y = 0; y < chunk_y_size; y++){
	  for(x = 0; x < chunk_x_size; x++){
	    flag_vals[((f * chunk_x_size * chunk_y_size) + (y * chunk_x_size) + x)] = -99;
	  }
	}
      }

      // read file containing mask into memory
      cout << "Reading source+noise flag file into memory . . . " << endl;
      fits_read_subset(mask_input,TINT,fits_read_start,fits_read_finish,fits_read_inc,NULL,flag_vals,NULL,&status);
      if(status > 0){ 
	
	fits_report_error(stderr,status); 
	return 1; 
	
      }

      // initialise data array
      cout << "Initialising data_vals array for this chunk . . . " << endl;
      for(f = 0; f < chunk_z_size; f++){
	for(y = 0; y < chunk_y_size; y++){
	  for(x = 0; x < chunk_x_size; x++){
	    data_vals[((f * chunk_x_size * chunk_y_size) + (y * chunk_x_size) + x)] = 0.0;
	  }
	}
      }
      
      // read source+noise cube into memory
      cout << "Reading source+noise file into memory . . . " << endl;
      fits_read_subset(fits_input,TFLOAT,fits_read_start,fits_read_finish,fits_read_inc,NULL,data_vals,NULL,&status);
      if(status > 0){ 
	
	fits_report_error(stderr,status); 
	return 1; 
	
      }
   
      // convert flux values into mJ/beam
      if(unit_scaling != 1.0){

	cout << "Converting fluxes to mJy/beam . . . " << endl;
	cout << "0 | |:| | : | |:| | 100% complete" << endl;
	progress = 0.0;
	for(f = 0; (f < chunk_z_size); f++){
	  for(y = 0; (y < chunk_y_size); y++){
	    for(x = 0; (x < chunk_x_size); x++){
	      
	      while(progress <= (((float) ((f * chunk_x_size * chunk_y_size) + (y * chunk_x_size) + x + 1)) / ((float) (chunk_x_size * chunk_y_size * NOf)))){ cout << "*"; cout.flush(); progress+=0.05; }
	      data_vals[((f * chunk_x_size * chunk_y_size) + (y * chunk_x_size) + x)]*=unit_scaling;
	      
	    }
	  }
	}
	cout << "* done." << endl;
      
      }

      // check that the flag_vals array consists of negative definite values
      cout << "Checking that flag_vals array consists of negative definite values . . . " << endl;
      cout << "0 | |:| | : | |:| | 100% complete" << endl;
      progress = 0.0;
      for(f = 0; (f < chunk_z_size); f++){
	for(y = 0; (y < chunk_y_size); y++){
	  for(x = 0; (x < chunk_x_size); x++){ 
	    
	    while(progress <= (((float) ((f * chunk_x_size * chunk_y_size) + (y * chunk_x_size) + x + 1)) / ((float) (chunk_x_size * chunk_y_size * NOf)))){ cout << "*"; cout.flush(); progress+=0.05; }

	    if(flag_vals[((f * chunk_x_size * chunk_y_size) + (y * chunk_x_size) + x)] == 0){ 

	      flag_vals[((f * chunk_x_size * chunk_y_size) + (y * chunk_x_size) + x)] = -99; 

	    } else if(flag_vals[((f * chunk_x_size * chunk_y_size) + (y * chunk_x_size) + x)] > 0){ 

	      // regular version
	      flag_vals[((f * chunk_x_size * chunk_y_size) + (y * chunk_x_size) + x)] = -1 * abs(flag_vals[((f * chunk_x_size * chunk_y_size) + (y * chunk_x_size) + x)]); 
	      // temporary version for processing Attilla's Duchamp code
	      //flag_vals[((f * chunk_x_size * chunk_y_size) + (y * chunk_x_size) + x)] = -1; 

	    }
	    	    
	  }
	}
      }
      cout << "* done." << endl;
      
      // create metric for accessing this data chunk in x,y,z order
      CreateMetric(data_metric,xyz_order,(fits_read_finish[0] - fits_read_start[0] + 1),(fits_read_finish[1] - fits_read_start[1] + 1),(fits_read_finish[2] - fits_read_start[2] + 1));

      // add existing objects to flag_vals array and check_obj_ids
      NO_check_obj_ids = AddObjsToChunk(flag_vals,detections,NOobj,obj_limit,chunk_x_start,chunk_y_start,chunk_z_start,(fits_read_finish[0] - fits_read_start[0] + 1),(fits_read_finish[1] - fits_read_start[1] + 1),(fits_read_finish[2] - fits_read_start[2] + 1),check_obj_ids,data_metric,xyz_order);

      // call function to create objects out of array of flagged voxels, and generate sparse representation of objects
      cout << "Creating objects using array of flag_vals and noise+source datacube . . . " << endl;
      NOobj = CreateObjects(data_vals,flag_vals,(fits_read_finish[0] - fits_read_start[0] + 1),(fits_read_finish[1] - fits_read_start[1] + 1),(fits_read_finish[2] - fits_read_start[2] + 1),chunk_x_start,chunk_y_start,chunk_z_start,merge_x,merge_y,merge_z,min_x_size,min_y_size,min_z_size,min_v_size,intens_thresh_min,intens_thresh_max,flag_val_em,NOobj,detections,obj_ids,check_obj_ids,obj_limit,NOx,NOy,NOf,ss_mode,data_metric,xyz_order);
 
      // count the number of objects remaining and display
      k = 0;
      for(obj = 0; obj < NOobj; obj++){

       	// calculate obj_batch number for this object
	obj_batch = floorf(((float) obj / (float) obj_limit));

	// increment count if not a re-initialised object
	if(detections[obj_batch][(obj - (obj_batch * obj_limit))].ShowVoxels() >= 1){ k++; }
	
      }   
      cout << "Maximum number of objects found so far is " << NOobj << ", " << k << " unique objects remain after merging then size thresholding." << endl;
   
      // for(i = 0; i < temp_chunk_x_size; i++)
    }

    // for(j = 0; j < temp_chunk_y_size; j++)
  }

  // close input files
  cout << "\nClosing fits input files . . . "; cout.flush();
  fits_close_file(fits_input,&status);
  cout << "data file closed, status = " << status << endl;
  fits_close_file(mask_input,&status);
  cout << "mask file closed, status = " << status << endl;
        
  // Apply final size threshold and intensity threshold
  cout << "\nApplying final size threshold and additional intensity threshold . . . " << endl;
  ThresholdObjs(detections,NOobj,obj_limit,min_x_size,min_y_size,min_z_size,min_v_size,intens_thresh_min,intens_thresh_max,min_LoS_count);

  // create output catalogue
  cout << "\nCreating output catalogues . . . " << endl;
  dummy1 = "";
  dummy1.clear();
  dummy1 = output_code+"_obj.cat";
  outputfile.open(dummy1.c_str(),ios::out);
  outputfile << "# Catalogue produced by Prototype_Catalog_v6.cpp.\n# " << endl;
  outputfile << "# Input parameters are . . . " << endl;
  outputfile << "# Output code: " << output_code << endl;
  outputfile << "# Mask file: " << maskfile << endl;
  outputfile << "# Source+noise cube: " << snfile << endl;
  outputfile << "# Minimum extent in voxels along RA: " << min_x_size << endl;
  outputfile << "# Minimum extent in voxels along Dec: " << min_y_size << endl;
  outputfile << "# Minimum extent in voxels along frequency: " << min_z_size << endl;
  outputfile << "# Sources must occupy at least " << min_LoS_count << " lines-of-sight through the datacube." << endl;
  outputfile << "# For filling factor of " << fill_factor << ", minimum size in voxels: " << min_v_size << endl;
  outputfile << "# Pseudo-total flux of sources must range from " << intens_thresh_min << " to " << intens_thresh_max << endl;
  outputfile << "# Treating values of " << flag_val_em << " in mask file as emission source voxels." << endl;
  outputfile << "# Using an overlap of: x = " << chunk_x_overlap << ", y = " << chunk_y_overlap << ", z = " << chunk_z_overlap << endl;
  if(ss_mode == 1){ 

    outputfile << "# Using rectangle to link spatial dimensions." << endl;

  } else {

    outputfile << "# Using ellipse to link spatial dimensions." << endl;

  }
  outputfile << "# Merging objects within a distance of: " << (merge_x + 1) << " (x), " << (merge_y + 1) << " (y), " << (merge_z + 1) << " (z)" << endl;
  outputfile << "# Mapping RA to data dimension: " << xyz_order[0] << ", Dec to data dimension: " << xyz_order[1] << " & Freq. to data dimension: " << xyz_order[2] << endl;
  outputfile << "# \n# " << endl;
  PrintCatalogueHeader(outputfile,cat_mode);
  m = CreateCatalogue(outputfile,detections,NOobj,obj_limit,cat_mode);
  outputfile.close();

  // create .fits file containing labelled mask
  if(mask_output_flag == 1){
    
    //CreateFitsMask(output_code,NOx,NOy,NOf,detections,NOobj,obj_limit);
    CreateFitsMask(output_code,NOx,NOy,NOf,detections,NOobj,obj_limit,flag_vals,chunk_x_size,chunk_y_size,temp_chunk_x_size,temp_chunk_y_size,data_metric,xyz_order);
  
    // if(mask_output_flag == 1)
  }
    
  // free up memory
  delete [] data_vals;
  delete [] flag_vals;
  FreeObjGen(detections,data_metric,xyz_order);

}





