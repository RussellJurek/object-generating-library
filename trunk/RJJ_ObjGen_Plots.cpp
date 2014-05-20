#include<iostream>
#include<RJJ_ObjGen.h>
#include<RJJ_ObjGen_Plots.h>
extern "C" {

#include<cpgplot.h>

}

using namespace std;

// functions using floats

void CreateObjPlots(int m, int plot_mode, std::string output_code, int * xyz_order, int NOx, int NOy, int NOf, vector<object_props *> & detections, int NOobj, int obj_limit, int ctype3, float crpix3, float crval3, float cdelt3, float restfreq){
  
  int g, i, j, k, obj, obj_batch, plot_points, temp_x, temp_y, temp_z;
  std::string dummy1;
  std::stringstream dummy2;
  float plot_min, plot_max, * plot_x, * plot_y, * plot_array, tr[6] = {-1.0,1.0,0.0,-1.0,0.0,1.0};
  float progress, vp_x_min, vp_x_max, vp_y_min, vp_y_max;

  std::cout << "Input limits are: " << NOx << " " << NOy << " " << NOf << std::endl;

  // reorder the datacube and subcube limits to be in x,y,z order
  temp_x = NOx;
  temp_y = NOy;
  temp_z = NOf;
  switch(xyz_order[0]){
  case 1:
    NOx = temp_x;
    break;
  case 2:
    NOx = temp_y;
    break;
  case 3:
    NOx = temp_z;
    break;
  default:
    NOx = temp_x;
    break;
  }
  switch(xyz_order[1]){
  case 1:
    NOy = temp_x;
    break;
  case 2:
    NOy = temp_y;
    break;
  case 3:
    NOy = temp_z;
    break;
  default:
    NOy = temp_y;
    break;
  }
  switch(xyz_order[2]){
  case 1:
    NOf = temp_x;
    break;
  case 2:
    NOf = temp_y;
    break;
  case 3:
    NOf = temp_z;
    break;
   default:
    NOf = temp_z;
    break;
  }

  std::cout << "Adjusting limits to: " << NOx << " " << NOy << " " << NOf << " because order is: " << xyz_order[0] << " " << xyz_order[1] << " " << xyz_order[2] << std::endl;

  // generate plots
  if((m > 1) && (plot_mode > 0)){

    if((plot_mode == 1) || (plot_mode == 3)){
      
      std::cout << "\nCreating global moment-0 and position-velocity (RA & Dec) plots . . . " << std::endl;

      // open plotting environment
      dummy1 = "";
      dummy1.clear();
      dummy1 = output_code+"_plots.ps/cps";
      std::cout << "Creating output file: " << dummy1 << std::endl;
      cpgopen(dummy1.c_str());
      cpgslw(3);
      
      // create boundary arrays
      plot_x = new float[(NOx * NOy)];
      plot_y = new float[(NOx * NOy)];
      
      // create plotting array
      plot_array = new float[(NOx * NOy)];
      
      // create array of moment 0 map values
      std::cout << "Creating moment 0 map of total intensity . . . " << std::endl;
      plot_max = CreateMoment0Map(plot_array,NOobj,detections,NOx,NOy,obj_limit);
      
      // re-scale plot_max
      plot_max = 0.3 * plot_max;

      // check if plot_max is 0.0
      if(plot_max == 0.0){ plot_max = 1.0; }

      // create plot of moment 0 map for sources with boundaries over-plotted
      std::cout << "Plotting moment 0 map of total intensity . . . " << std::endl;
      //cpgenv(-1,NOx,-1,NOy,0,0);
      cpgvstd();
      cpgqvp(1,&vp_x_min,&vp_x_max,&vp_y_min,&vp_y_max);
      vp_x_min = ((float) NOx) / (vp_x_max - vp_x_min);
      vp_y_min = ((float) NOy) / (vp_y_max - vp_y_min);
      cpgswin((-1.0 - (3.0 * vp_x_min / 200.0)),(NOx + (3.0 * vp_x_min / 200.0)),(-1.0 - (3.0 * vp_y_min / 200.0)),(NOy + (3.0 * vp_y_min / 200.0)));
      cpglab("RA","Dec","moment 0 map");
      cpgslw(2);
      cpggray(plot_array,NOx,NOy,1,NOx,1,NOy,plot_max,0.0,tr);
      cpgwedg("R",0.0,3.0,plot_max,0.0,"Integrated Intensity (mJy/beam)");
      cpgslw(3);
      
      // plot detection boundaries
      std::cout << "Creating boundaries of objects . . . " << std::endl;
      cpgsci(3);
      cpgslw(1.5);
      for(k = 0; k < NOobj; k++){
	
	obj_batch = floor(((float) k / (float) obj_limit));
	if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	plot_points = CreateMoment0Bounds(detections,(detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin() + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin() + 1),detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin(),detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin(),plot_x,plot_y,k,obj_limit);
	if(detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() < 0.0){ cpgsci(4); }
	cpgline(plot_points,plot_x,plot_y);
	cpgsci(3);
	
      } 
      cpgslw(3);
      
      // create list of object pixel centres
      plot_points = 0;
      for(k = 0; k < NOobj; k++){
	
	obj_batch = floorf(((float) k / (float) obj_limit));
	if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	plot_x[plot_points] = detections[obj_batch][(k - (obj_batch * obj_limit))].GetRA();
	plot_y[plot_points] = detections[obj_batch][(k - (obj_batch * obj_limit))].GetDEC();
	plot_points++;
	
	//for(k = 0; k < NOobj; k++) 
      }
      
      // plot object pixel centres
      cpgsci(2);
      cpgpt(plot_points,plot_x,plot_y,5);
      
      // plot axes
      cpgsci(1);
      cpgbox("BCNST",0,0,"BCNST",0,0);
      
      if(NOf > 1){

	// adjust size of plot_array if required
	if(NOy != NOf){
	
	  std::cout << "Dec and frequency dimensions are different. Adjusting size of plot_array pointer . . . " << std::endl;
	  
	  // free up memory used for plotting array
	  delete [] plot_array;
	  delete [] plot_x;
	  delete [] plot_y;
	  
	  // create plotting array
	  plot_array = new float[(NOx * NOf)];
	  plot_x = new float[(NOx * NOf)];
	  plot_y = new float[(NOx * NOf)];
	  
	}
	
	// de-bugging, create array of RA Position-Velocity map values
	std::cout << "Creating RA Position-Velocity map of total intensity . . . " << std::endl;
	plot_max = CreateRAPVPlot(plot_array,NOobj,detections,NOx,NOf,obj_limit);
	
	// re-scale plot_max
	plot_max = 0.3 * plot_max;
	
	// check if plot_max is 0.0
	if(plot_max == 0.0){ plot_max = 1.0; }
	
	// create RA position velocity intensity map
	std::cout << "Plotting RA position-velocity map of total intensity . . . " << std::endl;
	//cpgenv(-2,(NOx + 1),-2,(NOf + 1),0,0);
	cpgpage();
	cpgvstd();
	cpgqvp(1,&vp_x_min,&vp_x_max,&vp_y_min,&vp_y_max);
	vp_x_min = ((float) NOx) / (vp_x_max - vp_x_min);
	vp_y_min = ((float) NOf) / (vp_y_max - vp_y_min);
	cpgswin((-1.0 - (3.0 * vp_x_min / 200.0)),(NOx + (3.0 * vp_x_min / 200.0)),(-1.0 - (3.0 * vp_y_min / 200.0)),(NOf + (3.0 * vp_y_min / 200.0)));
	cpglab("RA","Velocity","RA Position-Velocity intensity map");
	cpgslw(2);
	cpggray(plot_array,NOx,NOf,1,NOx,1,NOf,plot_max,0.0,tr);
	cpgwedg("R",0.0,3.0,plot_max,0.0,"Integrated Intensity (mJy/beam)");
	cpgslw(3);
	
	// plot detection boundaries
	std::cout << "Creating boundaries of objects . . . " << std::endl;
	cpgsci(3);
	cpgslw(1.5);
	for(k = 0; k < NOobj; k++){
	  
	  obj_batch = floor(((float) k / (float) obj_limit));
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	  plot_points = CreateRAPVBounds(detections,(detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin() + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin() + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmin() + 1),detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin(),detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmin(),plot_x,plot_y,k,obj_limit);
	  
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() < 0.0){ cpgsci(4); }
	  cpgline(plot_points,plot_x,plot_y);
	  cpgsci(3);
	  
	}
	cpgslw(3);
	
	// create list of object pixel centres
	plot_points = 0;
	for(k = 0; k < NOobj; k++){
	  
	  obj_batch = floorf(((float) k / (float) obj_limit));   
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	  plot_x[plot_points] = detections[obj_batch][(k - (obj_batch * obj_limit))].GetRA();
	  plot_y[plot_points] = detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQ();
	  plot_points++;
	  
	  //for(k = 0; k < NOobj; k++) 
	}
	
	// plot object pixel centres
	cpgsci(2);
	cpgpt(plot_points,plot_x,plot_y,5);
	
	// plot axes
	cpgsci(1);
	cpgbox("BCNST",0,0,"BCNST",0,0);
	
	// adjust size of plot_array if required
	if(NOy != NOx){
	  
	  std::cout << "RA and Dec dimensions are different. Adjusting size of plot_array pointer . . . " << std::endl;
	  
	  // free up memory used for plotting array
	  delete [] plot_array;
	  delete [] plot_x;
	  delete [] plot_y;
	  
	  // create plotting array
	  plot_array = new float[(NOy * NOf)];
	  plot_x = new float[(NOy * NOf)];
	  plot_y = new float[(NOy * NOf)];
	  
	}
	
	// de-bugging, create array of Dec Position-Velocity map values
	std::cout << "Creating Dec Position-Velocity map of total intensity . . . " << std::endl;
	plot_max = CreateDecPVPlot(plot_array,NOobj,detections,NOy,NOf,obj_limit);
	
	// re-scale plot_max
	plot_max = 0.3 * plot_max;
	
	// check if plot_max is 0.0
	if(plot_max == 0.0){ plot_max = 1.0; }
	
	// create DEC position velocity intensity map
	std::cout << "Plotting Dec position-velocity map of total intensity . . . " << std::endl;
	//cpgenv(-2,(NOy + 1),-2,(NOf + 1),0,0);
	cpgpage();
	cpgvstd();
	cpgqvp(1,&vp_x_min,&vp_x_max,&vp_y_min,&vp_y_max);
	vp_x_min = ((float) NOy) / (vp_x_max - vp_x_min);
	vp_y_min = ((float) NOf) / (vp_y_max - vp_y_min);
	cpgswin((-1.0 - (3.0 * vp_x_min / 200.0)),(NOy + (3.0 * vp_x_min / 200.0)),(-1.0 - (3.0 * vp_y_min / 200.0)),(NOf + (3.0 * vp_y_min / 200.0)));
	cpglab("Dec","Velocity","Dec Position-Velocity intensity map");
	cpgslw(2);
	cpggray(plot_array,NOy,NOf,1,NOy,1,NOf,plot_max,0.0,tr);
	cpgwedg("R",0.0,3.0,plot_max,0.0,"Integrated Intensity (mJy/beam)");
	cpgslw(3);
	
	// plot detection boundaries
	std::cout << "Creating boundaries of objects . . . " << std::endl;
	cpgsci(3);
	cpgslw(1.5);
	for(k = 0; k < NOobj; k++){
	  
	  obj_batch = floor(((float) k / (float) obj_limit));
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	  plot_points = CreateDecPVBounds(detections,(detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin() + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin() + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmin() + 1),detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin(),detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmin(),plot_x,plot_y,k,obj_limit);
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() < 0.0){ cpgsci(4); }
	  cpgline(plot_points,plot_x,plot_y);
	  cpgsci(3);
	  
	}
	cpgslw(3);
	
	// create list of object pixel centres
	plot_points = 0;
	for(k = 0; k < NOobj; k++){
	  
	  obj_batch = floorf(((float) k / (float) obj_limit));
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	  plot_x[plot_points] = detections[obj_batch][(k - (obj_batch * obj_limit))].GetDEC();
	  plot_y[plot_points] = detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQ();
	  plot_points++;
	  
	  //for(k = 0; k < NOobj; k++) 
	}
	
	// plot object pixel centres
	cpgsci(2);
	cpgpt(plot_points,plot_x,plot_y,5);
	
	// plot axes
	cpgsci(1);
	cpgbox("BCNST",0,0,"BCNST",0,0);
	
	// if(NOf > 1)
      }
      
      // free up memory used for plotting array
      delete [] plot_array;
      delete [] plot_x;
      delete [] plot_y;
      
      // close plotting device
      cpgclos();
      
      std::cout << "Finished constructing moment-0 and position-velocity plots for entire datacube." << std::endl;
      
    }

    if((plot_mode == 2) || (plot_mode == 3)){
      
      // scale velocity field postage stamp images to km/s
      if(ctype3 >= 0){
	
	// z axis is already in velocity
	
	for(k = 0; ((k < NOobj) && (k < SPEC_PLOT_LIMIT)); k++){
	
	  // calculate obj_batch number of object
	  obj_batch = floorf(((float) k / (float) obj_limit));
	
	  // skip if this object has been re-initialised
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	 
	  for(g = 0; g < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); g++){

	    detections[obj_batch][(k - (obj_batch * obj_limit))].Set_vfield(g,(1.0 * (fabs(cdelt3)) * detections[obj_batch][(k - (obj_batch * obj_limit))].Get_vfield(g)));

	  }

	}

      } else {

	// z axis is in frequency

	for(k = 0; ((k < NOobj) && (k < SPEC_PLOT_LIMIT)); k++){
	
	  // calculate obj_batch number of object
	  obj_batch = floorf(((float) k / (float) obj_limit));
	  
	  // skip if this object has been re-initialised
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	  
	  for(g = 0; g < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); g++){
	    
	    detections[obj_batch][(k - (obj_batch * obj_limit))].Set_vfield(g,(299792.458 * restfreq * ((1.0/(crval3 + ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_vfield(g) - crpix3) * cdelt3))) - (1.0/(crval3 + ((detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQi() - crpix3) * cdelt3))))));

	  }

	}

      }

      // create plots of integrated spectra
      std::cout << "Plotting integrated spectra of individual objects . . . " << std::endl;
      
      // create plot arrays
      if(NOf > 20){

	plot_x = new float[(2 * NOf)];
	plot_y = new float[(2 * NOf)];

      } else {

	plot_x = new float[21];
	plot_y = new float[21];

      }

      // create plotting file for integrated spectra
      dummy1 = "";
      dummy1.clear();
      dummy1 = output_code+"_spectra.ps/vcps";
      std::cout << "Creating output file: " << dummy1 << std::endl;
      cpgopen(dummy1.c_str());
      cpgslw(3);
      
      progress = 0.0;
      std::cout << "0 | |:| | : | |:| | 100% complete" << std::endl;

      // plotting 10 integrated spectra on each page along with postage stamps
      i = -1;
      for(k = 0; ((k < NOobj) && (k < SPEC_PLOT_LIMIT)); k++){
	
	// calculate obj_batch number of object
	obj_batch = floorf(((float) k / (float) obj_limit));
	
	// skip if this object has been re-initialised
	if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	
	// create plotting array for object spectrum
	plot_array = new float[(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1)];

	// increment i - tracks the number of objects plotted on the current page
	i++;
	
	// move to a new page if this is a new block of 10 spectra to be plotted
	if((i % 10) == 0){ cpgpage(); }
	
	// determine range of spectra
	plot_min = 1E10;
	plot_max = -1E10;
	if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1) >= 10){
	  
	  for(g = 0; g < (2 * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); g++){
	    
	    if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g) >= plot_max){ plot_max = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g); }
	    if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g) <= plot_min){ plot_min = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g); }
	    
	  }
	  
	} else {
	  
	  for(g = 0; g < (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 11); g++){
	    
	    if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g) >= plot_max){ plot_max = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g); }
	    if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g) <= plot_min){ plot_min = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g); }
	    
	  }
	  
	}
	for(g = 0; g < (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1); g++){
	  
	  plot_array[g] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_ospec(g);
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_ospec(g) >= plot_max){ plot_max = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_ospec(g); }
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_ospec(g) <= plot_min){ plot_min = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_ospec(g); }
	  
	}
	
	// expand plot scale
	if(plot_min == plot_max){

	  plot_min-=0.025;
	  plot_max+=0.025;

	} else {

	  plot_min = plot_min - 0.025 * (plot_max - plot_min);
	  plot_max = plot_max + 0.025 * (plot_max - plot_min);
	
	}

	// set viewport for object
	g = (int) floorf(((float) i / 10.0));
	cpgsvp(0.05,0.46,(1.0 - 0.1 * ((float) i + 1 - ((float) g * 10.0)) + 0.02),(1.0 - 0.1 * ((float) i - ((float) g * 10.0)) - 0.01));
	
	// set window (taking into account datacube size) and plot spectra for object (object + reference spectra)
	if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1) >= 10){
	  
	  if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - ((int) floorf((0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))))) > 0){
	    
	    if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + ((int) floorf((0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))))) < NOf){
	      
	      cpgswin((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - ((int) floorf((0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))))),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + ((int) floorf((0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))))),plot_min,plot_max);
	      
	    } else {
	      
	      cpgswin((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - ((int) floorf((0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))))),NOf,plot_min,plot_max);
	      
	    }
	    
	  } else {
	    
	    if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + ((int) floorf((0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))))) < NOf){
	      
	      cpgswin(0,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + ((int) floorf((0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))))),plot_min,plot_max);
	      
	    } else {
	      
	      cpgswin(0,NOf,plot_min,plot_max);
	      
	    }
	    
	  }
	  
	  for(g = 0; g < (2 * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); g++){ 

	    plot_x[g] = g + (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - ((int) floorf((0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))))); 
	    plot_y[g] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g);

	  }
	  
	  cpgline((2 * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1)),plot_x,plot_y);
	  
	} else {
	  
	  if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - 5) > 0){
	    
	    if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + 5) < NOf){
	      
	      cpgswin((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - 5),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + 5),plot_min,plot_max);
	      
	    } else {
	      
	      cpgswin((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - 5),NOf,plot_min,plot_max);
	      
	    }
	    
	  } else {
	    
	    if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + 5) < NOf){
	      
	      cpgswin(0,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + 5),plot_min,plot_max);
	      
	    } else {
	      
	      cpgswin(0,NOf,plot_min,plot_max);
	      
	    }
	    
	  }
	  
	  for(g = 0; g < (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 11); g++){ 

	    plot_x[g] = g + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - 5;
	    plot_y[g] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g);

	  }
	  
	  cpgline((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 11),plot_x,plot_y);
	  
	}
	
	// plot positions of W_50 and W_20 estimates
	
	// a. plot conventional width estimates
	cpgsci(4);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w_max(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w_max(),plot_max);
	cpgsls(2);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w50_min(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w50_min(),plot_max);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w50_max(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w50_max(),plot_max);
	cpgsls(4);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w20_min(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w20_min(),plot_max);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w20_max(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w20_max(),plot_max);
	cpgsls(1);
	
	// b. plot c.d.f. based width estimates
	cpgsci(3);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw_max(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw_max(),plot_max);
	cpgsls(2);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw50_min(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw50_min(),plot_max);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw50_max(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw50_max(),plot_max);
	cpgsls(4);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw20_min(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw20_min(),plot_max);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw20_max(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw20_max(),plot_max);
	cpgsls(1);
	cpgsci(1);

	// c. plot frequency of pixel central moment, both unweighted and intensity weighted
	cpgsci(2);
	cpgsls(2);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQi(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQi(),plot_max);
	cpgsls(4);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQ(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQ(),plot_max);
	cpgsls(1);
	cpgsci(1);

	for(g = 0; g < (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1); g++){ 

	  plot_x[g] = g + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4); 
	  plot_array[g] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_ospec(g);
	  
	}
	cpgsci(2);
	if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1) > 1){
	  
	  cpgline((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1),plot_x,plot_array);
	  
	} else {
	  
	  cpgpt(1,plot_x,plot_array,-4);
	  
	}
	cpgsci(1);
	
	// add labels
	cpgsch(0.6);
	cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	cpgsch(0.7);
	dummy2.str("");
	dummy2.clear();
	dummy2 << (i + 1);
	dummy1 = "";
	dummy2 >> dummy1;
	cpglab("",dummy1.c_str(),"");
	cpgsch(1.0);
		
	// d. free plotting arrays
	delete [] plot_array;

	// plot postage stamp moment-0 map

	// 0. create plotting array
	plot_array = new float[((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1)*(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1))];
	
	// a. set viewport for postage stamp
	g = (int) floorf(((float) i / 10.0));
	cpgsvp(0.49,0.59,(1.0 - 0.1 * ((float) i + 1 - ((float) g * 10.0)) + 0.02),(1.0 - 0.1 * ((float) i - ((float) g * 10.0)) - 0.01));
	
	// b. set window for postage stamp
	if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0)) >= (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2))){

	  plot_max = 0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2));  
	  cpgswin(-1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),(-1.0 - plot_max),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1 - plot_max));
	  
	} else {
	  
	  plot_max = 0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0));
	  cpgswin((-1.0 - plot_max),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1.0 - plot_max),-1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1));	  
	    
	}
	
	// c. adjust plot_max
	plot_max = -1E10;
	plot_min = 1E10;
	for(g = 0; g < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); g++){ 
	  
	  plot_array[g] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_mom0(g);
	  if(plot_array[g] >= plot_max){ plot_max = plot_array[g]; } 
	  if(plot_array[g] <= plot_min){ plot_min = plot_array[g]; } 
	  
	}
	if(plot_min == plot_max){ 

	  plot_max+=1.0;
	  plot_min-=1.0; 

	}
	if(detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() < 0.0){ plot_max = plot_min; }
	if(plot_max == 0.0){ plot_max = 1.0; }

	// d. plot postage stamp
	cpgslw(2);
	cpggray(plot_array,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1),1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1),plot_max,0.0,tr);
	cpgslw(3);
	
	// e. plot borders
	cpgsch(0.7);
	cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	cpgsch(1.0);
	
	// f. free array --> redundant
	//delete [] plot_array;

	// plot postage stamp velocity field 

	// 0. create plotting array --> redundant
	//plot_array = new float[((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1)*(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1))];
	
	// a. set viewport for postage stamp
	g = (int) floorf(((float) i / 10.0));
	cpgsvp(0.62,0.72,(1.0 - 0.1 * ((float) i + 1 - ((float) g * 10.0)) + 0.02),(1.0 - 0.1 * ((float) i - ((float) g * 10.0)) - 0.01));
	
	// b. set window for postage stamp
	if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0)) >= (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2))){

	  plot_max = 0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2));  
	  cpgswin(-1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),(-1.0 - plot_max),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1 - plot_max));
	  
	} else {
	  
	  plot_max = 0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0));
	  cpgswin((-1.0 - plot_max),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1.0 - plot_max),-1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1));	  
	    
	}
	
	// c. adjust plot_max
	j = 0;
	for(g = 0; g < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); g++){ 
	  
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_mom0(g) != 0.0){

	    plot_array[j] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_vfield(g);
	    j++;

	  }

	}
	HeapSort(j,(plot_array - 1));
	g = ((int) floorf(((1.0/3.0) * ((float) j))));
	if(g >= 0){ plot_min = plot_array[g]; } else { plot_min = plot_array[0]; }
	g = 1 + ((int) floorf(((2.0/3.0) * ((float) j))));
	if(g < j){ plot_max = plot_array[g]; } else { plot_max = plot_array[(j - 1)]; }

	for(g = 0; g < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); g++){ 
	  
	  plot_array[g] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_vfield(g);

	}
	if((fabs(plot_min)) > (fabs(plot_max))){

	  plot_max = -1.0 * plot_min;

	} else if((fabs(plot_min)) < (fabs(plot_max))){

	  plot_min = -1.0 * plot_max;

	}
	if(plot_min == plot_max){ 

	  plot_max+=1.0;
	  plot_min-=1.0; 

	}
	if(plot_min == plot_max){

	  plot_max = 1.0;
	  plot_min = -1.0;

	}

	// d. plot postage stamp
	cpgslw(2);
	cpggray(plot_array,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1),1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1),plot_max,plot_min,tr);
	cpgslw(3);
	
	// e. plot borders
	cpgsch(0.7);
	cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	cpgsch(1.0);
	
	// f. free array
	delete [] plot_array;

	// plot postage stamp RA PV diagram

	// 0. create plotting array
	plot_array = new float[((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1)*(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))];
		
	// a. set viewport for postage stamp
	g = (int) floorf(((float) i / 10.0));
	cpgsvp(0.75,0.85,(1.0 - 0.1 * ((float) i + 1 - ((float) g * 10.0)) + 0.02),(1.0 - 0.1 * ((float) i - ((float) g * 10.0)) - 0.01));
	
	// b. set window for postage stamp
	if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0)) >= (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4))){

	  plot_max = 0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4));
	  cpgswin(-1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),(-1.0 - plot_max),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1 - plot_max));
	  
	} else {
	  
	  plot_max = 0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0));
	  cpgswin((-1.0 - plot_max),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1.0 - plot_max),-1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1));	  
	    
	}
	
	// c. adjust plot_max
	plot_max = -1E10;
	plot_min = 1E10;
	for(g = 0; g < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); g++){ 
	  
	  plot_array[g] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_RAPV(g);
	  if(plot_array[g] >= plot_max){ plot_max = plot_array[g]; } 
	  if(plot_array[g] <= plot_min){ plot_min = plot_array[g]; } 
	  
	}
	if(plot_min == plot_max){ 

	  plot_max+=1.0;
	  plot_min-=1.0; 

	}
	if(detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() < 0.0){ plot_max = plot_min; }
	if(plot_max == 0.0){ plot_max = 1.0; }

	// d. plot postage stamp
	cpgslw(2);
	cpggray(plot_array,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1),1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1),plot_max,0.0,tr);
	cpgslw(3);
	
	// e. plot borders
	cpgsch(0.7);
	cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	cpgsch(1.0);
	
	// f. free plotting array
	delete [] plot_array;

	// plot postage stamp Dec PV diagram

	// 0. create plotting array
	plot_array = new float[((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1)*(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))];
		
	// a. set viewport for postage stamp
	g = (int) floorf(((float) i / 10.0));
	cpgsvp(0.88,0.98,(1.0 - 0.1 * ((float) i + 1 - ((float) g * 10.0)) + 0.02),(1.0 - 0.1 * ((float) i - ((float) g * 10.0)) - 0.01));
	
	// b. set window for postage stamp
	if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2)) >= (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4))){

	  plot_max = 0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4));
	  cpgswin(-1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1),(-1.0 - plot_max),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1 - plot_max));
	  
	} else {
	  
	  plot_max = 0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2));
	  cpgswin((-1.0 - plot_max),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1.0 - plot_max),-1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1));	  
	    
	}
	
	// c. adjust plot_max
	plot_max = -1E10;
	plot_min = 1E10;
	for(g = 0; g < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); g++){ 
	  
	  plot_array[g] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_DECPV(g);
	  if(plot_array[g] >= plot_max){ plot_max = plot_array[g]; } 
	  if(plot_array[g] <= plot_min){ plot_min = plot_array[g]; } 
	  
	}
	if(plot_min == plot_max){ 

	  plot_max+=1.0;
	  plot_min-=1.0; 

	}
	if(detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() < 0.0){ plot_max = plot_min; }
	if(plot_max == 0.0){ plot_max = 1.0; }

	// d. plot postage stamp
	cpgslw(2);
	cpggray(plot_array,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1),1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1),1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1),plot_max,0.0,tr);
	cpgslw(3);
	
	// e. plot borders
	cpgsch(0.7);
	cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	cpgsch(1.0);

	// f. free plotting array
	delete [] plot_array;

	// display progress
	while(progress <= (((float) (k + 1)) / ((float) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
	
	//for(k = 0; k < NOobj; k++) 
      }  

      std::cout << "* done." << std::endl;

      // free up remaining plotting arrays
      delete [] plot_x;
      delete [] plot_y;

      // close plotting device
      cpgclos();
            
    }

  }

}

void CreateObjPlots(int m, int plot_mode, std::string output_code, int * xyz_order, int NOx, int NOy, int NOf, vector<object_props *> & detections, long int NOobj, int obj_limit, int ctype3, float crpix3, float crval3, float cdelt3, float restfreq){
  
  long int k,obj,obj_batch;
  int g, i, j, plot_points, temp_x, temp_y, temp_z;
  std::string dummy1;
  std::stringstream dummy2;
  float plot_min, plot_max, * plot_x, * plot_y, * plot_array, tr[6] = {-1.0,1.0,0.0,-1.0,0.0,1.0};
  float progress, vp_x_min, vp_x_max, vp_y_min, vp_y_max;

  std::cout << "Input limits are: " << NOx << " " << NOy << " " << NOf << std::endl;

  // reorder the datacube and subcube limits to be in x,y,z order
  temp_x = NOx;
  temp_y = NOy;
  temp_z = NOf;
  switch(xyz_order[0]){
  case 1:
    NOx = temp_x;
    break;
  case 2:
    NOx = temp_y;
    break;
  case 3:
    NOx = temp_z;
    break;
  default:
    NOx = temp_x;
    break;
  }
  switch(xyz_order[1]){
  case 1:
    NOy = temp_x;
    break;
  case 2:
    NOy = temp_y;
    break;
  case 3:
    NOy = temp_z;
    break;
  default:
    NOy = temp_y;
    break;
  }
  switch(xyz_order[2]){
  case 1:
    NOf = temp_x;
    break;
  case 2:
    NOf = temp_y;
    break;
  case 3:
    NOf = temp_z;
    break;
   default:
    NOf = temp_z;
    break;
  }

  std::cout << "Adjusting limits to: " << NOx << " " << NOy << " " << NOf << " because order is: " << xyz_order[0] << " " << xyz_order[1] << " " << xyz_order[2] << std::endl;

  // generate plots
  if((m > 1) && (plot_mode > 0)){

    if((plot_mode == 1) || (plot_mode == 3)){
      
      std::cout << "\nCreating global moment-0 and position-velocity (RA & Dec) plots . . . " << std::endl;

      // open plotting environment
      dummy1 = "";
      dummy1.clear();
      dummy1 = output_code+"_plots.ps/cps";
      std::cout << "Creating output file: " << dummy1 << std::endl;
      cpgopen(dummy1.c_str());
      cpgslw(3);
      
      // create boundary arrays
      plot_x = new float[(NOx * NOy)];
      plot_y = new float[(NOx * NOy)];
      
      // create plotting array
      plot_array = new float[(NOx * NOy)];
      
      // create array of moment 0 map values
      std::cout << "Creating moment 0 map of total intensity . . . " << std::endl;
      plot_max = CreateMoment0Map(plot_array,NOobj,detections,NOx,NOy,obj_limit);
      
      // re-scale plot_max
      plot_max = 0.3 * plot_max;

      // check if plot_max is 0.0
      if(plot_max == 0.0){ plot_max = 1.0; }

      // create plot of moment 0 map for sources with boundaries over-plotted
      std::cout << "Plotting moment 0 map of total intensity . . . " << std::endl;
      //cpgenv(-1,NOx,-1,NOy,0,0);
      cpgvstd();
      cpgqvp(1,&vp_x_min,&vp_x_max,&vp_y_min,&vp_y_max);
      vp_x_min = ((float) NOx) / (vp_x_max - vp_x_min);
      vp_y_min = ((float) NOy) / (vp_y_max - vp_y_min);
      cpgswin((-1.0 - (3.0 * vp_x_min / 200.0)),(NOx + (3.0 * vp_x_min / 200.0)),(-1.0 - (3.0 * vp_y_min / 200.0)),(NOy + (3.0 * vp_y_min / 200.0)));
      cpglab("RA","Dec","moment 0 map");
      cpgslw(2);
      cpggray(plot_array,NOx,NOy,1,NOx,1,NOy,plot_max,0.0,tr);
      cpgwedg("R",0.0,3.0,plot_max,0.0,"Integrated Intensity (mJy/beam)");
      cpgslw(3);
      
      // plot detection boundaries
      std::cout << "Creating boundaries of objects . . . " << std::endl;
      cpgsci(3);
      cpgslw(1.5);
      for(k = 0; k < NOobj; k++){
	
	obj_batch = (long int) floor(((double) k / (double) obj_limit));
	if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	plot_points = CreateMoment0Bounds(detections,(detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin() + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin() + 1),detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin(),detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin(),plot_x,plot_y,k,obj_limit);
	if(detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() < 0.0){ cpgsci(4); }
	cpgline(plot_points,plot_x,plot_y);
	cpgsci(3);
	
      } 
      cpgslw(3);
      
      // create list of object pixel centres
      plot_points = 0;
      for(k = 0; k < NOobj; k++){
	
	obj_batch = (long int) floor(((double) k / (double) obj_limit));
	if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	plot_x[plot_points] = detections[obj_batch][(k - (obj_batch * obj_limit))].GetRA();
	plot_y[plot_points] = detections[obj_batch][(k - (obj_batch * obj_limit))].GetDEC();
	plot_points++;
	
	//for(k = 0; k < NOobj; k++) 
      }
      
      // plot object pixel centres
      cpgsci(2);
      cpgpt(plot_points,plot_x,plot_y,5);
      
      // plot axes
      cpgsci(1);
      cpgbox("BCNST",0,0,"BCNST",0,0);
      
      if(NOf > 1){

	// adjust size of plot_array if required
	if(NOy != NOf){
	
	  std::cout << "Dec and frequency dimensions are different. Adjusting size of plot_array pointer . . . " << std::endl;
	  
	  // free up memory used for plotting array
	  delete [] plot_array;
	  delete [] plot_x;
	  delete [] plot_y;
	  
	  // create plotting array
	  plot_array = new float[(NOx * NOf)];
	  plot_x = new float[(NOx * NOf)];
	  plot_y = new float[(NOx * NOf)];
	  
	}
	
	// de-bugging, create array of RA Position-Velocity map values
	std::cout << "Creating RA Position-Velocity map of total intensity . . . " << std::endl;
	plot_max = CreateRAPVPlot(plot_array,NOobj,detections,NOx,NOf,obj_limit);
	
	// re-scale plot_max
	plot_max = 0.3 * plot_max;
	
	// check if plot_max is 0.0
	if(plot_max == 0.0){ plot_max = 1.0; }
	
	// create RA position velocity intensity map
	std::cout << "Plotting RA position-velocity map of total intensity . . . " << std::endl;
	//cpgenv(-2,(NOx + 1),-2,(NOf + 1),0,0);
	cpgpage();
	cpgvstd();
	cpgqvp(1,&vp_x_min,&vp_x_max,&vp_y_min,&vp_y_max);
	vp_x_min = ((float) NOx) / (vp_x_max - vp_x_min);
	vp_y_min = ((float) NOf) / (vp_y_max - vp_y_min);
	cpgswin((-1.0 - (3.0 * vp_x_min / 200.0)),(NOx + (3.0 * vp_x_min / 200.0)),(-1.0 - (3.0 * vp_y_min / 200.0)),(NOf + (3.0 * vp_y_min / 200.0)));
	cpglab("RA","Velocity","RA Position-Velocity intensity map");
	cpgslw(2);
	cpggray(plot_array,NOx,NOf,1,NOx,1,NOf,plot_max,0.0,tr);
	cpgwedg("R",0.0,3.0,plot_max,0.0,"Integrated Intensity (mJy/beam)");
	cpgslw(3);
	
	// plot detection boundaries
	std::cout << "Creating boundaries of objects . . . " << std::endl;
	cpgsci(3);
	cpgslw(1.5);
	for(k = 0; k < NOobj; k++){
	  
	  obj_batch = (long int) floor(((double) k / (double) obj_limit));
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	  plot_points = CreateRAPVBounds(detections,(detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin() + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin() + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmin() + 1),detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin(),detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmin(),plot_x,plot_y,k,obj_limit);
	  
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() < 0.0){ cpgsci(4); }
	  cpgline(plot_points,plot_x,plot_y);
	  cpgsci(3);
	  
	}
	cpgslw(3);
	
	// create list of object pixel centres
	plot_points = 0;
	for(k = 0; k < NOobj; k++){
	  
	  obj_batch = (long int) floor(((double) k / (double) obj_limit));   
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	  plot_x[plot_points] = detections[obj_batch][(k - (obj_batch * obj_limit))].GetRA();
	  plot_y[plot_points] = detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQ();
	  plot_points++;
	  
	  //for(k = 0; k < NOobj; k++) 
	}
	
	// plot object pixel centres
	cpgsci(2);
	cpgpt(plot_points,plot_x,plot_y,5);
	
	// plot axes
	cpgsci(1);
	cpgbox("BCNST",0,0,"BCNST",0,0);
	
	// adjust size of plot_array if required
	if(NOy != NOx){
	  
	  std::cout << "RA and Dec dimensions are different. Adjusting size of plot_array pointer . . . " << std::endl;
	  
	  // free up memory used for plotting array
	  delete [] plot_array;
	  delete [] plot_x;
	  delete [] plot_y;
	  
	  // create plotting array
	  plot_array = new float[(NOy * NOf)];
	  plot_x = new float[(NOy * NOf)];
	  plot_y = new float[(NOy * NOf)];
	  
	}
	
	// de-bugging, create array of Dec Position-Velocity map values
	std::cout << "Creating Dec Position-Velocity map of total intensity . . . " << std::endl;
	plot_max = CreateDecPVPlot(plot_array,NOobj,detections,NOy,NOf,obj_limit);
	
	// re-scale plot_max
	plot_max = 0.3 * plot_max;
	
	// check if plot_max is 0.0
	if(plot_max == 0.0){ plot_max = 1.0; }
	
	// create DEC position velocity intensity map
	std::cout << "Plotting Dec position-velocity map of total intensity . . . " << std::endl;
	//cpgenv(-2,(NOy + 1),-2,(NOf + 1),0,0);
	cpgpage();
	cpgvstd();
	cpgqvp(1,&vp_x_min,&vp_x_max,&vp_y_min,&vp_y_max);
	vp_x_min = ((float) NOy) / (vp_x_max - vp_x_min);
	vp_y_min = ((float) NOf) / (vp_y_max - vp_y_min);
	cpgswin((-1.0 - (3.0 * vp_x_min / 200.0)),(NOy + (3.0 * vp_x_min / 200.0)),(-1.0 - (3.0 * vp_y_min / 200.0)),(NOf + (3.0 * vp_y_min / 200.0)));
	cpglab("Dec","Velocity","Dec Position-Velocity intensity map");
	cpgslw(2);
	cpggray(plot_array,NOy,NOf,1,NOy,1,NOf,plot_max,0.0,tr);
	cpgwedg("R",0.0,3.0,plot_max,0.0,"Integrated Intensity (mJy/beam)");
	cpgslw(3);
	
	// plot detection boundaries
	std::cout << "Creating boundaries of objects . . . " << std::endl;
	cpgsci(3);
	cpgslw(1.5);
	for(k = 0; k < NOobj; k++){
	  
	  obj_batch = (long int) floor(((double) k / (double) obj_limit));
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	  plot_points = CreateDecPVBounds(detections,(detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin() + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin() + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmin() + 1),detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin(),detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmin(),plot_x,plot_y,k,obj_limit);
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() < 0.0){ cpgsci(4); }
	  cpgline(plot_points,plot_x,plot_y);
	  cpgsci(3);
	  
	}
	cpgslw(3);
	
	// create list of object pixel centres
	plot_points = 0;
	for(k = 0; k < NOobj; k++){
	  
	  obj_batch = (long int) floor(((double) k / (double) obj_limit));
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	  plot_x[plot_points] = detections[obj_batch][(k - (obj_batch * obj_limit))].GetDEC();
	  plot_y[plot_points] = detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQ();
	  plot_points++;
	  
	  //for(k = 0; k < NOobj; k++) 
	}
	
	// plot object pixel centres
	cpgsci(2);
	cpgpt(plot_points,plot_x,plot_y,5);
	
	// plot axes
	cpgsci(1);
	cpgbox("BCNST",0,0,"BCNST",0,0);
	
	// if(NOf > 1)
      }
      
      // free up memory used for plotting array
      delete [] plot_array;
      delete [] plot_x;
      delete [] plot_y;
      
      // close plotting device
      cpgclos();
      
      std::cout << "Finished constructing moment-0 and position-velocity plots for entire datacube." << std::endl;
      
    }

    if((plot_mode == 2) || (plot_mode == 3)){
      
      // scale velocity field postage stamp images to km/s
      if(ctype3 >= 0){
	
	// z axis is already in velocity
	
	for(k = 0; ((k < NOobj) && (k < SPEC_PLOT_LIMIT)); k++){
	
	  // calculate obj_batch number of object
	  obj_batch = (long int) floor(((double) k / (double) obj_limit));
	
	  // skip if this object has been re-initialised
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	 
	  for(g = 0; g < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); g++){

	    detections[obj_batch][(k - (obj_batch * obj_limit))].Set_vfield(g,(1.0 * (fabs(cdelt3)) * detections[obj_batch][(k - (obj_batch * obj_limit))].Get_vfield(g)));

	  }

	}

      } else {

	// z axis is in frequency

	for(k = 0; ((k < NOobj) && (k < SPEC_PLOT_LIMIT)); k++){
	
	  // calculate obj_batch number of object
	  obj_batch = (long int) floor(((double) k / (double) obj_limit));
	  
	  // skip if this object has been re-initialised
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	  
	  for(g = 0; g < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); g++){
	    
	    detections[obj_batch][(k - (obj_batch * obj_limit))].Set_vfield(g,(299792.458 * restfreq * ((1.0/(crval3 + ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_vfield(g) - crpix3) * cdelt3))) - (1.0/(crval3 + ((detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQi() - crpix3) * cdelt3))))));

	  }

	}

      }

      // create plots of integrated spectra
      std::cout << "Plotting integrated spectra of individual objects . . . " << std::endl;
      
      // create plot arrays
      if(NOf > 20){

	plot_x = new float[(2 * NOf)];
	plot_y = new float[(2 * NOf)];

      } else {

	plot_x = new float[21];
	plot_y = new float[21];

      }

      // create plotting file for integrated spectra
      dummy1 = "";
      dummy1.clear();
      dummy1 = output_code+"_spectra.ps/vcps";
      std::cout << "Creating output file: " << dummy1 << std::endl;
      cpgopen(dummy1.c_str());
      cpgslw(3);
      
      progress = 0.0;
      std::cout << "0 | |:| | : | |:| | 100% complete" << std::endl;

      // plotting 10 integrated spectra on each page along with postage stamps
      i = -1;
      for(k = 0; ((k < NOobj) && (k < SPEC_PLOT_LIMIT)); k++){
	
	// calculate obj_batch number of object
	obj_batch = (long int) floor(((double) k / (double) obj_limit));
	
	// skip if this object has been re-initialised
	if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	
	// create plotting array for object spectrum
	plot_array = new float[(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1)];

	// increment i - tracks the number of objects plotted on the current page
	i++;
	
	// move to a new page if this is a new block of 10 spectra to be plotted
	if((i % 10) == 0){ cpgpage(); }
	
	// determine range of spectra
	plot_min = 1E10;
	plot_max = -1E10;
	if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1) >= 10){
	  
	  for(g = 0; g < (2 * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); g++){
	    
	    if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g) >= plot_max){ plot_max = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g); }
	    if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g) <= plot_min){ plot_min = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g); }
	    
	  }
	  
	} else {
	  
	  for(g = 0; g < (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 11); g++){
	    
	    if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g) >= plot_max){ plot_max = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g); }
	    if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g) <= plot_min){ plot_min = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g); }
	    
	  }
	  
	}
	for(g = 0; g < (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1); g++){
	  
	  plot_array[g] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_ospec(g);
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_ospec(g) >= plot_max){ plot_max = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_ospec(g); }
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_ospec(g) <= plot_min){ plot_min = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_ospec(g); }
	  
	}
	
	// expand plot scale
	if(plot_min == plot_max){

	  plot_min-=0.025;
	  plot_max+=0.025;

	} else {

	  plot_min = plot_min - 0.025 * (plot_max - plot_min);
	  plot_max = plot_max + 0.025 * (plot_max - plot_min);
	
	}

	// set viewport for object
	g = (int) floorf(((float) i / 10.0));
	cpgsvp(0.05,0.46,(1.0 - 0.1 * ((float) i + 1 - ((float) g * 10.0)) + 0.02),(1.0 - 0.1 * ((float) i - ((float) g * 10.0)) - 0.01));
	
	// set window (taking into account datacube size) and plot spectra for object (object + reference spectra)
	if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1) >= 10){
	  
	  if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - ((int) floorf((0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))))) > 0){
	    
	    if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + ((int) floorf((0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))))) < NOf){
	      
	      cpgswin((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - ((int) floorf((0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))))),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + ((int) floorf((0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))))),plot_min,plot_max);
	      
	    } else {
	      
	      cpgswin((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - ((int) floorf((0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))))),NOf,plot_min,plot_max);
	      
	    }
	    
	  } else {
	    
	    if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + ((int) floorf((0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))))) < NOf){
	      
	      cpgswin(0,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + ((int) floorf((0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))))),plot_min,plot_max);
	      
	    } else {
	      
	      cpgswin(0,NOf,plot_min,plot_max);
	      
	    }
	    
	  }
	  
	  for(g = 0; g < (2 * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); g++){ 

	    plot_x[g] = g + (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - ((int) floorf((0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))))); 
	    plot_y[g] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g);

	  }
	  
	  cpgline((2 * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1)),plot_x,plot_y);
	  
	} else {
	  
	  if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - 5) > 0){
	    
	    if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + 5) < NOf){
	      
	      cpgswin((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - 5),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + 5),plot_min,plot_max);
	      
	    } else {
	      
	      cpgswin((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - 5),NOf,plot_min,plot_max);
	      
	    }
	    
	  } else {
	    
	    if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + 5) < NOf){
	      
	      cpgswin(0,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + 5),plot_min,plot_max);
	      
	    } else {
	      
	      cpgswin(0,NOf,plot_min,plot_max);
	      
	    }
	    
	  }
	  
	  for(g = 0; g < (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 11); g++){ 

	    plot_x[g] = g + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - 5;
	    plot_y[g] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g);

	  }
	  
	  cpgline((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 11),plot_x,plot_y);
	  
	}
	
	// plot positions of W_50 and W_20 estimates
	
	// a. plot conventional width estimates
	cpgsci(4);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w_max(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w_max(),plot_max);
	cpgsls(2);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w50_min(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w50_min(),plot_max);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w50_max(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w50_max(),plot_max);
	cpgsls(4);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w20_min(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w20_min(),plot_max);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w20_max(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w20_max(),plot_max);
	cpgsls(1);
	
	// b. plot c.d.f. based width estimates
	cpgsci(3);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw_max(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw_max(),plot_max);
	cpgsls(2);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw50_min(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw50_min(),plot_max);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw50_max(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw50_max(),plot_max);
	cpgsls(4);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw20_min(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw20_min(),plot_max);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw20_max(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw20_max(),plot_max);
	cpgsls(1);
	cpgsci(1);

	// c. plot frequency of pixel central moment, both unweighted and intensity weighted
	cpgsci(2);
	cpgsls(2);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQi(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQi(),plot_max);
	cpgsls(4);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQ(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQ(),plot_max);
	cpgsls(1);
	cpgsci(1);

	for(g = 0; g < (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1); g++){ 

	  plot_x[g] = g + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4); 
	  plot_array[g] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_ospec(g);
	  
	}
	cpgsci(2);
	if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1) > 1){
	  
	  cpgline((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1),plot_x,plot_array);
	  
	} else {
	  
	  cpgpt(1,plot_x,plot_array,-4);
	  
	}
	cpgsci(1);
	
	// add labels
	cpgsch(0.6);
	cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	cpgsch(0.7);
	dummy2.str("");
	dummy2.clear();
	dummy2 << (i + 1);
	dummy1 = "";
	dummy2 >> dummy1;
	cpglab("",dummy1.c_str(),"");
	cpgsch(1.0);
		
	// d. free plotting arrays
	delete [] plot_array;

	// plot postage stamp moment-0 map

	// 0. create plotting array
	plot_array = new float[((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1)*(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1))];
	
	// a. set viewport for postage stamp
	g = (int) floorf(((float) i / 10.0));
	cpgsvp(0.49,0.59,(1.0 - 0.1 * ((float) i + 1 - ((float) g * 10.0)) + 0.02),(1.0 - 0.1 * ((float) i - ((float) g * 10.0)) - 0.01));
	
	// b. set window for postage stamp
	if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0)) >= (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2))){

	  plot_max = 0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2));  
	  cpgswin(-1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),(-1.0 - plot_max),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1 - plot_max));
	  
	} else {
	  
	  plot_max = 0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0));
	  cpgswin((-1.0 - plot_max),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1.0 - plot_max),-1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1));	  
	    
	}
	
	// c. adjust plot_max
	plot_max = -1E10;
	plot_min = 1E10;
	for(g = 0; g < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); g++){ 
	  
	  plot_array[g] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_mom0(g);
	  if(plot_array[g] >= plot_max){ plot_max = plot_array[g]; } 
	  if(plot_array[g] <= plot_min){ plot_min = plot_array[g]; } 
	  
	}
	if(plot_min == plot_max){ 

	  plot_max+=1.0;
	  plot_min-=1.0; 

	}
	if(detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() < 0.0){ plot_max = plot_min; }
	if(plot_max == 0.0){ plot_max = 1.0; }

	// d. plot postage stamp
	cpgslw(2);
	cpggray(plot_array,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1),1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1),plot_max,0.0,tr);
	cpgslw(3);
	
	// e. plot borders
	cpgsch(0.7);
	cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	cpgsch(1.0);
	
	// f. free array --> redundant
	//delete [] plot_array;

	// plot postage stamp velocity field 

	// 0. create plotting array --> redundant
	//plot_array = new float[((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1)*(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1))];
	
	// a. set viewport for postage stamp
	g = (int) floorf(((float) i / 10.0));
	cpgsvp(0.62,0.72,(1.0 - 0.1 * ((float) i + 1 - ((float) g * 10.0)) + 0.02),(1.0 - 0.1 * ((float) i - ((float) g * 10.0)) - 0.01));
	
	// b. set window for postage stamp
	if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0)) >= (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2))){

	  plot_max = 0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2));  
	  cpgswin(-1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),(-1.0 - plot_max),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1 - plot_max));
	  
	} else {
	  
	  plot_max = 0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0));
	  cpgswin((-1.0 - plot_max),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1.0 - plot_max),-1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1));	  
	    
	}
	
	// c. adjust plot_max
	j = 0;
	for(g = 0; g < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); g++){ 
	  
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_mom0(g) != 0.0){

	    plot_array[j] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_vfield(g);
	    j++;

	  }

	}
	HeapSort(j,(plot_array - 1));
	g = ((int) floorf(((1.0/3.0) * ((float) j))));
	if(g >= 0){ plot_min = plot_array[g]; } else { plot_min = plot_array[0]; }
	g = 1 + ((int) floorf(((2.0/3.0) * ((float) j))));
	if(g < j){ plot_max = plot_array[g]; } else { plot_max = plot_array[(j - 1)]; }

	for(g = 0; g < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); g++){ 
	  
	  plot_array[g] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_vfield(g);

	}
	if((fabs(plot_min)) > (fabs(plot_max))){

	  plot_max = -1.0 * plot_min;

	} else if((fabs(plot_min)) < (fabs(plot_max))){

	  plot_min = -1.0 * plot_max;

	}
	if(plot_min == plot_max){ 

	  plot_max+=1.0;
	  plot_min-=1.0; 

	}
	if(plot_min == plot_max){

	  plot_max = 1.0;
	  plot_min = -1.0;

	}

	// d. plot postage stamp
	cpgslw(2);
	cpggray(plot_array,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1),1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1),plot_max,plot_min,tr);
	cpgslw(3);
	
	// e. plot borders
	cpgsch(0.7);
	cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	cpgsch(1.0);
	
	// f. free array
	delete [] plot_array;

	// plot postage stamp RA PV diagram

	// 0. create plotting array
	plot_array = new float[((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1)*(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))];
		
	// a. set viewport for postage stamp
	g = (int) floorf(((float) i / 10.0));
	cpgsvp(0.75,0.85,(1.0 - 0.1 * ((float) i + 1 - ((float) g * 10.0)) + 0.02),(1.0 - 0.1 * ((float) i - ((float) g * 10.0)) - 0.01));
	
	// b. set window for postage stamp
	if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0)) >= (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4))){

	  plot_max = 0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4));
	  cpgswin(-1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),(-1.0 - plot_max),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1 - plot_max));
	  
	} else {
	  
	  plot_max = 0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0));
	  cpgswin((-1.0 - plot_max),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1.0 - plot_max),-1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1));	  
	    
	}
	
	// c. adjust plot_max
	plot_max = -1E10;
	plot_min = 1E10;
	for(g = 0; g < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); g++){ 
	  
	  plot_array[g] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_RAPV(g);
	  if(plot_array[g] >= plot_max){ plot_max = plot_array[g]; } 
	  if(plot_array[g] <= plot_min){ plot_min = plot_array[g]; } 
	  
	}
	if(plot_min == plot_max){ 

	  plot_max+=1.0;
	  plot_min-=1.0; 

	}
	if(detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() < 0.0){ plot_max = plot_min; }
	if(plot_max == 0.0){ plot_max = 1.0; }

	// d. plot postage stamp
	cpgslw(2);
	cpggray(plot_array,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1),1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1),plot_max,0.0,tr);
	cpgslw(3);
	
	// e. plot borders
	cpgsch(0.7);
	cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	cpgsch(1.0);
	
	// f. free plotting array
	delete [] plot_array;

	// plot postage stamp Dec PV diagram

	// 0. create plotting array
	plot_array = new float[((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1)*(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))];
		
	// a. set viewport for postage stamp
	g = (int) floorf(((float) i / 10.0));
	cpgsvp(0.88,0.98,(1.0 - 0.1 * ((float) i + 1 - ((float) g * 10.0)) + 0.02),(1.0 - 0.1 * ((float) i - ((float) g * 10.0)) - 0.01));
	
	// b. set window for postage stamp
	if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2)) >= (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4))){

	  plot_max = 0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4));
	  cpgswin(-1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1),(-1.0 - plot_max),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1 - plot_max));
	  
	} else {
	  
	  plot_max = 0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2));
	  cpgswin((-1.0 - plot_max),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1.0 - plot_max),-1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1));	  
	    
	}
	
	// c. adjust plot_max
	plot_max = -1E10;
	plot_min = 1E10;
	for(g = 0; g < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); g++){ 
	  
	  plot_array[g] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_DECPV(g);
	  if(plot_array[g] >= plot_max){ plot_max = plot_array[g]; } 
	  if(plot_array[g] <= plot_min){ plot_min = plot_array[g]; } 
	  
	}
	if(plot_min == plot_max){ 

	  plot_max+=1.0;
	  plot_min-=1.0; 

	}
	if(detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() < 0.0){ plot_max = plot_min; }
	if(plot_max == 0.0){ plot_max = 1.0; }

	// d. plot postage stamp
	cpgslw(2);
	cpggray(plot_array,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1),1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1),1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1),plot_max,0.0,tr);
	cpgslw(3);
	
	// e. plot borders
	cpgsch(0.7);
	cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	cpgsch(1.0);

	// f. free plotting array
	delete [] plot_array;

	// display progress
	while(progress <= (((float) (k + 1)) / ((float) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
	
	//for(k = 0; k < NOobj; k++) 
      }  

      std::cout << "* done." << std::endl;

      // free up remaining plotting arrays
      delete [] plot_x;
      delete [] plot_y;

      // close plotting device
      cpgclos();
            
    }

  }

}

void HeapSort(int n, float ra[]){

  int i,ir,j,l;
  float rra;

  if (n<2) return;
  l = (n >> 1) + 1;
  ir = n;
  for(;;){

    if (l > 1) {
      rra=ra[--l];
    } else {
      rra=ra[ir];
      ra[ir]=ra[1];
      if (--ir == 1) {
	ra[1]=rra;
	break;
      }
    }
    i=l;
    j=l+l;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) j++;
      if (rra < ra[j]) {
	ra[i]=ra[j];
	i=j;
	j <<= 1;
      } else break;
    }
    ra[i]=rra;

  }

}

// functions using doubles

void CreateObjPlots(int m, int plot_mode, std::string output_code, int * xyz_order, int NOx, int NOy, int NOf, vector<object_props_dbl *> & detections, int NOobj, int obj_limit, int ctype3, float crpix3, float crval3, float cdelt3, float restfreq){
  
  int g, i, j, k, obj, obj_batch, plot_points, temp_x, temp_y, temp_z;
  std::string dummy1;
  std::stringstream dummy2;
  float plot_min, plot_max, * plot_x, * plot_y, * plot_array, tr[6] = {-1.0,1.0,0.0,-1.0,0.0,1.0};
  float progress, vp_x_min, vp_x_max, vp_y_min, vp_y_max;

  std::cout << "Input limits are: " << NOx << " " << NOy << " " << NOf << std::endl;

  // reorder the datacube and subcube limits to be in x,y,z order
  temp_x = NOx;
  temp_y = NOy;
  temp_z = NOf;
  switch(xyz_order[0]){
  case 1:
    NOx = temp_x;
    break;
  case 2:
    NOx = temp_y;
    break;
  case 3:
    NOx = temp_z;
    break;
  default:
    NOx = temp_x;
    break;
  }
  switch(xyz_order[1]){
  case 1:
    NOy = temp_x;
    break;
  case 2:
    NOy = temp_y;
    break;
  case 3:
    NOy = temp_z;
    break;
  default:
    NOy = temp_y;
    break;
  }
  switch(xyz_order[2]){
  case 1:
    NOf = temp_x;
    break;
  case 2:
    NOf = temp_y;
    break;
  case 3:
    NOf = temp_z;
    break;
   default:
    NOf = temp_z;
    break;
  }

  std::cout << "Adjusting limits to: " << NOx << " " << NOy << " " << NOf << " because order is: " << xyz_order[0] << " " << xyz_order[1] << " " << xyz_order[2] << std::endl;

  // generate plots
  if((m > 1) && (plot_mode > 0)){

    if((plot_mode == 1) || (plot_mode == 3)){
      
      std::cout << "\nCreating global moment-0 and position-velocity (RA & Dec) plots . . . " << std::endl;

      // open plotting environment
      dummy1 = "";
      dummy1.clear();
      dummy1 = output_code+"_plots.ps/cps";
      std::cout << "Creating output file: " << dummy1 << std::endl;
      cpgopen(dummy1.c_str());
      cpgslw(3);
      
      // create boundary arrays
      plot_x = new float[(NOx * NOy)];
      plot_y = new float[(NOx * NOy)];
      
      // create plotting array
      plot_array = new float[(NOx * NOy)];
      
      // create array of moment 0 map values
      std::cout << "Creating moment 0 map of total intensity . . . " << std::endl;
      plot_max = CreateMoment0Map(plot_array,NOobj,detections,NOx,NOy,obj_limit);
      
      // re-scale plot_max
      plot_max = 0.3 * plot_max;

      // check if plot_max is 0.0
      if(plot_max == 0.0){ plot_max = 1.0; }

      // create plot of moment 0 map for sources with boundaries over-plotted
      std::cout << "Plotting moment 0 map of total intensity . . . " << std::endl;
      //cpgenv(-1,NOx,-1,NOy,0,0);
      cpgvstd();
      cpgqvp(1,&vp_x_min,&vp_x_max,&vp_y_min,&vp_y_max);
      vp_x_min = ((float) NOx) / (vp_x_max - vp_x_min);
      vp_y_min = ((float) NOy) / (vp_y_max - vp_y_min);
      cpgswin((-1.0 - (3.0 * vp_x_min / 200.0)),(NOx + (3.0 * vp_x_min / 200.0)),(-1.0 - (3.0 * vp_y_min / 200.0)),(NOy + (3.0 * vp_y_min / 200.0)));
      cpglab("RA","Dec","moment 0 map");
      cpgslw(2);
      cpggray(plot_array,NOx,NOy,1,NOx,1,NOy,plot_max,0.0,tr);
      cpgwedg("R",0.0,3.0,plot_max,0.0,"Integrated Intensity (mJy/beam)");
      cpgslw(3);
      
      // plot detection boundaries
      std::cout << "Creating boundaries of objects . . . " << std::endl;
      cpgsci(3);
      cpgslw(1.5);
      for(k = 0; k < NOobj; k++){
	
	obj_batch = floor(((float) k / (float) obj_limit));
	if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	plot_points = CreateMoment0Bounds(detections,(detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin() + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin() + 1),detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin(),detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin(),plot_x,plot_y,k,obj_limit);
	if(detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() < 0.0){ cpgsci(4); }
	cpgline(plot_points,plot_x,plot_y);
	cpgsci(3);
	
      } 
      cpgslw(3);
      
      // create list of object pixel centres
      plot_points = 0;
      for(k = 0; k < NOobj; k++){
	
	obj_batch = floorf(((float) k / (float) obj_limit));
	if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	plot_x[plot_points] = detections[obj_batch][(k - (obj_batch * obj_limit))].GetRA();
	plot_y[plot_points] = detections[obj_batch][(k - (obj_batch * obj_limit))].GetDEC();
	plot_points++;
	
	//for(k = 0; k < NOobj; k++) 
      }
      
      // plot object pixel centres
      cpgsci(2);
      cpgpt(plot_points,plot_x,plot_y,5);
      
      // plot axes
      cpgsci(1);
      cpgbox("BCNST",0,0,"BCNST",0,0);
      
      if(NOf > 1){

	// adjust size of plot_array if required
	if(NOy != NOf){
	
	  std::cout << "Dec and frequency dimensions are different. Adjusting size of plot_array pointer . . . " << std::endl;
	  
	  // free up memory used for plotting array
	  delete [] plot_array;
	  delete [] plot_x;
	  delete [] plot_y;
	  
	  // create plotting array
	  plot_array = new float[(NOx * NOf)];
	  plot_x = new float[(NOx * NOf)];
	  plot_y = new float[(NOx * NOf)];
	  
	}
	
	// de-bugging, create array of RA Position-Velocity map values
	std::cout << "Creating RA Position-Velocity map of total intensity . . . " << std::endl;
	plot_max = CreateRAPVPlot(plot_array,NOobj,detections,NOx,NOf,obj_limit);
	
	// re-scale plot_max
	plot_max = 0.3 * plot_max;
	
	// check if plot_max is 0.0
	if(plot_max == 0.0){ plot_max = 1.0; }
	
	// create RA position velocity intensity map
	std::cout << "Plotting RA position-velocity map of total intensity . . . " << std::endl;
	//cpgenv(-2,(NOx + 1),-2,(NOf + 1),0,0);
	cpgpage();
	cpgvstd();
	cpgqvp(1,&vp_x_min,&vp_x_max,&vp_y_min,&vp_y_max);
	vp_x_min = ((float) NOx) / (vp_x_max - vp_x_min);
	vp_y_min = ((float) NOf) / (vp_y_max - vp_y_min);
	cpgswin((-1.0 - (3.0 * vp_x_min / 200.0)),(NOx + (3.0 * vp_x_min / 200.0)),(-1.0 - (3.0 * vp_y_min / 200.0)),(NOf + (3.0 * vp_y_min / 200.0)));
	cpglab("RA","Velocity","RA Position-Velocity intensity map");
	cpgslw(2);
	cpggray(plot_array,NOx,NOf,1,NOx,1,NOf,plot_max,0.0,tr);
	cpgwedg("R",0.0,3.0,plot_max,0.0,"Integrated Intensity (mJy/beam)");
	cpgslw(3);
	
	// plot detection boundaries
	std::cout << "Creating boundaries of objects . . . " << std::endl;
	cpgsci(3);
	cpgslw(1.5);
	for(k = 0; k < NOobj; k++){
	  
	  obj_batch = floor(((float) k / (float) obj_limit));
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	  plot_points = CreateRAPVBounds(detections,(detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin() + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin() + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmin() + 1),detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin(),detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmin(),plot_x,plot_y,k,obj_limit);
	  
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() < 0.0){ cpgsci(4); }
	  cpgline(plot_points,plot_x,plot_y);
	  cpgsci(3);
	  
	}
	cpgslw(3);
	
	// create list of object pixel centres
	plot_points = 0;
	for(k = 0; k < NOobj; k++){
	  
	  obj_batch = floorf(((float) k / (float) obj_limit));   
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	  plot_x[plot_points] = detections[obj_batch][(k - (obj_batch * obj_limit))].GetRA();
	  plot_y[plot_points] = detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQ();
	  plot_points++;
	  
	  //for(k = 0; k < NOobj; k++) 
	}
	
	// plot object pixel centres
	cpgsci(2);
	cpgpt(plot_points,plot_x,plot_y,5);
	
	// plot axes
	cpgsci(1);
	cpgbox("BCNST",0,0,"BCNST",0,0);
	
	// adjust size of plot_array if required
	if(NOy != NOx){
	  
	  std::cout << "RA and Dec dimensions are different. Adjusting size of plot_array pointer . . . " << std::endl;
	  
	  // free up memory used for plotting array
	  delete [] plot_array;
	  delete [] plot_x;
	  delete [] plot_y;
	  
	  // create plotting array
	  plot_array = new float[(NOy * NOf)];
	  plot_x = new float[(NOy * NOf)];
	  plot_y = new float[(NOy * NOf)];
	  
	}
	
	// de-bugging, create array of Dec Position-Velocity map values
	std::cout << "Creating Dec Position-Velocity map of total intensity . . . " << std::endl;
	plot_max = CreateDecPVPlot(plot_array,NOobj,detections,NOy,NOf,obj_limit);
	
	// re-scale plot_max
	plot_max = 0.3 * plot_max;
	
	// check if plot_max is 0.0
	if(plot_max == 0.0){ plot_max = 1.0; }
	
	// create DEC position velocity intensity map
	std::cout << "Plotting Dec position-velocity map of total intensity . . . " << std::endl;
	//cpgenv(-2,(NOy + 1),-2,(NOf + 1),0,0);
	cpgpage();
	cpgvstd();
	cpgqvp(1,&vp_x_min,&vp_x_max,&vp_y_min,&vp_y_max);
	vp_x_min = ((float) NOy) / (vp_x_max - vp_x_min);
	vp_y_min = ((float) NOf) / (vp_y_max - vp_y_min);
	cpgswin((-1.0 - (3.0 * vp_x_min / 200.0)),(NOy + (3.0 * vp_x_min / 200.0)),(-1.0 - (3.0 * vp_y_min / 200.0)),(NOf + (3.0 * vp_y_min / 200.0)));
	cpglab("Dec","Velocity","Dec Position-Velocity intensity map");
	cpgslw(2);
	cpggray(plot_array,NOy,NOf,1,NOy,1,NOf,plot_max,0.0,tr);
	cpgwedg("R",0.0,3.0,plot_max,0.0,"Integrated Intensity (mJy/beam)");
	cpgslw(3);
	
	// plot detection boundaries
	std::cout << "Creating boundaries of objects . . . " << std::endl;
	cpgsci(3);
	cpgslw(1.5);
	for(k = 0; k < NOobj; k++){
	  
	  obj_batch = floor(((float) k / (float) obj_limit));
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	  plot_points = CreateDecPVBounds(detections,(detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin() + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin() + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmin() + 1),detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin(),detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmin(),plot_x,plot_y,k,obj_limit);
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() < 0.0){ cpgsci(4); }
	  cpgline(plot_points,plot_x,plot_y);
	  cpgsci(3);
	  
	}
	cpgslw(3);
	
	// create list of object pixel centres
	plot_points = 0;
	for(k = 0; k < NOobj; k++){
	  
	  obj_batch = floorf(((float) k / (float) obj_limit));
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	  plot_x[plot_points] = detections[obj_batch][(k - (obj_batch * obj_limit))].GetDEC();
	  plot_y[plot_points] = detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQ();
	  plot_points++;
	  
	  //for(k = 0; k < NOobj; k++) 
	}
	
	// plot object pixel centres
	cpgsci(2);
	cpgpt(plot_points,plot_x,plot_y,5);
	
	// plot axes
	cpgsci(1);
	cpgbox("BCNST",0,0,"BCNST",0,0);
	
	// if(NOf > 1)
      }
      
      // free up memory used for plotting array
      delete [] plot_array;
      delete [] plot_x;
      delete [] plot_y;
      
      // close plotting device
      cpgclos();
      
      std::cout << "Finished constructing moment-0 and position-velocity plots for entire datacube." << std::endl;
      
    }

    if((plot_mode == 2) || (plot_mode == 3)){
      
      // scale velocity field postage stamp images to km/s
      if(ctype3 >= 0){
	
	// z axis is already in velocity
	
	for(k = 0; ((k < NOobj) && (k < SPEC_PLOT_LIMIT)); k++){
	
	  // calculate obj_batch number of object
	  obj_batch = floorf(((float) k / (float) obj_limit));
	
	  // skip if this object has been re-initialised
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	 
	  for(g = 0; g < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); g++){

	    detections[obj_batch][(k - (obj_batch * obj_limit))].Set_vfield(g,(1.0 * (fabs(cdelt3)) * detections[obj_batch][(k - (obj_batch * obj_limit))].Get_vfield(g)));

	  }

	}

      } else {

	// z axis is in frequency

	for(k = 0; ((k < NOobj) && (k < SPEC_PLOT_LIMIT)); k++){
	
	  // calculate obj_batch number of object
	  obj_batch = floorf(((float) k / (float) obj_limit));
	  
	  // skip if this object has been re-initialised
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	  
	  for(g = 0; g < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); g++){
	    
	    detections[obj_batch][(k - (obj_batch * obj_limit))].Set_vfield(g,(299792.458 * restfreq * ((1.0/(crval3 + ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_vfield(g) - crpix3) * cdelt3))) - (1.0/(crval3 + ((detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQi() - crpix3) * cdelt3))))));

	  }

	}

      }

      // create plots of integrated spectra
      std::cout << "Plotting integrated spectra of individual objects . . . " << std::endl;
      
      // create plot arrays
      if(NOf > 20){

	plot_x = new float[(2 * NOf)];
	plot_y = new float[(2 * NOf)];

      } else {

	plot_x = new float[21];
	plot_y = new float[21];

      }

      // create plotting file for integrated spectra
      dummy1 = "";
      dummy1.clear();
      dummy1 = output_code+"_spectra.ps/vcps";
      std::cout << "Creating output file: " << dummy1 << std::endl;
      cpgopen(dummy1.c_str());
      cpgslw(3);
      
      progress = 0.0;
      std::cout << "0 | |:| | : | |:| | 100% complete" << std::endl;

      // plotting 10 integrated spectra on each page along with postage stamps
      i = -1;
      for(k = 0; ((k < NOobj) && (k < SPEC_PLOT_LIMIT)); k++){
	
	// calculate obj_batch number of object
	obj_batch = floorf(((float) k / (float) obj_limit));
	
	// skip if this object has been re-initialised
	if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	
	// create plotting array for object spectrum
	plot_array = new float[(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1)];

	// increment i - tracks the number of objects plotted on the current page
	i++;
	
	// move to a new page if this is a new block of 10 spectra to be plotted
	if((i % 10) == 0){ cpgpage(); }
	
	// determine range of spectra
	plot_min = 1E10;
	plot_max = -1E10;
	if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1) >= 10){
	  
	  for(g = 0; g < (2 * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); g++){
	    
	    if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g) >= plot_max){ plot_max = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g); }
	    if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g) <= plot_min){ plot_min = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g); }
	    
	  }
	  
	} else {
	  
	  for(g = 0; g < (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 11); g++){
	    
	    if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g) >= plot_max){ plot_max = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g); }
	    if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g) <= plot_min){ plot_min = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g); }
	    
	  }
	  
	}
	for(g = 0; g < (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1); g++){
	  
	  plot_array[g] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_ospec(g);
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_ospec(g) >= plot_max){ plot_max = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_ospec(g); }
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_ospec(g) <= plot_min){ plot_min = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_ospec(g); }
	  
	}
	
	// expand plot scale
	if(plot_min == plot_max){

	  plot_min-=0.025;
	  plot_max+=0.025;

	} else {

	  plot_min = plot_min - 0.025 * (plot_max - plot_min);
	  plot_max = plot_max + 0.025 * (plot_max - plot_min);
	
	}

	// set viewport for object
	g = (int) floorf(((float) i / 10.0));
	cpgsvp(0.05,0.46,(1.0 - 0.1 * ((float) i + 1 - ((float) g * 10.0)) + 0.02),(1.0 - 0.1 * ((float) i - ((float) g * 10.0)) - 0.01));
	
	// set window (taking into account datacube size) and plot spectra for object (object + reference spectra)
	if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1) >= 10){
	  
	  if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - ((int) floorf((0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))))) > 0){
	    
	    if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + ((int) floorf((0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))))) < NOf){
	      
	      cpgswin((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - ((int) floorf((0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))))),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + ((int) floorf((0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))))),plot_min,plot_max);
	      
	    } else {
	      
	      cpgswin((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - ((int) floorf((0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))))),NOf,plot_min,plot_max);
	      
	    }
	    
	  } else {
	    
	    if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + ((int) floorf((0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))))) < NOf){
	      
	      cpgswin(0,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + ((int) floorf((0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))))),plot_min,plot_max);
	      
	    } else {
	      
	      cpgswin(0,NOf,plot_min,plot_max);
	      
	    }
	    
	  }
	  
	  for(g = 0; g < (2 * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); g++){ 

	    plot_x[g] = g + (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - ((int) floorf((0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))))); 
	    plot_y[g] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g);

	  }
	  
	  cpgline((2 * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1)),plot_x,plot_y);
	  
	} else {
	  
	  if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - 5) > 0){
	    
	    if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + 5) < NOf){
	      
	      cpgswin((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - 5),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + 5),plot_min,plot_max);
	      
	    } else {
	      
	      cpgswin((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - 5),NOf,plot_min,plot_max);
	      
	    }
	    
	  } else {
	    
	    if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + 5) < NOf){
	      
	      cpgswin(0,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + 5),plot_min,plot_max);
	      
	    } else {
	      
	      cpgswin(0,NOf,plot_min,plot_max);
	      
	    }
	    
	  }
	  
	  for(g = 0; g < (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 11); g++){ 

	    plot_x[g] = g + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - 5;
	    plot_y[g] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g);

	  }
	  
	  cpgline((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 11),plot_x,plot_y);
	  
	}
	
	// plot positions of W_50 and W_20 estimates
	
	// a. plot conventional width estimates
	cpgsci(4);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w_max(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w_max(),plot_max);
	cpgsls(2);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w50_min(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w50_min(),plot_max);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w50_max(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w50_max(),plot_max);
	cpgsls(4);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w20_min(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w20_min(),plot_max);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w20_max(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w20_max(),plot_max);
	cpgsls(1);
	
	// b. plot c.d.f. based width estimates
	cpgsci(3);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw_max(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw_max(),plot_max);
	cpgsls(2);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw50_min(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw50_min(),plot_max);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw50_max(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw50_max(),plot_max);
	cpgsls(4);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw20_min(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw20_min(),plot_max);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw20_max(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw20_max(),plot_max);
	cpgsls(1);
	cpgsci(1);

	// c. plot frequency of pixel central moment, both unweighted and intensity weighted
	cpgsci(2);
	cpgsls(2);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQi(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQi(),plot_max);
	cpgsls(4);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQ(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQ(),plot_max);
	cpgsls(1);
	cpgsci(1);

	for(g = 0; g < (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1); g++){ 

	  plot_x[g] = g + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4); 
	  plot_array[g] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_ospec(g);
	  
	}
	cpgsci(2);
	if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1) > 1){
	  
	  cpgline((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1),plot_x,plot_array);
	  
	} else {
	  
	  cpgpt(1,plot_x,plot_array,-4);
	  
	}
	cpgsci(1);
	
	// add labels
	cpgsch(0.6);
	cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	cpgsch(0.7);
	dummy2.str("");
	dummy2.clear();
	dummy2 << (i + 1);
	dummy1 = "";
	dummy2 >> dummy1;
	cpglab("",dummy1.c_str(),"");
	cpgsch(1.0);
		
	// d. free plotting arrays
	delete [] plot_array;

	// plot postage stamp moment-0 map

	// 0. create plotting array
	plot_array = new float[((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1)*(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1))];
	
	// a. set viewport for postage stamp
	g = (int) floorf(((float) i / 10.0));
	cpgsvp(0.49,0.59,(1.0 - 0.1 * ((float) i + 1 - ((float) g * 10.0)) + 0.02),(1.0 - 0.1 * ((float) i - ((float) g * 10.0)) - 0.01));
	
	// b. set window for postage stamp
	if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0)) >= (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2))){

	  plot_max = 0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2));  
	  cpgswin(-1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),(-1.0 - plot_max),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1 - plot_max));
	  
	} else {
	  
	  plot_max = 0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0));
	  cpgswin((-1.0 - plot_max),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1.0 - plot_max),-1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1));	  
	    
	}
	
	// c. adjust plot_max
	plot_max = -1E10;
	plot_min = 1E10;
	for(g = 0; g < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); g++){ 
	  
	  plot_array[g] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_mom0(g);
	  if(plot_array[g] >= plot_max){ plot_max = plot_array[g]; } 
	  if(plot_array[g] <= plot_min){ plot_min = plot_array[g]; } 
	  
	}
	if(plot_min == plot_max){ 

	  plot_max+=1.0;
	  plot_min-=1.0; 

	}
	if(detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() < 0.0){ plot_max = plot_min; }
	if(plot_max == 0.0){ plot_max = 1.0; }

	// d. plot postage stamp
	cpgslw(2);
	cpggray(plot_array,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1),1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1),plot_max,0.0,tr);
	cpgslw(3);
	
	// e. plot borders
	cpgsch(0.7);
	cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	cpgsch(1.0);
	
	// f. free array --> redundant
	//delete [] plot_array;

	// plot postage stamp velocity field 

	// 0. create plotting array --> redundant
	//plot_array = new float[((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1)*(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1))];
	
	// a. set viewport for postage stamp
	g = (int) floorf(((float) i / 10.0));
	cpgsvp(0.62,0.72,(1.0 - 0.1 * ((float) i + 1 - ((float) g * 10.0)) + 0.02),(1.0 - 0.1 * ((float) i - ((float) g * 10.0)) - 0.01));
	
	// b. set window for postage stamp
	if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0)) >= (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2))){

	  plot_max = 0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2));  
	  cpgswin(-1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),(-1.0 - plot_max),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1 - plot_max));
	  
	} else {
	  
	  plot_max = 0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0));
	  cpgswin((-1.0 - plot_max),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1.0 - plot_max),-1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1));	  
	    
	}
	
	// c. adjust plot_max
	j = 0;
	for(g = 0; g < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); g++){ 
	  
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_mom0(g) != 0.0){

	    plot_array[j] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_vfield(g);
	    j++;

	  }

	}
	HeapSort(j,(plot_array - 1));
	g = ((int) floorf(((1.0/3.0) * ((float) j))));
	if(g >= 0){ plot_min = plot_array[g]; } else { plot_min = plot_array[0]; }
	g = 1 + ((int) floorf(((2.0/3.0) * ((float) j))));
	if(g < j){ plot_max = plot_array[g]; } else { plot_max = plot_array[(j - 1)]; }

	for(g = 0; g < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); g++){ 
	  
	  plot_array[g] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_vfield(g);

	}
	if((fabs(plot_min)) > (fabs(plot_max))){

	  plot_max = -1.0 * plot_min;

	} else if((fabs(plot_min)) < (fabs(plot_max))){

	  plot_min = -1.0 * plot_max;

	}
	if(plot_min == plot_max){ 

	  plot_max+=1.0;
	  plot_min-=1.0; 

	}
	if(plot_min == plot_max){

	  plot_max = 1.0;
	  plot_min = -1.0;

	}

	// d. plot postage stamp
	cpgslw(2);
	cpggray(plot_array,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1),1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1),plot_max,plot_min,tr);
	cpgslw(3);
	
	// e. plot borders
	cpgsch(0.7);
	cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	cpgsch(1.0);
	
	// f. free array
	delete [] plot_array;

	// plot postage stamp RA PV diagram

	// 0. create plotting array
	plot_array = new float[((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1)*(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))];
		
	// a. set viewport for postage stamp
	g = (int) floorf(((float) i / 10.0));
	cpgsvp(0.75,0.85,(1.0 - 0.1 * ((float) i + 1 - ((float) g * 10.0)) + 0.02),(1.0 - 0.1 * ((float) i - ((float) g * 10.0)) - 0.01));
	
	// b. set window for postage stamp
	if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0)) >= (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4))){

	  plot_max = 0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4));
	  cpgswin(-1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),(-1.0 - plot_max),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1 - plot_max));
	  
	} else {
	  
	  plot_max = 0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0));
	  cpgswin((-1.0 - plot_max),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1.0 - plot_max),-1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1));	  
	    
	}
	
	// c. adjust plot_max
	plot_max = -1E10;
	plot_min = 1E10;
	for(g = 0; g < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); g++){ 
	  
	  plot_array[g] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_RAPV(g);
	  if(plot_array[g] >= plot_max){ plot_max = plot_array[g]; } 
	  if(plot_array[g] <= plot_min){ plot_min = plot_array[g]; } 
	  
	}
	if(plot_min == plot_max){ 

	  plot_max+=1.0;
	  plot_min-=1.0; 

	}
	if(detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() < 0.0){ plot_max = plot_min; }
	if(plot_max == 0.0){ plot_max = 1.0; }

	// d. plot postage stamp
	cpgslw(2);
	cpggray(plot_array,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1),1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1),plot_max,0.0,tr);
	cpgslw(3);
	
	// e. plot borders
	cpgsch(0.7);
	cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	cpgsch(1.0);
	
	// f. free plotting array
	delete [] plot_array;

	// plot postage stamp Dec PV diagram

	// 0. create plotting array
	plot_array = new float[((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1)*(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))];
		
	// a. set viewport for postage stamp
	g = (int) floorf(((float) i / 10.0));
	cpgsvp(0.88,0.98,(1.0 - 0.1 * ((float) i + 1 - ((float) g * 10.0)) + 0.02),(1.0 - 0.1 * ((float) i - ((float) g * 10.0)) - 0.01));
	
	// b. set window for postage stamp
	if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2)) >= (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4))){

	  plot_max = 0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4));
	  cpgswin(-1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1),(-1.0 - plot_max),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1 - plot_max));
	  
	} else {
	  
	  plot_max = 0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2));
	  cpgswin((-1.0 - plot_max),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1.0 - plot_max),-1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1));	  
	    
	}
	
	// c. adjust plot_max
	plot_max = -1E10;
	plot_min = 1E10;
	for(g = 0; g < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); g++){ 
	  
	  plot_array[g] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_DECPV(g);
	  if(plot_array[g] >= plot_max){ plot_max = plot_array[g]; } 
	  if(plot_array[g] <= plot_min){ plot_min = plot_array[g]; } 
	  
	}
	if(plot_min == plot_max){ 

	  plot_max+=1.0;
	  plot_min-=1.0; 

	}
	if(detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() < 0.0){ plot_max = plot_min; }
	if(plot_max == 0.0){ plot_max = 1.0; }

	// d. plot postage stamp
	cpgslw(2);
	cpggray(plot_array,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1),1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1),1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1),plot_max,0.0,tr);
	cpgslw(3);
	
	// e. plot borders
	cpgsch(0.7);
	cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	cpgsch(1.0);

	// f. free plotting array
	delete [] plot_array;

	// display progress
	while(progress <= (((float) (k + 1)) / ((float) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
	
	//for(k = 0; k < NOobj; k++) 
      }  

      std::cout << "* done." << std::endl;

      // free up remaining plotting arrays
      delete [] plot_x;
      delete [] plot_y;

      // close plotting device
      cpgclos();
            
    }

  }

}

void CreateObjPlots(int m, int plot_mode, std::string output_code, int * xyz_order, int NOx, int NOy, int NOf, vector<object_props_dbl *> & detections, long int NOobj, int obj_limit, int ctype3, float crpix3, float crval3, float cdelt3, float restfreq){
  
  long k, obj, obj_batch;
  int g, i, j, plot_points, temp_x, temp_y, temp_z;
  std::string dummy1;
  std::stringstream dummy2;
  float plot_min, plot_max, * plot_x, * plot_y, * plot_array, tr[6] = {-1.0,1.0,0.0,-1.0,0.0,1.0};
  float progress, vp_x_min, vp_x_max, vp_y_min, vp_y_max;

  std::cout << "Input limits are: " << NOx << " " << NOy << " " << NOf << std::endl;

  // reorder the datacube and subcube limits to be in x,y,z order
  temp_x = NOx;
  temp_y = NOy;
  temp_z = NOf;
  switch(xyz_order[0]){
  case 1:
    NOx = temp_x;
    break;
  case 2:
    NOx = temp_y;
    break;
  case 3:
    NOx = temp_z;
    break;
  default:
    NOx = temp_x;
    break;
  }
  switch(xyz_order[1]){
  case 1:
    NOy = temp_x;
    break;
  case 2:
    NOy = temp_y;
    break;
  case 3:
    NOy = temp_z;
    break;
  default:
    NOy = temp_y;
    break;
  }
  switch(xyz_order[2]){
  case 1:
    NOf = temp_x;
    break;
  case 2:
    NOf = temp_y;
    break;
  case 3:
    NOf = temp_z;
    break;
   default:
    NOf = temp_z;
    break;
  }

  std::cout << "Adjusting limits to: " << NOx << " " << NOy << " " << NOf << " because order is: " << xyz_order[0] << " " << xyz_order[1] << " " << xyz_order[2] << std::endl;

  // generate plots
  if((m > 1) && (plot_mode > 0)){

    if((plot_mode == 1) || (plot_mode == 3)){
      
      std::cout << "\nCreating global moment-0 and position-velocity (RA & Dec) plots . . . " << std::endl;

      // open plotting environment
      dummy1 = "";
      dummy1.clear();
      dummy1 = output_code+"_plots.ps/cps";
      std::cout << "Creating output file: " << dummy1 << std::endl;
      cpgopen(dummy1.c_str());
      cpgslw(3);
      
      // create boundary arrays
      plot_x = new float[(NOx * NOy)];
      plot_y = new float[(NOx * NOy)];
      
      // create plotting array
      plot_array = new float[(NOx * NOy)];
      
      // create array of moment 0 map values
      std::cout << "Creating moment 0 map of total intensity . . . " << std::endl;
      plot_max = CreateMoment0Map(plot_array,NOobj,detections,NOx,NOy,obj_limit);
      
      // re-scale plot_max
      plot_max = 0.3 * plot_max;

      // check if plot_max is 0.0
      if(plot_max == 0.0){ plot_max = 1.0; }

      // create plot of moment 0 map for sources with boundaries over-plotted
      std::cout << "Plotting moment 0 map of total intensity . . . " << std::endl;
      //cpgenv(-1,NOx,-1,NOy,0,0);
      cpgvstd();
      cpgqvp(1,&vp_x_min,&vp_x_max,&vp_y_min,&vp_y_max);
      vp_x_min = ((float) NOx) / (vp_x_max - vp_x_min);
      vp_y_min = ((float) NOy) / (vp_y_max - vp_y_min);
      cpgswin((-1.0 - (3.0 * vp_x_min / 200.0)),(NOx + (3.0 * vp_x_min / 200.0)),(-1.0 - (3.0 * vp_y_min / 200.0)),(NOy + (3.0 * vp_y_min / 200.0)));
      cpglab("RA","Dec","moment 0 map");
      cpgslw(2);
      cpggray(plot_array,NOx,NOy,1,NOx,1,NOy,plot_max,0.0,tr);
      cpgwedg("R",0.0,3.0,plot_max,0.0,"Integrated Intensity (mJy/beam)");
      cpgslw(3);
      
      // plot detection boundaries
      std::cout << "Creating boundaries of objects . . . " << std::endl;
      cpgsci(3);
      cpgslw(1.5);
      for(k = 0; k < NOobj; k++){
	
	obj_batch = (long int) floor(((double) k / (double) obj_limit));
	if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	plot_points = CreateMoment0Bounds(detections,(detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin() + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin() + 1),detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin(),detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin(),plot_x,plot_y,k,obj_limit);
	if(detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() < 0.0){ cpgsci(4); }
	cpgline(plot_points,plot_x,plot_y);
	cpgsci(3);
	
      } 
      cpgslw(3);
      
      // create list of object pixel centres
      plot_points = 0;
      for(k = 0; k < NOobj; k++){
	
	obj_batch = (long int) floor(((double) k / (double) obj_limit));
	if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	plot_x[plot_points] = detections[obj_batch][(k - (obj_batch * obj_limit))].GetRA();
	plot_y[plot_points] = detections[obj_batch][(k - (obj_batch * obj_limit))].GetDEC();
	plot_points++;
	
	//for(k = 0; k < NOobj; k++) 
      }
      
      // plot object pixel centres
      cpgsci(2);
      cpgpt(plot_points,plot_x,plot_y,5);
      
      // plot axes
      cpgsci(1);
      cpgbox("BCNST",0,0,"BCNST",0,0);
      
      if(NOf > 1){

	// adjust size of plot_array if required
	if(NOy != NOf){
	
	  std::cout << "Dec and frequency dimensions are different. Adjusting size of plot_array pointer . . . " << std::endl;
	  
	  // free up memory used for plotting array
	  delete [] plot_array;
	  delete [] plot_x;
	  delete [] plot_y;
	  
	  // create plotting array
	  plot_array = new float[(NOx * NOf)];
	  plot_x = new float[(NOx * NOf)];
	  plot_y = new float[(NOx * NOf)];
	  
	}
	
	// de-bugging, create array of RA Position-Velocity map values
	std::cout << "Creating RA Position-Velocity map of total intensity . . . " << std::endl;
	plot_max = CreateRAPVPlot(plot_array,NOobj,detections,NOx,NOf,obj_limit);
	
	// re-scale plot_max
	plot_max = 0.3 * plot_max;
	
	// check if plot_max is 0.0
	if(plot_max == 0.0){ plot_max = 1.0; }
	
	// create RA position velocity intensity map
	std::cout << "Plotting RA position-velocity map of total intensity . . . " << std::endl;
	//cpgenv(-2,(NOx + 1),-2,(NOf + 1),0,0);
	cpgpage();
	cpgvstd();
	cpgqvp(1,&vp_x_min,&vp_x_max,&vp_y_min,&vp_y_max);
	vp_x_min = ((float) NOx) / (vp_x_max - vp_x_min);
	vp_y_min = ((float) NOf) / (vp_y_max - vp_y_min);
	cpgswin((-1.0 - (3.0 * vp_x_min / 200.0)),(NOx + (3.0 * vp_x_min / 200.0)),(-1.0 - (3.0 * vp_y_min / 200.0)),(NOf + (3.0 * vp_y_min / 200.0)));
	cpglab("RA","Velocity","RA Position-Velocity intensity map");
	cpgslw(2);
	cpggray(plot_array,NOx,NOf,1,NOx,1,NOf,plot_max,0.0,tr);
	cpgwedg("R",0.0,3.0,plot_max,0.0,"Integrated Intensity (mJy/beam)");
	cpgslw(3);
	
	// plot detection boundaries
	std::cout << "Creating boundaries of objects . . . " << std::endl;
	cpgsci(3);
	cpgslw(1.5);
	for(k = 0; k < NOobj; k++){
	  
	  obj_batch = (long int) floor(((double) k / (double) obj_limit));
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	  plot_points = CreateRAPVBounds(detections,(detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin() + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin() + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmin() + 1),detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin(),detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmin(),plot_x,plot_y,k,obj_limit);
	  
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() < 0.0){ cpgsci(4); }
	  cpgline(plot_points,plot_x,plot_y);
	  cpgsci(3);
	  
	}
	cpgslw(3);
	
	// create list of object pixel centres
	plot_points = 0;
	for(k = 0; k < NOobj; k++){
	  
	  obj_batch = (long int) floor(((double) k / (double) obj_limit));   
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	  plot_x[plot_points] = detections[obj_batch][(k - (obj_batch * obj_limit))].GetRA();
	  plot_y[plot_points] = detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQ();
	  plot_points++;
	  
	  //for(k = 0; k < NOobj; k++) 
	}
	
	// plot object pixel centres
	cpgsci(2);
	cpgpt(plot_points,plot_x,plot_y,5);
	
	// plot axes
	cpgsci(1);
	cpgbox("BCNST",0,0,"BCNST",0,0);
	
	// adjust size of plot_array if required
	if(NOy != NOx){
	  
	  std::cout << "RA and Dec dimensions are different. Adjusting size of plot_array pointer . . . " << std::endl;
	  
	  // free up memory used for plotting array
	  delete [] plot_array;
	  delete [] plot_x;
	  delete [] plot_y;
	  
	  // create plotting array
	  plot_array = new float[(NOy * NOf)];
	  plot_x = new float[(NOy * NOf)];
	  plot_y = new float[(NOy * NOf)];
	  
	}
	
	// de-bugging, create array of Dec Position-Velocity map values
	std::cout << "Creating Dec Position-Velocity map of total intensity . . . " << std::endl;
	plot_max = CreateDecPVPlot(plot_array,NOobj,detections,NOy,NOf,obj_limit);
	
	// re-scale plot_max
	plot_max = 0.3 * plot_max;
	
	// check if plot_max is 0.0
	if(plot_max == 0.0){ plot_max = 1.0; }
	
	// create DEC position velocity intensity map
	std::cout << "Plotting Dec position-velocity map of total intensity . . . " << std::endl;
	//cpgenv(-2,(NOy + 1),-2,(NOf + 1),0,0);
	cpgpage();
	cpgvstd();
	cpgqvp(1,&vp_x_min,&vp_x_max,&vp_y_min,&vp_y_max);
	vp_x_min = ((float) NOy) / (vp_x_max - vp_x_min);
	vp_y_min = ((float) NOf) / (vp_y_max - vp_y_min);
	cpgswin((-1.0 - (3.0 * vp_x_min / 200.0)),(NOy + (3.0 * vp_x_min / 200.0)),(-1.0 - (3.0 * vp_y_min / 200.0)),(NOf + (3.0 * vp_y_min / 200.0)));
	cpglab("Dec","Velocity","Dec Position-Velocity intensity map");
	cpgslw(2);
	cpggray(plot_array,NOy,NOf,1,NOy,1,NOf,plot_max,0.0,tr);
	cpgwedg("R",0.0,3.0,plot_max,0.0,"Integrated Intensity (mJy/beam)");
	cpgslw(3);
	
	// plot detection boundaries
	std::cout << "Creating boundaries of objects . . . " << std::endl;
	cpgsci(3);
	cpgslw(1.5);
	for(k = 0; k < NOobj; k++){
	  
	  obj_batch = (long int) floor(((double) k / (double) obj_limit));
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	  plot_points = CreateDecPVBounds(detections,(detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetRAmin() + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin() + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmax() - detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmin() + 1),detections[obj_batch][(k - (obj_batch * obj_limit))].GetDECmin(),detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQmin(),plot_x,plot_y,k,obj_limit);
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() < 0.0){ cpgsci(4); }
	  cpgline(plot_points,plot_x,plot_y);
	  cpgsci(3);
	  
	}
	cpgslw(3);
	
	// create list of object pixel centres
	plot_points = 0;
	for(k = 0; k < NOobj; k++){
	  
	  obj_batch = (long int) floor(((double) k / (double) obj_limit));
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	  plot_x[plot_points] = detections[obj_batch][(k - (obj_batch * obj_limit))].GetDEC();
	  plot_y[plot_points] = detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQ();
	  plot_points++;
	  
	  //for(k = 0; k < NOobj; k++) 
	}
	
	// plot object pixel centres
	cpgsci(2);
	cpgpt(plot_points,plot_x,plot_y,5);
	
	// plot axes
	cpgsci(1);
	cpgbox("BCNST",0,0,"BCNST",0,0);
	
	// if(NOf > 1)
      }
      
      // free up memory used for plotting array
      delete [] plot_array;
      delete [] plot_x;
      delete [] plot_y;
      
      // close plotting device
      cpgclos();
      
      std::cout << "Finished constructing moment-0 and position-velocity plots for entire datacube." << std::endl;
      
    }

    if((plot_mode == 2) || (plot_mode == 3)){
      
      // scale velocity field postage stamp images to km/s
      if(ctype3 >= 0){
	
	// z axis is already in velocity
	
	for(k = 0; ((k < NOobj) && (k < SPEC_PLOT_LIMIT)); k++){
	
	  // calculate obj_batch number of object
	  obj_batch = (long int) floor(((double) k / (double) obj_limit));
	
	  // skip if this object has been re-initialised
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	 
	  for(g = 0; g < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); g++){

	    detections[obj_batch][(k - (obj_batch * obj_limit))].Set_vfield(g,(1.0 * (fabs(cdelt3)) * detections[obj_batch][(k - (obj_batch * obj_limit))].Get_vfield(g)));

	  }

	}

      } else {

	// z axis is in frequency

	for(k = 0; ((k < NOobj) && (k < SPEC_PLOT_LIMIT)); k++){
	
	  // calculate obj_batch number of object
	  obj_batch = (long int) floor(((double) k / (double) obj_limit));
	  
	  // skip if this object has been re-initialised
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	  
	  for(g = 0; g < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); g++){
	    
	    detections[obj_batch][(k - (obj_batch * obj_limit))].Set_vfield(g,(299792.458 * restfreq * ((1.0/(crval3 + ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_vfield(g) - crpix3) * cdelt3))) - (1.0/(crval3 + ((detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQi() - crpix3) * cdelt3))))));

	  }

	}

      }

      // create plots of integrated spectra
      std::cout << "Plotting integrated spectra of individual objects . . . " << std::endl;
      
      // create plot arrays
      if(NOf > 20){

	plot_x = new float[(2 * NOf)];
	plot_y = new float[(2 * NOf)];

      } else {

	plot_x = new float[21];
	plot_y = new float[21];

      }

      // create plotting file for integrated spectra
      dummy1 = "";
      dummy1.clear();
      dummy1 = output_code+"_spectra.ps/vcps";
      std::cout << "Creating output file: " << dummy1 << std::endl;
      cpgopen(dummy1.c_str());
      cpgslw(3);
      
      progress = 0.0;
      std::cout << "0 | |:| | : | |:| | 100% complete" << std::endl;

      // plotting 10 integrated spectra on each page along with postage stamps
      i = -1;
      for(k = 0; ((k < NOobj) && (k < SPEC_PLOT_LIMIT)); k++){
	
	// calculate obj_batch number of object
	obj_batch = (long int) floor(((double) k / (double) obj_limit));
	
	// skip if this object has been re-initialised
	if(detections[obj_batch][(k - (obj_batch * obj_limit))].ShowVoxels() < 1){ continue; }
	
	// create plotting array for object spectrum
	plot_array = new float[(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1)];

	// increment i - tracks the number of objects plotted on the current page
	i++;
	
	// move to a new page if this is a new block of 10 spectra to be plotted
	if((i % 10) == 0){ cpgpage(); }
	
	// determine range of spectra
	plot_min = 1E10;
	plot_max = -1E10;
	if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1) >= 10){
	  
	  for(g = 0; g < (2 * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); g++){
	    
	    if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g) >= plot_max){ plot_max = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g); }
	    if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g) <= plot_min){ plot_min = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g); }
	    
	  }
	  
	} else {
	  
	  for(g = 0; g < (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 11); g++){
	    
	    if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g) >= plot_max){ plot_max = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g); }
	    if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g) <= plot_min){ plot_min = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g); }
	    
	  }
	  
	}
	for(g = 0; g < (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1); g++){
	  
	  plot_array[g] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_ospec(g);
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_ospec(g) >= plot_max){ plot_max = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_ospec(g); }
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_ospec(g) <= plot_min){ plot_min = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_ospec(g); }
	  
	}
	
	// expand plot scale
	if(plot_min == plot_max){

	  plot_min-=0.025;
	  plot_max+=0.025;

	} else {

	  plot_min = plot_min - 0.025 * (plot_max - plot_min);
	  plot_max = plot_max + 0.025 * (plot_max - plot_min);
	
	}

	// set viewport for object
	g = (int) floorf(((float) i / 10.0));
	cpgsvp(0.05,0.46,(1.0 - 0.1 * ((float) i + 1 - ((float) g * 10.0)) + 0.02),(1.0 - 0.1 * ((float) i - ((float) g * 10.0)) - 0.01));
	
	// set window (taking into account datacube size) and plot spectra for object (object + reference spectra)
	if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1) >= 10){
	  
	  if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - ((int) floorf((0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))))) > 0){
	    
	    if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + ((int) floorf((0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))))) < NOf){
	      
	      cpgswin((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - ((int) floorf((0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))))),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + ((int) floorf((0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))))),plot_min,plot_max);
	      
	    } else {
	      
	      cpgswin((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - ((int) floorf((0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))))),NOf,plot_min,plot_max);
	      
	    }
	    
	  } else {
	    
	    if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + ((int) floorf((0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))))) < NOf){
	      
	      cpgswin(0,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + ((int) floorf((0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))))),plot_min,plot_max);
	      
	    } else {
	      
	      cpgswin(0,NOf,plot_min,plot_max);
	      
	    }
	    
	  }
	  
	  for(g = 0; g < (2 * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); g++){ 

	    plot_x[g] = g + (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - ((int) floorf((0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))))); 
	    plot_y[g] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g);

	  }
	  
	  cpgline((2 * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1)),plot_x,plot_y);
	  
	} else {
	  
	  if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - 5) > 0){
	    
	    if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + 5) < NOf){
	      
	      cpgswin((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - 5),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + 5),plot_min,plot_max);
	      
	    } else {
	      
	      cpgswin((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - 5),NOf,plot_min,plot_max);
	      
	    }
	    
	  } else {
	    
	    if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + 5) < NOf){
	      
	      cpgswin(0,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + 5),plot_min,plot_max);
	      
	    } else {
	      
	      cpgswin(0,NOf,plot_min,plot_max);
	      
	    }
	    
	  }
	  
	  for(g = 0; g < (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 11); g++){ 

	    plot_x[g] = g + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - 5;
	    plot_y[g] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_rspec(g);

	  }
	  
	  cpgline((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 11),plot_x,plot_y);
	  
	}
	
	// plot positions of W_50 and W_20 estimates
	
	// a. plot conventional width estimates
	cpgsci(4);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w_max(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w_max(),plot_max);
	cpgsls(2);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w50_min(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w50_min(),plot_max);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w50_max(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w50_max(),plot_max);
	cpgsls(4);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w20_min(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w20_min(),plot_max);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w20_max(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_w20_max(),plot_max);
	cpgsls(1);
	
	// b. plot c.d.f. based width estimates
	cpgsci(3);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw_max(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw_max(),plot_max);
	cpgsls(2);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw50_min(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw50_min(),plot_max);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw50_max(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw50_max(),plot_max);
	cpgsls(4);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw20_min(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw20_min(),plot_max);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw20_max(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_cw20_max(),plot_max);
	cpgsls(1);
	cpgsci(1);

	// c. plot frequency of pixel central moment, both unweighted and intensity weighted
	cpgsci(2);
	cpgsls(2);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQi(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQi(),plot_max);
	cpgsls(4);
	cpgmove(detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQ(),plot_min);
	cpgdraw(detections[obj_batch][(k - (obj_batch * obj_limit))].GetFREQ(),plot_max);
	cpgsls(1);
	cpgsci(1);

	for(g = 0; g < (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1); g++){ 

	  plot_x[g] = g + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4); 
	  plot_array[g] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_ospec(g);
	  
	}
	cpgsci(2);
	if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1) > 1){
	  
	  cpgline((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1),plot_x,plot_array);
	  
	} else {
	  
	  cpgpt(1,plot_x,plot_array,-4);
	  
	}
	cpgsci(1);
	
	// add labels
	cpgsch(0.6);
	cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	cpgsch(0.7);
	dummy2.str("");
	dummy2.clear();
	dummy2 << (i + 1);
	dummy1 = "";
	dummy2 >> dummy1;
	cpglab("",dummy1.c_str(),"");
	cpgsch(1.0);
		
	// d. free plotting arrays
	delete [] plot_array;

	// plot postage stamp moment-0 map

	// 0. create plotting array
	plot_array = new float[((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1)*(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1))];
	
	// a. set viewport for postage stamp
	g = (int) floorf(((float) i / 10.0));
	cpgsvp(0.49,0.59,(1.0 - 0.1 * ((float) i + 1 - ((float) g * 10.0)) + 0.02),(1.0 - 0.1 * ((float) i - ((float) g * 10.0)) - 0.01));
	
	// b. set window for postage stamp
	if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0)) >= (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2))){

	  plot_max = 0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2));  
	  cpgswin(-1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),(-1.0 - plot_max),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1 - plot_max));
	  
	} else {
	  
	  plot_max = 0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0));
	  cpgswin((-1.0 - plot_max),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1.0 - plot_max),-1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1));	  
	    
	}
	
	// c. adjust plot_max
	plot_max = -1E10;
	plot_min = 1E10;
	for(g = 0; g < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); g++){ 
	  
	  plot_array[g] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_mom0(g);
	  if(plot_array[g] >= plot_max){ plot_max = plot_array[g]; } 
	  if(plot_array[g] <= plot_min){ plot_min = plot_array[g]; } 
	  
	}
	if(plot_min == plot_max){ 

	  plot_max+=1.0;
	  plot_min-=1.0; 

	}
	if(detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() < 0.0){ plot_max = plot_min; }
	if(plot_max == 0.0){ plot_max = 1.0; }

	// d. plot postage stamp
	cpgslw(2);
	cpggray(plot_array,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1),1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1),plot_max,0.0,tr);
	cpgslw(3);
	
	// e. plot borders
	cpgsch(0.7);
	cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	cpgsch(1.0);
	
	// f. free array --> redundant
	//delete [] plot_array;

	// plot postage stamp velocity field 

	// 0. create plotting array --> redundant
	//plot_array = new float[((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1)*(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1))];
	
	// a. set viewport for postage stamp
	g = (int) floorf(((float) i / 10.0));
	cpgsvp(0.62,0.72,(1.0 - 0.1 * ((float) i + 1 - ((float) g * 10.0)) + 0.02),(1.0 - 0.1 * ((float) i - ((float) g * 10.0)) - 0.01));
	
	// b. set window for postage stamp
	if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0)) >= (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2))){

	  plot_max = 0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2));  
	  cpgswin(-1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),(-1.0 - plot_max),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1 - plot_max));
	  
	} else {
	  
	  plot_max = 0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0));
	  cpgswin((-1.0 - plot_max),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1.0 - plot_max),-1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1));	  
	    
	}
	
	// c. adjust plot_max
	j = 0;
	for(g = 0; g < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); g++){ 
	  
	  if(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_mom0(g) != 0.0){

	    plot_array[j] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_vfield(g);
	    j++;

	  }

	}
	HeapSort(j,(plot_array - 1));
	g = ((int) floorf(((1.0/3.0) * ((float) j))));
	if(g >= 0){ plot_min = plot_array[g]; } else { plot_min = plot_array[0]; }
	g = 1 + ((int) floorf(((2.0/3.0) * ((float) j))));
	if(g < j){ plot_max = plot_array[g]; } else { plot_max = plot_array[(j - 1)]; }

	for(g = 0; g < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1)); g++){ 
	  
	  plot_array[g] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_vfield(g);

	}
	if((fabs(plot_min)) > (fabs(plot_max))){

	  plot_max = -1.0 * plot_min;

	} else if((fabs(plot_min)) < (fabs(plot_max))){

	  plot_min = -1.0 * plot_max;

	}
	if(plot_min == plot_max){ 

	  plot_max+=1.0;
	  plot_min-=1.0; 

	}
	if(plot_min == plot_max){

	  plot_max = 1.0;
	  plot_min = -1.0;

	}

	// d. plot postage stamp
	cpgslw(2);
	cpggray(plot_array,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1),1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1),plot_max,plot_min,tr);
	cpgslw(3);
	
	// e. plot borders
	cpgsch(0.7);
	cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	cpgsch(1.0);
	
	// f. free array
	delete [] plot_array;

	// plot postage stamp RA PV diagram

	// 0. create plotting array
	plot_array = new float[((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1)*(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))];
		
	// a. set viewport for postage stamp
	g = (int) floorf(((float) i / 10.0));
	cpgsvp(0.75,0.85,(1.0 - 0.1 * ((float) i + 1 - ((float) g * 10.0)) + 0.02),(1.0 - 0.1 * ((float) i - ((float) g * 10.0)) - 0.01));
	
	// b. set window for postage stamp
	if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0)) >= (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4))){

	  plot_max = 0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4));
	  cpgswin(-1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),(-1.0 - plot_max),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1 - plot_max));
	  
	} else {
	  
	  plot_max = 0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0));
	  cpgswin((-1.0 - plot_max),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1.0 - plot_max),-1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1));	  
	    
	}
	
	// c. adjust plot_max
	plot_max = -1E10;
	plot_min = 1E10;
	for(g = 0; g < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); g++){ 
	  
	  plot_array[g] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_RAPV(g);
	  if(plot_array[g] >= plot_max){ plot_max = plot_array[g]; } 
	  if(plot_array[g] <= plot_min){ plot_min = plot_array[g]; } 
	  
	}
	if(plot_min == plot_max){ 

	  plot_max+=1.0;
	  plot_min-=1.0; 

	}
	if(detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() < 0.0){ plot_max = plot_min; }
	if(plot_max == 0.0){ plot_max = 1.0; }

	// d. plot postage stamp
	cpgslw(2);
	cpggray(plot_array,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1),1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(1) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(0) + 1),1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1),plot_max,0.0,tr);
	cpgslw(3);
	
	// e. plot borders
	cpgsch(0.7);
	cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	cpgsch(1.0);
	
	// f. free plotting array
	delete [] plot_array;

	// plot postage stamp Dec PV diagram

	// 0. create plotting array
	plot_array = new float[((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1)*(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1))];
		
	// a. set viewport for postage stamp
	g = (int) floorf(((float) i / 10.0));
	cpgsvp(0.88,0.98,(1.0 - 0.1 * ((float) i + 1 - ((float) g * 10.0)) + 0.02),(1.0 - 0.1 * ((float) i - ((float) g * 10.0)) - 0.01));
	
	// b. set window for postage stamp
	if((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2)) >= (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4))){

	  plot_max = 0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4));
	  cpgswin(-1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1),(-1.0 - plot_max),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1 - plot_max));
	  
	} else {
	  
	  plot_max = 0.5 * (float) (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) + detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2));
	  cpgswin((-1.0 - plot_max),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1.0 - plot_max),-1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1));	  
	    
	}
	
	// c. adjust plot_max
	plot_max = -1E10;
	plot_min = 1E10;
	for(g = 0; g < ((detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1) * (detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1)); g++){ 
	  
	  plot_array[g] = detections[obj_batch][(k - (obj_batch * obj_limit))].Get_DECPV(g);
	  if(plot_array[g] >= plot_max){ plot_max = plot_array[g]; } 
	  if(plot_array[g] <= plot_min){ plot_min = plot_array[g]; } 
	  
	}
	if(plot_min == plot_max){ 

	  plot_max+=1.0;
	  plot_min-=1.0; 

	}
	if(detections[obj_batch][(k - (obj_batch * obj_limit))].GetTI() < 0.0){ plot_max = plot_min; }
	if(plot_max == 0.0){ plot_max = 1.0; }

	// d. plot postage stamp
	cpgslw(2);
	cpggray(plot_array,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1),(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1),1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(3) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(2) + 1),1,(detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(5) - detections[obj_batch][(k - (obj_batch * obj_limit))].Get_srep_size(4) + 1),plot_max,0.0,tr);
	cpgslw(3);
	
	// e. plot borders
	cpgsch(0.7);
	cpgbox("BCNST",0.0,0,"BCNST",0.0,0);
	cpgsch(1.0);

	// f. free plotting array
	delete [] plot_array;

	// display progress
	while(progress <= (((double) (k + 1)) / ((double) NOobj))){ std::cout << "*"; std::cout.flush(); progress+=0.05; }
	
	//for(k = 0; k < NOobj; k++) 
      }  

      std::cout << "* done." << std::endl;

      // free up remaining plotting arrays
      delete [] plot_x;
      delete [] plot_y;

      // close plotting device
      cpgclos();
            
    }

  }

}

