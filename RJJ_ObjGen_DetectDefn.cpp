#include<iostream>
#include<RJJ_ObjGen.h>

object_props::object_props(){ 
    
  NOvox = 0; 
  ra = dec = freq = ra_i = dec_i = freq_i = tot_intens = avg_intens = sigma_intens = rms = 0.0;
  ra_min = dec_min = freq_min = min_intens = 1E10;
  ra_max = dec_max = freq_max = max_intens = -1E10;
  w_max = w20_min = w50_min = w20_max = w50_max = -1E10;
  cw_max = cw20_min = cw50_min = cw20_max = cw50_max = -1E10;
  srep_size[0] = srep_size[1] = srep_size[2] = srep_size[3] = srep_size[4] = srep_size[5] = -99;
  srep_update = 0;
  srep_grid = NULL;
  srep_strings = NULL;
  mini_mom0 = NULL;
  mini_RAPV = NULL;
  mini_DECPV = NULL;
  mini_obj_spec = NULL;
  mini_ref_spec = NULL;
  mini_vfield = NULL;

  // new parameters --- multiple central moment calculations
  p_tot_intens = n_tot_intens = 0.0;
  p_ra_i = p_dec_i = p_freq_i = 0.0;
  n_ra_i = n_dec_i = n_freq_i = 0.0;

}

object_props::~object_props(){ 

  if(srep_grid != NULL){ delete [] srep_grid; srep_grid = NULL; } 
  if(srep_strings != NULL){ delete [] srep_strings; srep_strings = NULL; }
  if(mini_mom0 != NULL){ delete [] mini_mom0; mini_mom0 = NULL; }
  if(mini_RAPV != NULL){ delete [] mini_RAPV; mini_RAPV = NULL; }
  if(mini_DECPV != NULL){ delete [] mini_DECPV; mini_DECPV = NULL; }
  if(mini_obj_spec != NULL){ delete [] mini_obj_spec; mini_obj_spec = NULL; }
  if(mini_ref_spec != NULL){ delete [] mini_ref_spec; mini_ref_spec = NULL; }
  if(mini_vfield != NULL){ delete [] mini_vfield; mini_vfield = NULL; }

}

void object_props::ReInit(){

  NOvox = 0;
  ra = dec = freq = ra_i = dec_i = freq_i = tot_intens = avg_intens = sigma_intens = rms = 0.0;
  ra_min = dec_min = freq_min = min_intens = 1E10;
  ra_max = dec_max = freq_max = max_intens = -1E10;
  w_max = w20_min = w50_min = w20_max = w50_max = -1E10;
  cw_max = cw20_min = cw50_min = cw20_max = cw50_max = -1E10;

}

void object_props::AddVoxel(int value){ NOvox+=value; }

void object_props::AddRa(float value){ ra+=value; }

void object_props::AddDec(float value){ dec+=value; }

void object_props::AddFreq(float value){ freq+=value; }

void object_props::AddRA_i(float pos, float value){ 
  ra_i+=(pos * value);
  if(value >= 0.0){ p_ra_i+=(pos * value); } else { n_ra_i+=(pos * value); }
}

void object_props::AddDec_i(float pos, float value){ 
  dec_i+=value; 
  if(value >= 0.0){ p_dec_i+=(pos * value); } else { n_dec_i+=(pos * value); }
}

void object_props::AddFreq_i(float pos, float value){ 
  freq_i+=value; 
  if(value >= 0.0){ p_freq_i+=(pos * value); } else { n_freq_i+=(pos * value); }
}

void object_props::AddTotIntens(float value){ 
  tot_intens+=value; 
  if(value >= 0.0){ p_tot_intens+=value; } else { n_tot_intens+=value; }
}

void object_props::AddAvgIntens(float value){ avg_intens+=value; }

void object_props::AddSigmaItens(float value){ sigma_intens+=(value * value); }

void object_props::AdjustRange(float value){
    
  if(value <= min_intens){ min_intens = value; }
  if(value >= max_intens){ max_intens = value; }

}

void object_props::CalcProps(){
  
  float dummy, flip;
  int g;

  // calculate basic properties
  ra = ra / NOvox; 
  dec = dec / NOvox; 
  freq = freq / NOvox;
  ra_i = ra_i / tot_intens;
  dec_i = dec_i / tot_intens;
  freq_i = freq_i / tot_intens;
  avg_intens = avg_intens / NOvox;
  rms = sqrtf((sigma_intens / NOvox));
  sigma_intens = sqrtf(((sigma_intens / NOvox) - (avg_intens * avg_intens)));
  
  // new parameters --- multiple central moment calculations
  p_ra_i = p_ra_i / p_tot_intens;
  p_dec_i = p_dec_i / p_tot_intens;
  p_freq_i = p_freq_i / p_tot_intens;
  n_ra_i = n_ra_i / n_tot_intens;
  n_dec_i = n_dec_i / n_tot_intens;
  n_freq_i = n_freq_i / n_tot_intens;
  if(tot_intens >= 0.0){
    ra_i = p_ra_i;
    dec_i = p_dec_i;
    freq_i = p_freq_i;
  } else {
    ra_i = n_ra_i;
    dec_i = n_dec_i;
    freq_i = n_freq_i;
  }

  // remove inf's and nan's from the various arrays
  for(g = 0; g < ((srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1)); g++){
    if((std::isinf(mini_mom0[g])) || (std::isnan(mini_mom0[g]))){ mini_mom0[g] = 0.0; }
    if((std::isinf(mini_vfield[g])) || (std::isnan(mini_vfield[g]))){ mini_vfield[g] = 0.0; }
  }
  for(g = 0; g < ((srep_size[1] - srep_size[0] + 1) * (srep_size[5] - srep_size[4] + 1)); g++){
    if((std::isinf(mini_RAPV[g])) || (std::isnan(mini_RAPV[g]))){ mini_RAPV[g] = 0.0; }
  }
  for(g = 0; g < ((srep_size[3] - srep_size[2] + 1) * (srep_size[5] - srep_size[4] + 1)); g++){
    if((std::isinf(mini_DECPV[g])) || (std::isnan(mini_DECPV[g]))){ mini_DECPV[g] = 0.0; }
  }
  for(g = 0; g < ((srep_size[5] - srep_size[4] + 1)); g++){
    if((std::isinf(mini_obj_spec[g])) || (std::isnan(mini_obj_spec[g]))){ mini_obj_spec[g] = 0.0; }
  }  

  // normalise velocity field array
  for(g = 0; g < ((srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1)); g++){ 

    if(mini_mom0[g] != 0.0){

      mini_vfield[g] = (mini_vfield[g] / mini_mom0[g]) - freq_i; 

    } else {

      mini_vfield[g] = 0.0;

    }

  }
 
  // calculate widths of source
  
  // determine maximum flux in integrated spectrum and W_50, W_20 limits
  if((srep_size[5] - srep_size[4]) > 0){
    
    flip = 1.0;
    if(tot_intens < 0.0){ flip = -1.0; }
      
    // a. calculate W_50 and W_20 in conventional manner
    dummy = -1E10;
    for(g = 0; g < (srep_size[5] - srep_size[4] + 1); g++){ if((flip * mini_obj_spec[g]) >= dummy){ dummy = (flip * mini_obj_spec[g]); w_max = (float) (g + srep_size[4]); } }
    w20_min = (float) (srep_size[4]);
    for(g = 0; g < (srep_size[5] - srep_size[4] + 1); g++){ 
      
      if((flip * mini_obj_spec[g]) >= (0.2 * dummy)){ 
	
	if(g > 0){
	  
	  w20_min = (float) (g + srep_size[4]) - 1.0 + (((0.2 * dummy) - (flip * mini_obj_spec[(g - 1)])) / (flip * (mini_obj_spec[g] - mini_obj_spec[(g - 1)]))); 
	  
	} else {
	  
	  w20_min = (float) (g + srep_size[4]); 
	  
	}
	
	break; 
	
      } 
      
    }
    w20_max = (float) (srep_size[5]);
    for(g = srep_size[5] - srep_size[4]; g >= 0; g--){ 
      
      if((flip * mini_obj_spec[g]) >= (0.2 * dummy)){ 
	
	if(g < (srep_size[5] - srep_size[4])){
	  
	  w20_max = (float) (g + srep_size[4]) + 1.0 - (((0.2 * dummy) - (flip * mini_obj_spec[(g + 1)])) / (flip * (mini_obj_spec[g] - mini_obj_spec[(g + 1)])));
	  
	} else {
	  
	  w20_max = (float) (g + srep_size[4]); 
	  
	}
	
	break; 
	
      } 
	
    }
    w50_min = (float) (srep_size[4]);
    for(g = (int) floorf((w20_min - (float) srep_size[4])); g < (srep_size[5] - srep_size[4] + 1); g++){ 
      
      if((flip * mini_obj_spec[g]) >= (0.5 * dummy)){ 
	
	if(g > 0){
	  
	  w50_min = (float) (g + srep_size[4]) - 1.0 + (((0.5 * dummy) - (flip * mini_obj_spec[(g - 1)])) / (flip * (mini_obj_spec[g] - mini_obj_spec[(g - 1)]))); 
	  
	} else {
	  
	  w50_min = (float) (g + srep_size[4]); 
	  
	}	
	
	break; 
	
      } 
      
    }
    w50_max = (float) (srep_size[5]);
    for(g = (int) floorf((w20_max - (float) srep_size[4])); g >= 0; g--){ 
      
      if((flip * mini_obj_spec[g]) >= (0.5 * dummy)){ 
	
	if(g < (srep_size[5] - srep_size[4])){
	  
	  w50_max = (float) (g + srep_size[4]) + 1.0 - (((0.5 * dummy) - (flip * mini_obj_spec[(g + 1)])) / (flip * (mini_obj_spec[g] - mini_obj_spec[(g + 1)])));
	  
	} else {
	  
	  w50_max = (float) (g + srep_size[4]); 
	  
	}
	
	break; 
	
      } 
      
    }
    
    // b. calculate W_50 and W_20 according to c.f.d. values that correspond to FWHM and 1/5th of peak for a gaussian profile
    cw_max = cw20_min = cw50_min = dummy = 0.0;
    for(g = 0; g < (srep_size[5] - srep_size[4] + 1); g++){ 
      
      if(((dummy + (mini_obj_spec[g]/tot_intens)) >= 0.5) && (dummy < 0.5)){ 
	
	if(g > 0){
	  
	  cw_max = (float) (g + srep_size[4]) - 1.0 + ((0.5 - dummy) * tot_intens / mini_obj_spec[g]);
	  
	} else {
	  
	  cw_max = (float) (g + srep_size[4]);
	  
	}
	
      } 
      
      if(((dummy + (mini_obj_spec[g]/tot_intens)) >= 0.036397) && (dummy < 0.036397)){ 
	
	if(g > 0){
	  
	  cw20_min = (float) (g + srep_size[4]) - 1.0 + ((0.036397 - dummy) * tot_intens / mini_obj_spec[g]);
	  
	} else {
	  
	  cw20_min = (float) (g + srep_size[4]);
	  
	}	  
	
      } 
      
      if(((dummy + (mini_obj_spec[g]/tot_intens)) >= 0.119516) && (dummy < 0.119516)){ 
	
	if(g > 0){
	  
	  cw50_min = (float) (g + srep_size[4]) - 1.0 + ((0.119516 - dummy) * tot_intens / mini_obj_spec[g]);
	  
	} else {
	  
	  cw50_min = (float) (g + srep_size[4]);
	  
	}	  
	
      } 
      
      dummy+=(mini_obj_spec[g]/tot_intens);
      
    }
    cw20_max = cw50_max = dummy = 1.0;
    for(g = srep_size[5] - srep_size[4]; ((g >= 0) && (g > (cw_max - 1 - srep_size[4]))); g--){ 
      
      if(((dummy - (mini_obj_spec[g]/tot_intens)) <= 0.880484) && (dummy > 0.880484)){ 
	
	if(g < (srep_size[5] - srep_size[4])){
	  
	  cw50_max = (float) (g + srep_size[4]) + 1.0 - ((dummy - 0.880484) * tot_intens / mini_obj_spec[g]);
	  
	} else {
	  
	  cw50_max = (float) (g + srep_size[4]);
	  
	}
	
      } 
      
      if(((dummy - (mini_obj_spec[g]/tot_intens)) <= 0.963603) && (dummy > 0.963603)){ 
	
	if(g < (srep_size[5] - srep_size[4])){
	  
	  cw20_max = (float) (g + srep_size[4]) + 1.0 - ((dummy - 0.963603) * tot_intens / mini_obj_spec[g]);
	  
	} else {
	  
	  cw20_max = (float) (g + srep_size[4]);
	  
	}
	
      } 
      
      dummy-=(mini_obj_spec[g]/tot_intens);
      
    }
    
  } else {
    
    w_max = w50_min = w50_max = w20_min = w20_max = cw_max = cw50_min = cw50_max = cw20_min = cw20_max = (float) srep_size[4];
    
  }

}

void object_props::ShowAll_file_WCS(int id, std::fstream& output_file, int cat_mode, double wcs_vals[6]){

  if(cat_mode > 0){
    
    object_props::ShowProps_file_WCS(id,output_file,wcs_vals);
    object_props::ShowSrep_file(output_file);
    output_file << std::endl;
    
  } else {
    
    object_props::ShowProps_file_WCS(id,output_file,wcs_vals); 
    output_file << std::endl;
    
  }

}

void object_props::ShowAll_file(int id, std::fstream& output_file, int cat_mode){

  if(cat_mode > 0){
    
    object_props::ShowProps_file(id,output_file);
    object_props::ShowSrep_file(output_file);
    output_file << std::endl;
    
  } else {
    
    object_props::ShowProps_file(id,output_file); 
    output_file << std::endl;
    
  }

}

void object_props::ShowProps_file_WCS(int id, std::fstream& output_file, double wcs_vals[6]){
    
  output_file << id << " " << NOvox << " " << wcs_vals[0] << " " << wcs_vals[1] << " " << wcs_vals[2] << " " << wcs_vals[3] << " " << wcs_vals[4] << " " << wcs_vals[5] << " " << ra << " " << dec << " " << freq << " " << ra_i << " " << dec_i << " " << freq_i << " " << tot_intens << " " << avg_intens << " " << sigma_intens << " " << rms << " " << min_intens << " " << max_intens << " " << ra_min << " " << ra_max << " " << dec_min << " " << dec_max << " " << freq_min << " " << freq_max << " " << w_max << " " << w50_min << " " << w50_max << " " << w20_min << " " << w20_max << " " << cw_max << " " << cw50_min << " " << cw50_max << " " << cw20_min << " " << cw20_max << " " << std::flush;
  
}

void object_props::ShowProps_file(int id, std::fstream& output_file){
    
  output_file << id << " " << NOvox << " " << ra << " " << dec << " " << freq << " " << ra_i << " " << dec_i << " " << freq_i << " " << tot_intens << " " << avg_intens << " " << sigma_intens << " " << rms << " " << min_intens << " " << max_intens << " " << ra_min << " " << ra_max << " " << dec_min << " " << dec_max << " " << freq_min << " " << freq_max << " " << w_max << " " << w50_min << " " << w50_max << " " << w20_min << " " << w20_max << " " << cw_max << " " << cw50_min << " " << cw50_max << " " << cw20_min << " " << cw20_max << " " << std::flush;
  
}

void object_props::ShowSrep_file(std::fstream& output_file){

  int g;
  
  output_file << ": " << std::flush;
  for(g = 0; g < (1 + ((srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1))); g++){ output_file << srep_grid[g] << " "; }
  output_file << ": " << std::flush;
  for(g = 0; g < (2 * srep_grid[((srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1))]); g++){ output_file << srep_strings[g] << " "; }
  output_file << std::flush;

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

int object_props::ShowArea(){ return (int) ((fabs(ra_max - ra_min) + 1) * (fabs(dec_max - dec_min) + 1)); }

int object_props::ShowVoxels(){ return NOvox; }

int object_props::ShowRArange(){ return (int) (fabs(ra_max - ra_min) + 1); }

int object_props::ShowDECrange(){ return (int) (fabs(dec_max - dec_min) + 1); }

int object_props::ShowFREQrange(){ return (int) (fabs(freq_max - freq_min) + 1); }

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

void object_props::Set_w_max(float value){ w_max = value;}

float object_props::Get_w_max(){ return w_max; }

void object_props::Set_w20_min(float value){ w20_min = value; }

float object_props::Get_w20_min(){ return w20_min; }

void object_props::Set_w20_max(float value){ w20_max = value; }

float object_props::Get_w20_max(){ return w20_max; }

void object_props::Set_w50_min(float value){ w50_min = value; }

float object_props::Get_w50_min(){ return w50_min; }

void object_props::Set_w50_max(float value){ w50_max = value; }

float object_props::Get_w50_max(){ return w50_max; }

void object_props::Set_cw_max(float value){ cw_max = value; }

float object_props::Get_cw_max(){ return cw_max; }

void object_props::Set_cw20_min(float value){ cw20_min = value; }

float object_props::Get_cw20_min(){ return cw20_min; }

void object_props::Set_cw20_max(float value){ cw20_max = value; }

float object_props::Get_cw20_max(){ return cw20_max; }

void object_props::Set_cw50_min(float value){ cw50_min = value; }

float object_props::Get_cw50_min(){ return cw50_min; }

void object_props::Set_cw50_max(float value){ cw50_max = value; }

float object_props::Get_cw50_max(){ return cw50_max; }

void object_props::Set_srep_update(int value){ srep_update = value; }

int object_props::Get_srep_update(){ return srep_update; }

void object_props::Set_srep_size(int index, int value){ srep_size[index] = value; }

int object_props::Get_srep_size(int index){ return srep_size[index]; }

void object_props::Create_srep_grid(int value){ while(srep_grid == NULL){ srep_grid = new int[value]; } }

void object_props::Set_srep_grid(int index, int value){ srep_grid[index] = value; }

int object_props::Get_srep_grid(int index){ return srep_grid[index]; }

void object_props::Free_srep_grid(){ if(srep_grid != NULL){ delete [] srep_grid; srep_grid = NULL; } }

void object_props::Create_srep_strings(int value){ while(srep_strings == NULL){ srep_strings = new int[value]; } }

void object_props::Set_srep_strings(int index, int value){ srep_strings[index] = value; }

int object_props::Get_srep_strings(int index){ return srep_strings[index]; }

void object_props::Free_srep_strings(){ if(srep_strings != NULL){ delete [] srep_strings; srep_strings = NULL; } }

void object_props::Create_mom0(int value){ while(mini_mom0 == NULL){ mini_mom0 = new float[value]; } }

void object_props::Set_mom0(int index, float value){ mini_mom0[index] = value; }

void object_props::Add_mom0(int index, float value){ mini_mom0[index]+=value; }

float object_props::Get_mom0(int index){ return mini_mom0[index]; }

void object_props::Free_mom0(){ if(mini_mom0 != NULL){ delete [] mini_mom0; mini_mom0 = NULL; } }

void object_props::Create_RAPV(int value){ while(mini_RAPV == NULL){ mini_RAPV = new float[value]; } }

void object_props::Set_RAPV(int index, float value){ mini_RAPV[index] = value; }

void object_props::Add_RAPV(int index, float value){ mini_RAPV[index]+=value; }

float object_props::Get_RAPV(int index){ return mini_RAPV[index]; }

void object_props::Free_RAPV(){ if(mini_RAPV != NULL){ delete [] mini_RAPV; mini_RAPV = NULL; } }

void object_props::Create_DECPV(int value){ while(mini_DECPV == NULL){ mini_DECPV = new float[value]; } }

void object_props::Set_DECPV(int index, float value){ mini_DECPV[index] = value; }

void object_props::Add_DECPV(int index, float value){ mini_DECPV[index]+=value; }

float object_props::Get_DECPV(int index){ return mini_DECPV[index]; }

void object_props::Free_DECPV(){ if(mini_DECPV != NULL){ delete [] mini_DECPV; mini_DECPV = NULL; } }

void object_props::Create_ospec(int value){ while(mini_obj_spec == NULL){ mini_obj_spec = new float[value]; } }

void object_props::Set_ospec(int index, float value){ mini_obj_spec[index] = value; }

void object_props::Add_ospec(int index, float value){ mini_obj_spec[index]+=value; }

float object_props::Get_ospec(int index){ return mini_obj_spec[index]; }

void object_props::Free_ospec(){ if(mini_obj_spec != NULL){ delete [] mini_obj_spec; mini_obj_spec = NULL; } }

void object_props::Create_rspec(int value){ while(mini_ref_spec == NULL){ mini_ref_spec = new float[value]; } }

void object_props::Set_rspec(int index, float value){ mini_ref_spec[index] = value; }

void object_props::Add_rspec(int index, float value){ mini_ref_spec[index]+=value; }

float object_props::Get_rspec(int index){ return mini_ref_spec[index]; }

void object_props::Free_rspec(){ if(mini_ref_spec != NULL){ delete [] mini_ref_spec; mini_ref_spec = NULL; } }

void object_props::Create_vfield(int value){ while(mini_vfield == NULL){ mini_vfield = new float[value]; } }

void object_props::Set_vfield(int index, float value){ mini_vfield[index] = value; }

void object_props::Add_vfield(int index, float value){ mini_vfield[index]+=value; }

float object_props::Get_vfield(int index){ return mini_vfield[index]; }

void object_props::Free_vfield(){ if(mini_vfield != NULL){ delete [] mini_vfield; mini_vfield = NULL; } }

void object_props::ReInit_srep(){

  if(srep_grid != NULL){ delete [] srep_grid; }
  srep_grid = NULL;
  if(srep_strings != NULL){ delete [] srep_strings; }
  srep_strings = NULL;
  
}

void object_props::ReInit_mini(){

  if(mini_mom0 != NULL){ delete [] mini_mom0; }
  mini_mom0 = NULL;
  if(mini_RAPV != NULL){ delete [] mini_RAPV; }
  mini_RAPV = NULL;
  if(mini_DECPV != NULL){ delete [] mini_DECPV; }
  mini_DECPV = NULL;
  if(mini_obj_spec != NULL){ delete [] mini_obj_spec; }
  mini_obj_spec = NULL;
  if(mini_ref_spec != NULL){ delete [] mini_ref_spec; }
  mini_ref_spec = NULL;
  if(mini_vfield != NULL){ delete [] mini_vfield; }
  mini_vfield = NULL;

}

void object_props::ReInit_size(){

  srep_size[0] = srep_size[1] = srep_size[2] = srep_size[3] = srep_size[4] = srep_size[5] = -99;

}


