#include<iostream>
#include<RJJ_ObjGen.h>

using namespace std;

// member definitions for float precision class definition

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

object_props::object_props(const object_props & copied){

  int i;

  NOvox = copied.NOvox; 
  ra = copied.ra;
  dec = copied.dec;
  freq = copied.freq;
  ra_i = copied.ra_i;
  dec_i = copied.dec_i;
  freq_i = copied.freq_i;
  tot_intens = copied.tot_intens;
  avg_intens = copied.avg_intens;
  sigma_intens = copied.sigma_intens;
  rms = copied.rms;
  ra_min = copied.ra_min;
  dec_min = copied.dec_min;
  freq_min = copied.freq_min;
  ra_max = copied.ra_max;
  dec_max = copied.dec_max;
  freq_max = copied.freq_max;
  min_intens = copied.min_intens;
  max_intens = copied.max_intens;

  w_max = copied.w_max;
  w20_min = copied.w20_min;
  w20_max = copied.w20_max;
  w50_min = copied.w50_min;
  w50_max = copied.w50_max;
  cw_max = copied.cw_max;
  cw20_min = copied.cw20_min;
  cw20_max = copied.cw20_max;
  cw50_min = copied.cw50_min;
  cw50_max = copied.cw50_max;

  srep_size[0] = copied.srep_size[0];
  srep_size[1] = copied.srep_size[1];
  srep_size[2] = copied.srep_size[2];
  srep_size[3] = copied.srep_size[3];
  srep_size[4] = copied.srep_size[4];
  srep_size[5] = copied.srep_size[5];
  srep_update = copied.srep_update;

  if(copied.srep_grid != NULL){

    srep_grid = new int[(srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1)];
    for(i = 0; i < ((srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1)); ++i){ srep_grid[i] = copied.srep_grid[i]; }    

  } else { srep_grid = NULL; }
  if(copied.srep_strings != NULL){

    srep_strings = new int[(2 * srep_grid[((srep_size[3] - srep_size[2] + 1) * (srep_size[1] - srep_size[0] + 1))])];
    for(i = 0; i < (2 * srep_grid[((srep_size[3] - srep_size[2] + 1) * (srep_size[1] - srep_size[0] + 1))]); ++i){ srep_strings[i] = copied.srep_strings[i]; }

  } else { srep_strings = NULL; }
  if(copied.mini_mom0 != NULL){

    mini_mom0 = new float[(srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1)];
    for(i = 0; i < ((srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1)); ++i){ mini_mom0[i] = copied.mini_mom0[i]; }

  } else { mini_mom0 = NULL; }
  if(copied.mini_RAPV != NULL){

    mini_RAPV = new float[(srep_size[1] - srep_size[0] + 1) * (srep_size[5] - srep_size[4] + 1)];
    for(i = 0; i < ((srep_size[1] - srep_size[0] + 1) * (srep_size[5] - srep_size[4] + 1)); ++i){ mini_RAPV[i] = copied.mini_RAPV[i]; }

  } else { mini_RAPV = NULL; }
  if(copied.mini_DECPV != NULL){

    mini_DECPV = new float[(srep_size[3] - srep_size[2] + 1) * (srep_size[5] - srep_size[4] + 1)];
    for(i = 0; i < ((srep_size[3] - srep_size[2] + 1) * (srep_size[5] - srep_size[4] + 1)); ++i){ mini_DECPV[i] = copied.mini_DECPV[i]; }
    
  } else { mini_DECPV = NULL; }
  if(mini_obj_spec != NULL){

    mini_obj_spec = new float[(srep_size[5] - srep_size[4] + 1)];
    for(i = 0; i < ((srep_size[5] - srep_size[4] + 1)); ++i){ mini_obj_spec[i] = copied.mini_obj_spec[i]; }

  } else { mini_obj_spec = NULL; }
  if(mini_ref_spec != NULL){

    if((copied.srep_size[5] - copied.srep_size[4] + 1) >= 10){

      mini_ref_spec = new float[(srep_size[5] - srep_size[4] + 1)];
      for(i = 0; i < ((srep_size[5] - srep_size[4] + 1)); ++i){ mini_ref_spec[i] = copied.mini_ref_spec[i]; }    

    } else {

      mini_ref_spec = new float[(srep_size[5] - srep_size[4] + 11)];
      for(i = 0; i < ((srep_size[5] - srep_size[4] + 11)); ++i){ mini_ref_spec[i] = copied.mini_ref_spec[i]; }    

    }

  } else { mini_ref_spec = NULL; }
  if(copied.mini_vfield != NULL){

    mini_vfield = new float[(srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1)];
    for(i = 0; i < ((srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1)); ++i){ mini_vfield[i] = copied.mini_vfield[i]; }

  } else { mini_vfield = NULL; }

  // new parameters --- multiple central moment calculations
  p_tot_intens = copied.p_tot_intens;
  n_tot_intens = copied.n_tot_intens;
  p_ra_i = copied.p_tot_intens;
  p_dec_i = copied.p_dec_i;
  p_freq_i = copied.p_freq_i;
  n_ra_i = copied.n_ra_i;
  n_dec_i = copied.n_dec_i;
  n_freq_i = copied.n_freq_i;

}

object_props & object_props::operator = (const object_props & copied){

  int i;

  if(this != &copied){

    NOvox = copied.NOvox; 
    ra = copied.ra;
    dec = copied.dec;
    freq = copied.freq;
    ra_i = copied.ra_i;
    dec_i = copied.dec_i;
    freq_i = copied.freq_i;
    tot_intens = copied.tot_intens;
    avg_intens = copied.avg_intens;
    sigma_intens = copied.sigma_intens;
    rms = copied.rms;
    ra_min = copied.ra_min;
    dec_min = copied.dec_min;
    freq_min = copied.freq_min;
    ra_max = copied.ra_max;
    dec_max = copied.dec_max;
    freq_max = copied.freq_max;
    min_intens = copied.min_intens;
    max_intens = copied.max_intens;
    
    w_max = copied.w_max;
    w20_min = copied.w20_min;
    w20_max = copied.w20_max;
    w50_min = copied.w50_min;
    w50_max = copied.w50_max;
    cw_max = copied.cw_max;
    cw20_min = copied.cw20_min;
    cw20_max = copied.cw20_max;
    cw50_min = copied.cw50_min;
    cw50_max = copied.cw50_max;
    
    srep_size[0] = copied.srep_size[0];
    srep_size[1] = copied.srep_size[1];
    srep_size[2] = copied.srep_size[2];
    srep_size[3] = copied.srep_size[3];
    srep_size[4] = copied.srep_size[4];
    srep_size[5] = copied.srep_size[5];
    srep_update = copied.srep_update;
    
    if(copied.srep_grid != NULL){
      
      srep_grid = new int[(srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1)];
      for(i = 0; i < ((srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1)); ++i){ srep_grid[i] = copied.srep_grid[i]; }    
      
    } else { srep_grid = NULL; }
    if(copied.srep_strings != NULL){
      
      srep_strings = new int[(2 * srep_grid[((srep_size[3] - srep_size[2] + 1) * (srep_size[1] - srep_size[0] + 1))])];
      for(i = 0; i < (2 * srep_grid[((srep_size[3] - srep_size[2] + 1) * (srep_size[1] - srep_size[0] + 1))]); ++i){ srep_strings[i] = copied.srep_strings[i]; }
      
    } else { srep_strings = NULL; }
    if(copied.mini_mom0 != NULL){
      
      mini_mom0 = new float[(srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1)];
      for(i = 0; i < ((srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1)); ++i){ mini_mom0[i] = copied.mini_mom0[i]; }
      
    } else { mini_mom0 = NULL; }
    if(copied.mini_RAPV != NULL){
      
      mini_RAPV = new float[(srep_size[1] - srep_size[0] + 1) * (srep_size[5] - srep_size[4] + 1)];
      for(i = 0; i < ((srep_size[1] - srep_size[0] + 1) * (srep_size[5] - srep_size[4] + 1)); ++i){ mini_RAPV[i] = copied.mini_RAPV[i]; }
      
    } else { mini_RAPV = NULL; }
    if(copied.mini_DECPV != NULL){
      
      mini_DECPV = new float[(srep_size[3] - srep_size[2] + 1) * (srep_size[5] - srep_size[4] + 1)];
      for(i = 0; i < ((srep_size[3] - srep_size[2] + 1) * (srep_size[5] - srep_size[4] + 1)); ++i){ mini_DECPV[i] = copied.mini_DECPV[i]; }
      
    } else { mini_DECPV = NULL; }
    if(mini_obj_spec != NULL){
      
      mini_obj_spec = new float[(srep_size[5] - srep_size[4] + 1)];
      for(i = 0; i < ((srep_size[5] - srep_size[4] + 1)); ++i){ mini_obj_spec[i] = copied.mini_obj_spec[i]; }
      
    } else { mini_obj_spec = NULL; }
    if(mini_ref_spec != NULL){
      
      if((copied.srep_size[5] - copied.srep_size[4] + 1) >= 10){
	
	mini_ref_spec = new float[(srep_size[5] - srep_size[4] + 1)];
	for(i = 0; i < ((srep_size[5] - srep_size[4] + 1)); ++i){ mini_ref_spec[i] = copied.mini_ref_spec[i]; }    
	
      } else {
	
	mini_ref_spec = new float[(srep_size[5] - srep_size[4] + 11)];
	for(i = 0; i < ((srep_size[5] - srep_size[4] + 11)); ++i){ mini_ref_spec[i] = copied.mini_ref_spec[i]; }    
	
      }
      
    } else { mini_ref_spec = NULL; }
    if(copied.mini_vfield != NULL){
      
      mini_vfield = new float[(srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1)];
      for(i = 0; i < ((srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1)); ++i){ mini_vfield[i] = copied.mini_vfield[i]; }
      
    } else { mini_vfield = NULL; }
    
    // new parameters --- multiple central moment calculations
    p_tot_intens = copied.p_tot_intens;
    n_tot_intens = copied.n_tot_intens;
    p_ra_i = copied.p_tot_intens;
    p_dec_i = copied.p_dec_i;
    p_freq_i = copied.p_freq_i;
    n_ra_i = copied.n_ra_i;
    n_dec_i = copied.n_dec_i;
    n_freq_i = copied.n_freq_i;

  }
  return *this;
    
}
  
object_props::object_props(object_props && moved){

  int i;

  NOvox = moved.NOvox; 
  ra = moved.ra;
  dec = moved.dec;
  freq = moved.freq;
  ra_i = moved.ra_i;
  dec_i = moved.dec_i;
  freq_i = moved.freq_i;
  tot_intens = moved.tot_intens;
  avg_intens = moved.avg_intens;
  sigma_intens = moved.sigma_intens;
  rms = moved.rms;
  ra_min = moved.ra_min;
  dec_min = moved.dec_min;
  freq_min = moved.freq_min;
  ra_max = moved.ra_max;
  dec_max = moved.dec_max;
  freq_max = moved.freq_max;
  min_intens = moved.min_intens;
  max_intens = moved.max_intens;

  w_max = moved.w_max;
  w20_min = moved.w20_min;
  w20_max = moved.w20_max;
  w50_min = moved.w50_min;
  w50_max = moved.w50_max;
  cw_max = moved.cw_max;
  cw20_min = moved.cw20_min;
  cw20_max = moved.cw20_max;
  cw50_min = moved.cw50_min;
  cw50_max = moved.cw50_max;

  srep_size[0] = moved.srep_size[0];
  srep_size[1] = moved.srep_size[1];
  srep_size[2] = moved.srep_size[2];
  srep_size[3] = moved.srep_size[3];
  srep_size[4] = moved.srep_size[4];
  srep_size[5] = moved.srep_size[5];
  srep_update = moved.srep_update;

  if(moved.srep_grid != NULL){

    srep_grid = moved.srep_grid;
    moved.srep_grid = NULL;

  } else { srep_grid = NULL; }
  if(moved.srep_strings != NULL){

    srep_strings = moved.srep_strings;
    moved.srep_strings = NULL;

  } else { srep_strings = NULL; }
  if(moved.mini_mom0 != NULL){

    mini_mom0 = moved.mini_mom0;
    moved.mini_mom0 = NULL;

  } else { mini_mom0 = NULL; }
  if(moved.mini_RAPV != NULL){

    mini_RAPV = moved.mini_RAPV;
    moved.mini_RAPV = NULL;

  } else { moved.mini_RAPV = NULL; }
  if(moved.mini_DECPV != NULL){

    mini_DECPV = moved.mini_DECPV;
    moved.mini_DECPV = NULL;
    
  } else { mini_DECPV = NULL; }
  if(mini_obj_spec != NULL){

    mini_obj_spec = moved.mini_obj_spec;
    moved.mini_obj_spec = NULL;

  } else { mini_obj_spec = NULL; }
  if(mini_ref_spec != NULL){

    mini_ref_spec = moved.mini_ref_spec;
    moved.mini_ref_spec = NULL;

  } else { mini_ref_spec = NULL; }
  if(moved.mini_vfield != NULL){

    mini_vfield = moved.mini_vfield;
    moved.mini_vfield = NULL;

  } else { mini_vfield = NULL; }

  // new parameters --- multiple central moment calculations
  p_tot_intens = moved.p_tot_intens;
  n_tot_intens = moved.n_tot_intens;
  p_ra_i = moved.p_tot_intens;
  p_dec_i = moved.p_dec_i;
  p_freq_i = moved.p_freq_i;
  n_ra_i = moved.n_ra_i;
  n_dec_i = moved.n_dec_i;
  n_freq_i = moved.n_freq_i;

}

object_props & object_props::operator = (object_props && moved){

  int i;

  if(this != &moved){

    NOvox = moved.NOvox; 
    ra = moved.ra;
    dec = moved.dec;
    freq = moved.freq;
    ra_i = moved.ra_i;
    dec_i = moved.dec_i;
    freq_i = moved.freq_i;
    tot_intens = moved.tot_intens;
    avg_intens = moved.avg_intens;
    sigma_intens = moved.sigma_intens;
    rms = moved.rms;
    ra_min = moved.ra_min;
    dec_min = moved.dec_min;
    freq_min = moved.freq_min;
    ra_max = moved.ra_max;
    dec_max = moved.dec_max;
    freq_max = moved.freq_max;
    min_intens = moved.min_intens;
    max_intens = moved.max_intens;
    
    w_max = moved.w_max;
    w20_min = moved.w20_min;
    w20_max = moved.w20_max;
    w50_min = moved.w50_min;
    w50_max = moved.w50_max;
    cw_max = moved.cw_max;
    cw20_min = moved.cw20_min;
    cw20_max = moved.cw20_max;
    cw50_min = moved.cw50_min;
    cw50_max = moved.cw50_max;
    
    srep_size[0] = moved.srep_size[0];
    srep_size[1] = moved.srep_size[1];
    srep_size[2] = moved.srep_size[2];
    srep_size[3] = moved.srep_size[3];
    srep_size[4] = moved.srep_size[4];
    srep_size[5] = moved.srep_size[5];
    srep_update = moved.srep_update;
    
    if(moved.srep_grid != NULL){
      
      srep_grid = moved.srep_grid;
      moved.srep_grid = NULL;
      
    } else { srep_grid = NULL; }
    if(moved.srep_strings != NULL){
      
      srep_strings = moved.srep_strings;
      moved.srep_strings = NULL;
      
    } else { srep_strings = NULL; }
    if(moved.mini_mom0 != NULL){
      
      mini_mom0 = moved.mini_mom0;
      moved.mini_mom0 = NULL;
      
    } else { mini_mom0 = NULL; }
    if(moved.mini_RAPV != NULL){
      
      mini_RAPV = moved.mini_RAPV;
      moved.mini_RAPV = NULL;
      
    } else { moved.mini_RAPV = NULL; }
    if(moved.mini_DECPV != NULL){
      
      mini_DECPV = moved.mini_DECPV;
      moved.mini_DECPV = NULL;
      
    } else { mini_DECPV = NULL; }
    if(mini_obj_spec != NULL){
      
      mini_obj_spec = moved.mini_obj_spec;
      moved.mini_obj_spec = NULL;
      
    } else { mini_obj_spec = NULL; }
    if(mini_ref_spec != NULL){
      
      mini_ref_spec = moved.mini_ref_spec;
      moved.mini_ref_spec = NULL;
      
    } else { mini_ref_spec = NULL; }
    if(moved.mini_vfield != NULL){
      
      mini_vfield = moved.mini_vfield;
      moved.mini_vfield = NULL;
      
    } else { mini_vfield = NULL; }
    
    // new parameters --- multiple central moment calculations
    p_tot_intens = moved.p_tot_intens;
    n_tot_intens = moved.n_tot_intens;
    p_ra_i = moved.p_tot_intens;
    p_dec_i = moved.p_dec_i;
    p_freq_i = moved.p_freq_i;
    n_ra_i = moved.n_ra_i;
    n_dec_i = moved.n_dec_i;
    n_freq_i = moved.n_freq_i;

  }
  return *this;

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

bool object_props::operator == (const object_props & compareTo){

  int i;

  if(this == &compareTo){ return true; }

  if(NOvox != compareTo.NOvox){ return false; } 
  if(ra != compareTo.ra){ return false; }
  if(dec != compareTo.dec){ return false; }
  if(freq != compareTo.freq){ return false; }
  if(ra_i != compareTo.ra_i){ return false; }
  if(dec_i != compareTo.dec_i){ return false; }
  if(freq_i != compareTo.freq_i){ return false; }
  if(tot_intens != compareTo.tot_intens){ return false; }
  if(avg_intens != compareTo.avg_intens){ return false; }
  if(sigma_intens != compareTo.sigma_intens){ return false; }
  if(rms != compareTo.rms){ return false; }
  if(ra_min != compareTo.ra_min){ return false; }
  if(dec_min != compareTo.dec_min){ return false; }
  if(freq_min != compareTo.freq_min){ return false; }
  if(ra_max != compareTo.ra_max){ return false; }
  if(dec_max != compareTo.dec_max){ return false; }
  if(freq_max != compareTo.freq_max){ return false; }
  if(min_intens != compareTo.min_intens){ return false; }
  if(max_intens != compareTo.max_intens){ return false; }
  
  if(w_max != compareTo.w_max){ return false; }
  if(w20_min != compareTo.w20_min){ return false; }
  if(w20_max != compareTo.w20_max){ return false; }
  if(w50_min != compareTo.w50_min){ return false; }
  if(w50_max != compareTo.w50_max){ return false; }
  if(cw_max != compareTo.cw_max){ return false; }
  if(cw20_min != compareTo.cw20_min){ return false; }
  if(cw20_max != compareTo.cw20_max){ return false; }
  if(cw50_min != compareTo.cw50_min){ return false; }
  if(cw50_max != compareTo.cw50_max){ return false; }
  
  if(srep_size[0] != compareTo.srep_size[0]){ return false; }
  if(srep_size[1] != compareTo.srep_size[1]){ return false; }
  if(srep_size[2] != compareTo.srep_size[2]){ return false; }
  if(srep_size[3] != compareTo.srep_size[3]){ return false; }
  if(srep_size[4] != compareTo.srep_size[4]){ return false; }
  if(srep_size[5] != compareTo.srep_size[5]){ return false; }
  if(srep_update != compareTo.srep_update){ return false; }
  
  // new parameters --- multiple central moment calculations
  if(p_tot_intens != compareTo.p_tot_intens){ return false; }
  if(n_tot_intens != compareTo.n_tot_intens){ return false; }
  if(p_ra_i != compareTo.p_tot_intens){ return false; }
  if(p_dec_i != compareTo.p_dec_i){ return false; }
  if(p_freq_i != compareTo.p_freq_i){ return false; }
  if(n_ra_i != compareTo.n_ra_i){ return false; }
  if(n_dec_i != compareTo.n_dec_i){ return false; }
  if(n_freq_i != compareTo.n_freq_i){ return false; }

  if(compareTo.srep_grid != NULL){
    
    for(i = 0; i < ((srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1)); ++i){ if(srep_grid[i] != compareTo.srep_grid[i]){ return false; } }    
    
  } else { if(srep_grid != NULL){ return false; } }
  if(compareTo.srep_strings != NULL){
    
    for(i = 0; i < (2 * srep_grid[((srep_size[3] - srep_size[2] + 1) * (srep_size[1] - srep_size[0] + 1))]); ++i){ if(srep_strings[i] != compareTo.srep_strings[i]){ return false; } }
    
  } else { if(srep_strings != NULL){ return false; } }
  if(compareTo.mini_mom0 != NULL){
    
    for(i = 0; i < ((srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1)); ++i){ if(mini_mom0[i] != compareTo.mini_mom0[i]){ return false; } }
    
  } else { if(mini_mom0 != NULL){ return false; } }
  if(compareTo.mini_RAPV != NULL){
    
    for(i = 0; i < ((srep_size[1] - srep_size[0] + 1) * (srep_size[5] - srep_size[4] + 1)); ++i){ if(mini_RAPV[i] != compareTo.mini_RAPV[i]){ return false; } }
    
  } else { if(compareTo.mini_RAPV != NULL){ return false; } }
  if(compareTo.mini_DECPV != NULL){
    
    for(i = 0; i < ((srep_size[3] - srep_size[2] + 1) * (srep_size[5] - srep_size[4] + 1)); ++i){ if(mini_DECPV[i] != compareTo.mini_DECPV[i]){ return false; } }
    
  } else { if(mini_DECPV != NULL){ return false; } }
  if(mini_obj_spec != NULL){
    
    for(i = 0; i < ((srep_size[5] - srep_size[4] + 1)); ++i){ if(mini_obj_spec[i] != compareTo.mini_obj_spec[i]){ return false; } }
    
  } else { if(mini_obj_spec != NULL){ return false; } }
  if(mini_ref_spec != NULL){
    
    if((compareTo.srep_size[5] - compareTo.srep_size[4] + 1) >= 10){
      
      for(i = 0; i < ((srep_size[5] - srep_size[4] + 1)); ++i){ if(mini_ref_spec[i] != compareTo.mini_ref_spec[i]){ return false; } }    
      
    } else {
      
      for(i = 0; i < ((srep_size[5] - srep_size[4] + 11)); ++i){ if(mini_ref_spec[i] != compareTo.mini_ref_spec[i]){ return false; } }    
      
    }
    
  } else { if(mini_ref_spec != NULL){ return false; } }
  if(compareTo.mini_vfield != NULL){
    
    for(i = 0; i < ((srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1)); ++i){ if(mini_vfield[i] != compareTo.mini_vfield[i]){ return false; } }

  } else { if(mini_vfield != NULL){ return false; } }
  
  return true;

}

bool object_props::operator != (const object_props & compareTo){

  if(this == &compareTo){ 

    return false; 

  } else {

    return (!(*this == compareTo));

  } 

}

void object_props::operator += (object_props & added){

  this->AddObject(added);

}

object_props object_props::operator + (object_props & added){

  object_props temp_result(*this);
  temp_result.AddObject(added);
  return temp_result;

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

float object_props::GetRAi_p(){ return p_ra_i; }

float object_props::GetDECi_p(){ return p_dec_i; }

float object_props::GetFREQi_p(){ return p_freq_i; }

float object_props::GetTI_p(){ return p_tot_intens; }

float object_props::GetRAi_n(){ return n_ra_i; }

float object_props::GetDECi_n(){ return n_dec_i; }

float object_props::GetFREQi_n(){ return n_freq_i; }

float object_props::GetTI_n(){ return n_tot_intens; }

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

void object_props::AddPoint(float x_pos, float y_pos, float z_pos, float value){

  // update bounding box
  this->AdjustRArange(x_pos);
  this->AdjustDECrange(y_pos);
  this->AdjustFREQrange(z_pos);
	    
  // update voxel count
  this->AddVoxel(1);
	    
  // call member functions for object_props to adjust appropriate values
  this->AddRa(x_pos);
  this->AddDec(y_pos);
  this->AddFreq(z_pos);
  this->AddRA_i(x_pos,value);
  this->AddDec_i(y_pos,value);
  this->AddFreq_i(z_pos,value);
  this->AddTotIntens(value);
  this->AddAvgIntens(value);
  this->AddSigmaItens(value);
  this->AdjustRange(value);
  
  // change sparse_reps_update
  this->Set_srep_update(1);


}

void object_props::AddObject(object_props & merged){

  int j,k,g,sx,sy,sz;
  vector<int> temp_sparse_reps_grid, temp_sparse_reps_strings;
  vector<float> temp_mom0, temp_RAPV, temp_DECPV, temp_ref_spec, temp_obj_spec, temp_vfield;
    
  if(this != &merged){ 

    // update the properties of this object with the object being merged in
    this->AddVoxel(merged.ShowVoxels());
    this->AdjustRArange(merged.GetRAmin());
    this->AdjustRArange(merged.GetRAmax());
    this->AdjustDECrange(merged.GetDECmin());
    this->AdjustDECrange(merged.GetDECmax());
    this->AdjustFREQrange(merged.GetFREQmin());
    this->AdjustFREQrange(merged.GetFREQmax());
    this->AddRa(merged.GetRA());
    this->AddDec(merged.GetDEC());
    this->AddFreq(merged.GetFREQ());
    this->AddRA_i(1.0,merged.GetRAi());
    this->AddDec_i(1.0,merged.GetDECi());
    this->AddFreq_i(1.0,merged.GetFREQi());
    this->AddTotIntens(merged.GetTI());
    this->AddAvgIntens(merged.GetAvgI());
    this->AddSigmaItens(merged.GetSigmaI());
    this->AdjustRange(merged.GetMinI());		
    this->AdjustRange(merged.GetMaxI());	

    // initialise temporary working arrays
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

    // combine the sparse representations using the temporary arrays --- if they exist
    if((merged.Get_srep_update() != 0) && (merged.Get_srep_size(0) >= 0)){
		  
      // write the existing object's sparse representations into temporary arrays and initialise temporary arrays at the same time, provided it exists
      
      // a. grid
      temp_sparse_reps_grid.resize(0);
      for(g = 0; g < (1 + ((this->GetRAmax() - this->GetRAmin() + 1) * (this->GetDECmax() - this->GetDECmin() + 1))); g++){ temp_sparse_reps_grid.push_back(0); }
      
      if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	
	for(sx = 0; sx < (this->Get_srep_size(1) - this->Get_srep_size(0) + 1); sx++){
	  for(sy = 0; sy < (this->Get_srep_size(3) - this->Get_srep_size(2) + 1); sy++){
	    
	    temp_sparse_reps_grid[(((sy + this->Get_srep_size(2) - this->GetDECmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx + this->Get_srep_size(0) - this->GetRAmin())] = this->Get_srep_grid(((sy * (this->Get_srep_size(1) - this->Get_srep_size(0) + 1)) + sx + 1)) - this->Get_srep_grid(((sy * (this->Get_srep_size(1) - this->Get_srep_size(0) + 1)) + sx));
	    
	  }
	  
	}
	
      }
      
      // b. mom-0
      temp_mom0.resize(0);
      for(g = 0; g < ((this->GetRAmax() - this->GetRAmin() + 1) * (this->GetDECmax() - this->GetDECmin() + 1)); g++){ temp_mom0.push_back(0.0); }
      
      if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	
	for(sx = 0; sx < (this->Get_srep_size(1) - this->Get_srep_size(0) + 1); sx++){
	  
	  for(sy = 0; sy < (this->Get_srep_size(3) - this->Get_srep_size(2) + 1); sy++){
	    
	    temp_mom0[(((sy + this->Get_srep_size(2) - this->GetDECmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx + this->Get_srep_size(0) - this->GetRAmin())]+=this->Get_mom0(((sy * (this->Get_srep_size(1) - this->Get_srep_size(0) + 1)) + sx));
	    
	  }
	  
	}
	
      }
      
      // c. RAPV
      temp_RAPV.resize(0);
      for(g = 0; g < ((this->GetRAmax() - this->GetRAmin() + 1) * (this->GetFREQmax() - this->GetFREQmin() + 1)); g++){ temp_RAPV.push_back(0.0); }
      
      if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	
	for(sx = 0; sx < (this->Get_srep_size(1) - this->Get_srep_size(0) + 1); sx++){
	  
	  for(sz = 0; sz < (this->Get_srep_size(5) - this->Get_srep_size(4) + 1); sz++){
	    
	    temp_RAPV[(((sz + this->Get_srep_size(4) - this->GetFREQmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx + this->Get_srep_size(0) - this->GetRAmin())]+=this->Get_RAPV(((sz * (this->Get_srep_size(1) - this->Get_srep_size(0) + 1)) + sx));
	    
	  }
	  
	}
	
      }
      
      // d. DECPV
      temp_DECPV.resize(0);
      for(g = 0; g < ((this->GetDECmax() - this->GetDECmin() + 1) * (this->GetFREQmax() - this->GetFREQmin() + 1)); g++){ temp_DECPV.push_back(0.0); }
      
      if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	
	for(sy = 0; sy < (this->Get_srep_size(3) - this->Get_srep_size(2) + 1); sy++){
	  
	  for(sz = 0; sz < (this->Get_srep_size(5) - this->Get_srep_size(4) + 1); sz++){
	    
	    temp_DECPV[(((sz + this->Get_srep_size(4) - this->GetFREQmin()) * (this->GetDECmax() - this->GetDECmin() + 1)) + sy + this->Get_srep_size(2) - this->GetDECmin())]+=this->Get_DECPV(((sz * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)) + sy));
	    
	  }
	  
	}
	
      }
      
      // e. ref_spec
      if((this->Get_srep_size(5) - this->Get_srep_size(4) + 1) >= 10){
	
	temp_ref_spec.resize(0);
	for(g = 0; g < (2 * (this->GetFREQmax() - this->GetFREQmin() + 1)); g++){ temp_ref_spec.push_back(0.0); }
	
	if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	  
	  for(sz = 0; sz < (2 * (this->Get_srep_size(5) - this->Get_srep_size(4) + 1)); sz++){ temp_ref_spec[(sz + this->Get_srep_size(4) - this->GetFREQmin())]+=this->Get_rspec(sz); } 
	  
	}
	
      } else {
	
	temp_ref_spec.resize(0);
	for(g = 0; g < (this->GetFREQmax() - this->GetFREQmin() + 11); g++){ temp_ref_spec.push_back(0.0); }
	
	if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	  
	  for(sz = 0; sz < (this->Get_srep_size(5) - this->Get_srep_size(4) + 11); sz++){ temp_ref_spec[(sz + this->Get_srep_size(4) - this->GetFREQmin())]+=this->Get_rspec(sz); } 		    
	  
	}
	
      }
      
      // f. obj_spec
      temp_obj_spec.resize(0);
      for(g = 0; g < (this->Get_srep_size(5) - this->Get_srep_size(4) + 1); g++){ temp_obj_spec.push_back(0.0); }		  
      if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	
	for(sz = 0; sz < (this->Get_srep_size(5) - this->Get_srep_size(4) + 1); sz++){
	  
	  temp_obj_spec[(sz + this->Get_srep_size(4) - this->GetFREQmin())]+=this->Get_ospec(sz);
	  
	}
	
      }
      
      // g. vfield
      temp_vfield.resize(0);
      for(g = 0; g < ((this->GetRAmax() - this->GetRAmin() + 1) * (this->GetDECmax() - this->GetDECmin() + 1)); g++){ temp_vfield.push_back(0.0); }
      
      if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	
	for(sx = 0; sx < (this->Get_srep_size(1) - this->Get_srep_size(0) + 1); sx++){
	  
	  for(sy = 0; sy < (this->Get_srep_size(3) - this->Get_srep_size(2) + 1); sy++){
	    
	    temp_vfield[(((sy + this->Get_srep_size(2) - this->GetDECmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx + this->Get_srep_size(0) - this->GetRAmin())]+=this->Get_vfield(((sy * (this->Get_srep_size(1) - this->Get_srep_size(0) + 1)) + sx));
	    
	  }
	  
	}
	
      }
      
      // write the merged object's sparse representations into temporary arrays
      
      // a. grid
      for(sx = 0; sx < (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1); sx++){
	
	for(sy = 0; sy < (merged.Get_srep_size(3) - merged.Get_srep_size(2) + 1); sy++){
	  
	  temp_sparse_reps_grid[(((sy + merged.Get_srep_size(2) - this->GetDECmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx + merged.Get_srep_size(0) - this->GetRAmin())]+=(merged.Get_srep_grid(((sy * (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1)) + sx + 1)) - merged.Get_srep_grid(((sy * (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1)) + sx)));
	  
	}
	
      }
      
      // b. mom-0
      for(sx = 0; sx < (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1); sx++){
	
	for(sy = 0; sy < (merged.Get_srep_size(3) - merged.Get_srep_size(2) + 1); sy++){
	  
	  temp_mom0[(((sy + merged.Get_srep_size(2) - this->GetDECmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx + merged.Get_srep_size(0) - this->GetRAmin())]+=merged.Get_mom0(((sy * (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1)) + sx));
	  
	}
	
      }
      
      // c. RAPV
      for(sx = 0; sx < (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1); sx++){
	
	for(sz = 0; sz < (merged.Get_srep_size(5) - merged.Get_srep_size(4) + 1); sz++){
	  
	  temp_RAPV[(((sz + merged.Get_srep_size(4) - this->GetFREQmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx + merged.Get_srep_size(0) - this->GetRAmin())]+=merged.Get_RAPV(((sz * (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1)) + sx));
	  
	}
	
      }
      
      // d. DECPV
      for(sy = 0; sy < (merged.Get_srep_size(3) - merged.Get_srep_size(2) + 1); sy++){
	
	for(sz = 0; sz < (merged.Get_srep_size(5) - merged.Get_srep_size(4) + 1); sz++){
	  
	  temp_DECPV[(((sz + merged.Get_srep_size(4) - this->GetFREQmin()) * (this->GetDECmax() - this->GetDECmin() + 1)) + sy + merged.Get_srep_size(2) - this->GetDECmin())]+=merged.Get_DECPV(((sz * (merged.Get_srep_size(3) - merged.Get_srep_size(2) + 1)) + sy));
	  
	}
	
      }
      
      // e. ref_spec
      if((merged.Get_srep_size(5) - merged.Get_srep_size(4) + 1) >= 10){
	
	for(sz = 0; sz < (2 * (merged.Get_srep_size(5) - merged.Get_srep_size(4) + 1)); sz++){ 
	  
	  if((sz + merged.Get_srep_size(4) - merged.GetFREQmin()) >= 0){
	    
	    temp_ref_spec[(sz + merged.Get_srep_size(4) - this->GetFREQmin())]+=merged.Get_rspec(sz); 
	    
	  }
	  
	} 
	
      } else {
	
	for(sz = 0; sz < (merged.Get_srep_size(5) - merged.Get_srep_size(4) + 11); sz++){ 
	  
	  if((sz + merged.Get_srep_size(4) - this->GetFREQmin()) >= 0){
	    
	    temp_ref_spec[(sz + merged.Get_srep_size(4) - this->GetFREQmin())]+=merged.Get_rspec(sz); 
	    
	  }
	  
	} 		    
	
      }
      
      // f. obj_spec
      for(sz = 0; sz < (merged.Get_srep_size(5) - merged.Get_srep_size(4) + 1); sz++){
	
	temp_obj_spec[(sz + merged.Get_srep_size(4) - this->GetFREQmin())]+=merged.Get_ospec(sz);
	
      }
      
      // g. vfield
      for(sx = 0; sx < (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1); sx++){
	
	for(sy = 0; sy < (merged.Get_srep_size(3) - merged.Get_srep_size(2) + 1); sy++){
	  
	  temp_vfield[(((sy + merged.Get_srep_size(2) - this->GetDECmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx + merged.Get_srep_size(0) - this->GetRAmin())]+=merged.Get_vfield(((sy * (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1)) + sx));
	  
	}
	
      }
      
      // convert temp_sparse_reps_grid from differential to cumulative counts using temp_sparse_reps_string as an intermediary
      temp_sparse_reps_strings.resize(1);
      temp_sparse_reps_strings[0] = 0;
      for(g = 1; g < (((this->GetRAmax() - this->GetRAmin() + 1) * (this->GetDECmax() - this->GetDECmin() + 1)) + 1); g++){ temp_sparse_reps_strings.push_back(temp_sparse_reps_strings[(g - 1)] + temp_sparse_reps_grid[(g - 1)]); }
      for(g = 0; g < (((this->GetRAmax() - this->GetRAmin() + 1) * (this->GetDECmax() - this->GetDECmin() + 1)) + 1); g++){ temp_sparse_reps_grid[g] = temp_sparse_reps_strings[g]; }
      
      // write new sparse_reps_strings value to temp_sparse_reps_strings array, using various grids to achieve indexing
      
      // a. initialise temp_sparse_reps_strings
      temp_sparse_reps_strings.resize(0);
      for(g = 0; g < (2 * temp_sparse_reps_grid[((this->GetDECmax() - this->GetDECmin() + 1) * (this->GetRAmax() - this->GetRAmin() + 1))]); g++){ temp_sparse_reps_strings.push_back(0); }
      
      
      // b. for each line of sight through the new existing object bounding box, retrieve the channel range of each object string along this LoS
      for(sy = this->GetDECmin(); sy <= this->GetDECmax(); sy++){
	
	for(sx = this->GetRAmin(); sx <= this->GetRAmax(); sx++){
	  
	  // initialise the number of object strings written to this LoS
	  k = 0;
	  
	  // retrieve the starting index for object strings along this LoS
	  j = 2 * temp_sparse_reps_grid[(((sy - this->GetDECmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx - this->GetRAmin())];
	  
	  // write existing object's object strings to temp_strings_array
	  if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	    
	    if((sx >= this->Get_srep_size(0)) && (sx <= this->Get_srep_size(1)) && (sy >= this->Get_srep_size(2)) && (sy <= this->Get_srep_size(3))){
	      
	      for(g = this->Get_srep_grid((((sy - this->Get_srep_size(2)) * (this->Get_srep_size(1) - this->Get_srep_size(0) + 1)) + sx - this->Get_srep_size(0))); g < this->Get_srep_grid((((sy - this->Get_srep_size(2)) * (this->Get_srep_size(1) - this->Get_srep_size(0) + 1)) + sx - this->Get_srep_size(0) + 1)); g++){
		
		temp_sparse_reps_strings[(j + k)] = this->Get_srep_strings((2 * g));
		k++;
		temp_sparse_reps_strings[(j + k)] = this->Get_srep_strings(((2 * g) + 1));
		k++;
		
	      }
	      
	      // if((sx >= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0]) && (sx <= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][1]) && (sy >= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][2]) && (sy <= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][3]))
	    }
	    
	    // if((sparse_reps_update[obj_batch][(existing - (obj_batch * obj_limit))] != 0) && (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] >= 0))
	  }
	  
	  // write match_init[i] object strings to temp_strings_array
	  if((sx >= merged.Get_srep_size(0)) && (sx <= merged.Get_srep_size(1)) && (sy >= merged.Get_srep_size(2)) && (sy <= merged.Get_srep_size(3))){
	    
	    for(g = merged.Get_srep_grid((((sy - merged.Get_srep_size(2)) * (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1)) + sx - merged.Get_srep_size(0))); g < merged.Get_srep_grid((((sy - merged.Get_srep_size(2)) * (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1)) + sx - merged.Get_srep_size(0) + 1)); g++){
	      
	      temp_sparse_reps_strings[(j + k)] = merged.Get_srep_strings((2 * g));
	      k++;
	      temp_sparse_reps_strings[(j + k)] = merged.Get_srep_strings(((2 * g) + 1));
	      k++;
	      
	    }
	    
	    // if((sx >= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0]) && (sx <= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][1]) && (sy >= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][2]) && (sy <= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][3]))
	  }
	  
	  // for(sx = 0; sx < (); sx++)
	}
	
	// for(sy = 0; sy < (); sy++)
      }
      
      // over-write the existing object's sparse representations with the existing+merged sparse representations
      this->Set_srep_size(0,this->GetRAmin());
      this->Set_srep_size(1,this->GetRAmax());
      this->Set_srep_size(2,this->GetDECmin());
      this->Set_srep_size(3,this->GetDECmax());
      this->Set_srep_size(4,this->GetFREQmin());
      this->Set_srep_size(5,this->GetFREQmax());
      
      // a. grid
      this->Free_srep_grid();
      this->Create_srep_grid((((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)) + 1));
      for(g = 0; g < (((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)) + 1); g++){ this->Set_srep_grid(g,temp_sparse_reps_grid[g]); }
      
      // b. mini_mom0
      this->Free_mom0();
      this->Create_mom0(((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)));
      for(g = 0; g < ((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)); g++){ this->Set_mom0(g,temp_mom0[g]); }
      
      // c. mini_RAPV
      this->Free_RAPV();
      this->Create_RAPV(((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(5) - this->Get_srep_size(4) + 1)));
      for(g = 0; g < ((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(5) - this->Get_srep_size(4) + 1)); g++){ this->Set_RAPV(g,temp_RAPV[g]); }
      
      // d. mini_DECPV
      this->Free_DECPV();
      this->Create_DECPV(((this->Get_srep_size(5) - this->Get_srep_size(4) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)));
      for(g = 0; g < ((this->Get_srep_size(3) - this->Get_srep_size(2) + 1) * (this->Get_srep_size(5) - this->Get_srep_size(4) + 1)); g++){ this->Set_DECPV(g,temp_DECPV[g]); }
      
      // e. mini_obj_spec
      this->Free_ospec();
      this->Create_ospec((this->Get_srep_size(5) - this->Get_srep_size(4) + 1));
      for(g = 0; g < (this->Get_srep_size(5) - this->Get_srep_size(4) + 1); g++){ this->Set_ospec(g,temp_obj_spec[g]); }
      
      // f. mini_ref_spec
      this->Free_rspec();
      if((this->Get_srep_size(5) - this->Get_srep_size(4) + 1) >= 10){
	
	this->Create_rspec((2 * (this->Get_srep_size(5) - this->Get_srep_size(4) + 1)));
	for(g = 0; g < (2 * (this->Get_srep_size(5) - this->Get_srep_size(4) + 1)); g++){ this->Set_rspec(g,temp_ref_spec[g]); }
	
      } else {
	
	this->Create_rspec((this->Get_srep_size(5) - this->Get_srep_size(4) + 11));
	for(g = 0; g < (this->Get_srep_size(5) - this->Get_srep_size(4) + 11); g++){ this->Set_rspec(g,temp_ref_spec[g]); }
	
      }
      
      // g. mini_vfield
      this->Free_vfield();
      this->Create_vfield(((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)));
      for(g = 0; g < ((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)); g++){ this->Set_vfield(g,temp_vfield[g]); }
      
      // h. sparse_reps_strings
      this->Free_srep_strings();
      this->Create_srep_strings((2 * this->Get_srep_grid(((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)))));
      for(g = 0; g < (2 * this->Get_srep_grid(((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)))); g++){ this->Set_srep_strings(g,temp_sparse_reps_strings[g]); }
      
    }
  
    // if(this != &merged)
  }

}

void object_props::AddObject(object_props & merged,vector<int> & temp_sparse_reps_grid,vector<int> & temp_sparse_reps_strings,vector<float> & temp_mom0,vector<float> & temp_RAPV,vector<float> & temp_DECPV,vector<float> & temp_ref_spec,vector<float> & temp_obj_spec,vector<float> & temp_vfield){

  int j,k,g,sx,sy,sz;
    
  if(this != &merged){ 

    // update the properties of this object with the object being merged in
    this->AddVoxel(merged.ShowVoxels());
    this->AdjustRArange(merged.GetRAmin());
    this->AdjustRArange(merged.GetRAmax());
    this->AdjustDECrange(merged.GetDECmin());
    this->AdjustDECrange(merged.GetDECmax());
    this->AdjustFREQrange(merged.GetFREQmin());
    this->AdjustFREQrange(merged.GetFREQmax());
    this->AddRa(merged.GetRA());
    this->AddDec(merged.GetDEC());
    this->AddFreq(merged.GetFREQ());
    this->AddRA_i(1.0,merged.GetRAi());
    this->AddDec_i(1.0,merged.GetDECi());
    this->AddFreq_i(1.0,merged.GetFREQi());
    this->AddTotIntens(merged.GetTI());
    this->AddAvgIntens(merged.GetAvgI());
    this->AddSigmaItens(merged.GetSigmaI());
    this->AdjustRange(merged.GetMinI());		
    this->AdjustRange(merged.GetMaxI());	

    // combine the sparse representations using the temporary arrays --- if they exist
    if((merged.Get_srep_update() != 0) && (merged.Get_srep_size(0) >= 0)){
		  
      // write the existing object's sparse representations into temporary arrays and initialise temporary arrays at the same time, provided it exists
      
      // a. grid
      temp_sparse_reps_grid.resize(0);
      for(g = 0; g < (1 + ((this->GetRAmax() - this->GetRAmin() + 1) * (this->GetDECmax() - this->GetDECmin() + 1))); g++){ temp_sparse_reps_grid.push_back(0); }
      
      if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	
	for(sx = 0; sx < (this->Get_srep_size(1) - this->Get_srep_size(0) + 1); sx++){
	  for(sy = 0; sy < (this->Get_srep_size(3) - this->Get_srep_size(2) + 1); sy++){
	    
	    temp_sparse_reps_grid[(((sy + this->Get_srep_size(2) - this->GetDECmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx + this->Get_srep_size(0) - this->GetRAmin())] = this->Get_srep_grid(((sy * (this->Get_srep_size(1) - this->Get_srep_size(0) + 1)) + sx + 1)) - this->Get_srep_grid(((sy * (this->Get_srep_size(1) - this->Get_srep_size(0) + 1)) + sx));
	    
	  }
	  
	}
	
      }
      
      // b. mom-0
      temp_mom0.resize(0);
      for(g = 0; g < ((this->GetRAmax() - this->GetRAmin() + 1) * (this->GetDECmax() - this->GetDECmin() + 1)); g++){ temp_mom0.push_back(0.0); }
      
      if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	
	for(sx = 0; sx < (this->Get_srep_size(1) - this->Get_srep_size(0) + 1); sx++){
	  
	  for(sy = 0; sy < (this->Get_srep_size(3) - this->Get_srep_size(2) + 1); sy++){
	    
	    temp_mom0[(((sy + this->Get_srep_size(2) - this->GetDECmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx + this->Get_srep_size(0) - this->GetRAmin())]+=this->Get_mom0(((sy * (this->Get_srep_size(1) - this->Get_srep_size(0) + 1)) + sx));
	    
	  }
	  
	}
	
      }
      
      // c. RAPV
      temp_RAPV.resize(0);
      for(g = 0; g < ((this->GetRAmax() - this->GetRAmin() + 1) * (this->GetFREQmax() - this->GetFREQmin() + 1)); g++){ temp_RAPV.push_back(0.0); }
      
      if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	
	for(sx = 0; sx < (this->Get_srep_size(1) - this->Get_srep_size(0) + 1); sx++){
	  
	  for(sz = 0; sz < (this->Get_srep_size(5) - this->Get_srep_size(4) + 1); sz++){
	    
	    temp_RAPV[(((sz + this->Get_srep_size(4) - this->GetFREQmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx + this->Get_srep_size(0) - this->GetRAmin())]+=this->Get_RAPV(((sz * (this->Get_srep_size(1) - this->Get_srep_size(0) + 1)) + sx));
	    
	  }
	  
	}
	
      }
      
      // d. DECPV
      temp_DECPV.resize(0);
      for(g = 0; g < ((this->GetDECmax() - this->GetDECmin() + 1) * (this->GetFREQmax() - this->GetFREQmin() + 1)); g++){ temp_DECPV.push_back(0.0); }
      
      if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	
	for(sy = 0; sy < (this->Get_srep_size(3) - this->Get_srep_size(2) + 1); sy++){
	  
	  for(sz = 0; sz < (this->Get_srep_size(5) - this->Get_srep_size(4) + 1); sz++){
	    
	    temp_DECPV[(((sz + this->Get_srep_size(4) - this->GetFREQmin()) * (this->GetDECmax() - this->GetDECmin() + 1)) + sy + this->Get_srep_size(2) - this->GetDECmin())]+=this->Get_DECPV(((sz * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)) + sy));
	    
	  }
	  
	}
	
      }
      
      // e. ref_spec
      if((this->Get_srep_size(5) - this->Get_srep_size(4) + 1) >= 10){
	
	temp_ref_spec.resize(0);
	for(g = 0; g < (2 * (this->GetFREQmax() - this->GetFREQmin() + 1)); g++){ temp_ref_spec.push_back(0.0); }
	
	if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	  
	  for(sz = 0; sz < (2 * (this->Get_srep_size(5) - this->Get_srep_size(4) + 1)); sz++){ temp_ref_spec[(sz + this->Get_srep_size(4) - this->GetFREQmin())]+=this->Get_rspec(sz); } 
	  
	}
	
      } else {
	
	temp_ref_spec.resize(0);
	for(g = 0; g < (this->GetFREQmax() - this->GetFREQmin() + 11); g++){ temp_ref_spec.push_back(0.0); }
	
	if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	  
	  for(sz = 0; sz < (this->Get_srep_size(5) - this->Get_srep_size(4) + 11); sz++){ temp_ref_spec[(sz + this->Get_srep_size(4) - this->GetFREQmin())]+=this->Get_rspec(sz); } 		    
	  
	}
	
      }
      
      // f. obj_spec
      temp_obj_spec.resize(0);
      for(g = 0; g < (this->Get_srep_size(5) - this->Get_srep_size(4) + 1); g++){ temp_obj_spec.push_back(0.0); }		  
      if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	
	for(sz = 0; sz < (this->Get_srep_size(5) - this->Get_srep_size(4) + 1); sz++){
	  
	  temp_obj_spec[(sz + this->Get_srep_size(4) - this->GetFREQmin())]+=this->Get_ospec(sz);
	  
	}
	
      }
      
      // g. vfield
      temp_vfield.resize(0);
      for(g = 0; g < ((this->GetRAmax() - this->GetRAmin() + 1) * (this->GetDECmax() - this->GetDECmin() + 1)); g++){ temp_vfield.push_back(0.0); }
      
      if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	
	for(sx = 0; sx < (this->Get_srep_size(1) - this->Get_srep_size(0) + 1); sx++){
	  
	  for(sy = 0; sy < (this->Get_srep_size(3) - this->Get_srep_size(2) + 1); sy++){
	    
	    temp_vfield[(((sy + this->Get_srep_size(2) - this->GetDECmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx + this->Get_srep_size(0) - this->GetRAmin())]+=this->Get_vfield(((sy * (this->Get_srep_size(1) - this->Get_srep_size(0) + 1)) + sx));
	    
	  }
	  
	}
	
      }
      
      // write the merged object's sparse representations into temporary arrays
      
      // a. grid
      for(sx = 0; sx < (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1); sx++){
	
	for(sy = 0; sy < (merged.Get_srep_size(3) - merged.Get_srep_size(2) + 1); sy++){
	  
	  temp_sparse_reps_grid[(((sy + merged.Get_srep_size(2) - this->GetDECmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx + merged.Get_srep_size(0) - this->GetRAmin())]+=(merged.Get_srep_grid(((sy * (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1)) + sx + 1)) - merged.Get_srep_grid(((sy * (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1)) + sx)));
	  
	}
	
      }
      
      // b. mom-0
      for(sx = 0; sx < (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1); sx++){
	
	for(sy = 0; sy < (merged.Get_srep_size(3) - merged.Get_srep_size(2) + 1); sy++){
	  
	  temp_mom0[(((sy + merged.Get_srep_size(2) - this->GetDECmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx + merged.Get_srep_size(0) - this->GetRAmin())]+=merged.Get_mom0(((sy * (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1)) + sx));
	  
	}
	
      }
      
      // c. RAPV
      for(sx = 0; sx < (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1); sx++){
	
	for(sz = 0; sz < (merged.Get_srep_size(5) - merged.Get_srep_size(4) + 1); sz++){
	  
	  temp_RAPV[(((sz + merged.Get_srep_size(4) - this->GetFREQmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx + merged.Get_srep_size(0) - this->GetRAmin())]+=merged.Get_RAPV(((sz * (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1)) + sx));
	  
	}
	
      }
      
      // d. DECPV
      for(sy = 0; sy < (merged.Get_srep_size(3) - merged.Get_srep_size(2) + 1); sy++){
	
	for(sz = 0; sz < (merged.Get_srep_size(5) - merged.Get_srep_size(4) + 1); sz++){
	  
	  temp_DECPV[(((sz + merged.Get_srep_size(4) - this->GetFREQmin()) * (this->GetDECmax() - this->GetDECmin() + 1)) + sy + merged.Get_srep_size(2) - this->GetDECmin())]+=merged.Get_DECPV(((sz * (merged.Get_srep_size(3) - merged.Get_srep_size(2) + 1)) + sy));
	  
	}
	
      }
      
      // e. ref_spec
      if((merged.Get_srep_size(5) - merged.Get_srep_size(4) + 1) >= 10){
	
	for(sz = 0; sz < (2 * (merged.Get_srep_size(5) - merged.Get_srep_size(4) + 1)); sz++){ 
	  
	  if((sz + merged.Get_srep_size(4) - merged.GetFREQmin()) >= 0){
	    
	    temp_ref_spec[(sz + merged.Get_srep_size(4) - this->GetFREQmin())]+=merged.Get_rspec(sz); 
	    
	  }
	  
	} 
	
      } else {
	
	for(sz = 0; sz < (merged.Get_srep_size(5) - merged.Get_srep_size(4) + 11); sz++){ 
	  
	  if((sz + merged.Get_srep_size(4) - this->GetFREQmin()) >= 0){
	    
	    temp_ref_spec[(sz + merged.Get_srep_size(4) - this->GetFREQmin())]+=merged.Get_rspec(sz); 
	    
	  }
	  
	} 		    
	
      }
      
      // f. obj_spec
      for(sz = 0; sz < (merged.Get_srep_size(5) - merged.Get_srep_size(4) + 1); sz++){
	
	temp_obj_spec[(sz + merged.Get_srep_size(4) - this->GetFREQmin())]+=merged.Get_ospec(sz);
	
      }
      
      // g. vfield
      for(sx = 0; sx < (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1); sx++){
	
	for(sy = 0; sy < (merged.Get_srep_size(3) - merged.Get_srep_size(2) + 1); sy++){
	  
	  temp_vfield[(((sy + merged.Get_srep_size(2) - this->GetDECmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx + merged.Get_srep_size(0) - this->GetRAmin())]+=merged.Get_vfield(((sy * (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1)) + sx));
	  
	}
	
      }
      
      // convert temp_sparse_reps_grid from differential to cumulative counts using temp_sparse_reps_string as an intermediary
      temp_sparse_reps_strings.resize(1);
      temp_sparse_reps_strings[0] = 0;
      for(g = 1; g < (((this->GetRAmax() - this->GetRAmin() + 1) * (this->GetDECmax() - this->GetDECmin() + 1)) + 1); g++){ temp_sparse_reps_strings.push_back(temp_sparse_reps_strings[(g - 1)] + temp_sparse_reps_grid[(g - 1)]); }
      for(g = 0; g < (((this->GetRAmax() - this->GetRAmin() + 1) * (this->GetDECmax() - this->GetDECmin() + 1)) + 1); g++){ temp_sparse_reps_grid[g] = temp_sparse_reps_strings[g]; }
      
      // write new sparse_reps_strings value to temp_sparse_reps_strings array, using various grids to achieve indexing
      
      // a. initialise temp_sparse_reps_strings
      temp_sparse_reps_strings.resize(0);
      for(g = 0; g < (2 * temp_sparse_reps_grid[((this->GetDECmax() - this->GetDECmin() + 1) * (this->GetRAmax() - this->GetRAmin() + 1))]); g++){ temp_sparse_reps_strings.push_back(0); }
      
      
      // b. for each line of sight through the new existing object bounding box, retrieve the channel range of each object string along this LoS
      for(sy = this->GetDECmin(); sy <= this->GetDECmax(); sy++){
	
	for(sx = this->GetRAmin(); sx <= this->GetRAmax(); sx++){
	  
	  // initialise the number of object strings written to this LoS
	  k = 0;
	  
	  // retrieve the starting index for object strings along this LoS
	  j = 2 * temp_sparse_reps_grid[(((sy - this->GetDECmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx - this->GetRAmin())];
	  
	  // write existing object's object strings to temp_strings_array
	  if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	    
	    if((sx >= this->Get_srep_size(0)) && (sx <= this->Get_srep_size(1)) && (sy >= this->Get_srep_size(2)) && (sy <= this->Get_srep_size(3))){
	      
	      for(g = this->Get_srep_grid((((sy - this->Get_srep_size(2)) * (this->Get_srep_size(1) - this->Get_srep_size(0) + 1)) + sx - this->Get_srep_size(0))); g < this->Get_srep_grid((((sy - this->Get_srep_size(2)) * (this->Get_srep_size(1) - this->Get_srep_size(0) + 1)) + sx - this->Get_srep_size(0) + 1)); g++){
		
		temp_sparse_reps_strings[(j + k)] = this->Get_srep_strings((2 * g));
		k++;
		temp_sparse_reps_strings[(j + k)] = this->Get_srep_strings(((2 * g) + 1));
		k++;
		
	      }
	      
	      // if((sx >= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0]) && (sx <= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][1]) && (sy >= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][2]) && (sy <= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][3]))
	    }
	    
	    // if((sparse_reps_update[obj_batch][(existing - (obj_batch * obj_limit))] != 0) && (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] >= 0))
	  }
	  
	  // write match_init[i] object strings to temp_strings_array
	  if((sx >= merged.Get_srep_size(0)) && (sx <= merged.Get_srep_size(1)) && (sy >= merged.Get_srep_size(2)) && (sy <= merged.Get_srep_size(3))){
	    
	    for(g = merged.Get_srep_grid((((sy - merged.Get_srep_size(2)) * (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1)) + sx - merged.Get_srep_size(0))); g < merged.Get_srep_grid((((sy - merged.Get_srep_size(2)) * (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1)) + sx - merged.Get_srep_size(0) + 1)); g++){
	      
	      temp_sparse_reps_strings[(j + k)] = merged.Get_srep_strings((2 * g));
	      k++;
	      temp_sparse_reps_strings[(j + k)] = merged.Get_srep_strings(((2 * g) + 1));
	      k++;
	      
	    }
	    
	    // if((sx >= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0]) && (sx <= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][1]) && (sy >= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][2]) && (sy <= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][3]))
	  }
	  
	  // for(sx = 0; sx < (); sx++)
	}
	
	// for(sy = 0; sy < (); sy++)
      }
      
      // over-write the existing object's sparse representations with the existing+merged sparse representations
      this->Set_srep_size(0,this->GetRAmin());
      this->Set_srep_size(1,this->GetRAmax());
      this->Set_srep_size(2,this->GetDECmin());
      this->Set_srep_size(3,this->GetDECmax());
      this->Set_srep_size(4,this->GetFREQmin());
      this->Set_srep_size(5,this->GetFREQmax());
      
      // a. grid
      this->Free_srep_grid();
      this->Create_srep_grid((((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)) + 1));
      for(g = 0; g < (((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)) + 1); g++){ this->Set_srep_grid(g,temp_sparse_reps_grid[g]); }
      
      // b. mini_mom0
      this->Free_mom0();
      this->Create_mom0(((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)));
      for(g = 0; g < ((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)); g++){ this->Set_mom0(g,temp_mom0[g]); }
      
      // c. mini_RAPV
      this->Free_RAPV();
      this->Create_RAPV(((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(5) - this->Get_srep_size(4) + 1)));
      for(g = 0; g < ((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(5) - this->Get_srep_size(4) + 1)); g++){ this->Set_RAPV(g,temp_RAPV[g]); }
      
      // d. mini_DECPV
      this->Free_DECPV();
      this->Create_DECPV(((this->Get_srep_size(5) - this->Get_srep_size(4) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)));
      for(g = 0; g < ((this->Get_srep_size(3) - this->Get_srep_size(2) + 1) * (this->Get_srep_size(5) - this->Get_srep_size(4) + 1)); g++){ this->Set_DECPV(g,temp_DECPV[g]); }
      
      // e. mini_obj_spec
      this->Free_ospec();
      this->Create_ospec((this->Get_srep_size(5) - this->Get_srep_size(4) + 1));
      for(g = 0; g < (this->Get_srep_size(5) - this->Get_srep_size(4) + 1); g++){ this->Set_ospec(g,temp_obj_spec[g]); }
      
      // f. mini_ref_spec
      this->Free_rspec();
      if((this->Get_srep_size(5) - this->Get_srep_size(4) + 1) >= 10){
	
	this->Create_rspec((2 * (this->Get_srep_size(5) - this->Get_srep_size(4) + 1)));
	for(g = 0; g < (2 * (this->Get_srep_size(5) - this->Get_srep_size(4) + 1)); g++){ this->Set_rspec(g,temp_ref_spec[g]); }
	
      } else {
	
	this->Create_rspec((this->Get_srep_size(5) - this->Get_srep_size(4) + 11));
	for(g = 0; g < (this->Get_srep_size(5) - this->Get_srep_size(4) + 11); g++){ this->Set_rspec(g,temp_ref_spec[g]); }
	
      }
      
      // g. mini_vfield
      this->Free_vfield();
      this->Create_vfield(((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)));
      for(g = 0; g < ((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)); g++){ this->Set_vfield(g,temp_vfield[g]); }
      
      // h. sparse_reps_strings
      this->Free_srep_strings();
      this->Create_srep_strings((2 * this->Get_srep_grid(((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)))));
      for(g = 0; g < (2 * this->Get_srep_grid(((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)))); g++){ this->Set_srep_strings(g,temp_sparse_reps_strings[g]); }
      
    }
  
    // if(this != &merged)
  }

}

// member definitions for double precision class definition

object_props_dbl::object_props_dbl(){ 
    
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

object_props_dbl::object_props_dbl(const object_props_dbl & copied){

  int i;

  NOvox = copied.NOvox; 
  ra = copied.ra;
  dec = copied.dec;
  freq = copied.freq;
  ra_i = copied.ra_i;
  dec_i = copied.dec_i;
  freq_i = copied.freq_i;
  tot_intens = copied.tot_intens;
  avg_intens = copied.avg_intens;
  sigma_intens = copied.sigma_intens;
  rms = copied.rms;
  ra_min = copied.ra_min;
  dec_min = copied.dec_min;
  freq_min = copied.freq_min;
  ra_max = copied.ra_max;
  dec_max = copied.dec_max;
  freq_max = copied.freq_max;
  min_intens = copied.min_intens;
  max_intens = copied.max_intens;

  w_max = copied.w_max;
  w20_min = copied.w20_min;
  w20_max = copied.w20_max;
  w50_min = copied.w50_min;
  w50_max = copied.w50_max;
  cw_max = copied.cw_max;
  cw20_min = copied.cw20_min;
  cw20_max = copied.cw20_max;
  cw50_min = copied.cw50_min;
  cw50_max = copied.cw50_max;

  srep_size[0] = copied.srep_size[0];
  srep_size[1] = copied.srep_size[1];
  srep_size[2] = copied.srep_size[2];
  srep_size[3] = copied.srep_size[3];
  srep_size[4] = copied.srep_size[4];
  srep_size[5] = copied.srep_size[5];
  srep_update = copied.srep_update;

  if(copied.srep_grid != NULL){

    srep_grid = new int[(srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1)];
    for(i = 0; i < ((srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1)); ++i){ srep_grid[i] = copied.srep_grid[i]; }    

  } else { srep_grid = NULL; }
  if(copied.srep_strings != NULL){

    srep_strings = new int[(2 * srep_grid[((srep_size[3] - srep_size[2] + 1) * (srep_size[1] - srep_size[0] + 1))])];
    for(i = 0; i < (2 * srep_grid[((srep_size[3] - srep_size[2] + 1) * (srep_size[1] - srep_size[0] + 1))]); ++i){ srep_strings[i] = copied.srep_strings[i]; }

  } else { srep_strings = NULL; }
  if(copied.mini_mom0 != NULL){

    mini_mom0 = new double[(srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1)];
    for(i = 0; i < ((srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1)); ++i){ mini_mom0[i] = copied.mini_mom0[i]; }

  } else { mini_mom0 = NULL; }
  if(copied.mini_RAPV != NULL){

    mini_RAPV = new double[(srep_size[1] - srep_size[0] + 1) * (srep_size[5] - srep_size[4] + 1)];
    for(i = 0; i < ((srep_size[1] - srep_size[0] + 1) * (srep_size[5] - srep_size[4] + 1)); ++i){ mini_RAPV[i] = copied.mini_RAPV[i]; }

  } else { mini_RAPV = NULL; }
  if(copied.mini_DECPV != NULL){

    mini_DECPV = new double[(srep_size[3] - srep_size[2] + 1) * (srep_size[5] - srep_size[4] + 1)];
    for(i = 0; i < ((srep_size[3] - srep_size[2] + 1) * (srep_size[5] - srep_size[4] + 1)); ++i){ mini_DECPV[i] = copied.mini_DECPV[i]; }
    
  } else { mini_DECPV = NULL; }
  if(mini_obj_spec != NULL){

    mini_obj_spec = new double[(srep_size[5] - srep_size[4] + 1)];
    for(i = 0; i < ((srep_size[5] - srep_size[4] + 1)); ++i){ mini_obj_spec[i] = copied.mini_obj_spec[i]; }

  } else { mini_obj_spec = NULL; }
  if(mini_ref_spec != NULL){

    if((copied.srep_size[5] - copied.srep_size[4] + 1) >= 10){

      mini_ref_spec = new double[(srep_size[5] - srep_size[4] + 1)];
      for(i = 0; i < ((srep_size[5] - srep_size[4] + 1)); ++i){ mini_ref_spec[i] = copied.mini_ref_spec[i]; }    

    } else {

      mini_ref_spec = new double[(srep_size[5] - srep_size[4] + 11)];
      for(i = 0; i < ((srep_size[5] - srep_size[4] + 11)); ++i){ mini_ref_spec[i] = copied.mini_ref_spec[i]; }    

    }

  } else { mini_ref_spec = NULL; }
  if(copied.mini_vfield != NULL){

    mini_vfield = new double[(srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1)];
    for(i = 0; i < ((srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1)); ++i){ mini_vfield[i] = copied.mini_vfield[i]; }

  } else { mini_vfield = NULL; }

  // new parameters --- multiple central moment calculations
  p_tot_intens = copied.p_tot_intens;
  n_tot_intens = copied.n_tot_intens;
  p_ra_i = copied.p_tot_intens;
  p_dec_i = copied.p_dec_i;
  p_freq_i = copied.p_freq_i;
  n_ra_i = copied.n_ra_i;
  n_dec_i = copied.n_dec_i;
  n_freq_i = copied.n_freq_i;

}

object_props_dbl & object_props_dbl::operator = (const object_props_dbl & copied){

  int i;

  if(this != &copied){

    NOvox = copied.NOvox; 
    ra = copied.ra;
    dec = copied.dec;
    freq = copied.freq;
    ra_i = copied.ra_i;
    dec_i = copied.dec_i;
    freq_i = copied.freq_i;
    tot_intens = copied.tot_intens;
    avg_intens = copied.avg_intens;
    sigma_intens = copied.sigma_intens;
    rms = copied.rms;
    ra_min = copied.ra_min;
    dec_min = copied.dec_min;
    freq_min = copied.freq_min;
    ra_max = copied.ra_max;
    dec_max = copied.dec_max;
    freq_max = copied.freq_max;
    min_intens = copied.min_intens;
    max_intens = copied.max_intens;
    
    w_max = copied.w_max;
    w20_min = copied.w20_min;
    w20_max = copied.w20_max;
    w50_min = copied.w50_min;
    w50_max = copied.w50_max;
    cw_max = copied.cw_max;
    cw20_min = copied.cw20_min;
    cw20_max = copied.cw20_max;
    cw50_min = copied.cw50_min;
    cw50_max = copied.cw50_max;
    
    srep_size[0] = copied.srep_size[0];
    srep_size[1] = copied.srep_size[1];
    srep_size[2] = copied.srep_size[2];
    srep_size[3] = copied.srep_size[3];
    srep_size[4] = copied.srep_size[4];
    srep_size[5] = copied.srep_size[5];
    srep_update = copied.srep_update;
    
    if(copied.srep_grid != NULL){
      
      srep_grid = new int[(srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1)];
      for(i = 0; i < ((srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1)); ++i){ srep_grid[i] = copied.srep_grid[i]; }    
      
    } else { srep_grid = NULL; }
    if(copied.srep_strings != NULL){
      
      srep_strings = new int[(2 * srep_grid[((srep_size[3] - srep_size[2] + 1) * (srep_size[1] - srep_size[0] + 1))])];
      for(i = 0; i < (2 * srep_grid[((srep_size[3] - srep_size[2] + 1) * (srep_size[1] - srep_size[0] + 1))]); ++i){ srep_strings[i] = copied.srep_strings[i]; }
      
    } else { srep_strings = NULL; }
    if(copied.mini_mom0 != NULL){
      
      mini_mom0 = new double[(srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1)];
      for(i = 0; i < ((srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1)); ++i){ mini_mom0[i] = copied.mini_mom0[i]; }
      
    } else { mini_mom0 = NULL; }
    if(copied.mini_RAPV != NULL){
      
      mini_RAPV = new double[(srep_size[1] - srep_size[0] + 1) * (srep_size[5] - srep_size[4] + 1)];
      for(i = 0; i < ((srep_size[1] - srep_size[0] + 1) * (srep_size[5] - srep_size[4] + 1)); ++i){ mini_RAPV[i] = copied.mini_RAPV[i]; }
      
    } else { mini_RAPV = NULL; }
    if(copied.mini_DECPV != NULL){
      
      mini_DECPV = new double[(srep_size[3] - srep_size[2] + 1) * (srep_size[5] - srep_size[4] + 1)];
      for(i = 0; i < ((srep_size[3] - srep_size[2] + 1) * (srep_size[5] - srep_size[4] + 1)); ++i){ mini_DECPV[i] = copied.mini_DECPV[i]; }
      
    } else { mini_DECPV = NULL; }
    if(mini_obj_spec != NULL){
      
      mini_obj_spec = new double[(srep_size[5] - srep_size[4] + 1)];
      for(i = 0; i < ((srep_size[5] - srep_size[4] + 1)); ++i){ mini_obj_spec[i] = copied.mini_obj_spec[i]; }
      
    } else { mini_obj_spec = NULL; }
    if(mini_ref_spec != NULL){
      
      if((copied.srep_size[5] - copied.srep_size[4] + 1) >= 10){
	
	mini_ref_spec = new double[(srep_size[5] - srep_size[4] + 1)];
	for(i = 0; i < ((srep_size[5] - srep_size[4] + 1)); ++i){ mini_ref_spec[i] = copied.mini_ref_spec[i]; }    
	
      } else {
	
	mini_ref_spec = new double[(srep_size[5] - srep_size[4] + 11)];
	for(i = 0; i < ((srep_size[5] - srep_size[4] + 11)); ++i){ mini_ref_spec[i] = copied.mini_ref_spec[i]; }    
	
      }
      
    } else { mini_ref_spec = NULL; }
    if(copied.mini_vfield != NULL){
      
      mini_vfield = new double[(srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1)];
      for(i = 0; i < ((srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1)); ++i){ mini_vfield[i] = copied.mini_vfield[i]; }
      
    } else { mini_vfield = NULL; }
    
    // new parameters --- multiple central moment calculations
    p_tot_intens = copied.p_tot_intens;
    n_tot_intens = copied.n_tot_intens;
    p_ra_i = copied.p_tot_intens;
    p_dec_i = copied.p_dec_i;
    p_freq_i = copied.p_freq_i;
    n_ra_i = copied.n_ra_i;
    n_dec_i = copied.n_dec_i;
    n_freq_i = copied.n_freq_i;

  }
  return *this;
    
}
  
object_props_dbl::object_props_dbl(object_props_dbl && moved){

  int i;

  NOvox = moved.NOvox; 
  ra = moved.ra;
  dec = moved.dec;
  freq = moved.freq;
  ra_i = moved.ra_i;
  dec_i = moved.dec_i;
  freq_i = moved.freq_i;
  tot_intens = moved.tot_intens;
  avg_intens = moved.avg_intens;
  sigma_intens = moved.sigma_intens;
  rms = moved.rms;
  ra_min = moved.ra_min;
  dec_min = moved.dec_min;
  freq_min = moved.freq_min;
  ra_max = moved.ra_max;
  dec_max = moved.dec_max;
  freq_max = moved.freq_max;
  min_intens = moved.min_intens;
  max_intens = moved.max_intens;

  w_max = moved.w_max;
  w20_min = moved.w20_min;
  w20_max = moved.w20_max;
  w50_min = moved.w50_min;
  w50_max = moved.w50_max;
  cw_max = moved.cw_max;
  cw20_min = moved.cw20_min;
  cw20_max = moved.cw20_max;
  cw50_min = moved.cw50_min;
  cw50_max = moved.cw50_max;

  srep_size[0] = moved.srep_size[0];
  srep_size[1] = moved.srep_size[1];
  srep_size[2] = moved.srep_size[2];
  srep_size[3] = moved.srep_size[3];
  srep_size[4] = moved.srep_size[4];
  srep_size[5] = moved.srep_size[5];
  srep_update = moved.srep_update;

  if(moved.srep_grid != NULL){

    srep_grid = moved.srep_grid;
    moved.srep_grid = NULL;

  } else { srep_grid = NULL; }
  if(moved.srep_strings != NULL){

    srep_strings = moved.srep_strings;
    moved.srep_strings = NULL;

  } else { srep_strings = NULL; }
  if(moved.mini_mom0 != NULL){

    mini_mom0 = moved.mini_mom0;
    moved.mini_mom0 = NULL;

  } else { mini_mom0 = NULL; }
  if(moved.mini_RAPV != NULL){

    mini_RAPV = moved.mini_RAPV;
    moved.mini_RAPV = NULL;

  } else { moved.mini_RAPV = NULL; }
  if(moved.mini_DECPV != NULL){

    mini_DECPV = moved.mini_DECPV;
    moved.mini_DECPV = NULL;
    
  } else { mini_DECPV = NULL; }
  if(mini_obj_spec != NULL){

    mini_obj_spec = moved.mini_obj_spec;
    moved.mini_obj_spec = NULL;

  } else { mini_obj_spec = NULL; }
  if(mini_ref_spec != NULL){

    mini_ref_spec = moved.mini_ref_spec;
    moved.mini_ref_spec = NULL;

  } else { mini_ref_spec = NULL; }
  if(moved.mini_vfield != NULL){

    mini_vfield = moved.mini_ref_spec;
    moved.mini_ref_spec = NULL;

  } else { mini_vfield = NULL; }

  // new parameters --- multiple central moment calculations
  p_tot_intens = moved.p_tot_intens;
  n_tot_intens = moved.n_tot_intens;
  p_ra_i = moved.p_tot_intens;
  p_dec_i = moved.p_dec_i;
  p_freq_i = moved.p_freq_i;
  n_ra_i = moved.n_ra_i;
  n_dec_i = moved.n_dec_i;
  n_freq_i = moved.n_freq_i;

}

object_props_dbl & object_props_dbl::operator = (object_props_dbl && moved){

  int i;

  if(this != &moved){

    NOvox = moved.NOvox; 
    ra = moved.ra;
    dec = moved.dec;
    freq = moved.freq;
    ra_i = moved.ra_i;
    dec_i = moved.dec_i;
    freq_i = moved.freq_i;
    tot_intens = moved.tot_intens;
    avg_intens = moved.avg_intens;
    sigma_intens = moved.sigma_intens;
    rms = moved.rms;
    ra_min = moved.ra_min;
    dec_min = moved.dec_min;
    freq_min = moved.freq_min;
    ra_max = moved.ra_max;
    dec_max = moved.dec_max;
    freq_max = moved.freq_max;
    min_intens = moved.min_intens;
    max_intens = moved.max_intens;
    
    w_max = moved.w_max;
    w20_min = moved.w20_min;
    w20_max = moved.w20_max;
    w50_min = moved.w50_min;
    w50_max = moved.w50_max;
    cw_max = moved.cw_max;
    cw20_min = moved.cw20_min;
    cw20_max = moved.cw20_max;
    cw50_min = moved.cw50_min;
    cw50_max = moved.cw50_max;
    
    srep_size[0] = moved.srep_size[0];
    srep_size[1] = moved.srep_size[1];
    srep_size[2] = moved.srep_size[2];
    srep_size[3] = moved.srep_size[3];
    srep_size[4] = moved.srep_size[4];
    srep_size[5] = moved.srep_size[5];
    srep_update = moved.srep_update;
    
    if(moved.srep_grid != NULL){
      
      srep_grid = moved.srep_grid;
      moved.srep_grid = NULL;
      
    } else { srep_grid = NULL; }
    if(moved.srep_strings != NULL){
      
      srep_strings = moved.srep_strings;
      moved.srep_strings = NULL;
      
    } else { srep_strings = NULL; }
    if(moved.mini_mom0 != NULL){
      
      mini_mom0 = moved.mini_mom0;
      moved.mini_mom0 = NULL;
      
    } else { mini_mom0 = NULL; }
    if(moved.mini_RAPV != NULL){
      
      mini_RAPV = moved.mini_RAPV;
      moved.mini_RAPV = NULL;
      
    } else { moved.mini_RAPV = NULL; }
    if(moved.mini_DECPV != NULL){
      
      mini_DECPV = moved.mini_DECPV;
      moved.mini_DECPV = NULL;
      
    } else { mini_DECPV = NULL; }
    if(mini_obj_spec != NULL){
      
      mini_obj_spec = moved.mini_obj_spec;
      moved.mini_obj_spec = NULL;
      
    } else { mini_obj_spec = NULL; }
    if(mini_ref_spec != NULL){
      
      mini_ref_spec = moved.mini_ref_spec;
      moved.mini_ref_spec = NULL;
      
    } else { mini_ref_spec = NULL; }
    if(moved.mini_vfield != NULL){
      
      mini_vfield = moved.mini_vfield;
      moved.mini_vfield = NULL;
      
    } else { mini_vfield = NULL; }
    
    // new parameters --- multiple central moment calculations
    p_tot_intens = moved.p_tot_intens;
    n_tot_intens = moved.n_tot_intens;
    p_ra_i = moved.p_tot_intens;
    p_dec_i = moved.p_dec_i;
    p_freq_i = moved.p_freq_i;
    n_ra_i = moved.n_ra_i;
    n_dec_i = moved.n_dec_i;
    n_freq_i = moved.n_freq_i;

  }
  return *this;

}

object_props_dbl::~object_props_dbl(){ 

  if(srep_grid != NULL){ delete [] srep_grid; srep_grid = NULL; } 
  if(srep_strings != NULL){ delete [] srep_strings; srep_strings = NULL; }
  if(mini_mom0 != NULL){ delete [] mini_mom0; mini_mom0 = NULL; }
  if(mini_RAPV != NULL){ delete [] mini_RAPV; mini_RAPV = NULL; }
  if(mini_DECPV != NULL){ delete [] mini_DECPV; mini_DECPV = NULL; }
  if(mini_obj_spec != NULL){ delete [] mini_obj_spec; mini_obj_spec = NULL; }
  if(mini_ref_spec != NULL){ delete [] mini_ref_spec; mini_ref_spec = NULL; }
  if(mini_vfield != NULL){ delete [] mini_vfield; mini_vfield = NULL; }

}

bool object_props_dbl::operator == (const object_props_dbl & compareTo){

  int i;

  if(this == &compareTo){ return true; }

  if(NOvox != compareTo.NOvox){ return false; } 
  if(ra != compareTo.ra){ return false; }
  if(dec != compareTo.dec){ return false; }
  if(freq != compareTo.freq){ return false; }
  if(ra_i != compareTo.ra_i){ return false; }
  if(dec_i != compareTo.dec_i){ return false; }
  if(freq_i != compareTo.freq_i){ return false; }
  if(tot_intens != compareTo.tot_intens){ return false; }
  if(avg_intens != compareTo.avg_intens){ return false; }
  if(sigma_intens != compareTo.sigma_intens){ return false; }
  if(rms != compareTo.rms){ return false; }
  if(ra_min != compareTo.ra_min){ return false; }
  if(dec_min != compareTo.dec_min){ return false; }
  if(freq_min != compareTo.freq_min){ return false; }
  if(ra_max != compareTo.ra_max){ return false; }
  if(dec_max != compareTo.dec_max){ return false; }
  if(freq_max != compareTo.freq_max){ return false; }
  if(min_intens != compareTo.min_intens){ return false; }
  if(max_intens != compareTo.max_intens){ return false; }
  
  if(w_max != compareTo.w_max){ return false; }
  if(w20_min != compareTo.w20_min){ return false; }
  if(w20_max != compareTo.w20_max){ return false; }
  if(w50_min != compareTo.w50_min){ return false; }
  if(w50_max != compareTo.w50_max){ return false; }
  if(cw_max != compareTo.cw_max){ return false; }
  if(cw20_min != compareTo.cw20_min){ return false; }
  if(cw20_max != compareTo.cw20_max){ return false; }
  if(cw50_min != compareTo.cw50_min){ return false; }
  if(cw50_max != compareTo.cw50_max){ return false; }
  
  if(srep_size[0] != compareTo.srep_size[0]){ return false; }
  if(srep_size[1] != compareTo.srep_size[1]){ return false; }
  if(srep_size[2] != compareTo.srep_size[2]){ return false; }
  if(srep_size[3] != compareTo.srep_size[3]){ return false; }
  if(srep_size[4] != compareTo.srep_size[4]){ return false; }
  if(srep_size[5] != compareTo.srep_size[5]){ return false; }
  if(srep_update != compareTo.srep_update){ return false; }
  
  // new parameters --- multiple central moment calculations
  if(p_tot_intens != compareTo.p_tot_intens){ return false; }
  if(n_tot_intens != compareTo.n_tot_intens){ return false; }
  if(p_ra_i != compareTo.p_tot_intens){ return false; }
  if(p_dec_i != compareTo.p_dec_i){ return false; }
  if(p_freq_i != compareTo.p_freq_i){ return false; }
  if(n_ra_i != compareTo.n_ra_i){ return false; }
  if(n_dec_i != compareTo.n_dec_i){ return false; }
  if(n_freq_i != compareTo.n_freq_i){ return false; }

  if(compareTo.srep_grid != NULL){
    
    for(i = 0; i < ((srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1)); ++i){ if(srep_grid[i] != compareTo.srep_grid[i]){ return false; } }    
    
  } else { if(srep_grid != NULL){ return false; } }
  if(compareTo.srep_strings != NULL){
    
    for(i = 0; i < (2 * srep_grid[((srep_size[3] - srep_size[2] + 1) * (srep_size[1] - srep_size[0] + 1))]); ++i){ if(srep_strings[i] != compareTo.srep_strings[i]){ return false; } }
    
  } else { if(srep_strings != NULL){ return false; } }
  if(compareTo.mini_mom0 != NULL){
    
    for(i = 0; i < ((srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1)); ++i){ if(mini_mom0[i] != compareTo.mini_mom0[i]){ return false; } }
    
  } else { if(mini_mom0 != NULL){ return false; } }
  if(compareTo.mini_RAPV != NULL){
    
    for(i = 0; i < ((srep_size[1] - srep_size[0] + 1) * (srep_size[5] - srep_size[4] + 1)); ++i){ if(mini_RAPV[i] != compareTo.mini_RAPV[i]){ return false; } }
    
  } else { if(compareTo.mini_RAPV != NULL){ return false; } }
  if(compareTo.mini_DECPV != NULL){
    
    for(i = 0; i < ((srep_size[3] - srep_size[2] + 1) * (srep_size[5] - srep_size[4] + 1)); ++i){ if(mini_DECPV[i] != compareTo.mini_DECPV[i]){ return false; } }
    
  } else { if(mini_DECPV != NULL){ return false; } }
  if(mini_obj_spec != NULL){
    
    for(i = 0; i < ((srep_size[5] - srep_size[4] + 1)); ++i){ if(mini_obj_spec[i] != compareTo.mini_obj_spec[i]){ return false; } }
    
  } else { if(mini_obj_spec != NULL){ return false; } }
  if(mini_ref_spec != NULL){
    
    if((compareTo.srep_size[5] - compareTo.srep_size[4] + 1) >= 10){
      
      for(i = 0; i < ((srep_size[5] - srep_size[4] + 1)); ++i){ if(mini_ref_spec[i] != compareTo.mini_ref_spec[i]){ return false; } }    
      
    } else {
      
      for(i = 0; i < ((srep_size[5] - srep_size[4] + 11)); ++i){ if(mini_ref_spec[i] != compareTo.mini_ref_spec[i]){ return false; } }    
      
    }
    
  } else { if(mini_ref_spec != NULL){ return false; } }
  if(compareTo.mini_vfield != NULL){
    
    for(i = 0; i < ((srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1)); ++i){ if(mini_vfield[i] != compareTo.mini_vfield[i]){ return false; } }

  } else { if(mini_vfield != NULL){ return false; } }
  
  return true;

}

bool object_props_dbl::operator != (const object_props_dbl & compareTo){

  if(this == &compareTo){ 

    return false; 

  } else {

    return (!(*this == compareTo));

  } 

}

void object_props_dbl::operator += (object_props_dbl & added){

  this->AddObject(added);

}

object_props_dbl object_props_dbl::operator + (object_props_dbl & added){

  object_props_dbl temp_result(*this);
  temp_result.AddObject(added);
  return temp_result;

}

void object_props_dbl::ReInit(){

  NOvox = 0;
  ra = dec = freq = ra_i = dec_i = freq_i = tot_intens = avg_intens = sigma_intens = rms = 0.0;
  ra_min = dec_min = freq_min = min_intens = 1E10;
  ra_max = dec_max = freq_max = max_intens = -1E10;
  w_max = w20_min = w50_min = w20_max = w50_max = -1E10;
  cw_max = cw20_min = cw50_min = cw20_max = cw50_max = -1E10;

}

void object_props_dbl::AddVoxel(int value){ NOvox+=value; }

void object_props_dbl::AddRa(double value){ ra+=value; }

void object_props_dbl::AddDec(double value){ dec+=value; }

void object_props_dbl::AddFreq(double value){ freq+=value; }

void object_props_dbl::AddRA_i(double pos, double value){ 
  ra_i+=(pos * value);
  if(value >= 0.0){ p_ra_i+=(pos * value); } else { n_ra_i+=(pos * value); }
}

void object_props_dbl::AddDec_i(double pos, double value){ 
  dec_i+=value; 
  if(value >= 0.0){ p_dec_i+=(pos * value); } else { n_dec_i+=(pos * value); }
}

void object_props_dbl::AddFreq_i(double pos, double value){ 
  freq_i+=value; 
  if(value >= 0.0){ p_freq_i+=(pos * value); } else { n_freq_i+=(pos * value); }
}

void object_props_dbl::AddTotIntens(double value){ 
  tot_intens+=value; 
  if(value >= 0.0){ p_tot_intens+=value; } else { n_tot_intens+=value; }
}

void object_props_dbl::AddAvgIntens(double value){ avg_intens+=value; }

void object_props_dbl::AddSigmaItens(double value){ sigma_intens+=(value * value); }

void object_props_dbl::AdjustRange(double value){
    
  if(value <= min_intens){ min_intens = value; }
  if(value >= max_intens){ max_intens = value; }

}

void object_props_dbl::CalcProps(){
  
  double dummy, flip;
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
    for(g = 0; g < (srep_size[5] - srep_size[4] + 1); g++){ if((flip * mini_obj_spec[g]) >= dummy){ dummy = (flip * mini_obj_spec[g]); w_max = (double) (g + srep_size[4]); } }
    w20_min = (double) (srep_size[4]);
    for(g = 0; g < (srep_size[5] - srep_size[4] + 1); g++){ 
      
      if((flip * mini_obj_spec[g]) >= (0.2 * dummy)){ 
	
	if(g > 0){
	  
	  w20_min = (double) (g + srep_size[4]) - 1.0 + (((0.2 * dummy) - (flip * mini_obj_spec[(g - 1)])) / (flip * (mini_obj_spec[g] - mini_obj_spec[(g - 1)]))); 
	  
	} else {
	  
	  w20_min = (double) (g + srep_size[4]); 
	  
	}
	
	break; 
	
      } 
      
    }
    w20_max = (double) (srep_size[5]);
    for(g = srep_size[5] - srep_size[4]; g >= 0; g--){ 
      
      if((flip * mini_obj_spec[g]) >= (0.2 * dummy)){ 
	
	if(g < (srep_size[5] - srep_size[4])){
	  
	  w20_max = (double) (g + srep_size[4]) + 1.0 - (((0.2 * dummy) - (flip * mini_obj_spec[(g + 1)])) / (flip * (mini_obj_spec[g] - mini_obj_spec[(g + 1)])));
	  
	} else {
	  
	  w20_max = (double) (g + srep_size[4]); 
	  
	}
	
	break; 
	
      } 
	
    }
    w50_min = (double) (srep_size[4]);
    for(g = (int) floorf((w20_min - (double) srep_size[4])); g < (srep_size[5] - srep_size[4] + 1); g++){ 
      
      if((flip * mini_obj_spec[g]) >= (0.5 * dummy)){ 
	
	if(g > 0){
	  
	  w50_min = (double) (g + srep_size[4]) - 1.0 + (((0.5 * dummy) - (flip * mini_obj_spec[(g - 1)])) / (flip * (mini_obj_spec[g] - mini_obj_spec[(g - 1)]))); 
	  
	} else {
	  
	  w50_min = (double) (g + srep_size[4]); 
	  
	}	
	
	break; 
	
      } 
      
    }
    w50_max = (double) (srep_size[5]);
    for(g = (int) floorf((w20_max - (double) srep_size[4])); g >= 0; g--){ 
      
      if((flip * mini_obj_spec[g]) >= (0.5 * dummy)){ 
	
	if(g < (srep_size[5] - srep_size[4])){
	  
	  w50_max = (double) (g + srep_size[4]) + 1.0 - (((0.5 * dummy) - (flip * mini_obj_spec[(g + 1)])) / (flip * (mini_obj_spec[g] - mini_obj_spec[(g + 1)])));
	  
	} else {
	  
	  w50_max = (double) (g + srep_size[4]); 
	  
	}
	
	break; 
	
      } 
      
    }
    
    // b. calculate W_50 and W_20 according to c.f.d. values that correspond to FWHM and 1/5th of peak for a gaussian profile
    cw_max = cw20_min = cw50_min = dummy = 0.0;
    for(g = 0; g < (srep_size[5] - srep_size[4] + 1); g++){ 
      
      if(((dummy + (mini_obj_spec[g]/tot_intens)) >= 0.5) && (dummy < 0.5)){ 
	
	if(g > 0){
	  
	  cw_max = (double) (g + srep_size[4]) - 1.0 + ((0.5 - dummy) * tot_intens / mini_obj_spec[g]);
	  
	} else {
	  
	  cw_max = (double) (g + srep_size[4]);
	  
	}
	
      } 
      
      if(((dummy + (mini_obj_spec[g]/tot_intens)) >= 0.036397) && (dummy < 0.036397)){ 
	
	if(g > 0){
	  
	  cw20_min = (double) (g + srep_size[4]) - 1.0 + ((0.036397 - dummy) * tot_intens / mini_obj_spec[g]);
	  
	} else {
	  
	  cw20_min = (double) (g + srep_size[4]);
	  
	}	  
	
      } 
      
      if(((dummy + (mini_obj_spec[g]/tot_intens)) >= 0.119516) && (dummy < 0.119516)){ 
	
	if(g > 0){
	  
	  cw50_min = (double) (g + srep_size[4]) - 1.0 + ((0.119516 - dummy) * tot_intens / mini_obj_spec[g]);
	  
	} else {
	  
	  cw50_min = (double) (g + srep_size[4]);
	  
	}	  
	
      } 
      
      dummy+=(mini_obj_spec[g]/tot_intens);
      
    }
    cw20_max = cw50_max = dummy = 1.0;
    for(g = srep_size[5] - srep_size[4]; ((g >= 0) && (g > (cw_max - 1 - srep_size[4]))); g--){ 
      
      if(((dummy - (mini_obj_spec[g]/tot_intens)) <= 0.880484) && (dummy > 0.880484)){ 
	
	if(g < (srep_size[5] - srep_size[4])){
	  
	  cw50_max = (double) (g + srep_size[4]) + 1.0 - ((dummy - 0.880484) * tot_intens / mini_obj_spec[g]);
	  
	} else {
	  
	  cw50_max = (double) (g + srep_size[4]);
	  
	}
	
      } 
      
      if(((dummy - (mini_obj_spec[g]/tot_intens)) <= 0.963603) && (dummy > 0.963603)){ 
	
	if(g < (srep_size[5] - srep_size[4])){
	  
	  cw20_max = (double) (g + srep_size[4]) + 1.0 - ((dummy - 0.963603) * tot_intens / mini_obj_spec[g]);
	  
	} else {
	  
	  cw20_max = (double) (g + srep_size[4]);
	  
	}
	
      } 
      
      dummy-=(mini_obj_spec[g]/tot_intens);
      
    }
    
  } else {
    
    w_max = w50_min = w50_max = w20_min = w20_max = cw_max = cw50_min = cw50_max = cw20_min = cw20_max = (double) srep_size[4];
    
  }

}

void object_props_dbl::ShowAll_file_WCS(int id, std::fstream& output_file, int cat_mode, double wcs_vals[6]){

  if(cat_mode > 0){
    
    object_props_dbl::ShowProps_file_WCS(id,output_file,wcs_vals);
    object_props_dbl::ShowSrep_file(output_file);
    output_file << std::endl;
    
  } else {
    
    object_props_dbl::ShowProps_file_WCS(id,output_file,wcs_vals); 
    output_file << std::endl;
    
  }

}

void object_props_dbl::ShowAll_file(int id, std::fstream& output_file, int cat_mode){

  if(cat_mode > 0){
    
    object_props_dbl::ShowProps_file(id,output_file);
    object_props_dbl::ShowSrep_file(output_file);
    output_file << std::endl;
    
  } else {
    
    object_props_dbl::ShowProps_file(id,output_file); 
    output_file << std::endl;
    
  }

}

void object_props_dbl::ShowProps_file_WCS(int id, std::fstream& output_file, double wcs_vals[6]){
    
  output_file << id << " " << NOvox << " " << wcs_vals[0] << " " << wcs_vals[1] << " " << wcs_vals[2] << " " << wcs_vals[3] << " " << wcs_vals[4] << " " << wcs_vals[5] << " " << ra << " " << dec << " " << freq << " " << ra_i << " " << dec_i << " " << freq_i << " " << tot_intens << " " << avg_intens << " " << sigma_intens << " " << rms << " " << min_intens << " " << max_intens << " " << ra_min << " " << ra_max << " " << dec_min << " " << dec_max << " " << freq_min << " " << freq_max << " " << w_max << " " << w50_min << " " << w50_max << " " << w20_min << " " << w20_max << " " << cw_max << " " << cw50_min << " " << cw50_max << " " << cw20_min << " " << cw20_max << " " << std::flush;
  
}

void object_props_dbl::ShowProps_file(int id, std::fstream& output_file){
    
  output_file << id << " " << NOvox << " " << ra << " " << dec << " " << freq << " " << ra_i << " " << dec_i << " " << freq_i << " " << tot_intens << " " << avg_intens << " " << sigma_intens << " " << rms << " " << min_intens << " " << max_intens << " " << ra_min << " " << ra_max << " " << dec_min << " " << dec_max << " " << freq_min << " " << freq_max << " " << w_max << " " << w50_min << " " << w50_max << " " << w20_min << " " << w20_max << " " << cw_max << " " << cw50_min << " " << cw50_max << " " << cw20_min << " " << cw20_max << " " << std::flush;
  
}

void object_props_dbl::ShowSrep_file(std::fstream& output_file){

  int g;
  
  output_file << ": " << std::flush;
  for(g = 0; g < (1 + ((srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1))); g++){ output_file << srep_grid[g] << " "; }
  output_file << ": " << std::flush;
  for(g = 0; g < (2 * srep_grid[((srep_size[1] - srep_size[0] + 1) * (srep_size[3] - srep_size[2] + 1))]); g++){ output_file << srep_strings[g] << " "; }
  output_file << std::flush;

}

void object_props_dbl::AdjustRArange(double value){

  if(value <= ra_min){ ra_min = value; }
  if(value >= ra_max){ ra_max = value; }

}

void object_props_dbl::AdjustDECrange(double value){

  if(value <= dec_min){ dec_min = value; }
  if(value >= dec_max){ dec_max = value; }

}

void object_props_dbl::AdjustFREQrange(double value){

  if(value <= freq_min){ freq_min = value; }
  if(value >= freq_max){ freq_max = value; }

}

int object_props_dbl::ShowArea(){ return (int) ((fabs(ra_max - ra_min) + 1) * (fabs(dec_max - dec_min) + 1)); }

int object_props_dbl::ShowVoxels(){ return NOvox; }

int object_props_dbl::ShowRArange(){ return (int) (fabs(ra_max - ra_min) + 1); }

int object_props_dbl::ShowDECrange(){ return (int) (fabs(dec_max - dec_min) + 1); }

int object_props_dbl::ShowFREQrange(){ return (int) (fabs(freq_max - freq_min) + 1); }

int object_props_dbl::GetRAmin(){ return ra_min; }

int object_props_dbl::GetRAmax(){ return ra_max; }

int object_props_dbl::GetDECmin(){ return dec_min; }

int object_props_dbl::GetDECmax(){ return dec_max; }

int object_props_dbl::GetFREQmin(){ return freq_min; }

int object_props_dbl::GetFREQmax(){ return freq_max; }

double object_props_dbl::GetRA(){ return ra; }

double object_props_dbl::GetDEC(){ return dec; }

double object_props_dbl::GetFREQ(){ return freq; }

double object_props_dbl::GetRAi(){ return ra_i; }

double object_props_dbl::GetDECi(){ return dec_i; }

double object_props_dbl::GetFREQi(){ return freq_i; }

double object_props_dbl::GetTI(){ return tot_intens; }

float object_props_dbl::GetRAi_p(){ return p_ra_i; }

float object_props_dbl::GetDECi_p(){ return p_dec_i; }

float object_props_dbl::GetFREQi_p(){ return p_freq_i; }

float object_props_dbl::GetTI_p(){ return p_tot_intens; }

float object_props_dbl::GetRAi_n(){ return n_ra_i; }

float object_props_dbl::GetDECi_n(){ return n_dec_i; }

float object_props_dbl::GetFREQi_n(){ return n_freq_i; }

float object_props_dbl::GetTI_n(){ return n_tot_intens; }

double object_props_dbl::GetSigmaI(){ return sigma_intens; }

double object_props_dbl::GetAvgI(){ return avg_intens; }

double object_props_dbl::GetMinI(){ return min_intens; }

double object_props_dbl::GetMaxI(){ return max_intens; }

void object_props_dbl::Set_w_max(double value){ w_max = value;}

double object_props_dbl::Get_w_max(){ return w_max; }

void object_props_dbl::Set_w20_min(double value){ w20_min = value; }

double object_props_dbl::Get_w20_min(){ return w20_min; }

void object_props_dbl::Set_w20_max(double value){ w20_max = value; }

double object_props_dbl::Get_w20_max(){ return w20_max; }

void object_props_dbl::Set_w50_min(double value){ w50_min = value; }

double object_props_dbl::Get_w50_min(){ return w50_min; }

void object_props_dbl::Set_w50_max(double value){ w50_max = value; }

double object_props_dbl::Get_w50_max(){ return w50_max; }

void object_props_dbl::Set_cw_max(double value){ cw_max = value; }

double object_props_dbl::Get_cw_max(){ return cw_max; }

void object_props_dbl::Set_cw20_min(double value){ cw20_min = value; }

double object_props_dbl::Get_cw20_min(){ return cw20_min; }

void object_props_dbl::Set_cw20_max(double value){ cw20_max = value; }

double object_props_dbl::Get_cw20_max(){ return cw20_max; }

void object_props_dbl::Set_cw50_min(double value){ cw50_min = value; }

double object_props_dbl::Get_cw50_min(){ return cw50_min; }

void object_props_dbl::Set_cw50_max(double value){ cw50_max = value; }

double object_props_dbl::Get_cw50_max(){ return cw50_max; }

void object_props_dbl::Set_srep_update(int value){ srep_update = value; }

int object_props_dbl::Get_srep_update(){ return srep_update; }

void object_props_dbl::Set_srep_size(int index, int value){ srep_size[index] = value; }

int object_props_dbl::Get_srep_size(int index){ return srep_size[index]; }

void object_props_dbl::Create_srep_grid(int value){ while(srep_grid == NULL){ srep_grid = new int[value]; } }

void object_props_dbl::Set_srep_grid(int index, int value){ srep_grid[index] = value; }

int object_props_dbl::Get_srep_grid(int index){ return srep_grid[index]; }

void object_props_dbl::Free_srep_grid(){ if(srep_grid != NULL){ delete [] srep_grid; srep_grid = NULL; } }

void object_props_dbl::Create_srep_strings(int value){ while(srep_strings == NULL){ srep_strings = new int[value]; } }

void object_props_dbl::Set_srep_strings(int index, int value){ srep_strings[index] = value; }

int object_props_dbl::Get_srep_strings(int index){ return srep_strings[index]; }

void object_props_dbl::Free_srep_strings(){ if(srep_strings != NULL){ delete [] srep_strings; srep_strings = NULL; } }

void object_props_dbl::Create_mom0(int value){ while(mini_mom0 == NULL){ mini_mom0 = new double[value]; } }

void object_props_dbl::Set_mom0(int index, double value){ mini_mom0[index] = value; }

void object_props_dbl::Add_mom0(int index, double value){ mini_mom0[index]+=value; }

double object_props_dbl::Get_mom0(int index){ return mini_mom0[index]; }

void object_props_dbl::Free_mom0(){ if(mini_mom0 != NULL){ delete [] mini_mom0; mini_mom0 = NULL; } }

void object_props_dbl::Create_RAPV(int value){ while(mini_RAPV == NULL){ mini_RAPV = new double[value]; } }

void object_props_dbl::Set_RAPV(int index, double value){ mini_RAPV[index] = value; }

void object_props_dbl::Add_RAPV(int index, double value){ mini_RAPV[index]+=value; }

double object_props_dbl::Get_RAPV(int index){ return mini_RAPV[index]; }

void object_props_dbl::Free_RAPV(){ if(mini_RAPV != NULL){ delete [] mini_RAPV; mini_RAPV = NULL; } }

void object_props_dbl::Create_DECPV(int value){ while(mini_DECPV == NULL){ mini_DECPV = new double[value]; } }

void object_props_dbl::Set_DECPV(int index, double value){ mini_DECPV[index] = value; }

void object_props_dbl::Add_DECPV(int index, double value){ mini_DECPV[index]+=value; }

double object_props_dbl::Get_DECPV(int index){ return mini_DECPV[index]; }

void object_props_dbl::Free_DECPV(){ if(mini_DECPV != NULL){ delete [] mini_DECPV; mini_DECPV = NULL; } }

void object_props_dbl::Create_ospec(int value){ while(mini_obj_spec == NULL){ mini_obj_spec = new double[value]; } }

void object_props_dbl::Set_ospec(int index, double value){ mini_obj_spec[index] = value; }

void object_props_dbl::Add_ospec(int index, double value){ mini_obj_spec[index]+=value; }

double object_props_dbl::Get_ospec(int index){ return mini_obj_spec[index]; }

void object_props_dbl::Free_ospec(){ if(mini_obj_spec != NULL){ delete [] mini_obj_spec; mini_obj_spec = NULL; } }

void object_props_dbl::Create_rspec(int value){ while(mini_ref_spec == NULL){ mini_ref_spec = new double[value]; } }

void object_props_dbl::Set_rspec(int index, double value){ mini_ref_spec[index] = value; }

void object_props_dbl::Add_rspec(int index, double value){ mini_ref_spec[index]+=value; }

double object_props_dbl::Get_rspec(int index){ return mini_ref_spec[index]; }

void object_props_dbl::Free_rspec(){ if(mini_ref_spec != NULL){ delete [] mini_ref_spec; mini_ref_spec = NULL; } }

void object_props_dbl::Create_vfield(int value){ while(mini_vfield == NULL){ mini_vfield = new double[value]; } }

void object_props_dbl::Set_vfield(int index, double value){ mini_vfield[index] = value; }

void object_props_dbl::Add_vfield(int index, double value){ mini_vfield[index]+=value; }

double object_props_dbl::Get_vfield(int index){ return mini_vfield[index]; }

void object_props_dbl::Free_vfield(){ if(mini_vfield != NULL){ delete [] mini_vfield; mini_vfield = NULL; } }

void object_props_dbl::ReInit_srep(){

  if(srep_grid != NULL){ delete [] srep_grid; }
  srep_grid = NULL;
  if(srep_strings != NULL){ delete [] srep_strings; }
  srep_strings = NULL;
  
}

void object_props_dbl::ReInit_mini(){

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

void object_props_dbl::ReInit_size(){

  srep_size[0] = srep_size[1] = srep_size[2] = srep_size[3] = srep_size[4] = srep_size[5] = -99;

}

void object_props_dbl::AddPoint(double x_pos, double y_pos, double z_pos, double value){

  // update bounding box
  this->AdjustRArange(x_pos);
  this->AdjustDECrange(y_pos);
  this->AdjustFREQrange(z_pos);
	    
  // update voxel count
  this->AddVoxel(1);
	    
  // call member functions for object_props_dbl to adjust appropriate values
  this->AddRa(x_pos);
  this->AddDec(y_pos);
  this->AddFreq(z_pos);
  this->AddRA_i(x_pos,value);
  this->AddDec_i(y_pos,value);
  this->AddFreq_i(z_pos,value);
  this->AddTotIntens(value);
  this->AddAvgIntens(value);
  this->AddSigmaItens(value);
  this->AdjustRange(value);
  
  // change sparse_reps_update
  this->Set_srep_update(1);


}

void object_props_dbl::AddObject(object_props_dbl & merged){

  int j,k,g,sx,sy,sz;
  vector<int> temp_sparse_reps_grid, temp_sparse_reps_strings;
  vector<double> temp_mom0, temp_RAPV, temp_DECPV, temp_ref_spec, temp_obj_spec, temp_vfield;
    
  if(this != &merged){ 

    // update the properties of this object with the object being merged in
    this->AddVoxel(merged.ShowVoxels());
    this->AdjustRArange(merged.GetRAmin());
    this->AdjustRArange(merged.GetRAmax());
    this->AdjustDECrange(merged.GetDECmin());
    this->AdjustDECrange(merged.GetDECmax());
    this->AdjustFREQrange(merged.GetFREQmin());
    this->AdjustFREQrange(merged.GetFREQmax());
    this->AddRa(merged.GetRA());
    this->AddDec(merged.GetDEC());
    this->AddFreq(merged.GetFREQ());
    this->AddRA_i(1.0,merged.GetRAi());
    this->AddDec_i(1.0,merged.GetDECi());
    this->AddFreq_i(1.0,merged.GetFREQi());
    this->AddTotIntens(merged.GetTI());
    this->AddAvgIntens(merged.GetAvgI());
    this->AddSigmaItens(merged.GetSigmaI());
    this->AdjustRange(merged.GetMinI());		
    this->AdjustRange(merged.GetMaxI());	

    // initialise temporary working arrays
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

    // combine the sparse representations using the temporary arrays --- if they exist
    if((merged.Get_srep_update() != 0) && (merged.Get_srep_size(0) >= 0)){
		  
      // write the existing object's sparse representations into temporary arrays and initialise temporary arrays at the same time, provided it exists
      
      // a. grid
      temp_sparse_reps_grid.resize(0);
      for(g = 0; g < (1 + ((this->GetRAmax() - this->GetRAmin() + 1) * (this->GetDECmax() - this->GetDECmin() + 1))); g++){ temp_sparse_reps_grid.push_back(0); }
      
      if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	
	for(sx = 0; sx < (this->Get_srep_size(1) - this->Get_srep_size(0) + 1); sx++){
	  for(sy = 0; sy < (this->Get_srep_size(3) - this->Get_srep_size(2) + 1); sy++){
	    
	    temp_sparse_reps_grid[(((sy + this->Get_srep_size(2) - this->GetDECmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx + this->Get_srep_size(0) - this->GetRAmin())] = this->Get_srep_grid(((sy * (this->Get_srep_size(1) - this->Get_srep_size(0) + 1)) + sx + 1)) - this->Get_srep_grid(((sy * (this->Get_srep_size(1) - this->Get_srep_size(0) + 1)) + sx));
	    
	  }
	  
	}
	
      }
      
      // b. mom-0
      temp_mom0.resize(0);
      for(g = 0; g < ((this->GetRAmax() - this->GetRAmin() + 1) * (this->GetDECmax() - this->GetDECmin() + 1)); g++){ temp_mom0.push_back(0.0); }
      
      if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	
	for(sx = 0; sx < (this->Get_srep_size(1) - this->Get_srep_size(0) + 1); sx++){
	  
	  for(sy = 0; sy < (this->Get_srep_size(3) - this->Get_srep_size(2) + 1); sy++){
	    
	    temp_mom0[(((sy + this->Get_srep_size(2) - this->GetDECmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx + this->Get_srep_size(0) - this->GetRAmin())]+=this->Get_mom0(((sy * (this->Get_srep_size(1) - this->Get_srep_size(0) + 1)) + sx));
	    
	  }
	  
	}
	
      }
      
      // c. RAPV
      temp_RAPV.resize(0);
      for(g = 0; g < ((this->GetRAmax() - this->GetRAmin() + 1) * (this->GetFREQmax() - this->GetFREQmin() + 1)); g++){ temp_RAPV.push_back(0.0); }
      
      if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	
	for(sx = 0; sx < (this->Get_srep_size(1) - this->Get_srep_size(0) + 1); sx++){
	  
	  for(sz = 0; sz < (this->Get_srep_size(5) - this->Get_srep_size(4) + 1); sz++){
	    
	    temp_RAPV[(((sz + this->Get_srep_size(4) - this->GetFREQmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx + this->Get_srep_size(0) - this->GetRAmin())]+=this->Get_RAPV(((sz * (this->Get_srep_size(1) - this->Get_srep_size(0) + 1)) + sx));
	    
	  }
	  
	}
	
      }
      
      // d. DECPV
      temp_DECPV.resize(0);
      for(g = 0; g < ((this->GetDECmax() - this->GetDECmin() + 1) * (this->GetFREQmax() - this->GetFREQmin() + 1)); g++){ temp_DECPV.push_back(0.0); }
      
      if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	
	for(sy = 0; sy < (this->Get_srep_size(3) - this->Get_srep_size(2) + 1); sy++){
	  
	  for(sz = 0; sz < (this->Get_srep_size(5) - this->Get_srep_size(4) + 1); sz++){
	    
	    temp_DECPV[(((sz + this->Get_srep_size(4) - this->GetFREQmin()) * (this->GetDECmax() - this->GetDECmin() + 1)) + sy + this->Get_srep_size(2) - this->GetDECmin())]+=this->Get_DECPV(((sz * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)) + sy));
	    
	  }
	  
	}
	
      }
      
      // e. ref_spec
      if((this->Get_srep_size(5) - this->Get_srep_size(4) + 1) >= 10){
	
	temp_ref_spec.resize(0);
	for(g = 0; g < (2 * (this->GetFREQmax() - this->GetFREQmin() + 1)); g++){ temp_ref_spec.push_back(0.0); }
	
	if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	  
	  for(sz = 0; sz < (2 * (this->Get_srep_size(5) - this->Get_srep_size(4) + 1)); sz++){ temp_ref_spec[(sz + this->Get_srep_size(4) - this->GetFREQmin())]+=this->Get_rspec(sz); } 
	  
	}
	
      } else {
	
	temp_ref_spec.resize(0);
	for(g = 0; g < (this->GetFREQmax() - this->GetFREQmin() + 11); g++){ temp_ref_spec.push_back(0.0); }
	
	if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	  
	  for(sz = 0; sz < (this->Get_srep_size(5) - this->Get_srep_size(4) + 11); sz++){ temp_ref_spec[(sz + this->Get_srep_size(4) - this->GetFREQmin())]+=this->Get_rspec(sz); } 		    
	  
	}
	
      }
      
      // f. obj_spec
      temp_obj_spec.resize(0);
      for(g = 0; g < (this->Get_srep_size(5) - this->Get_srep_size(4) + 1); g++){ temp_obj_spec.push_back(0.0); }		  
      if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	
	for(sz = 0; sz < (this->Get_srep_size(5) - this->Get_srep_size(4) + 1); sz++){
	  
	  temp_obj_spec[(sz + this->Get_srep_size(4) - this->GetFREQmin())]+=this->Get_ospec(sz);
	  
	}
	
      }
      
      // g. vfield
      temp_vfield.resize(0);
      for(g = 0; g < ((this->GetRAmax() - this->GetRAmin() + 1) * (this->GetDECmax() - this->GetDECmin() + 1)); g++){ temp_vfield.push_back(0.0); }
      
      if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	
	for(sx = 0; sx < (this->Get_srep_size(1) - this->Get_srep_size(0) + 1); sx++){
	  
	  for(sy = 0; sy < (this->Get_srep_size(3) - this->Get_srep_size(2) + 1); sy++){
	    
	    temp_vfield[(((sy + this->Get_srep_size(2) - this->GetDECmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx + this->Get_srep_size(0) - this->GetRAmin())]+=this->Get_vfield(((sy * (this->Get_srep_size(1) - this->Get_srep_size(0) + 1)) + sx));
	    
	  }
	  
	}
	
      }
      
      // write the merged object's sparse representations into temporary arrays
      
      // a. grid
      for(sx = 0; sx < (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1); sx++){
	
	for(sy = 0; sy < (merged.Get_srep_size(3) - merged.Get_srep_size(2) + 1); sy++){
	  
	  temp_sparse_reps_grid[(((sy + merged.Get_srep_size(2) - this->GetDECmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx + merged.Get_srep_size(0) - this->GetRAmin())]+=(merged.Get_srep_grid(((sy * (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1)) + sx + 1)) - merged.Get_srep_grid(((sy * (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1)) + sx)));
	  
	}
	
      }
      
      // b. mom-0
      for(sx = 0; sx < (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1); sx++){
	
	for(sy = 0; sy < (merged.Get_srep_size(3) - merged.Get_srep_size(2) + 1); sy++){
	  
	  temp_mom0[(((sy + merged.Get_srep_size(2) - this->GetDECmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx + merged.Get_srep_size(0) - this->GetRAmin())]+=merged.Get_mom0(((sy * (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1)) + sx));
	  
	}
	
      }
      
      // c. RAPV
      for(sx = 0; sx < (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1); sx++){
	
	for(sz = 0; sz < (merged.Get_srep_size(5) - merged.Get_srep_size(4) + 1); sz++){
	  
	  temp_RAPV[(((sz + merged.Get_srep_size(4) - this->GetFREQmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx + merged.Get_srep_size(0) - this->GetRAmin())]+=merged.Get_RAPV(((sz * (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1)) + sx));
	  
	}
	
      }
      
      // d. DECPV
      for(sy = 0; sy < (merged.Get_srep_size(3) - merged.Get_srep_size(2) + 1); sy++){
	
	for(sz = 0; sz < (merged.Get_srep_size(5) - merged.Get_srep_size(4) + 1); sz++){
	  
	  temp_DECPV[(((sz + merged.Get_srep_size(4) - this->GetFREQmin()) * (this->GetDECmax() - this->GetDECmin() + 1)) + sy + merged.Get_srep_size(2) - this->GetDECmin())]+=merged.Get_DECPV(((sz * (merged.Get_srep_size(3) - merged.Get_srep_size(2) + 1)) + sy));
	  
	}
	
      }
      
      // e. ref_spec
      if((merged.Get_srep_size(5) - merged.Get_srep_size(4) + 1) >= 10){
	
	for(sz = 0; sz < (2 * (merged.Get_srep_size(5) - merged.Get_srep_size(4) + 1)); sz++){ 
	  
	  if((sz + merged.Get_srep_size(4) - merged.GetFREQmin()) >= 0){
	    
	    temp_ref_spec[(sz + merged.Get_srep_size(4) - this->GetFREQmin())]+=merged.Get_rspec(sz); 
	    
	  }
	  
	} 
	
      } else {
	
	for(sz = 0; sz < (merged.Get_srep_size(5) - merged.Get_srep_size(4) + 11); sz++){ 
	  
	  if((sz + merged.Get_srep_size(4) - this->GetFREQmin()) >= 0){
	    
	    temp_ref_spec[(sz + merged.Get_srep_size(4) - this->GetFREQmin())]+=merged.Get_rspec(sz); 
	    
	  }
	  
	} 		    
	
      }
      
      // f. obj_spec
      for(sz = 0; sz < (merged.Get_srep_size(5) - merged.Get_srep_size(4) + 1); sz++){
	
	temp_obj_spec[(sz + merged.Get_srep_size(4) - this->GetFREQmin())]+=merged.Get_ospec(sz);
	
      }
      
      // g. vfield
      for(sx = 0; sx < (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1); sx++){
	
	for(sy = 0; sy < (merged.Get_srep_size(3) - merged.Get_srep_size(2) + 1); sy++){
	  
	  temp_vfield[(((sy + merged.Get_srep_size(2) - this->GetDECmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx + merged.Get_srep_size(0) - this->GetRAmin())]+=merged.Get_vfield(((sy * (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1)) + sx));
	  
	}
	
      }
      
      // convert temp_sparse_reps_grid from differential to cumulative counts using temp_sparse_reps_string as an intermediary
      temp_sparse_reps_strings.resize(1);
      temp_sparse_reps_strings[0] = 0;
      for(g = 1; g < (((this->GetRAmax() - this->GetRAmin() + 1) * (this->GetDECmax() - this->GetDECmin() + 1)) + 1); g++){ temp_sparse_reps_strings.push_back(temp_sparse_reps_strings[(g - 1)] + temp_sparse_reps_grid[(g - 1)]); }
      for(g = 0; g < (((this->GetRAmax() - this->GetRAmin() + 1) * (this->GetDECmax() - this->GetDECmin() + 1)) + 1); g++){ temp_sparse_reps_grid[g] = temp_sparse_reps_strings[g]; }
      
      // write new sparse_reps_strings value to temp_sparse_reps_strings array, using various grids to achieve indexing
      
      // a. initialise temp_sparse_reps_strings
      temp_sparse_reps_strings.resize(0);
      for(g = 0; g < (2 * temp_sparse_reps_grid[((this->GetDECmax() - this->GetDECmin() + 1) * (this->GetRAmax() - this->GetRAmin() + 1))]); g++){ temp_sparse_reps_strings.push_back(0); }
      
      
      // b. for each line of sight through the new existing object bounding box, retrieve the channel range of each object string along this LoS
      for(sy = this->GetDECmin(); sy <= this->GetDECmax(); sy++){
	
	for(sx = this->GetRAmin(); sx <= this->GetRAmax(); sx++){
	  
	  // initialise the number of object strings written to this LoS
	  k = 0;
	  
	  // retrieve the starting index for object strings along this LoS
	  j = 2 * temp_sparse_reps_grid[(((sy - this->GetDECmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx - this->GetRAmin())];
	  
	  // write existing object's object strings to temp_strings_array
	  if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	    
	    if((sx >= this->Get_srep_size(0)) && (sx <= this->Get_srep_size(1)) && (sy >= this->Get_srep_size(2)) && (sy <= this->Get_srep_size(3))){
	      
	      for(g = this->Get_srep_grid((((sy - this->Get_srep_size(2)) * (this->Get_srep_size(1) - this->Get_srep_size(0) + 1)) + sx - this->Get_srep_size(0))); g < this->Get_srep_grid((((sy - this->Get_srep_size(2)) * (this->Get_srep_size(1) - this->Get_srep_size(0) + 1)) + sx - this->Get_srep_size(0) + 1)); g++){
		
		temp_sparse_reps_strings[(j + k)] = this->Get_srep_strings((2 * g));
		k++;
		temp_sparse_reps_strings[(j + k)] = this->Get_srep_strings(((2 * g) + 1));
		k++;
		
	      }
	      
	      // if((sx >= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0]) && (sx <= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][1]) && (sy >= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][2]) && (sy <= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][3]))
	    }
	    
	    // if((sparse_reps_update[obj_batch][(existing - (obj_batch * obj_limit))] != 0) && (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] >= 0))
	  }
	  
	  // write match_init[i] object strings to temp_strings_array
	  if((sx >= merged.Get_srep_size(0)) && (sx <= merged.Get_srep_size(1)) && (sy >= merged.Get_srep_size(2)) && (sy <= merged.Get_srep_size(3))){
	    
	    for(g = merged.Get_srep_grid((((sy - merged.Get_srep_size(2)) * (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1)) + sx - merged.Get_srep_size(0))); g < merged.Get_srep_grid((((sy - merged.Get_srep_size(2)) * (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1)) + sx - merged.Get_srep_size(0) + 1)); g++){
	      
	      temp_sparse_reps_strings[(j + k)] = merged.Get_srep_strings((2 * g));
	      k++;
	      temp_sparse_reps_strings[(j + k)] = merged.Get_srep_strings(((2 * g) + 1));
	      k++;
	      
	    }
	    
	    // if((sx >= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0]) && (sx <= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][1]) && (sy >= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][2]) && (sy <= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][3]))
	  }
	  
	  // for(sx = 0; sx < (); sx++)
	}
	
	// for(sy = 0; sy < (); sy++)
      }
      
      // over-write the existing object's sparse representations with the existing+merged sparse representations
      this->Set_srep_size(0,this->GetRAmin());
      this->Set_srep_size(1,this->GetRAmax());
      this->Set_srep_size(2,this->GetDECmin());
      this->Set_srep_size(3,this->GetDECmax());
      this->Set_srep_size(4,this->GetFREQmin());
      this->Set_srep_size(5,this->GetFREQmax());
      
      // a. grid
      this->Free_srep_grid();
      this->Create_srep_grid((((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)) + 1));
      for(g = 0; g < (((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)) + 1); g++){ this->Set_srep_grid(g,temp_sparse_reps_grid[g]); }
      
      // b. mini_mom0
      this->Free_mom0();
      this->Create_mom0(((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)));
      for(g = 0; g < ((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)); g++){ this->Set_mom0(g,temp_mom0[g]); }
      
      // c. mini_RAPV
      this->Free_RAPV();
      this->Create_RAPV(((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(5) - this->Get_srep_size(4) + 1)));
      for(g = 0; g < ((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(5) - this->Get_srep_size(4) + 1)); g++){ this->Set_RAPV(g,temp_RAPV[g]); }
      
      // d. mini_DECPV
      this->Free_DECPV();
      this->Create_DECPV(((this->Get_srep_size(5) - this->Get_srep_size(4) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)));
      for(g = 0; g < ((this->Get_srep_size(3) - this->Get_srep_size(2) + 1) * (this->Get_srep_size(5) - this->Get_srep_size(4) + 1)); g++){ this->Set_DECPV(g,temp_DECPV[g]); }
      
      // e. mini_obj_spec
      this->Free_ospec();
      this->Create_ospec((this->Get_srep_size(5) - this->Get_srep_size(4) + 1));
      for(g = 0; g < (this->Get_srep_size(5) - this->Get_srep_size(4) + 1); g++){ this->Set_ospec(g,temp_obj_spec[g]); }
      
      // f. mini_ref_spec
      this->Free_rspec();
      if((this->Get_srep_size(5) - this->Get_srep_size(4) + 1) >= 10){
	
	this->Create_rspec((2 * (this->Get_srep_size(5) - this->Get_srep_size(4) + 1)));
	for(g = 0; g < (2 * (this->Get_srep_size(5) - this->Get_srep_size(4) + 1)); g++){ this->Set_rspec(g,temp_ref_spec[g]); }
	
      } else {
	
	this->Create_rspec((this->Get_srep_size(5) - this->Get_srep_size(4) + 11));
	for(g = 0; g < (this->Get_srep_size(5) - this->Get_srep_size(4) + 11); g++){ this->Set_rspec(g,temp_ref_spec[g]); }
	
      }
      
      // g. mini_vfield
      this->Free_vfield();
      this->Create_vfield(((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)));
      for(g = 0; g < ((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)); g++){ this->Set_vfield(g,temp_vfield[g]); }
      
      // h. sparse_reps_strings
      this->Free_srep_strings();
      this->Create_srep_strings((2 * this->Get_srep_grid(((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)))));
      for(g = 0; g < (2 * this->Get_srep_grid(((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)))); g++){ this->Set_srep_strings(g,temp_sparse_reps_strings[g]); }
      
    }
  
    // if(this != &merged)
  }

}

void object_props_dbl::AddObject(object_props_dbl & merged,vector<int> & temp_sparse_reps_grid,vector<int> & temp_sparse_reps_strings,vector<double> & temp_mom0,vector<double> & temp_RAPV,vector<double> & temp_DECPV,vector<double> & temp_ref_spec,vector<double> & temp_obj_spec,vector<double> & temp_vfield){

  int j,k,g,sx,sy,sz;
    
  if(this != &merged){ 

    // update the properties of this object with the object being merged in
    this->AddVoxel(merged.ShowVoxels());
    this->AdjustRArange(merged.GetRAmin());
    this->AdjustRArange(merged.GetRAmax());
    this->AdjustDECrange(merged.GetDECmin());
    this->AdjustDECrange(merged.GetDECmax());
    this->AdjustFREQrange(merged.GetFREQmin());
    this->AdjustFREQrange(merged.GetFREQmax());
    this->AddRa(merged.GetRA());
    this->AddDec(merged.GetDEC());
    this->AddFreq(merged.GetFREQ());
    this->AddRA_i(1.0,merged.GetRAi());
    this->AddDec_i(1.0,merged.GetDECi());
    this->AddFreq_i(1.0,merged.GetFREQi());
    this->AddTotIntens(merged.GetTI());
    this->AddAvgIntens(merged.GetAvgI());
    this->AddSigmaItens(merged.GetSigmaI());
    this->AdjustRange(merged.GetMinI());		
    this->AdjustRange(merged.GetMaxI());	

    // combine the sparse representations using the temporary arrays --- if they exist
    if((merged.Get_srep_update() != 0) && (merged.Get_srep_size(0) >= 0)){
		  
      // write the existing object's sparse representations into temporary arrays and initialise temporary arrays at the same time, provided it exists
      
      // a. grid
      temp_sparse_reps_grid.resize(0);
      for(g = 0; g < (1 + ((this->GetRAmax() - this->GetRAmin() + 1) * (this->GetDECmax() - this->GetDECmin() + 1))); g++){ temp_sparse_reps_grid.push_back(0); }
      
      if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	
	for(sx = 0; sx < (this->Get_srep_size(1) - this->Get_srep_size(0) + 1); sx++){
	  for(sy = 0; sy < (this->Get_srep_size(3) - this->Get_srep_size(2) + 1); sy++){
	    
	    temp_sparse_reps_grid[(((sy + this->Get_srep_size(2) - this->GetDECmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx + this->Get_srep_size(0) - this->GetRAmin())] = this->Get_srep_grid(((sy * (this->Get_srep_size(1) - this->Get_srep_size(0) + 1)) + sx + 1)) - this->Get_srep_grid(((sy * (this->Get_srep_size(1) - this->Get_srep_size(0) + 1)) + sx));
	    
	  }
	  
	}
	
      }
      
      // b. mom-0
      temp_mom0.resize(0);
      for(g = 0; g < ((this->GetRAmax() - this->GetRAmin() + 1) * (this->GetDECmax() - this->GetDECmin() + 1)); g++){ temp_mom0.push_back(0.0); }
      
      if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	
	for(sx = 0; sx < (this->Get_srep_size(1) - this->Get_srep_size(0) + 1); sx++){
	  
	  for(sy = 0; sy < (this->Get_srep_size(3) - this->Get_srep_size(2) + 1); sy++){
	    
	    temp_mom0[(((sy + this->Get_srep_size(2) - this->GetDECmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx + this->Get_srep_size(0) - this->GetRAmin())]+=this->Get_mom0(((sy * (this->Get_srep_size(1) - this->Get_srep_size(0) + 1)) + sx));
	    
	  }
	  
	}
	
      }
      
      // c. RAPV
      temp_RAPV.resize(0);
      for(g = 0; g < ((this->GetRAmax() - this->GetRAmin() + 1) * (this->GetFREQmax() - this->GetFREQmin() + 1)); g++){ temp_RAPV.push_back(0.0); }
      
      if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	
	for(sx = 0; sx < (this->Get_srep_size(1) - this->Get_srep_size(0) + 1); sx++){
	  
	  for(sz = 0; sz < (this->Get_srep_size(5) - this->Get_srep_size(4) + 1); sz++){
	    
	    temp_RAPV[(((sz + this->Get_srep_size(4) - this->GetFREQmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx + this->Get_srep_size(0) - this->GetRAmin())]+=this->Get_RAPV(((sz * (this->Get_srep_size(1) - this->Get_srep_size(0) + 1)) + sx));
	    
	  }
	  
	}
	
      }
      
      // d. DECPV
      temp_DECPV.resize(0);
      for(g = 0; g < ((this->GetDECmax() - this->GetDECmin() + 1) * (this->GetFREQmax() - this->GetFREQmin() + 1)); g++){ temp_DECPV.push_back(0.0); }
      
      if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	
	for(sy = 0; sy < (this->Get_srep_size(3) - this->Get_srep_size(2) + 1); sy++){
	  
	  for(sz = 0; sz < (this->Get_srep_size(5) - this->Get_srep_size(4) + 1); sz++){
	    
	    temp_DECPV[(((sz + this->Get_srep_size(4) - this->GetFREQmin()) * (this->GetDECmax() - this->GetDECmin() + 1)) + sy + this->Get_srep_size(2) - this->GetDECmin())]+=this->Get_DECPV(((sz * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)) + sy));
	    
	  }
	  
	}
	
      }
      
      // e. ref_spec
      if((this->Get_srep_size(5) - this->Get_srep_size(4) + 1) >= 10){
	
	temp_ref_spec.resize(0);
	for(g = 0; g < (2 * (this->GetFREQmax() - this->GetFREQmin() + 1)); g++){ temp_ref_spec.push_back(0.0); }
	
	if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	  
	  for(sz = 0; sz < (2 * (this->Get_srep_size(5) - this->Get_srep_size(4) + 1)); sz++){ temp_ref_spec[(sz + this->Get_srep_size(4) - this->GetFREQmin())]+=this->Get_rspec(sz); } 
	  
	}
	
      } else {
	
	temp_ref_spec.resize(0);
	for(g = 0; g < (this->GetFREQmax() - this->GetFREQmin() + 11); g++){ temp_ref_spec.push_back(0.0); }
	
	if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	  
	  for(sz = 0; sz < (this->Get_srep_size(5) - this->Get_srep_size(4) + 11); sz++){ temp_ref_spec[(sz + this->Get_srep_size(4) - this->GetFREQmin())]+=this->Get_rspec(sz); } 		    
	  
	}
	
      }
      
      // f. obj_spec
      temp_obj_spec.resize(0);
      for(g = 0; g < (this->Get_srep_size(5) - this->Get_srep_size(4) + 1); g++){ temp_obj_spec.push_back(0.0); }		  
      if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	
	for(sz = 0; sz < (this->Get_srep_size(5) - this->Get_srep_size(4) + 1); sz++){
	  
	  temp_obj_spec[(sz + this->Get_srep_size(4) - this->GetFREQmin())]+=this->Get_ospec(sz);
	  
	}
	
      }
      
      // g. vfield
      temp_vfield.resize(0);
      for(g = 0; g < ((this->GetRAmax() - this->GetRAmin() + 1) * (this->GetDECmax() - this->GetDECmin() + 1)); g++){ temp_vfield.push_back(0.0); }
      
      if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	
	for(sx = 0; sx < (this->Get_srep_size(1) - this->Get_srep_size(0) + 1); sx++){
	  
	  for(sy = 0; sy < (this->Get_srep_size(3) - this->Get_srep_size(2) + 1); sy++){
	    
	    temp_vfield[(((sy + this->Get_srep_size(2) - this->GetDECmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx + this->Get_srep_size(0) - this->GetRAmin())]+=this->Get_vfield(((sy * (this->Get_srep_size(1) - this->Get_srep_size(0) + 1)) + sx));
	    
	  }
	  
	}
	
      }
      
      // write the merged object's sparse representations into temporary arrays
      
      // a. grid
      for(sx = 0; sx < (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1); sx++){
	
	for(sy = 0; sy < (merged.Get_srep_size(3) - merged.Get_srep_size(2) + 1); sy++){
	  
	  temp_sparse_reps_grid[(((sy + merged.Get_srep_size(2) - this->GetDECmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx + merged.Get_srep_size(0) - this->GetRAmin())]+=(merged.Get_srep_grid(((sy * (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1)) + sx + 1)) - merged.Get_srep_grid(((sy * (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1)) + sx)));
	  
	}
	
      }
      
      // b. mom-0
      for(sx = 0; sx < (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1); sx++){
	
	for(sy = 0; sy < (merged.Get_srep_size(3) - merged.Get_srep_size(2) + 1); sy++){
	  
	  temp_mom0[(((sy + merged.Get_srep_size(2) - this->GetDECmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx + merged.Get_srep_size(0) - this->GetRAmin())]+=merged.Get_mom0(((sy * (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1)) + sx));
	  
	}
	
      }
      
      // c. RAPV
      for(sx = 0; sx < (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1); sx++){
	
	for(sz = 0; sz < (merged.Get_srep_size(5) - merged.Get_srep_size(4) + 1); sz++){
	  
	  temp_RAPV[(((sz + merged.Get_srep_size(4) - this->GetFREQmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx + merged.Get_srep_size(0) - this->GetRAmin())]+=merged.Get_RAPV(((sz * (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1)) + sx));
	  
	}
	
      }
      
      // d. DECPV
      for(sy = 0; sy < (merged.Get_srep_size(3) - merged.Get_srep_size(2) + 1); sy++){
	
	for(sz = 0; sz < (merged.Get_srep_size(5) - merged.Get_srep_size(4) + 1); sz++){
	  
	  temp_DECPV[(((sz + merged.Get_srep_size(4) - this->GetFREQmin()) * (this->GetDECmax() - this->GetDECmin() + 1)) + sy + merged.Get_srep_size(2) - this->GetDECmin())]+=merged.Get_DECPV(((sz * (merged.Get_srep_size(3) - merged.Get_srep_size(2) + 1)) + sy));
	  
	}
	
      }
      
      // e. ref_spec
      if((merged.Get_srep_size(5) - merged.Get_srep_size(4) + 1) >= 10){
	
	for(sz = 0; sz < (2 * (merged.Get_srep_size(5) - merged.Get_srep_size(4) + 1)); sz++){ 
	  
	  if((sz + merged.Get_srep_size(4) - merged.GetFREQmin()) >= 0){
	    
	    temp_ref_spec[(sz + merged.Get_srep_size(4) - this->GetFREQmin())]+=merged.Get_rspec(sz); 
	    
	  }
	  
	} 
	
      } else {
	
	for(sz = 0; sz < (merged.Get_srep_size(5) - merged.Get_srep_size(4) + 11); sz++){ 
	  
	  if((sz + merged.Get_srep_size(4) - this->GetFREQmin()) >= 0){
	    
	    temp_ref_spec[(sz + merged.Get_srep_size(4) - this->GetFREQmin())]+=merged.Get_rspec(sz); 
	    
	  }
	  
	} 		    
	
      }
      
      // f. obj_spec
      for(sz = 0; sz < (merged.Get_srep_size(5) - merged.Get_srep_size(4) + 1); sz++){
	
	temp_obj_spec[(sz + merged.Get_srep_size(4) - this->GetFREQmin())]+=merged.Get_ospec(sz);
	
      }
      
      // g. vfield
      for(sx = 0; sx < (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1); sx++){
	
	for(sy = 0; sy < (merged.Get_srep_size(3) - merged.Get_srep_size(2) + 1); sy++){
	  
	  temp_vfield[(((sy + merged.Get_srep_size(2) - this->GetDECmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx + merged.Get_srep_size(0) - this->GetRAmin())]+=merged.Get_vfield(((sy * (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1)) + sx));
	  
	}
	
      }
      
      // convert temp_sparse_reps_grid from differential to cumulative counts using temp_sparse_reps_string as an intermediary
      temp_sparse_reps_strings.resize(1);
      temp_sparse_reps_strings[0] = 0;
      for(g = 1; g < (((this->GetRAmax() - this->GetRAmin() + 1) * (this->GetDECmax() - this->GetDECmin() + 1)) + 1); g++){ temp_sparse_reps_strings.push_back(temp_sparse_reps_strings[(g - 1)] + temp_sparse_reps_grid[(g - 1)]); }
      for(g = 0; g < (((this->GetRAmax() - this->GetRAmin() + 1) * (this->GetDECmax() - this->GetDECmin() + 1)) + 1); g++){ temp_sparse_reps_grid[g] = temp_sparse_reps_strings[g]; }
      
      // write new sparse_reps_strings value to temp_sparse_reps_strings array, using various grids to achieve indexing
      
      // a. initialise temp_sparse_reps_strings
      temp_sparse_reps_strings.resize(0);
      for(g = 0; g < (2 * temp_sparse_reps_grid[((this->GetDECmax() - this->GetDECmin() + 1) * (this->GetRAmax() - this->GetRAmin() + 1))]); g++){ temp_sparse_reps_strings.push_back(0); }
      
      
      // b. for each line of sight through the new existing object bounding box, retrieve the channel range of each object string along this LoS
      for(sy = this->GetDECmin(); sy <= this->GetDECmax(); sy++){
	
	for(sx = this->GetRAmin(); sx <= this->GetRAmax(); sx++){
	  
	  // initialise the number of object strings written to this LoS
	  k = 0;
	  
	  // retrieve the starting index for object strings along this LoS
	  j = 2 * temp_sparse_reps_grid[(((sy - this->GetDECmin()) * (this->GetRAmax() - this->GetRAmin() + 1)) + sx - this->GetRAmin())];
	  
	  // write existing object's object strings to temp_strings_array
	  if((this->Get_srep_update() != 0) && (this->Get_srep_size(0) >= 0)){
	    
	    if((sx >= this->Get_srep_size(0)) && (sx <= this->Get_srep_size(1)) && (sy >= this->Get_srep_size(2)) && (sy <= this->Get_srep_size(3))){
	      
	      for(g = this->Get_srep_grid((((sy - this->Get_srep_size(2)) * (this->Get_srep_size(1) - this->Get_srep_size(0) + 1)) + sx - this->Get_srep_size(0))); g < this->Get_srep_grid((((sy - this->Get_srep_size(2)) * (this->Get_srep_size(1) - this->Get_srep_size(0) + 1)) + sx - this->Get_srep_size(0) + 1)); g++){
		
		temp_sparse_reps_strings[(j + k)] = this->Get_srep_strings((2 * g));
		k++;
		temp_sparse_reps_strings[(j + k)] = this->Get_srep_strings(((2 * g) + 1));
		k++;
		
	      }
	      
	      // if((sx >= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0]) && (sx <= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][1]) && (sy >= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][2]) && (sy <= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][3]))
	    }
	    
	    // if((sparse_reps_update[obj_batch][(existing - (obj_batch * obj_limit))] != 0) && (sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0] >= 0))
	  }
	  
	  // write match_init[i] object strings to temp_strings_array
	  if((sx >= merged.Get_srep_size(0)) && (sx <= merged.Get_srep_size(1)) && (sy >= merged.Get_srep_size(2)) && (sy <= merged.Get_srep_size(3))){
	    
	    for(g = merged.Get_srep_grid((((sy - merged.Get_srep_size(2)) * (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1)) + sx - merged.Get_srep_size(0))); g < merged.Get_srep_grid((((sy - merged.Get_srep_size(2)) * (merged.Get_srep_size(1) - merged.Get_srep_size(0) + 1)) + sx - merged.Get_srep_size(0) + 1)); g++){
	      
	      temp_sparse_reps_strings[(j + k)] = merged.Get_srep_strings((2 * g));
	      k++;
	      temp_sparse_reps_strings[(j + k)] = merged.Get_srep_strings(((2 * g) + 1));
	      k++;
	      
	    }
	    
	    // if((sx >= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][0]) && (sx <= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][1]) && (sy >= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][2]) && (sy <= sparse_reps_size[obj_batch][(existing - (obj_batch * obj_limit))][3]))
	  }
	  
	  // for(sx = 0; sx < (); sx++)
	}
	
	// for(sy = 0; sy < (); sy++)
      }
      
      // over-write the existing object's sparse representations with the existing+merged sparse representations
      this->Set_srep_size(0,this->GetRAmin());
      this->Set_srep_size(1,this->GetRAmax());
      this->Set_srep_size(2,this->GetDECmin());
      this->Set_srep_size(3,this->GetDECmax());
      this->Set_srep_size(4,this->GetFREQmin());
      this->Set_srep_size(5,this->GetFREQmax());
      
      // a. grid
      this->Free_srep_grid();
      this->Create_srep_grid((((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)) + 1));
      for(g = 0; g < (((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)) + 1); g++){ this->Set_srep_grid(g,temp_sparse_reps_grid[g]); }
      
      // b. mini_mom0
      this->Free_mom0();
      this->Create_mom0(((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)));
      for(g = 0; g < ((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)); g++){ this->Set_mom0(g,temp_mom0[g]); }
      
      // c. mini_RAPV
      this->Free_RAPV();
      this->Create_RAPV(((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(5) - this->Get_srep_size(4) + 1)));
      for(g = 0; g < ((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(5) - this->Get_srep_size(4) + 1)); g++){ this->Set_RAPV(g,temp_RAPV[g]); }
      
      // d. mini_DECPV
      this->Free_DECPV();
      this->Create_DECPV(((this->Get_srep_size(5) - this->Get_srep_size(4) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)));
      for(g = 0; g < ((this->Get_srep_size(3) - this->Get_srep_size(2) + 1) * (this->Get_srep_size(5) - this->Get_srep_size(4) + 1)); g++){ this->Set_DECPV(g,temp_DECPV[g]); }
      
      // e. mini_obj_spec
      this->Free_ospec();
      this->Create_ospec((this->Get_srep_size(5) - this->Get_srep_size(4) + 1));
      for(g = 0; g < (this->Get_srep_size(5) - this->Get_srep_size(4) + 1); g++){ this->Set_ospec(g,temp_obj_spec[g]); }
      
      // f. mini_ref_spec
      this->Free_rspec();
      if((this->Get_srep_size(5) - this->Get_srep_size(4) + 1) >= 10){
	
	this->Create_rspec((2 * (this->Get_srep_size(5) - this->Get_srep_size(4) + 1)));
	for(g = 0; g < (2 * (this->Get_srep_size(5) - this->Get_srep_size(4) + 1)); g++){ this->Set_rspec(g,temp_ref_spec[g]); }
	
      } else {
	
	this->Create_rspec((this->Get_srep_size(5) - this->Get_srep_size(4) + 11));
	for(g = 0; g < (this->Get_srep_size(5) - this->Get_srep_size(4) + 11); g++){ this->Set_rspec(g,temp_ref_spec[g]); }
	
      }
      
      // g. mini_vfield
      this->Free_vfield();
      this->Create_vfield(((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)));
      for(g = 0; g < ((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)); g++){ this->Set_vfield(g,temp_vfield[g]); }
      
      // h. sparse_reps_strings
      this->Free_srep_strings();
      this->Create_srep_strings((2 * this->Get_srep_grid(((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)))));
      for(g = 0; g < (2 * this->Get_srep_grid(((this->Get_srep_size(1) - this->Get_srep_size(0) + 1) * (this->Get_srep_size(3) - this->Get_srep_size(2) + 1)))); g++){ this->Set_srep_strings(g,temp_sparse_reps_strings[g]); }
      
    }
  
    // if(this != &merged)
  }

}

