#include<fstream>
#include<iostream>
#include<string>
#include<sstream>
#include<iomanip>
#include<cmath>
#include<ctime>
#include<ctype.h>
#include<vector>
#include<algorithm>

#ifndef RJJ_ObjGen 
#define RJJ_ObjGen

using namespace std;

// float precision class definition

class object_props {

 private:

  int NOvox;
  float ra, ra_i, ra_min, ra_max;
  float dec, dec_i, dec_min, dec_max;
  float freq, freq_i, freq_min, freq_max;
  float tot_intens, avg_intens, min_intens, max_intens, sigma_intens, rms;
  int srep_size[6], * srep_grid, * srep_strings, srep_update;
  float * mini_mom0, * mini_RAPV, * mini_DECPV, * mini_obj_spec, * mini_ref_spec, * mini_vfield;
  float w_max, w20_min, w20_max, w50_min, w50_max;
  float cw_max, cw20_min, cw20_max, cw50_min, cw50_max;

  // new parameters --- multiple central moment calculations
  float p_tot_intens, n_tot_intens;
  float p_ra_i, p_dec_i, p_freq_i;
  float n_ra_i, n_dec_i, n_freq_i;

  // define operators that don't make sense for this type of object
  object_props operator - (const object_props & subtracted);
  void operator -= (const object_props & subtracted);
  object_props operator * (const object_props & multiplier);
  void operator *= (const object_props & multiplier);
  object_props operator / (const object_props & divisor);
  void operator /= (const object_props & divisor);

  object_props& operator ++ ();
  object_props operator ++ (int);
  object_props& operator -- ();
  object_props operator -- (int);

  bool operator < (const object_props & compareTo);
  bool operator > (const object_props & compareTo);
  bool operator <= (const object_props & compareTo);
  bool operator >= (const object_props & compareTo); 

 public:
  
  object_props();
  object_props(const object_props & copied);
  object_props(object_props && moved);
  ~object_props();

  void AddVoxel(int value);
  void AddRa(float value);
  void AddDec(float value);
  void AddFreq(float value);
  void AddRA_i(float pos, float value);
  void AddDec_i(float pos, float value);
  void AddFreq_i(float pos, float value);
  void AddTotIntens(float value);
  void AddAvgIntens(float value);
  void AddSigmaItens(float value);
  void AdjustRange(float value);
  void CalcProps();  
  void ShowAll_file_WCS(int id, std::fstream& output_file, int cat_mode, double wcs_vals[6]);
  void ShowAll_file(int id, std::fstream& output_file, int cat_mode);
  void ShowProps_file_WCS(int id, std::fstream& output_file, double wcs_vals[6]);
  void ShowProps_file(int id, std::fstream& output_file);
  void ShowSrep_file(std::fstream& output_file);
  void AdjustRArange(float value);
  void AdjustDECrange(float value);
  void AdjustFREQrange(float value);
  int ShowArea();
  int ShowVoxels();
  int ShowRArange();
  int ShowDECrange();
  int ShowFREQrange();
  int GetRAmin();
  int GetRAmax();
  int GetDECmin();
  int GetDECmax();
  int GetFREQmin();
  int GetFREQmax();
  float GetRA();
  float GetDEC();
  float GetFREQ();
  float GetRAi();
  float GetDECi();
  float GetFREQi();
  float GetTI();
  float GetSigmaI();
  float GetAvgI();
  float GetMinI();
  float GetMaxI();
  void ReInit();
  
  float GetRAi_p();
  float GetDECi_p();
  float GetFREQi_p();
  float GetTI_p();
  float GetRAi_n();
  float GetDECi_n();
  float GetFREQi_n();
  float GetTI_n();

  void AddPoint(float x_pos, float y_pos, float z_pos, float value);
  void AddObject(object_props & merged);
  void AddObject(object_props & merged,vector<int> & temp_sparse_reps_grid,vector<int> & temp_sparse_reps_strings,vector<float> & temp_mom0,vector<float> & temp_RAPV,vector<float> & temp_DECPV,vector<float> & temp_ref_spec,vector<float> & temp_obj_spec,vector<float> & temp_vfield);

  void Set_w_max(float value);
  float Get_w_max();
  void Set_w20_min(float value);
  float Get_w20_min();
  void Set_w20_max(float value);
  float Get_w20_max();
  void Set_w50_min(float value);
  float Get_w50_min();
  void Set_w50_max(float value);
  float Get_w50_max();
  
  void Set_cw_max(float value);
  float Get_cw_max();
  void Set_cw20_min(float value);
  float Get_cw20_min();
  void Set_cw20_max(float value);
  float Get_cw20_max();
  void Set_cw50_min(float value);
  float Get_cw50_min();
  void Set_cw50_max(float value);
  float Get_cw50_max();
  
  void Set_srep_update(int value);
  int Get_srep_update();

  void Set_srep_size(int index, int value);
  int Get_srep_size(int index);
  
  void Create_srep_grid(int value);
  void Set_srep_grid(int index, int value);
  int Get_srep_grid(int index);
  void Free_srep_grid();

  void Create_srep_strings(int value);
  void Set_srep_strings(int index, int value);
  int Get_srep_strings(int index);
  void Free_srep_strings();

  void Create_mom0(int value);
  void Set_mom0(int index, float value);
  void Add_mom0(int index, float value);
  float Get_mom0(int index);
  void Free_mom0();

  void Create_RAPV(int value);
  void Set_RAPV(int index, float value);
  void Add_RAPV(int index, float value);
  float Get_RAPV(int index);
  void Free_RAPV();

  void Create_DECPV(int value);
  void Set_DECPV(int index, float value);
  void Add_DECPV(int index, float value);
  float Get_DECPV(int index);
  void Free_DECPV();

  void Create_ospec(int value);
  void Set_ospec(int index, float value);
  void Add_ospec(int index, float value);
  float Get_ospec(int index);
  void Free_ospec();

  void Create_rspec(int value);
  void Set_rspec(int index, float value);
  void Add_rspec(int index, float value);
  float Get_rspec(int index);
  void Free_rspec();

  void Create_vfield(int value);
  void Set_vfield(int index, float value);
  void Add_vfield(int index, float value);
  float Get_vfield(int index);
  void Free_vfield();

  void ReInit_srep();
  void ReInit_mini();
  void ReInit_size();

  // define operators that make sense for this type of object
  object_props & operator = (const object_props & copied);
  object_props & operator = (object_props && moved);

  bool operator == (const object_props & compareTo);
  bool operator != (const object_props & compareTo);
  object_props operator + (object_props & added);
  void operator += (object_props & added);

  // end of class definition
};

// double precision class definition

class object_props_dbl {

 private:

  int NOvox;
  double ra, ra_i, ra_min, ra_max;
  double dec, dec_i, dec_min, dec_max;
  double freq, freq_i, freq_min, freq_max;
  double tot_intens, avg_intens, min_intens, max_intens, sigma_intens, rms;
  int srep_size[6], * srep_grid, * srep_strings, srep_update;
  double * mini_mom0, * mini_RAPV, * mini_DECPV, * mini_obj_spec, * mini_ref_spec, * mini_vfield;
  double w_max, w20_min, w20_max, w50_min, w50_max;
  double cw_max, cw20_min, cw20_max, cw50_min, cw50_max;

  // new parameters --- multiple central moment calculations
  double p_tot_intens, n_tot_intens;
  double p_ra_i, p_dec_i, p_freq_i;
  double n_ra_i, n_dec_i, n_freq_i;

  // define operators that don't make sense for this type of object
  object_props_dbl operator - (const object_props_dbl & subtracted);
  void operator -= (const object_props_dbl & subtracted);
  object_props_dbl operator * (const object_props_dbl & multiplier);
  void operator *= (const object_props_dbl & multiplier);
  object_props_dbl operator / (const object_props_dbl & divisor);
  void operator /= (const object_props_dbl & divisor);

  object_props_dbl & operator ++ ();
  object_props_dbl operator ++ (int);
  object_props_dbl & operator -- ();
  object_props_dbl operator -- (int);

  bool operator < (const object_props_dbl & compareTo);
  bool operator > (const object_props_dbl & compareTo);
  bool operator <= (const object_props_dbl & compareTo);
  bool operator >= (const object_props_dbl & compareTo);

public:
  
  object_props_dbl();
  object_props_dbl(const object_props_dbl & copied);
  object_props_dbl(object_props_dbl && moved);
  ~object_props_dbl();

  void AddVoxel(int value);
  void AddRa(double value);
  void AddDec(double value);
  void AddFreq(double value);
  void AddRA_i(double pos, double value);
  void AddDec_i(double pos, double value);
  void AddFreq_i(double pos, double value);
  void AddTotIntens(double value);
  void AddAvgIntens(double value);
  void AddSigmaItens(double value);
  void AdjustRange(double value);
  void CalcProps();  
  void ShowAll_file_WCS(int id, std::fstream& output_file, int cat_mode, double wcs_vals[6]);
  void ShowAll_file(int id, std::fstream& output_file, int cat_mode);
  void ShowProps_file_WCS(int id, std::fstream& output_file, double wcs_vals[6]);
  void ShowProps_file(int id, std::fstream& output_file);
  void ShowSrep_file(std::fstream& output_file);
  void AdjustRArange(double value);
  void AdjustDECrange(double value);
  void AdjustFREQrange(double value);
  int ShowArea();
  int ShowVoxels();
  int ShowRArange();
  int ShowDECrange();
  int ShowFREQrange();
  int GetRAmin();
  int GetRAmax();
  int GetDECmin();
  int GetDECmax();
  int GetFREQmin();
  int GetFREQmax();
  double GetRA();
  double GetDEC();
  double GetFREQ();
  double GetRAi();
  double GetDECi();
  double GetFREQi();
  double GetTI();
  double GetSigmaI();
  double GetAvgI();
  double GetMinI();
  double GetMaxI();
  void ReInit();
  
  float GetRAi_p();
  float GetDECi_p();
  float GetFREQi_p();
  float GetTI_p();
  float GetRAi_n();
  float GetDECi_n();
  float GetFREQi_n();
  float GetTI_n();

  void AddPoint(double x_pos, double y_pos, double z_pos, double value);
  void AddObject(object_props_dbl & merged);
  void AddObject(object_props_dbl & merged,vector<int> & temp_sparse_reps_grid,vector<int> & temp_sparse_reps_strings,vector<double> & temp_mom0,vector<double> & temp_RAPV,vector<double> & temp_DECPV,vector<double> & temp_ref_spec,vector<double> & temp_obj_spec,vector<double> & temp_vfield);

  void Set_w_max(double value);
  double Get_w_max();
  void Set_w20_min(double value);
  double Get_w20_min();
  void Set_w20_max(double value);
  double Get_w20_max();
  void Set_w50_min(double value);
  double Get_w50_min();
  void Set_w50_max(double value);
  double Get_w50_max();
  
  void Set_cw_max(double value);
  double Get_cw_max();
  void Set_cw20_min(double value);
  double Get_cw20_min();
  void Set_cw20_max(double value);
  double Get_cw20_max();
  void Set_cw50_min(double value);
  double Get_cw50_min();
  void Set_cw50_max(double value);
  double Get_cw50_max();
  
  void Set_srep_update(int value);
  int Get_srep_update();

  void Set_srep_size(int index, int value);
  int Get_srep_size(int index);
  
  void Create_srep_grid(int value);
  void Set_srep_grid(int index, int value);
  int Get_srep_grid(int index);
  void Free_srep_grid();

  void Create_srep_strings(int value);
  void Set_srep_strings(int index, int value);
  int Get_srep_strings(int index);
  void Free_srep_strings();

  void Create_mom0(int value);
  void Set_mom0(int index, double value);
  void Add_mom0(int index, double value);
  double Get_mom0(int index);
  void Free_mom0();

  void Create_RAPV(int value);
  void Set_RAPV(int index, double value);
  void Add_RAPV(int index, double value);
  double Get_RAPV(int index);
  void Free_RAPV();

  void Create_DECPV(int value);
  void Set_DECPV(int index, double value);
  void Add_DECPV(int index, double value);
  double Get_DECPV(int index);
  void Free_DECPV();

  void Create_ospec(int value);
  void Set_ospec(int index, double value);
  void Add_ospec(int index, double value);
  double Get_ospec(int index);
  void Free_ospec();

  void Create_rspec(int value);
  void Set_rspec(int index, double value);
  void Add_rspec(int index, double value);
  double Get_rspec(int index);
  void Free_rspec();

  void Create_vfield(int value);
  void Set_vfield(int index, double value);
  void Add_vfield(int index, double value);
  double Get_vfield(int index);
  void Free_vfield();

  void ReInit_srep();
  void ReInit_mini();
  void ReInit_size();

  // define operators that make sense for this type of object
  object_props_dbl & operator = (const object_props_dbl & copied);
  object_props_dbl & operator = (object_props_dbl && moved);

  bool operator == (const object_props_dbl & compareTo);
  bool operator != (const object_props_dbl & compareTo);
  object_props_dbl operator + (object_props_dbl & added);
  void operator += (object_props_dbl & added);

  // end of class definition
};

extern void PrintCatalogueHeader(std::fstream& output_file, int cat_mode);

extern void CreateMetric(int * data_metric, int * xyz_order, int NOx, int NOy, int NOz);

// functions using floats

extern int CreateObjects(float * data_vals, int * flag_vals, int size_x, int size_y, int size_z, int chunk_x_start, int chunk_y_start, int chunk_z_start, int merge_x, int merge_y, int merge_z, int min_x_size, int min_y_size, int min_z_size, int min_v_size, float intens_thresh_min, float intens_thresh_max, int flag_value, int start_obj, vector<object_props *> & detections, vector<int> & obj_ids, vector<int> & check_obj_ids, int obj_limit, int max_x_val, int max_y_val, int max_z_val, int ss_mode, int * data_metric, int * xyz_order);
extern long int CreateObjects(float * data_vals, long int * flag_vals, int size_x, int size_y, int size_z, int chunk_x_start, int chunk_y_start, int chunk_z_start, int merge_x, int merge_y, int merge_z, int min_x_size, int min_y_size, int min_z_size, int min_v_size, float intens_thresh_min, float intens_thresh_max, long int flag_value, long int start_obj, vector<object_props *> & detections, vector<long int> & obj_ids, vector<long int> & check_obj_ids, int obj_limit, int max_x_val, int max_y_val, int max_z_val, int ss_mode, int * data_metric, int * xyz_order);

extern int AddObjsToChunk(int * flag_vals, vector<object_props *> & detections, int NOobj, int obj_limit, int chunk_x_start, int chunk_y_start, int chunk_z_start, int chunk_x_size, int chunk_y_size, int chunk_z_size, vector<int> & check_obj_ids, int * data_metric, int * xyz_order);
extern long int AddObjsToChunk(long int * flag_vals, vector<object_props *> & detections, long int NOobj, int obj_limit, int chunk_x_start, int chunk_y_start, int chunk_z_start, int chunk_x_size, int chunk_y_size, int chunk_z_size, vector<long int> & check_obj_ids, int * data_metric, int * xyz_order);

extern void ThresholdObjs(vector<object_props *> & detections, int NOobj, int obj_limit, int min_x_size, int min_y_size, int min_z_size, int min_v_size, float intens_thresh_min, float intens_thresh_max, int min_LoS_count);
extern void ThresholdObjs(vector<object_props *> & detections, long int NOobj, int obj_limit, int min_x_size, int min_y_size, int min_z_size, int min_v_size, float intens_thresh_min, float intens_thresh_max, int min_LoS_count);

extern float CreateMoment0Map(float * plot_array, int NOobj, vector<object_props *> & detections, int NOx, int NOy, int obj_limit);
extern float CreateMoment0Map(float * plot_array, long int NOobj, vector<object_props *> & detections, int NOx, int NOy, int obj_limit);

extern float CreateRAPVPlot(float * plot_array, int NOobj, vector<object_props *> & detections, int NOx, int NOz, int obj_limit);
extern float CreateRAPVPlot(float * plot_array, long int NOobj, vector<object_props *> & detections, int NOx, int NOz, int obj_limit);

extern float CreateDecPVPlot(float * plot_array, int NOobj, vector<object_props *> & detections, int NOy, int NOz, int obj_limit);
extern float CreateDecPVPlot(float * plot_array, long int NOobj, vector<object_props *> & detections, int NOy, int NOz, int obj_limit);

extern int CreateMoment0Bounds(vector<object_props *> & detections, int size_x, int size_y, int min_x, int min_y, float * plot_x, float * plot_y, int obj, int obj_limit);
extern int CreateMoment0Bounds(vector<object_props *> & detections, int size_x, int size_y, int min_x, int min_y, float * plot_x, float * plot_y, long int obj, int obj_limit);

extern int CreateRAPVBounds(vector<object_props *> & detections, int size_x, int size_y, int size_z, int min_x, int min_z, float * plot_x, float * plot_y, int obj, int obj_limit);
extern int CreateRAPVBounds(vector<object_props *> & detections, int size_x, int size_y, int size_z, int min_x, int min_z, float * plot_x, float * plot_y, long int obj, int obj_limit);

extern int CreateDecPVBounds(vector<object_props *> & detections, int size_x, int size_y, int size_z, int min_y, int min_z, float * plot_x, float * plot_y, int obj, int obj_limit);
extern int CreateDecPVBounds(vector<object_props *> & detections, int size_x, int size_y, int size_z, int min_y, int min_z, float * plot_x, float * plot_y, long int obj, int obj_limit);

extern int CreateCatalogue(std::fstream& output_file, vector<object_props *> & detections, int NOobj, int obj_limit, int cat_mode);
extern long int CreateCatalogue(std::fstream& output_file, vector<object_props *> & detections, long int NOobj, int obj_limit, int cat_mode);

extern void InitObjGen(vector<object_props *> & detections, int & NOobj, int obj_limit, vector<int> & obj_ids, vector<int> & check_obj_ids, int *& data_metric, int *& xyz_order);
extern void InitObjGen(vector<object_props *> & detections, long int & NOobj, int obj_limit, vector<long int> & obj_ids, vector<long int> & check_obj_ids, int *& data_metric, int *& xyz_order);

extern void FreeObjGen(vector<object_props *> & detections, int * & data_metric, int * & xyz_order);

extern void CreateFitsMask(std::string output_code, int NOx, int NOy, int NOf, vector<object_props *> & detections, int NOobj, int obj_limit, int * flag_vals, int chunk_x_size, int chunk_y_size, int temp_chunk_x_size, int temp_chunk_y_size, int * data_metric, int * xyz_order);
extern void CreateFitsMask(std::string output_code, int NOx, int NOy, int NOf, vector<object_props *> & detections, long int NOobj, int obj_limit, long int * flag_vals, int chunk_x_size, int chunk_y_size, int temp_chunk_x_size, int temp_chunk_y_size, int * data_metric, int * xyz_order);

// functions using doubles

extern int CreateObjects(double * data_vals, int * flag_vals, int size_x, int size_y, int size_z, int chunk_x_start, int chunk_y_start, int chunk_z_start, int merge_x, int merge_y, int merge_z, int min_x_size, int min_y_size, int min_z_size, int min_v_size, double intens_thresh_min, double intens_thresh_max, int flag_value, int start_obj, vector<object_props_dbl *> & detections, vector<int> & obj_ids, vector<int> & check_obj_ids, int obj_limit, int max_x_val, int max_y_val, int max_z_val, int ss_mode, int * data_metric, int * xyz_order);
extern long int CreateObjects(double * data_vals, long int * flag_vals, int size_x, int size_y, int size_z, int chunk_x_start, int chunk_y_start, int chunk_z_start, int merge_x, int merge_y, int merge_z, int min_x_size, int min_y_size, int min_z_size, int min_v_size, double intens_thresh_min, double intens_thresh_max, long int flag_value, long int start_obj, vector<object_props_dbl *> & detections, vector<long int> & obj_ids, vector<long int> & check_obj_ids, int obj_limit, int max_x_val, int max_y_val, int max_z_val, int ss_mode, int * data_metric, int * xyz_order);

extern int AddObjsToChunk(int * flag_vals, vector<object_props_dbl *> & detections, int NOobj, int obj_limit, int chunk_x_start, int chunk_y_start, int chunk_z_start, int chunk_x_size, int chunk_y_size, int chunk_z_size, vector<int> & check_obj_ids, int * data_metric, int * xyz_order);
extern long int AddObjsToChunk(long int * flag_vals, vector<object_props_dbl *> & detections, long int NOobj, int obj_limit, int chunk_x_start, int chunk_y_start, int chunk_z_start, int chunk_x_size, int chunk_y_size, int chunk_z_size, vector<long int> & check_obj_ids, int * data_metric, int * xyz_order);

extern void ThresholdObjs(vector<object_props_dbl *> & detections, int NOobj, int obj_limit, int min_x_size, int min_y_size, int min_z_size, int min_v_size, double intens_thresh_min, double intens_thresh_max, int min_LoS_count);
extern void ThresholdObjs(vector<object_props_dbl *> & detections, long int NOobj, int obj_limit, int min_x_size, int min_y_size, int min_z_size, int min_v_size, double intens_thresh_min, double intens_thresh_max, int min_LoS_count);

extern float CreateMoment0Map(float * plot_array, int NOobj, vector<object_props_dbl *> & detections, int NOx, int NOy, int obj_limit);
extern float CreateMoment0Map(float * plot_array, long int NOobj, vector<object_props_dbl *> & detections, int NOx, int NOy, int obj_limit);

extern float CreateRAPVPlot(float * plot_array, int NOobj, vector<object_props_dbl *> & detections, int NOx, int NOz, int obj_limit);
extern float CreateRAPVPlot(float * plot_array, long int NOobj, vector<object_props_dbl *> & detections, int NOx, int NOz, int obj_limit);

extern float CreateDecPVPlot(float * plot_array, int NOobj, vector<object_props_dbl *> & detections, int NOy, int NOz, int obj_limit);
extern float CreateDecPVPlot(float * plot_array, long int NOobj, vector<object_props_dbl *> & detections, int NOy, int NOz, int obj_limit);

extern int CreateMoment0Bounds(vector<object_props_dbl *> & detections, int size_x, int size_y, int min_x, int min_y, float * plot_x, float * plot_y, int obj, int obj_limit);
extern int CreateMoment0Bounds(vector<object_props_dbl *> & detections, int size_x, int size_y, int min_x, int min_y, float * plot_x, float * plot_y, long int obj, int obj_limit);

extern int CreateRAPVBounds(vector<object_props_dbl *> & detections, int size_x, int size_y, int size_z, int min_x, int min_z, float * plot_x, float * plot_y, int obj, int obj_limit);
extern int CreateRAPVBounds(vector<object_props_dbl *> & detections, int size_x, int size_y, int size_z, int min_x, int min_z, float * plot_x, float * plot_y, long int obj, int obj_limit);

extern int CreateDecPVBounds(vector<object_props_dbl *> & detections, int size_x, int size_y, int size_z, int min_y, int min_z, float * plot_x, float * plot_y, int obj, int obj_limit);
extern int CreateDecPVBounds(vector<object_props_dbl *> & detections, int size_x, int size_y, int size_z, int min_y, int min_z, float * plot_x, float * plot_y, long int obj, int obj_limit);

extern int CreateCatalogue(std::fstream& output_file, vector<object_props_dbl *> & detections, int NOobj, int obj_limit, int cat_mode);
extern long int CreateCatalogue(std::fstream& output_file, vector<object_props_dbl *> & detections, long int NOobj, int obj_limit, int cat_mode);

extern void InitObjGen(vector<object_props_dbl *> & detections, int & NOobj, int obj_limit, vector<int> & obj_ids, vector<int> & check_obj_ids, int *& data_metric, int *& xyz_order);
extern void InitObjGen(vector<object_props_dbl *> & detections, long int & NOobj, int obj_limit, vector<long int> & obj_ids, vector<long int> & check_obj_ids, int *& data_metric, int *& xyz_order);

extern void FreeObjGen(vector<object_props_dbl *> & detections, int * & data_metric, int * & xyz_order);

extern void CreateFitsMask(std::string output_code, int NOx, int NOy, int NOf, vector<object_props_dbl *> & detections, int NOobj, int obj_limit, int * flag_vals, int chunk_x_size, int chunk_y_size, int temp_chunk_x_size, int temp_chunk_y_size, int * data_metric, int * xyz_order);
extern void CreateFitsMask(std::string output_code, int NOx, int NOy, int NOf, vector<object_props_dbl *> & detections, long int NOobj, int obj_limit, long int * flag_vals, int chunk_x_size, int chunk_y_size, int temp_chunk_x_size, int temp_chunk_y_size, int * data_metric, int * xyz_order);

#endif



