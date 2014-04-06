#include<fstream>
#include<iostream>
#include<string>
#include<sstream>
#include<iomanip>
#include<cmath>
#include<ctime>
#include<ctype.h>
#include<vector>
#define SIZE_OBJ_IDS 1e6
#define SIZE_CHECK_OBJ_IDS 1e6

class object_props {

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

public:
  
  object_props();
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

  // end of class definition
};

extern void PrintCatalogueHeader(std::fstream& output_file, int cat_mode);

extern void HeapSort_ids(int n, int ra[]);

extern int CreateObjects(float * data_vals, int * flag_vals, int size_x, int size_y, int size_z, int chunk_x_start, int chunk_y_start, int chunk_z_start, int merge_x, int merge_y, int merge_z, int min_x_size, int min_y_size, int min_z_size, int min_v_size, float intens_thresh_min, float intens_thresh_max, int flag_value, int start_obj, object_props ** detections, int * obj_ids, int& NO_obj_ids, int * check_obj_ids, int& NO_check_obj_ids, int obj_limit, int obj_batch_limit, int max_x_val, int max_y_val, int max_z_val, int ss_mode, int * data_metric, int * xyz_order);

extern int AddObjsToChunk(int * flag_vals, object_props ** detections, int NOobj, int obj_limit, int chunk_x_start, int chunk_y_start, int chunk_z_start, int chunk_x_size, int chunk_y_size, int chunk_z_size, int * check_obj_ids, int * data_metric, int * xyz_order);

extern void ThresholdObjs(object_props ** detections, int NOobj, int obj_limit, int min_x_size, int min_y_size, int min_z_size, int min_v_size, float intens_thresh_min, float intens_thresh_max, int min_LoS_count);

extern float CreateMoment0Map(float * plot_array, int NOobj, object_props ** detections, int NOx, int NOy, int obj_limit, int obj_batch_limit);

extern float CreateRAPVPlot(float * plot_array, int NOobj, object_props ** detections, int NOx, int NOz, int obj_limit, int obj_batch_limit);

extern float CreateDecPVPlot(float * plot_array, int NOobj, object_props ** detections, int NOy, int NOz, int obj_limit, int obj_batch_limit);

extern int CreateMoment0Bounds(object_props ** detections, int size_x, int size_y, int min_x, int min_y, float * plot_x, float * plot_y, int obj, int obj_limit);

extern int CreateRAPVBounds(object_props ** detections, int size_x, int size_y, int size_z, int min_x, int min_z, float * plot_x, float * plot_y, int obj, int obj_limit);

extern int CreateDecPVBounds(object_props ** detections, int size_x, int size_y, int size_z, int min_y, int min_z, float * plot_x, float * plot_y, int obj, int obj_limit);

extern int CreateCatalogue(std::fstream& output_file, object_props ** detections, int NOobj, int obj_limit, int cat_mode);

extern void InitObjGen(object_props **& detections, int& NOobj, int obj_limit, int obj_batch_limit, int *& obj_ids, int& NO_obj_ids, int *& check_obj_ids, int *& data_metric, int *& xyz_order);

extern void FreeObjGen(object_props **& detections, int NOobj, int obj_batch_limit, int *& obj_ids, int *& check_obj_ids, int *& data_metric, int *& xyz_order);

extern void CreateFitsMask(std::string output_code, int NOx, int NOy, int NOf, object_props ** detections, int NOobj, int obj_limit, int * flag_vals, int chunk_x_size, int chunk_y_size, int temp_chunk_x_size, int temp_chunk_y_size, int * data_metric, int * xyz_order);
    
extern void CreateMetric(int * data_metric, int * xyz_order, int NOx, int NOy, int NOz);




