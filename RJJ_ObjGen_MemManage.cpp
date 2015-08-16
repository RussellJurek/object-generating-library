#include<iostream>
#include<RJJ_ObjGen.h>

using namespace std;

// functions using floats

void InitObjGen(vector <object_props *> & detections, int & NOobj, int obj_limit, vector<int> & obj_ids, vector<int> & check_obj_ids, size_t * & data_metric, int * & xyz_order){

  detections.reserve(1000);
  detections.resize(0);
  detections.push_back(new object_props[obj_limit]);
  obj_ids.reserve(static_cast<int>(1.0E6));
  obj_ids.resize(0);
  obj_ids.push_back(0);
  check_obj_ids.reserve(static_cast<int>(1.0E6));
  check_obj_ids.resize(0);
  NOobj = 0;
  data_metric = new size_t[3];
  xyz_order = new int[3];
  
}

void InitObjGen(vector <object_props *> & detections, long int & NOobj, int obj_limit, vector<long int> & obj_ids, vector<long int> & check_obj_ids, size_t * & data_metric, int * & xyz_order){

  detections.reserve(1000);
  detections.resize(0);
  detections.push_back(new object_props[obj_limit]);
  obj_ids.reserve(static_cast<int>(1.0E6));
  obj_ids.resize(0);
  obj_ids.push_back(0);
  check_obj_ids.reserve(static_cast<int>(1.0E6));
  check_obj_ids.resize(0);
  NOobj = 0;
  data_metric = new size_t[3];
  xyz_order = new int[3];
  
}

void FreeObjGen(vector <object_props *> & detections, size_t * & data_metric, int * & xyz_order){

  unsigned int i;

  for(i = 0; i < detections.size(); ++i){

    delete [] detections[i];

  }
  delete [] data_metric;
  delete [] xyz_order;

}

// functions using doubles

void InitObjGen(vector <object_props_dbl *> & detections, int & NOobj, int obj_limit, vector<int> & obj_ids, vector<int> & check_obj_ids, size_t * & data_metric, int * & xyz_order){

  detections.reserve(1000);
  detections.resize(0);
  detections.push_back(new object_props_dbl[obj_limit]);
  obj_ids.reserve(static_cast<int>(1.0E6));
  obj_ids.resize(0);
  obj_ids.push_back(0);
  check_obj_ids.reserve(static_cast<int>(1.0E6));
  check_obj_ids.resize(0);
  NOobj = 0;
  data_metric = new size_t[3];
  xyz_order = new int[3];
  
}

void InitObjGen(vector <object_props_dbl *> & detections, long int & NOobj, int obj_limit, vector<long int> & obj_ids, vector<long int> & check_obj_ids, size_t * & data_metric, int * & xyz_order){

  detections.reserve(1000);
  detections.resize(0);
  detections.push_back(new object_props_dbl[obj_limit]);
  obj_ids.reserve(static_cast<int>(1.0E6));
  obj_ids.resize(0);
  obj_ids.push_back(0);
  check_obj_ids.reserve(static_cast<int>(1.0E6));
  check_obj_ids.resize(0);
  NOobj = 0;
  data_metric = new size_t[3];
  xyz_order = new int[3];
  
}

void FreeObjGen(vector <object_props_dbl *> & detections, size_t * & data_metric, int * & xyz_order){

  unsigned int i;

  for(i = 0; i < detections.size(); ++i){

    delete [] detections[i];

  }
  delete [] data_metric;
  delete [] xyz_order;

}


