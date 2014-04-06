#include<iostream>
#include<RJJ_ObjGen.h>

void InitObjGen(object_props **& detections, int& NOobj, int obj_limit, int obj_batch_limit, int *& obj_ids, int& NO_obj_ids, int *& check_obj_ids, int *& data_metric, int *& xyz_order){

  detections = new object_props * [obj_batch_limit];
  detections[0] = new object_props[obj_limit];
  obj_ids = new int[((int) SIZE_OBJ_IDS)];
  obj_ids[0] = 0;
  NO_obj_ids = 1;
  check_obj_ids = new int[((int) SIZE_CHECK_OBJ_IDS)];
  NOobj = 0;
  data_metric = new int[3];
  xyz_order = new int[3];
  
}

void FreeObjGen(object_props **& detections, int NOobj, int obj_batch_limit, int *& obj_ids, int *& check_obj_ids, int *& data_metric, int *& xyz_order){

  int i;

  delete [] obj_ids;
  delete [] check_obj_ids;
  for(i = 0; i < ceilf((NOobj / obj_batch_limit)); i++){

    delete [] detections[i];

  }
  delete [] detections;
  delete [] data_metric;
  delete [] xyz_order;

}


