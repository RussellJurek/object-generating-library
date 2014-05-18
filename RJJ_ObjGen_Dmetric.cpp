#include<iostream>
#include<RJJ_ObjGen.h>

using namespace std;

void CreateMetric(int * data_metric, int * xyz_order, int NOx, int NOy, int NOz){

  // ensure that the x_order, y_order and z_order values are 1,2,3
  if((xyz_order[0] <= xyz_order[1]) && (xyz_order[0] <= xyz_order[2])){
    
    xyz_order[0] = 1;
    if(xyz_order[1] <= xyz_order[2]){ 

      xyz_order[1] = 2; 
      xyz_order[2] = 3;

    } else {

      xyz_order[1] = 3;
      xyz_order[2] = 2;

    }

  } else if((xyz_order[1] <= xyz_order[0]) && (xyz_order[1] <= xyz_order[2])){

    xyz_order[1] = 1;
    if(xyz_order[0] <= xyz_order[2]){ 

      xyz_order[0] = 2; 
      xyz_order[2] = 3;

    } else {

      xyz_order[0] = 3;
      xyz_order[2] = 2;

    }

  } else if((xyz_order[2] <= xyz_order[0]) && (xyz_order[2] <= xyz_order[1])){

    xyz_order[2] = 1;
    if(xyz_order[0] <= xyz_order[1]){ 

      xyz_order[0] = 2; 
      xyz_order[1] = 3;

    } else {

      xyz_order[0] = 3;
      xyz_order[1] = 2;

    }

  }

  // assign values to the data_metric array
  
  // 1. x value
  switch(xyz_order[0]){
  case 1:
    data_metric[0] = 1;
    break;
  case 2:
    if(xyz_order[1] == 1){ data_metric[0] = NOy; } else { data_metric[0] = NOz; }
    break;
  case 3:
    data_metric[0] = NOy * NOz;
    break;
  default:
    data_metric[0] = 1;
    break;
  }
    
  // 2. y value
  switch(xyz_order[1]){
  case 1:
    data_metric[1] = 1;
    break;
  case 2:
    if(xyz_order[0] == 1){ data_metric[1] = NOx; } else { data_metric[1] = NOz; }
    break;
  case 3:
    data_metric[1] = NOx * NOz;
    break;
  default:
    data_metric[1] = 1;
    break;
  }

  // 3. z value
  switch(xyz_order[2]){
  case 1:
    data_metric[2] = 1;
    break;
  case 2:
    if(xyz_order[0] == 1){ data_metric[2] = NOx; } else { data_metric[2] = NOy; }
    break;
  case 3:
    data_metric[2] = NOx * NOy;
    break;
  default:
    data_metric[2] = 1;
    break;
  }

}





