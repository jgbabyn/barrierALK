#define TMB_LIB_INIT R_init_spatialBarrier
#include <TMB.hpp>
#include "R_inla_barrier.hpp"
#include "prediction.hpp"
#include "residuals.hpp"
#include "penal.hpp"
#include "spatialo.hpp"
#include "spatialspde.hpp"
#include "spatialAR.hpp"
#include "nonspatial.hpp"
#include "spatialARNoB.hpp"

bool isNAINT(int x){
  return NA_INTEGER==x;
}


enum valid_model {
		  not_spatial = 0,
		  spatial_barrier = 1,
		  spatial_spde = 2,
		  spatial_AR = 3,
		  spatial_ARnoBar = 4
};






template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_INTEGER(model);

  switch(model){
  case not_spatial:
    return nonSpatial(this);
    break;
  case spatial_barrier:
    return spatialo(this);
    break;
  case spatial_spde:
    return spatialspde(this);
    break;
  case spatial_AR:
    return spatialAR(this);
    break;
  case spatial_ARnoBar:
    return spatialARnB(this);
    break;
  default:
    error("Incorrect model choice");
    break;
  }

  //If we get here, there's a problem!
  return 1;
    
}
