#define TMB_LIB_INIT R_init_spatBarr
#include <TMB.hpp>
#include "R_inla_barrier.hpp"

template<class Type>
Type rhoTrans(Type x){
  return Type(2)/(Type(1) + exp(-Type(2)*x))-Type(1);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla_barrier;
  using namespace density;

  DATA_VECTOR(y);
  DATA_MATRIX(X);
  DATA_SPARSE_MATRIX(A);
  DATA_FACTOR(cohort);
  DATA_STRUCT(fem,fem_barrier_t);
  DATA_INTEGER(ages);
  ages = ages-1;

  PARAMETER_VECTOR(beta);
  PARAMETER_MATRIX(ranX);
  PARAMETER(log_range);
  PARAMETER(log_sigma_u);
  PARAMETER(rhoT);

  vector<Type> ranges(2);
  ranges(0) = exp(log_range);
  ranges(1) = ranges(0)*0.1;
  Type sigma_u = exp(log_sigma_u);
  Type rho = rhoTrans(rhoT);

  
  Type nll = 0;
  
  SparseMatrix<Type> Q = Q_barrier(fem,ranges,sigma_u);

  for(int i = 0; i < ranX.cols();i++){
    nll += GMRF(Q)(ranX.col(i));
  }
  
  
  matrix<Type> Xr = ranX;
  for(int i = 1; i < ranX.cols();i++){
    Xr.col(i) = rho*Xr.col(i-1)+sqrt(1-pow(rho,2))*ranX.col(i);
  }

  vector<Type> eta = X*beta;
  vector< vector <Type> > AXs(ages);
  for(int i = 0; i < AXs.size();i++){
    AXs(i) = A*Xr.col(i);
  }

  for(int i = 0;i < eta.size();i++){
    eta(i) = eta(i) + AXs(cohort(i))(i);
    nll -= dbinom_robust(y(i),Type(1.0),eta(i),true);
  }

  return nll;


}

  
