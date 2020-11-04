
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template <class Type>
Type spatialARnB(objective_function<Type>* obj){
  using namespace R_inla;
  using namespace density;

  DATA_VECTOR(y);
  DATA_MATRIX(X);
  DATA_SPARSE_MATRIX(A);
  DATA_FACTOR(cohort);
  DATA_STRUCT(fem,spde_t);
  DATA_SCALAR(range_fraction); //Just ignore
  DATA_INTEGER(ages);
  DATA_VECTOR(size);
  ages = ages-1;
  DATA_FACTOR(year);

  DATA_IMATRIX(ageTimeKey);

  DATA_VECTOR(lambda);
  DATA_MATRIX(P);
  
  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(rhoT);

  PARAMETER(log_kappa);
  PARAMETER(log_tau);
  
  
  PARAMETER_MATRIX(ranX);

  Type kappa = exp(log_kappa);
  Type sigma_u = 1/exp(log_tau);
  vector<Type> rho = rhoTrans(rhoT);

  Type nll = 0.0;

  SparseMatrix<Type> Q = Q_spde(fem,kappa);

  for(int i = 0; i < ranX.cols();i++){
     nll += GMRF(Q)(ranX.col(i));
   }

  matrix<Type> Xr = makeXr(rho,sigma_u,ranX,ageTimeKey);
  vector<Type> eta = X*beta;

  vector< vector <Type> > AXs(ageTimeKey.rows()*ageTimeKey.cols());
  for(int i = 0; i < AXs.size();i++){
    AXs(i) = A*Xr.col(i);
  }

  for(int i = 0;i < eta.size();i++){
    int ya = ageTimeKey(year(i),cohort(i));
    eta(i) = eta(i) + AXs(ya)(i);
    nll -= dbinom_robust(y(i),size(i),eta(i),true);
  }

  REPORT(eta);
  REPORT(beta);
  REPORT(Xr);
  REPORT(ranX);
  REPORT(cohort);
  REPORT(year);
  REPORT(ageTimeKey);
  REPORT(rho);
  REPORT(rhoT);
  REPORT(log_kappa);
  REPORT(log_tau);
  
  nll += Ridge(lambda,P,beta);
  

  return nll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
