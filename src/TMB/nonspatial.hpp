#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template <class Type>
Type nonSpatial(objective_function<Type>* obj) {
  DATA_VECTOR(y);
  DATA_MATRIX(X);
  DATA_INTEGER(ages);
  ages = ages-1;
  DATA_VECTOR(size);

  //Penalization stuff
  DATA_VECTOR(lambda);
  DATA_MATRIX(P);

  PARAMETER_VECTOR(beta);

  Type nll = 0.0;

  vector<Type> eta = X*beta;

  vector<Type> deviance(eta.size());
  for(int i =0;i < eta.size();i++){
    deviance(i) = dbinom_robust(y(i),size(i),eta(i),true);
    nll -= deviance(i);
    deviance(i) = Type(-2.0)*deviance(i);
  }


  nll += Ridge(lambda,P,beta);
  
  vector<Type> devres = deviance_resids(deviance,y);
  vector<Type> pearres = pearson_resids(eta,y);

  REPORT(eta);
  REPORT(beta);

  return nll;

}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

  
