#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template <class Type>
Type spatialspde(objective_function<Type>* obj){
  using namespace R_inla;
  using namespace density;
  using namespace barrier_pred;

  DATA_VECTOR(y);
  DATA_MATRIX(X);
  DATA_INTEGER(ages);
  ages = ages-1;
  DATA_VECTOR(size);

  DATA_VECTOR(lambda);
  DATA_MATRIX(P);

  PARAMETER_VECTOR(beta);

  DATA_SPARSE_MATRIX(A);
  DATA_FACTOR(cohort);
  DATA_FACTOR(ranXind);
  DATA_IMATRIX(Qlinks);
  PARAMETER_MATRIX(ranX);

  DATA_STRUCT(fem,spde_t);

  PARAMETER_VECTOR(log_kappas);
  PARAMETER_VECTOR(log_taus);

  vector<vector <Type> > AXs(Qlinks.rows());
  vector<SparseMatrix <Type> > Qs(ranX.cols());

  Type nll = 0.0;
  vector<Type> eta = X*beta;

  for(int i = 0;i < log_kappas.size();i++){
    Type kappa = exp(log_kappas(i));

    Qs(i) = Q_spde(fem,kappa);
  }

  for(int i = 0;i < Qlinks.rows();i++){
    Type tau = exp(log_taus(Qlinks(i,1)));

    SparseMatrix<Type> Q = Qs(Qlinks(i,0));
    nll += GMRF(Q)(ranX.col(i));
    AXs(i) = A*ranX.col(i)/tau;
  }

  vector<int> inds = ranXind(cohort);

  vector<Type> deviance(eta.size());
  for(int i = 0;i < eta.size();i++){
    eta(i) = eta(i) + AXs(inds(i))(i);
    deviance(i) = dbinom_robust(y(i),size(i),eta(i),true);
    nll -= deviance(i);
    deviance(i) = Type(-2.0)*deviance(i);
  }

  vector<Type> devres = deviance_resids(deviance,y);
  vector<Type> pearres = pearson_resids(eta,y);
  REPORT(devres);
  REPORT(pearres);

  vector<Type> fullbeta(beta.size()+log_kappas.size()+log_taus.size());
  fullbeta.segment(0,beta.size()) = beta;
  fullbeta.segment(beta.size(),log_kappas.size()) = log_kappas;
  fullbeta.segment(beta.size()+log_kappas.size(),log_taus.size());
  
  nll += Ridge(lambda,P,fullbeta);

  REPORT(eta);
  REPORT(beta);
  REPORT(P);
  REPORT(ranX);
  REPORT(log_kappas);
  REPORT(log_taus);
  REPORT_SMV(Qs);

  return nll;

}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
