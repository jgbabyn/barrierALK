#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template <class Type>
Type spatialo(objective_function<Type>* obj){
  using namespace R_inla_barrier;
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

  DATA_STRUCT(fem,fem_barrier_t);
  DATA_SCALAR(range_fraction);

  PARAMETER_VECTOR(log_ranges);
  PARAMETER_VECTOR(log_sigma_us);

  //vector<vector <Type> > AXs(Qlinks.rows());
  vector<SparseMatrix <Type> > Qs(ranX.cols());

  Type nll = 0.0;
  vector<Type> eta = X*beta;
  
  for(int i = 0;i < log_ranges.size();i++){
    vector<Type> ranges(2);
    ranges(0) = exp(log_ranges(i));
    ranges(1) = ranges(0)*range_fraction;

    Qs(i) = Q_barrier(fem,ranges,Type(1.0));
  }

  for(int i = 0;i < Qlinks.rows();i++){
    Type sigma_u = exp(log_sigma_us(Qlinks(i,1)));
    Type sig2inv = 1/(pow(sigma_u,2));

    SparseMatrix<Type> Q = Qs(Qlinks(i,0))*sig2inv;
    nll += GMRF(Q)(ranX.col(Qlinks(i,0)));
    //AXs(i) = A*ranX.col(Qlinks(i,0));
  }

  vector<int> inds = ranXind(cohort);

  vector<Type> deviance(eta.size());
  for(int i = 0;i < eta.size();i++){
    //eta(i) = eta(i) + AXs(inds(i))(i);
    Type ax = A.row(i)*ranX.col(inds(i));
    eta(i) = eta(i) + ax;
    deviance(i) = dbinom_robust(y(i),size(i),eta(i),true);
    nll -= deviance(i);
    deviance(i) = Type(-2.0)*deviance(i);
  }

  vector<Type> devres = deviance_resids(deviance,y);
  vector<Type> pearres = pearson_resids(eta,y);
  REPORT(devres);
  REPORT(pearres);

  vector<Type> fullbeta(beta.size()+log_ranges.size()+log_sigma_us.size());
  fullbeta.segment(0,beta.size()) = beta;
  fullbeta.segment(beta.size(),log_ranges.size()) = log_ranges;
  fullbeta.segment(beta.size()+log_ranges.size(),log_sigma_us.size()) = log_sigma_us;
  
  nll += Ridge(lambda,P,fullbeta);

  REPORT(eta);
  REPORT(beta);
  REPORT(P);
  REPORT(ranX);
  REPORT(log_ranges);
  REPORT(log_sigma_us);
  REPORT_SMV(Qs);

  return nll;


}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
