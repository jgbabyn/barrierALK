template<class Type>
matrix<Type> makeXr(vector<Type> rho,Type sd, matrix<Type> ranX,matrix<int> ageTimeKey){
  matrix<Type> Xr = ranX;
  Type pc_age = sqrt(1-pow(rho(0),2));
  Type pc_year = sqrt(1-pow(rho(1),2));

  int A = ageTimeKey.cols();
  int Y = ageTimeKey.rows();

  Type s = 0;
  int i = 0;
  matrix<Type> m(Xr.rows(),1);
  for(int y = 0; y < Y; y++){
    for(int a = 0; a < A; a++){
      int i = ageTimeKey(y,a);
      if((y == 0) & (a == 0)){
	s = sd/(pc_age*pc_year);
	Xr.col(i) = Xr.col(i)*s;
      }else if((a > 0) & (y == 0)){
	int aol = ageTimeKey(y,a-1);
	m.col(0) = rho(0)*Xr.col(aol);
	s = sd/pc_year;
	Xr.col(i) = Xr.col(i)*s+m.col(0);
      }else if((a == 0) & (y > 0)){
	int yol = ageTimeKey(y-1,a);
	m.col(0) = rho(1)*Xr.col(yol);
	s = sd/pc_age;
	Xr.col(i) = Xr.col(i)*s+m.col(0);
      }else{
	int aol = ageTimeKey(y,a-1);
	int yol = ageTimeKey(y-1,a);
	int ayol = ageTimeKey(y-1,a-1);
	m.col(0) = rho(1)*Xr.col(yol)+rho(0)*(Xr.col(aol)-rho(1)*Xr.col(ayol));
	s = sd;
	Xr.col(i) = Xr.col(i)*s+m.col(0);
      }
    }
  }
	

  return Xr;
}

  
  

template<class Type>
Type rhoTrans(Type x){
  return Type(2)/(Type(1) + exp(-Type(2)*x))-Type(1);
}

VECTORIZE1_t(rhoTrans);

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template <class Type>
Type spatialAR(objective_function<Type>* obj){
  using namespace R_inla_barrier;
  using namespace density;

  DATA_VECTOR(y);
  DATA_MATRIX(X);
  DATA_SPARSE_MATRIX(A);
  DATA_FACTOR(cohort);
  DATA_STRUCT(fem,fem_barrier_t);
  DATA_SCALAR(range_fraction);
  DATA_INTEGER(ages);
  DATA_VECTOR(size);
  ages = ages-1;
  DATA_FACTOR(year);

  DATA_IMATRIX(ageTimeKey);

  DATA_VECTOR(lambda);
  DATA_MATRIX(P);
  

  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(rhoT);
  PARAMETER(log_range);
  PARAMETER(log_sigma_u);
  PARAMETER_MATRIX(ranX);

  vector<Type> ranges(2);
  ranges(0) = exp(log_range);
  ranges(1) = ranges(0)*range_fraction;
  Type sigma_u = exp(log_sigma_u);
  vector<Type> rho = rhoTrans(rhoT);

  Type nll = 0.0;

  SparseMatrix<Type> Q = Q_barrier(fem,ranges,Type(1.0));

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
  REPORT(log_range);
  REPORT(log_sigma_u);
  
  nll += Ridge(lambda,P,beta);
  

  return nll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
