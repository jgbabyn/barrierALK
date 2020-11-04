#include <Rcpp.h>
#include <RcppEigen.h>


// [[Rcpp::depends(RcppEigen)]]

//' quickPredict
//'
//' Quickly predict the output of a spatial ALK model
//'
//' @param X the design matrix
//' @param beta the vector fixed effect parameters
//' @param A the A prediction matrix
//' @param ranX the predicted random effects
//' @param cohort group cohort of CRL
//' @param ranXind convert cohort to ranX cols
// [[Rcpp::export]]
Eigen::VectorXd quickPredict(const Eigen::Map<Eigen::MatrixXd> X,const Eigen::Map<Eigen::VectorXd> beta,const Eigen::Map<Eigen::MatrixXd> ranX,const Eigen::MappedSparseMatrix<double> A,const Rcpp::IntegerVector cohort, const Rcpp::IntegerVector ranXind){
  Eigen::VectorXd eta = X*beta; 
  std::vector<Eigen::VectorXd> AXs(ranX.cols());

  for(int i = 0; i < ranX.cols();i++){
    AXs[i] = A*ranX.col(i);
  }

   Eigen::VectorXd AX(A.rows());
   AX = AXs[0];

   Rcpp::IntegerVector inds = ranXind[cohort];

   for(int i = 0; i < eta.size();i++){
     eta(i) = eta(i) + AXs[inds[i]](i);
   }
   
  return eta;

}

// [[Rcpp::export]]
Eigen::VectorXd quickPredictAR(const Eigen::Map<Eigen::MatrixXd> X,const Eigen::Map<Eigen::VectorXd> beta,const Eigen::Map<Eigen::MatrixXd> Xr,const Eigen::MappedSparseMatrix<double> A, const Rcpp::IntegerVector cohort,const Rcpp::IntegerVector year,
			       const Rcpp::IntegerMatrix ageTimeKey){
  Eigen::VectorXd eta = X*beta; 
  std::vector<Eigen::VectorXd> AXs(Xr.cols());

  for(int i = 0; i < Xr.cols();i++){
    AXs[i] = A*Xr.col(i);
  }

  for(int i = 0;i < eta.size();i++){
    int ya = ageTimeKey(year[i],cohort[i]);
    eta(i) = eta(i) + AXs[ya](i);
  }

  return eta;
}

  

// [[Rcpp::export]]
Eigen::VectorXd quickPredictB(const Eigen::Map<Eigen::MatrixXd> X,const Eigen::Map<Eigen::VectorXd> beta){
  Eigen::VectorXd eta = X*beta;
  return eta;
}

enum valid_ptype {
		  link = 0,
		  response = 1
};


// [[Rcpp::export]]
Eigen::MatrixXd predict_probs(Eigen::VectorXd etas,int ages, int p_type){
  int rows = etas.size()/ages;
  Eigen::MatrixXd etaMat(rows,ages);

  int j = -1;
  for(int i = 0; i < etas.size();i++){
    if((i % ages) == 0){
      j++;
    }
    etaMat(j,i % ages) = etas(i);
  }

Eigen::MatrixXd probs(etaMat.rows(),etaMat.cols()+1);

 for(int i = 0;i < probs.rows();i++){
   probs(i,0) = R::plogis(etaMat(i,0),0.0,1.0,1,0);
 }

 for(int i = 0; i < probs.cols()-1;i++){
   for(int j = 0;j < probs.rows();j++){
     double pi = R::plogis(etaMat(j,i),0.0,1.0,1,0);
     Eigen::VectorXd pis = probs.row(j).segment(0,(i-1)+1);
     double summing = pis.sum();
     probs(j,i) = pi*(1-summing);
   }
 }

for(int i = 0;i < probs.rows();i++){
  Eigen::VectorXd pis = probs.row(i).segment(0,probs.cols()-1);
  probs(i,probs.cols()-1) = 1-pis.sum();
 }
Eigen::MatrixXd ret;

if(p_type == link){
  ret = etaMat;
 }else if(p_type == response){
  ret = probs;
 }
 
return ret;

}

// [[Rcpp::export]]
Eigen::MatrixXd predict_full(const Eigen::Map<Eigen::MatrixXd> X,const Eigen::Map<Eigen::VectorXd> beta,
			     const Eigen::Map<Eigen::MatrixXd> ranX, const Eigen::MappedSparseMatrix<double> A, const Rcpp::IntegerVector cohort,const Rcpp::IntegerVector ranXind,int ages, int p_type){
  Eigen::VectorXd etas = quickPredict(X,beta,ranX,A,cohort,ranXind);
  Eigen::MatrixXd ret = predict_probs(etas, ages, p_type);
  return ret;
}

// [[Rcpp::export]]
Eigen::MatrixXd predict_full_AR(const Eigen::Map<Eigen::MatrixXd> X,const Eigen::Map<Eigen::VectorXd> beta,const Eigen::Map<Eigen::MatrixXd> Xr,const Eigen::MappedSparseMatrix<double> A, const Rcpp::IntegerVector cohort,const Rcpp::IntegerVector year,const Rcpp::IntegerMatrix ageTimeKey,int ages, int p_type){
  Eigen::VectorXd etas = quickPredictAR(X,beta,Xr,A,cohort,year,ageTimeKey);
  Eigen::MatrixXd ret = predict_probs(etas,ages,p_type);
  return ret;
}


// [[Rcpp::export]]
Eigen::MatrixXd predict_fullB(const Eigen::Map<Eigen::MatrixXd> X,const Eigen::Map<Eigen::VectorXd> beta,int ages,int p_type){
  Eigen::VectorXd etas = quickPredictB(X, beta);
  Eigen::MatrixXd ret = predict_probs(etas, ages, p_type);
  return ret;
}

// // [[Rcpp::export]]
// Eigen::MatrixXd age_fish(const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::VectorXd> beta,
// 			 const Eigen::Map<Eigen::MatrixXd> ranX, const Eigen::MappedSparseMatrix<double> A,
// 			 const Eigen::Map<Eigen::VectorXi> cohort, const Eigen::Map<Eigen::VectorXi> ranXind,
// 			 int ages,const Eigen::Map<Eigen::MatrixXd> lf,int lblen){
//   Eigen::VectorXd etas = quickPredict(X, beta, ranX, A, cohort, ranXind);
//   Eigen::MatrixXd probs = predict_probs(etas,ages, 1);

//   Eigen::MatrixXd aged(lf.rows(),probs.cols());

//   for(int i = 0;i < lf.rows();i++){
//     int ALKcount = i*lblen;
//     Eigen::MatrixXd ALK = probs.block(ALKcount,0,lblen,probs.cols());
//     aged.row(i) = lf.row(i)*ALK;
//   }

//   return aged;
// }


// [[Rcpp::export]]
Eigen::MatrixXd age_fishB(const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::VectorXd> beta,
			  int ages, const Eigen::Map<Eigen::MatrixXd> lf, int lblen){
  Eigen::VectorXd etas = quickPredictB(X, beta);
  Eigen::MatrixXd probs = predict_probs(etas,ages, 1);

  Eigen::MatrixXd aged(lf.rows(),probs.cols());

  for(int i = 0;i < lf.rows();i++){
    int ALKcount = i*lblen;
    Eigen::MatrixXd ALK = probs.block(ALKcount,0,lblen,probs.cols());
    aged.row(i) = lf.row(i)*ALK;
  }

  return aged;
}
  
  

  
