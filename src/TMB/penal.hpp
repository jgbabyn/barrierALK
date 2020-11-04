//' Penalize the likelihood using ridge regression
//'
//' Penalize the likelihood using ridge regression of the form
//' \frac{1}{2}\lambda\beta P \beta^T.
//'
//' @param lambda vector of penalty factors
//' @param P the penalty matrix
//' @param beta the vector of parameters to use
template<class Type>
Type Ridge(vector<Type> lambda,matrix<Type> P,vector<Type> beta){
  matrix<Type> pf = lambda.sqrt().matrix().asDiagonal();
  P = pf*P*pf;
  matrix<Type> betaT = beta.segment(1,beta.size()-1).matrix().transpose();
  matrix<Type> Pbeta = P*betaT.transpose(); //ugh
  matrix<Type> Penalty = 0.5*betaT*Pbeta;
  Type penalty = Penalty(0,0);
  return penalty;
}
  
