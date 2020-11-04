namespace barrier_pred{
  using namespace Eigen;
  using namespace tmbutils;
  
template<class Type>
matrix<Type> predictProbs(matrix<Type> X,vector<Type> beta, SparseMatrix<Type> A, matrix<Type> ranX)
{
  //add one for the plus group
  matrix<Type> etaMat(X.rows(),ranX.cols());

  for(int i = 0; i < etaMat.cols();i++){
    etaMat.col(i) = X*beta;
    etaMat.col(i) = etaMat.col(i) + A*ranX.col(i);
  }

  matrix<Type> probs(etaMat.rows(),etaMat.cols()+1);
  probs.col(0) = etaMat.col(0);
  
  for(int i = 1;i < probs.cols();i++){
    for(int j = 0;j < etaMat.rows();j++){
      Type pi = etaMat(j,i);
      for(int k = i-1;k >= 0;k--){
	pi *= 1 - etaMat(j,k);
      }
      probs(j,i) = pi;
    }
  }

  return probs;
}

}
