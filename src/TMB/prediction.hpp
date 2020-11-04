#include "logis.hpp"

namespace barrier_pred{

  enum valid_ptype {
		   link = 0,
		   response = 1,
  };

  template<class Type>
  matrix<Type> predict_probs(vector<Type> etas,int ages,int p_type){

    int rows = etas.size()/ages;
    matrix<Type> etaMat(rows,ages);
    int j = -1;
    for(int i = 0; i < etas.size();i++){
      if((i % ages) == 0){
	j++;
       }
      etaMat(j,i % ages) = etas(i);
    }
    
    matrix<Type> probs(etaMat.rows(),etaMat.cols()+1);

    for(int i = 0; i < probs.rows();i++){
      probs(i,0) = plogis(etaMat(i,0),Type(0.0),Type(1.0),true,false);
    }
  
    for(int i = 1;i < probs.cols()-1;i++){
      for(int j = 0;j < probs.rows();j++){
	Type pi = plogis(etaMat(j,i),Type(0.0),Type(1.0),true,false);
	vector<Type> pis = probs.row(j).segment(0,(i-1)+1);
	Type summing = pis.sum();
	probs(j,i) = pi*(1-summing);
      }
    }

    for(int i = 0; i < probs.rows();i++){
      vector<Type> pis = probs.row(i).segment(0,probs.cols()-1);
      probs(i,probs.cols()-1) = 1-pis.sum();
    }
    matrix<Type> ret;
    
    if(p_type == link){
      ret = etaMat;
    }else if(p_type == response){
      ret = probs;
    }else{
      error("Incorrect Prediction Type");
    }
    
    return ret;

  }

  
}
