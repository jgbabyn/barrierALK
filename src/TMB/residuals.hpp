#ifndef LOGIS_H
#include "logis.hpp"
#endif

template<class Type>
vector<Type> pearson_resids(vector<Type>& etas, vector<Type>& y_obs){
  vector<Type> peares(etas.size());
  for(int i = 0; i < etas.size();i++){
    Type predp = plogis(etas(i),Type(0.0),Type(1.0),1,0);
    Type num = y_obs(i) - predp;
    Type denom = sqrt(predp*(1-predp));
    peares(i) = num/denom;   
  }
  return peares;
}

template<class Type>
vector<Type> deviance_resids(vector<Type>& dev, vector<Type>& y_obs){
  vector<Type> devres(dev.size());
  for(int i = 0; i < dev.size();i++){
    Type signofy = (y_obs(i) == 1) ? 1 : -1;
    devres(i) = signofy*sqrt(dev(i));
  }
  return devres;
}

  
