#ifndef LOGIS_H
#define LOGIS_H

template <class Type>
Type plogis(Type q, Type location = Type(0.0), Type scale = Type(1.0), bool lower_tail= true,bool log_p  =false){
  Type invlogit = Type(1.0)/(Type(1.0)+exp(-(q-location)/scale));
  if(lower_tail == false){
    invlogit = Type(1.0) - invlogit;
  }
  if(log_p == true){
    invlogit = log(invlogit);
  }
  return invlogit;
}

template <class Type>
Type qlogis(Type p, Type location = Type(0.0), Type scale = Type(1.0), bool lower_tail = true, bool log_p=false){
  if(lower_tail == false){
    p = 1-p;
  }
  if(log_p == true){
    p = exp(p);
  }
  Type logit = location + scale*log(p/(Type(1.0)-p));
  return logit;
}

#endif /*LOGIS_H*/


  
  
			   
