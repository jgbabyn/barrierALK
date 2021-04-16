//Report a vector of sparse matrices as a list
#define REPORT_SMV(name)						\
  if(isDouble<Type>::value &&						\
     TMB_OBJECTIVE_PTR -> current_parallel_region < 0)			\
    {									\
      SEXP _TMB_temporary_sexp_ = PROTECT(allocVector(VECSXP,name.size())); \
      for(int i = 0; i < name.size(); i++){				\
	SEXP temp = asSEXP(name(i));					\
	SET_VECTOR_ELT(_TMB_temporary_sexp_,i,temp);			\
      }									\
      Rf_defineVar(Rf_install(#name),_TMB_temporary_sexp_,		\
			      TMB_OBJECTIVE_PTR -> report);		\
      UNPROTECT(1);							\
    }									\
  

  
namespace R_inla_barrier {
  using namespace Eigen;
  using namespace tmbutils;



  //Probably a more elegant way to do this using other structs
  template<class Type>
  struct fem_barrier_t{
    SparseMatrix<Type> I; // I matrix from inla.barrier.fem (J in Bakka 2018)
    SEXP CList; // C matrix list from inla.barrier.fem
    vector<vector<Type> > C; //The C++ vector of the C matrix
    SEXP DList; // D matrix list from inla.barrier.fem
    vector<SparseMatrix<Type> > D; //C++ vector of D matrix
    fem_barrier_t(SEXP x){
      I = asSparseMatrix<Type>(getListElement(x,"I"));
      CList = getListElement(x,"C");
      C.resize(LENGTH(CList));
       for(int i = 0;i < LENGTH(CList); i++){ 
        	SEXP v = VECTOR_ELT(CList,i);
		C(i) = asVector<Type>(v);
       }

      DList = getListElement(x,"D");
      D.resize(LENGTH(DList));
      for(int i = 0;i < LENGTH(DList); i++){
	SEXP sm = VECTOR_ELT(DList,i);
	if(!isValidSparseMatrix(sm))
	  error("Not a sparse matrix");
	D(i) = asSparseMatrix<Type>(sm);
      }
  }
  };

  //The actual function to handle the generation of Q for the GMRF...
  //Basically a C++ version of inla.barrier.q
  template<class Type>
  SparseMatrix<Type> Q_barrier(fem_barrier_t<Type>& fem, vector<Type> ranges, Type sigma){
    int xi = ranges.size();

    vector<Type> Cdiag = pow(ranges[0],2)*fem.C(0);
    if(xi > 1){
      for(int k = 1; k < xi;k++){
	Cdiag = Cdiag + pow(ranges[k],2)*fem.C(k);
      }
    }

    int N = Cdiag.size();
    vector<Type> CdiagInv = 1/Cdiag;
    SparseMatrix<Type> Cinv(N,N);
    Cinv.reserve(VectorXi::Constant(N,1));
    for(int i = 0;i < N;i++){
      Cinv.insert(i,i) = CdiagInv(i);
    }
    Cinv.makeCompressed();
    
    SparseMatrix<Type> A = fem.I;
    for(int k = 0;k < xi;k++){
      Type raiseRange = pow(ranges[k],2)/8;
      A +=  raiseRange*fem.D(k);
    }
    //M_PI is the pi constant defined in math.h which is probably loaded somewhere along the line...
    Type restOfStory = (1/pow(sigma,2))/(M_PI/2)*3; //Change a 4 to a 3!
    SparseMatrix<Type> Q = A.transpose()*Cinv*A*restOfStory;
    return Q;
  }
  

  
}




    
