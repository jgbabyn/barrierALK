#'Function to simulate data from a mesh given certain hyper-parameters
#'
#' This simulates spatial data from a mesh from a specified number of precision matrices. The number of
#' precision matrices is specified by the length of the hyperparameter vectors, which must be equal.
#' 
#' @param mesh The mesh to simulate from
#' @param size the number of locations to sample
#' @param sp_type whether spde or barrier type
#' @param hyperp1 vector of kappas for spde, vector of ranges for barrier
#' @param hyperp2 vector of taus for spde, vector of sigma_us for barrier
#' @param barrier.triangles the triangles defining the barrier for the barrier model
#' @param barrier.polygon the barrier polygon to avoid sampling on the barrier
#' @param seed The seed for the random number generation
#'
#' @export
#'
simFromMesh <- function(mesh,size=5000,sp_type="barrier",hyperp1,hyperp2,barrier.triangles=NULL,barrier.polygon=NULL,seed=NULL){
    if(length(hyperp1) != length(hyperp2))
        stop("hyperp1 and hyperp2 must be the same length")

    loc = data.frame(x=mesh$loc[,1],y=mesh$loc[,2])
    loc = sf::st_as_sf(loc,coords=c("x","y"))
 
    Qs <- list()
    if(sp_type == "spde"){
        spde.model = INLA::inla.spde2.matern(mesh)
        for(i in 1:length(hyperp1)){
            Qs[[i]] = INLA::inla.rgeneric.q(spde.model,cmd="Q",theta=c(log(hyperp1[i]),log(hyperp2[i])))
        }
    loc = sf::st_convex_hull(sf::st_union(loc))
    }else if(sp_type == "barrier"){
        ##Priors don't matter
        barrier.model = INLA::inla.barrier.pcmatern(mesh,barrier.triangles=barrier.triangles,prior.range =c(1.44,0.5),
                                                    prior.sigma = c(0.7,0.5),range.fraction = 0.1)
        for(i in 1:length(hyperp1)){
            Qs[[i]] = INLA::inla.rgeneric.q(barrier.model,cmd="Q",
                                            theta = c(log(hyperp2[i]),log(hyperp1[i])))
        }
        bbox = sf::st_bbox(loc)
        loc = sf::st_as_sfc(bbox)
        loc = sf::st_difference(loc,sf::st_as_sf(barrier.polygon))
    }else{
        stop("Spatial type doesn't exist!")
    }

    uMat = matrix(NA,nrow=nrow(Qs[[1]]),ncol=length(hyperp1))
    for(i in 1:length(hyperp1)){
        if(!is.null(seed)){
            u = INLA::inla.qsample(n=1,Q=Qs[[i]],seed=seed)
        }else{
            u = INLA::inla.qsample(n=1,Q=Qs[[i]])
        }
        
        u = u[,1]
        uMat[,i] = u
    }

    ##Sample from the locations the requested number of times
    samps = sf::st_sample(loc,size,exact=TRUE)
    A.data = INLA::inla.spde.make.A(mesh,sf::as_Spatial(samps))

    aMat = list()
    for(i in 1:ncol(uMat)){
        aMat[[i]] = A.data%*%uMat[,i]
    }
    aMat = do.call(cbind,aMat)
    ret <- list()
    ret$aMat <- aMat
    ret$samps <- samps

    ret
    
    
    }

#'Simulate CRL data with the specified parameters
#'
#' Simulate CRL ratio data with the specified parameters with the option to add
#' spatial data based off of a mesh using simFromMesh.
#'
#' @return A data frame containing the simulated response, any sampled coordinates from the mesh 
 #'
#' @param fParamMatrix the matrix of fixed parameters, rows are classes-1,columns are number of parameters, first column is
#' intercept if allowed
#' @param intercept use an intercept term?
#' @param size number of samples to simulate
#' @param hyperp1 hyperp1 for the spatial simulation
#' @param hyperp2 hyperp2 for the spatial simulation
#' @param mesh mesh for the spatial simulation
#' @param barrier.triangles the triangles defining the barrier for the barrier model
#' @param barrier.polygon the barrier polygon to avoid sampling on the barrier
#' @param sizeL number of spatial locations to sim, if not given, then unique location for each size
#' @param seed the seed to use
#' @param distro the distribution to use for the parameters
#'
#' @export
simCRL <- function(fParamMatrix,intercept=TRUE,size=1000,hyperp1=NULL,hyperp2=NULL,mesh=NULL,barrier.triangles=NULL,
                   barrier.polygon=NULL,sizeL=size,seed=NULL,distro=rnorm){
    
    mus = matrix(NA,nrow=size,ncol=nrow(fParamMatrix))    
    rand = matrix(NA,nrow=size,ncol=length(1:ncol(fParamMatrix)))
    rand = apply(rand,2,distro)
    randD = as.data.frame(rand)
    
    if(intercept == TRUE){
        rand[,1] = 1
        randD$V1 <- NULL
    }

    for(i in 1:nrow(fParamMatrix)){
        swept = sweep(rand,2,fParamMatrix[i,],"*")
        mu = apply(swept,1,sum)
        mus[,i] = mu
    }

    if(!is.null(hyperp1)){
        if(!is.null(barrier.triangles)){
            spat = simFromMesh(mesh,sizeL,sp_type="barrier",hyperp1=hyperp1,hyperp2=hyperp2,barrier.triangles=barrier.triangles,
                               barrier.polygon=barrier.polygon,seed=seed)
        }else{
            spat = simFromMesh(mesh,sizeL,sp_type="spde",hyperp1=hyperp1,hyperp2=hyperp2,barrier.triangles=NULL,
                               barrier.polygon=NULL,seed=seed)
        }

        if(sizeL != size){
            ind <- sample(1:nrow(spat$aMat),size,replace=TRUE) ##I could probably do this a better way to have more realistic hauls
        }else{
            ind <- 1:nrow(spat$aMat)
        }
        
        mus = mus+spat$aMat[ind,]
        coords <- as.data.frame(sf::st_coordinates(spat$samps))[ind,]
        names(coords) <- c("x","y")
        randD$x <- coords$x
        randD$y <- coords$y
            
        
    }

    
        
    pis = apply(mus,2,plogis)
    probs = matrix(NA,nrow=size,ncol=nrow(fParamMatrix)+1)
    probs[,1] = pis[,1]

    for(i in 2:(ncol(probs)-1)){
        if(i > 2){
            probs[,i] = pis[,i]*(1-apply(probs[,1:i-1],1,sum))
        }
        else{
            probs[,i] = pis[,i]*(1-probs[,1])
        }
        
    }
    
    probs[,ncol(probs)] = 1-apply(probs[,1:ncol(probs)-1],1,sum)
    response <- apply(probs,1,function(x) sample(1:ncol(probs),1,prob=x))

    ##The constructed simulated data frame to return
    ret <- cbind(response,randD)
    
    
    
}

#'Residuals for the barrierALK object
#'
#' @param object the barrierALK object to get the residuals of
#' @param type the type of residual to return, either deviance, pearson or RQR,PSR
#' @param ... unused
#'
#' Residuals of the barrierALK object including support for randomized quantile residuals (RQR) and probability scale residuals (PRS).
#' Note that for the deviance and pearson residuals there is one or more residuals for each observation from each level of the CRL they pass through. The replicated response that the pearson or deviance residuals actually refer to can be retrived from the model fit along with the indices of how observations were replicated from the original data using fit$yrep and fit$repind respectively.
#'
#' Randomized Quantile Residuals try to get a single residual for each observation that should be intereptable the same as residuals for a linear model. 
#' 
#' @export
#'
residuals.barrierALK <- function(object,type="deviance",...){
    if(type == "deviance"){
        resids = object$rep$devres
    }else if (type == "pearson"){
        resids = object$rep$pearres
    }else if (type == "RQR"){
        probs = predict(object,type="response")
        probs = cbind(rep(0,nrow(probs)),probs) ##Add a column of zeros for walking along the probs
        cdf = t(apply(probs,1,cumsum))
        age = object$age
        a = cdf[cbind(seq_along(age),age)]
        b = cdf[cbind(seq_along(age),age+1)]
        u = runif(n=length(age),min=a,max=b)
        resids = qnorm(u)
    }else if (type == "PSR"){
        probs = predict(object,type="response")
        cdf = t(apply(probs,1,cumsum))
        age = object$age
        resids = 2*cdf[cbind(seq_along(age),age)]*probs[cbind(seq_along(age),age)]-1
    }
        
    
    resids
}

