#'Barrier ALK, fit a barrier ALK model to age-growth data to generate ALKs
#'
#' Fit smooth age-length keys with the possibility of incorporating spatial error terms represented by
#' approximations of Guassian Fields using SPDEs
#'
#' @param formula formula for the barrier alk fit
#' @param lengthVar the name of the length variable in the data
#' @param equalSlopes TRUE or FALSE for equal slopes assumption, or formla specifiying which are
#' @param data optionally data frame containing data
#' @param silent whether to run the trace from TMB
#' @param runSA Apply Symbolic Analysis from TMB package in effort to up model
#' @param control control options for optimr optimizer
#' @param optMethod the optimizer within optimr to use
#' @param lambda the value or vector for penalization factors, default is 0, no penalization
#' @param penalty.matrix manually set the penalty matrix to use
#' @param penalize.hyper penalize the hyper parameters for the gf approximations?
#' @param starting.param named list to replace default starting parameter values with
#'
#' @useDynLib spatialBarrier
#' @useDynLib barrierALK
#' @importFrom Rcpp sourceCpp
#' 
#' @export
barrierALK <- function(formula,minAge,maxAge,lengthVar="length",equalSlopes=FALSE,data=NULL,silent=TRUE,runSA=TRUE,control=NULL,...,optMethod="nlminb",cores=parallel::detectCores(),lambda=0,penalty.matrix=NULL,penalize.hyper=FALSE,starting.param=NULL,scale=FALSE,edf.calc=TRUE){
    
    call <- match.call()
    formula <- as.formula(formula)
    tmbD <- interpret.salk(formula,equalSlopes,data,minAge,maxAge,pred=FALSE,lambda=lambda,penalty.matrix=penalty.matrix,penrand=penalize.hyper,starting.param=starting.param,scale=scale)

    if(tmbD$data$model == 0){
        obj <- TMB::MakeADFun(tmbD$data,tmbD$parm,DLL="spatialBarrier",silent=silent,...)
    }else{
        if(!is.null(tmbD$map)){
            obj <- TMB::MakeADFun(tmbD$data,tmbD$parm,random="ranX",DLL="spatialBarrier",silent=silent,map=tmbD$map,...)
        }else{
            obj <- TMB::MakeADFun(tmbD$data,tmbD$parm,random="ranX",DLL="spatialBarrier",silent=silent,...)
        }
                if(runSA == TRUE){
            
            if(exists('runSymbolicAnalysis',where=asNamespace('TMB'),mode='function')){
                TMB::runSymbolicAnalysis(obj)
            }else{
                warning("runSymbolicAnalysis not installed, see TMB install documentation for details.")
            }
        }

    }

    if(optMethod == "parallel"){
        cl <- parallel::makeCluster(cores,type="FORK")
        parallel::setDefaultCluster(cl=cl)
        opt = optimParallel::optimParallel(obj$par,obj$fn,obj$gr,control=control)
        parallel::setDefaultCluster(cl=NULL)
        parallel::stopCluster(cl)
    }else{
        opt <- optimx::optimr(obj$par,obj$fn,obj$gr,control=control,method=optMethod)
    }

    if(edf.calc==TRUE){    
        EDFl <- edf(tmbD,opt$par,obj$env$last.par.best)
        EDF <- EDFl$edf
        LR <- tmbD$other$nullDev + EDFl$dev0
        mAIC <- LR - 2*EDF

        
    }

        
    ret <- list()
    ret$opt <- opt
    ret$obj <- obj
    ret$rep <- obj$report()
    ret$yrep <- tmbD$other$repped$y
    ret$repind <- tmbD$other$repped$subs
    ret$response <- tmbD$other$response
    ret$age <- tmbD$other$age
    class(ret) <- "barrierALK"
    attr(ret,"formula") <- formula
    attr(ret,"data") <- data
    attr(ret,"nobs") <- nrow(data) ##Right now doesn't support weights
    attr(ret,"ageV") <- attr(tmbD,"ageV")
    attr(ret,"equalSlopes") <- equalSlopes
    attr(ret,"lengthVar") <- lengthVar
    attr(ret,"call") <- call
    attr(ret,"Xnames") <- colnames(tmbD$data$X)
    attr(ret,"lambda") <- lambda
    attr(ret,"penalty.matrix") <- penalty.matrix
    attr(ret,"stdobj") <- tmbD$other$stdobj
    attr(ret,"scale") <- scale
    attr(ret,"tmbD") <- tmbD
    attr(ret,"edf") <- ifelse(edf.calc,EDF,length(opt$par))
    attr(ret,"nullDev") <- tmbD$other$nullDev
    attr(ret,"LR") <- LR
    attr(ret,"mAIC") <- mAIC
    print(opt$convergence)
    print(opt$message)
    ret
}

#'Get the log likelihood of a barrierALK fit
#'
#' @param object a fit barrierALK object
#'
#' @export
logLik.barrierALK <- function(object,...)
{
    val = -object$opt$value
    attr(val,"df") = attr(object,"edf")
    attr(val,"nobs") <- attr(object,"nobs")
    class(val) = "logLik"
    val

}


Rpred <- function(X,beta,ranX,A,cohort,ranXind){
    eta = X%*%beta
    AXs = apply(ranX,2,function(x) A%*%x)

    ind = ranXind[cohort+1]+1
    
    eta = sapply(AXs,function(x) as.matrix(eta + x))
    eta[cbind(seq_along(ind),ind)]

}

#'Function to predict from a barrierALK fit
#'
#' This will predict from a barrierALK object with three possible types,
#' the default response will return the unconditional probability of an
#' observation belonging to a class, link the etas and class the predicted
#' class of that observation.
#'
#' @param object the barrierALK fit object
#' @param newdata optionally the newdata to predict on
#' @param type the type of prediction
#'
#' @export
predict.barrierALK <- function(object,newdata=NULL,type="response",...)
{
    if(is.null(newdata)){
        newdata = attr(object,"data")
    }
    form <- attr(object,"formula")
    ages <- attr(object,"ageV")
    
    if(data.table::is.data.table(newdata)){
        newdata = as.data.frame(newdata)
    }

    scale = attr(object,"scale")
    stdobj = attr(object,"stdobj")

    ##Dealing with the the special
    tf <- terms.formula(form,specials="gf")
    terms <- attr(tf,"term.labels")
    gf <- attr(tf,"specials")$gf
    gfInd <- grep("gf\\(",terms)

    if(length(gf) != 0){
        gfFun <- terms[gfInd]
        gfFun <- parse(text=gfFun)
        mcall <- match.call(barrierALK::gf,gfFun[[1]])
        #mcall$loc_x[2] <- quote(newdata())
        #mcall$loc_y[2] <- quote(newdata())
        pgf = mcall
    }else{
        pgf = NULL
    }
    
    
    
    equalSlopes <- attr(object,"equalSlopes")
    predy <- rep(max(ages),nrow(newdata))
    predy <- factor(predy,levels=sort(ages))
    ls = formula.tools::lhs.vars(form)
    newdata[,ls] = predy

    MINAGE = min(ages)
    MAXAGE = max(ages)

    lambda = attr(object,"lambda")
    penalty.matrix = attr(object,"penalty.matrix")

    pred = interpret.salk(form,equalSlopes,newdata,minAge=MINAGE,maxAge=MAXAGE,pgf=pgf,pred=TRUE,lambda=lambda,
                          penalty.matrix=penalty.matrix,scale=scale,stdobj=stdobj)

    if(type == "response" | type == "class"){
        pred$data$p_type = 1
    }else if(type == "link"){
        pred$data$p_type = 0
    }

    pred$data$ages <- length(unique(ages))
    X = pred$data$X
    ages = pred$data$ages-1
    rep <- object$obj$report()
    beta <- rep$beta
    
    if(pred$data$model == 0){
        ##etas = X%*%beta
        ret = predict_fullB(X,beta,ages,pred$data$p_type)
    }else if(pred$data$model == 1 || pred$data$model == 2){
        ranX = rep$ranX
        A = pred$data$A
        cohort = as.integer(pred$data$cohort)-1L
        ranXind = as.integer(pred$data$ranXind)-1L
        ##etas = Rpred(X,beta,ranX,A,cohort,ranXind)
        ret = predict_full(X,beta,ranX,A,cohort,ranXind,ages,pred$data$p_type)
    }else if (pred$data$model == 3 | pred$data$model == 4){
        Xr = rep$Xr
        A = pred$data$A
        cohort = as.integer(pred$data$cohort)-1L
        year = as.integer(pred$data$year)-1L
        ageTimeKey = pred$data$ageTimeKey
        ret = predict_full_AR(X,beta,Xr,A,cohort,year,ageTimeKey,ages,pred$data$p_type)
    }
    

    ##ret = predict_probs(etas,ages,pred$data$p_type)
    

    if(type == "class"){
        ret = apply(ret,1,which.max)
    }

    ret
    
}




#'Print a barrierALK object
#'
#'Print a barrierALK object
#'
#' @param object the barrierALK object to print
#' @param digits the number of digits to round to
#'
#' @export
print.barrierALK <- function(object,digits = 3L,...){
    call <- attr(object,"call")
    cat("\nCall:  ",
        paste(deparse(call), sep = "\n", collapse = "\n"), "\n\n",sep="")
    if(length(object$opt$par)) {
        cat("Coefficients:\n")
        coeff <- object$opt$par
        betN <- attr(object,"Xnames")
        betL <- grep("beta",names(coeff))
        rhoL <- grep("rhoT",names(coeff))
        if(!is.null(rhoL)){
            names(coeff)[rhoL] <- "rho"
            coeff[rhoL] <- object$rep$rho
        }
        names(coeff)[betL] <- betN
        print.default(format(coeff,digits=digits),
                      print.gap = 2, quote=FALSE)
    }else cat("No coefficients\n\n")
    lik <- logLik(object)
    df <- attr(lik,"df")
    nobs <- attr(lik,"nobs")
    aic <- AIC(object)
    bic <- BIC(object)
    cat("\nDegrees of Freedom:",df,"\nNumber of Observations:",nobs,"\nAIC:",format(signif(aic,digits)),"\nBIC:",format(signif(bic,digits)))
    cat("\n")
    invisible(object)
}


#'Summary of a barrierALK object fit
#'
#'Generate a summary of a barrierALK object. This uses and stores a copy of running
#' TMB::sdreport() on the fit object. This can be a little slow for for some models on the first run.
#' Gives a summary only of the fixed effects.
#'
#'@param object the barrierALK object to generate a summary of
#'
#' @export
#'
summary.barrierALK <- function(object){
    sdr <- attr(object,"sdreport")
    if(is.null(sdr)){
        sdr <- TMB::sdreport(object$obj)
        attr(object,"sdreport") <- sdr
    }
    summ <- summary(sdr,select="fixed",p.value=TRUE)
    betL <- grep("beta",rownames(summ))
    betN <- attr(object,"Xnames")
    rownames(summ)[betL] <- betN

    ret <- list()
    ret$coeff <- summ
    ret$AIC <- AIC(object)
    ret$BIC <- BIC(object)
    ret$lik <- logLik(object)
    ret$df <- attr(ret$lik,"df")
    ret$nobs <- attr(ret$lik,"nobs")
    ret$call <- attr(object,"call")
    class(ret) <- "summary.barrierALK"
    ret

}


#'Print a summary of a barrierALK object fit
#'
#' Print a summary of a barrierALK object fit
#'
#' @param object the summary.barrierALK object to print
#'
#' @export
#'
print.summary.barrierALK <- function(object,digits=3L){
    call <- object$call
    cat("\nCall:\n",
        paste(deparse(call),sep="\n",collapse="\n"),"\n\n",sep="")
    if(nrow(object$coeff)){
        cat("Coefficients:\n")
        printCoefmat(object$coeff,digits)
    }else{
        cat("No Coefficients\n")
    }
    cat("\nDegrees of Freedom:",object$df,"\nNumber of Observations:",object$nobs,"\nAIC:",format(signif(object$AIC,digits)),"\nBIC:",format(signif(object$BIC,digits)),"\nLog Likelihood:",format(signif(object$lik,digits)))
    cat("\n")
    invisible(object)
}



#'LRT testing for barrierALK objects
#'
#' @param object barrierALK object to test
#' @param ... other barrierALK objects to test
#'
#' @export
#'
anova.barrierALK <- function(object,...,dispersion=NULL){
    dotargs <- list(...)
    named <- if(is.null(names(dotargs)))
                 rep(FALSE,length(dotargs))
             else(names(dotargs) != "")
    if(any(named))
        warning("The following arguments are invalid and dropped: ",
                paste(dparse(dotargs[named]),collapse=","))
    dotargs <- dotargs[!named]
    is.barrierALK <- unlist(lapply(dotargs,function(x) inherits(x,"barrierALK")))
    dotargs <- dotargs[is.barrierALK]

    ##Dispersion is assumed to be 1 since this is a binomial model
    if(is.null(dispersion))
        dispersion = 1

    if(length(dotargs) > 0){
        n <- attr(logLik(object),"nobs")
        objs <- c(list(object),dotargs)
        ns <- sapply(objs,function(x) attr(logLik(x),"nobs"))
        if(any(ns != n))
            stop("Models do not have the same number of observations")
        resdf <- sapply(objs,function(x) attr(logLik(x),"nobs")-attr(logLik(x),"df"))
        resdev <- sapply(objs,function(x) -2*as.numeric(logLik(x)))

        table <- data.frame(resdf,resdev,c(NA,-diff(resdf)),
                            c(NA,-diff(resdev)))
        variables <- lapply(objs,function(x)
            paste(deparse(attr(x,"formula")),collapse="\n"))
        dimnames(table) <- list(1:length(objs),c("Resid. Df", "Resid. Dev", "Df",
                                                 "Deviance"))
        title <- "Analysis of Deviance Table\n"
        topnote <- paste0("Model ", format(1:length(objs)),": ",
                          variables,collapse="\n")
        df.dispersion <- ifelse(dispersion==1,Inf,min(resdf))
        table <- stat.anova(table=table,test="Chisq",scale=dispersion,df.scale=df.dispersion,n=NULL)

        structure(table,heading = c(title,topnote),
                  class = c("anova","data.frame"))

    }
    else{
        print(summary(object))
    }
}


#'Age fish from long length frequencies
#'
#' Age fish into a long set of age frequencies instead
#'
#' @param object the fit barrierALK object
#' @param longLF the long length frquencies with covariates
#' @param lfname name of the length frequency column
#'
#' @export
#'
ageFish <- function(object,longLF,lfname="length_freq"){
    probs = predict(object,newdata=longLF)

    longLFdf = as.data.frame(longLF)
    
    aged = longLFdf[,lfname]*probs
    aged = as.data.frame(aged)
    names(aged) = paste0("age",sort(attr(object,"ageV")))
    wide = cbind(longLFdf,aged)

    long = tidyr::gather(wide,"age","n",paste0("age",sort(attr(object,"ageV"))))
    long$age = as.numeric(gsub("age","",long$age))
    long
}

        
#'Function to create a smooth ALK from a barrierALK object
#'
#' @param object The barrierALK fit object
#' @param lengthBins vector containing the length bins you want to generate an age-length key for
#' @param setData data you want to apply an ALK to (must contain covariates with same name as in in barrierALK fit)
#'
#' @export
createALKs <- function(object,lengthBins,setData){
    newDat <- reshape::expand.grid.df(as.data.frame(lengthBins),setData)
    lvar <- attr(object,"lengthVar")
    names(newDat)[1] <- lvar
    ALK <- predict(object,newdata=newDat)
    ALK <- split.data.frame(ALK,rep(1:nrow(setData),each=length(lengthBins)))

    ALKs <- list()

    for(i in 1:nrow(setData)){
        ALKs[[i]] = ALK[[i]]
        rownames(ALKs[[i]]) <- lengthBins
        colnames(ALKs[[i]]) <- sort(attr(object,"ageV"))
        class(ALKs[[i]]) <- "predALK"
    }
    ALKs
}

#'Perform k-fold cross validation to try and determine the best value of lambda
#'
#' @param formula the barrierALK formula to use
#' @param k the number of cross validation folds to use, at least 10 recommended
#' @param minAge the minimum age to be included in the model
#' @param maxAge the maximum age to be included in the model as a plus group
#' @param data the data frame containing the data
#' @param lengthVar the name of the length variable in the data
#' @param equalSlopes TRUE or FALSE for equal slopes assumption, or formla specifiying which are
#' @param silent display trace from TMB?
#' @param runSA Apply Symbolic Analysis from TMB package in effort to up model
#' @param control control options for optimr optimizer
#' @param optMethod the optimizer within optimr to use
#' @param lambda vector of lambda values to try and find the optimal value of
#' @param penalty.matrix manually set the penalty matrix to use
#'
#' @export
#'
barrierALKCV <- function(formula,k=10,minAge,maxAge,data,equalSlopes=FALSE,lengthVar="length",silent=TRUE,runSA=FALSE,control,optMethod="nlminb",lambda,penalty.matrix=NULL,scale=FALSE){
    resp = formula.tools::lhs.vars(as.formula(formula))
 
    folds = data
    folds$.folds = balanced.cv.fold(data[,resp],num.cv=k)

    acc <- list()
    olp <- list()
    for(i in 1:k){
        test = folds[folds$.folds == i,]
        train = folds[folds$.folds != i,]

        acc[[i]] <- list()
        for(j in 1:length(lambda)){
            if(i == 1){
                mod = barrierALK(formula,minAge,maxAge,lengthVar,equalSlopes,train,silent,runSA,control,optMethod=optMethod,lambda=lambda[j],penalty.matrix = penalty.matrix,scale=scale)
                pred = predict(mod,newdata=test,type="class")
                acc[[i]][[j]] = sum(test[,resp] == pred)/length(pred)
                olp[[j]] = list(beta=mod$rep$beta,log_ranges=mod$rep$log_ranges,log_sigma_us=mod$rep$log_sigma_us,log_kappas=mod$rep$log_kappas,log_taus=mod$rep$log_taus,rhoT=mod$rep$rhoT,log_range=mod$rep$log_range,log_sigma_u=mod$rep$log_sigma_u)
            }else{
                mod = barrierALK(formula,minAge,maxAge,lengthVar,equalSlopes,train,silent,runSA,control,optMethod=optMethod,lambda=lambda[j],penalty.matrix = penalty.matrix,starting.param=olp[[j]],scale=scale)
                pred = predict(mod,newdata=test,type="class")
                acc[[i]][[j]] = sum(test[,resp] == pred)/length(pred)
                }
                
        }
    }
    ret = data.frame(acc=unlist(acc),lamb=rep(lambda,k))
    ret = dplyr::group_by(ret,lamb)
    ret = dplyr::summarize(ret,acc=mean(acc))
    ret
    }

##Grabbed from bmrm R package
balanced.cv.fold <- function(y,num.cv=10) {
  fold <- factor(rep_len(sample(seq_len(num.cv)),length.out=length(y)))
  i <- sample(seq_along(y))
  o <- order(y[i])
  fold[i[o]] <- fold
  fold
}

#'Calculate the effective degrees of freedom 
edf <- function(tmbD,parP,fullPar){

    tmbD2 = tmbD
    tmbD2$data$lambda = rep(0,length(tmbD$data$lambda))
    
    if(tmbD$data$model == 0){
        objP <- TMB::MakeADFun(tmbD$data,tmbD$parm,DLL="spatialBarrier",silent=TRUE)
        obj0 <- TMB::MakeADFun(tmbD2$data,tmbD$parm,DLL="spatialBarrier",silent=TRUE)
        
        
    }else{
        if(!is.null(tmbD$map)){
            objP <- TMB::MakeADFun(tmbD$data,tmbD$parm,random="ranX",DLL="spatialBarrier",silent=TRUE,map=tmbD$map)
            obj0 <- TMB::MakeADFun(tmbD2$data,tmbD2$parm,random="ranX",DLL="spatialBarrier",silent=TRUE,map=tmbD$map)
            
        }else{
            objP <- TMB::MakeADFun(tmbD$data,tmbD$parm,random="ranX",DLL="spatialBarrier",silent=TRUE)
            obj0 <- TMB::MakeADFun(tmbD2$data,tmbD2$parm,random="ranX",DLL="spatialBarrier",silent=TRUE)
        }
  
    }

    ##Different way of thinking about it
    obj0$env$last.par.best = fullPar
    objP$env$last.par.best = fullPar
    sdr0 <- TMB::sdreport(obj0,parP)
    sdrP <- TMB::sdreport(objP,parP)

    h1 = solve(sdr0$cov.fixed)
    V = sdrP$cov.fixed
    
    
    ## h1 = obj0$he(parP)
    ## h2 = objP$he(parP)
    ## V = solve(h2)

    dev0 = -2*obj0$env$f(fullPar)


    
    edf = sum(diag(h1%*%V))

    ret = list(edf=edf,sdr=sdrP,dev0=dev0)
    
    ret

    
}

#'Find the optimal lambda penalty value using modified AIC
#'
#' @param model the fit barrierALK model to test
#' @param lambda vector of lambda values to test
#'
#' @export
#'
penaltyOptim <- function(model,lambda){
    call = attr(model,"call")
    call$edf.calc=TRUE
    ret = data.frame(lambda=lambda)
    ret$df = NA
    ret$mAIC = NA
    models = list()
    for(i in 1:length(lambda)){
        call$lambda = lambda[i]
        run = eval(call)
        ret$df[i] = attr(run,"edf")
        ret$mAIC[i] = attr(run,"mAIC")
        models[[i]] = run
    }
    bestModel = which.max(ret$mAIC)
    bmod = models[[bestModel]]
    rret = list(table=ret,best.model=bmod)
    class(rret) <- "penopt"
    rret
    
}

#'Print function to print the optimal value table and model
#'
#' @param object the object to print
print.penopt <- function(object,...){
    call <- attr(object$best.model,"call")
    cat("\nBest Model Penalty Factor:\n",
        object$table$lambda[which.max(object$table$lambda)])
    cat("\nBest Model Call:\n",
        paste(deparse(call),sep="\n",collapse="\n"),"\n\n",sep="")
    object$best.model
    cat("\nResults Table:\n")
    object$table
}

    
