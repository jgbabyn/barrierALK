#'Interpret a spatial ALK model formula
#'
#'Reads in a spatial ALK formula and sets up the data for loading into the TMB model.
#'
#' @param sf a spatial ALK model formula
#' @param equalSlopes formula or boolean indicating whether to have equal
#' @param data the data
#' @param minAge the minimum age group to fit
#' @param maxAge the maximum age group to fit
#' @param pgf the predicted gaussian field workaround
#' @param pred are we predicting?
#' @param atilde customize the length factor too?
#' @param lambda for penalized smoothing
#' @param scale whether to rescale the fixed effects parameters
interpret.salk <- function(sf,equalSlopes=FALSE,data=NULL,minAge,maxAge,pgf=NULL,pred=FALSE,atilde=NULL,lambda=0,penalty.matrix=NULL,
                           penrand=FALSE,starting.param=NULL,scale,stdobj=NULL){
    if(minAge > maxAge){
        stop("minAge should not be greater than the maxAge")
    }

    ##data.tables don't work exactly like data.frames in all cases and are ugh
    if(data.table::is.data.table(data)){
        data <- as.data.frame(data)
    }
    

    p.env <- environment(sf)
    tf <- terms.formula(sf,specials="gf")
    terms <- attr(tf,"term.labels")
    nt <- length(terms)
    if(attr(tf,"response") == 0){
        stop("Need a response of ages!")
    }else{
        response <- as.character(attr(tf,"variables")[2])
    }

    gf <- attr(tf,"specials")$gf

    if(length(gf) > 1){
        stop("Only one Gaussian Field supported!")
    }
    off <-attr(tf,"offset")
    inter <- attr(tf,"intercept")

    gfInd <- grep("gf\\(",terms)
    if(!is.null(off)){
        terms[nt+1] <- as.character(attr(tf,"variables"))[1+off]
    }
    if(length(gf) == 0){
        linearCovar <- paste(terms,collapse="+")
        linearPred <- paste(response,"~",linearCovar)
        gfList <- NULL
    }else{
        linearCovar <- paste(terms[-gfInd],collapse="+")
        linearPred <- paste(response,"~",linearCovar)
        if(is.null(pgf)){
            gfList <- eval(parse(text = paste0("barrierALK::",terms[gfInd])),envir=p.env)
        }else{
            gfList <- eval(pgf,envir=parent.frame())
        }

    }

    linearPred <- as.formula(linearPred,p.env)
    mf <- model.frame(linearPred,data)
    vars <- all.vars(linearPred)
    missingVars <- !(names(mf) %in% vars)
    if(any(missingVars)){
        mf <- cbind(mf,data[,vars[missingVars],drop=FALSE])
    }
    ret <- list()
    ret$model <- 0

    ##Duplicate the data correctly to do
    age <- mf[,response]

    if(!is.factor(age)){
        age[age > maxAge & !is.na(age)] = maxAge
        age[age < minAge & !is.na(age)] = minAge
    }

    if(!is.null(gfList)){
        ##mf$loc_x <- gfList$loc_x
        ##mf$loc_y <- gfList$loc_y
        mf$loc_x <- data[,gfList$locxn]
        mf$loc_y <- data[,gfList$locyn]
        ret$model = gfList$model
        ret$fem <- gfList$fem
        ret$mesh <- gfList$mesh
        ret$range_fraction = 0.1

        if(!is.null(gfList$year)){
            mf$year <- data[,gfList$year]
        }
        
        if(gfList$model != 3){
            ##Setup the rangeKey, if it doesn't exist we just do it once per age
            if(!is.null(gfList$rangeKey)){
                ret$rangeKey = gfList$rangeKey
                if(length(gfList$rangeKey) != length(minAge:(maxAge-1))){
                    stop("Incorrect rangeKey setup, must be equal to number of ages-1!")
                }

            }else{
                ret$rangeKey = minAge:(maxAge-1)
            }

            ##Ditto for the sdUkey
            if(!is.null(gfList$sdUkey)){
                ret$sdUkey = gfList$sdUkey
                if(length(gfList$sdUkey) != length(minAge:(maxAge-1))){
                    stop("Incorrect sdUkey setup,length must be equal to number of ages-1!")
                }


            }else{
                ret$sdUkey = minAge:(maxAge-1)
            }

            ret$ranXind = paste(ret$rangeKey,ret$sdUkey,sep="-")
            ret$ranXind = factor(ret$ranXind)
            ret$Qlinks = as.matrix(
                unique(cbind(ret$rangeKey,ret$sdUkey))-1)
        }else{
            if(!is.null(gfList$ageKey)){
                ret$ageKey = gfList$ageKey
                if(length(gfList$ageKey) != length(minAge:(maxAge-1))){
                    stop("Incorrect ageKey setup, must be equal to the number of ages-1!")
                }
            }
            else{
                ret$ageKey = minAge:(maxAge-1)
            }
        } 

    }

    size = rep(1,nrow(mf))
   
    
    u <- rms::cr.setup(age)
    newMf <- data.frame(y = u$y,aGroup=u$cohort)
    newMf[,3:(ncol(mf)+2)] <- mf[u$subs,]
    ##NOPE
    ##Yes?

    ##Calcualte Null deviance for use in trying to avoid cross validation

    getNullDev <- function(weights,y){
        sumwts <- tapply(weights,y,sum)
        dev <- -2*sum(sumwts*logb(sumwts/sum(sumwts)))
        dev
    }
    nullDev <- getNullDev(rep(1,length(u$y)),u$y)

    
    if(pred==FALSE){
      newMf = dplyr::group_by_at(newMf,2:ncol(newMf))
      newMf = dplyr::summarize(newMf,y=sum(y),trials=n())
      size = newMf$trials
      newMf$trials = NULL
    }else{
        size <- size[u$subs]
    }
    
    ##Create a subset factor of the cohorts
    if(!is.null(atilde)){
        subcohort = u$cohort
        levels(subcohort) = paste0("tildeGroup",atilde)
        newMf$tildeGroup = subcohort
    }
    
        

    ##NEEDS to be after we rep observations
    if(!is.null(gfList)){
        ret$A <- INLA::inla.spde.make.A(ret$mesh,as.matrix(cbind(newMf$loc_x,newMf$loc_y)))
        if(sum(ret$A > 0) == 0){
            warning("Projection matrix A is empty! Check spatial coordinates?")
        }

        ret$cohort <- newMf$aGroup
    }


    ##This is to adjust for equal/unequal slope assumptions
    if(equalSlopes == TRUE){
        newForm <- paste("y ~ aGroup +",linearCovar)
    }else if(equalSlopes == FALSE){
        newForm <- paste("y ~ aGroup*(",linearCovar,")")
    }else if(plyr::is.formula(equalSlopes)){
        newTF <- terms.formula(equalSlopes)
        newTerms <- attr(newTF,"term.labels")
        notIn <- terms[!(newTerms %in% terms)]
        if(is.null(gfList)){
            newForm <- paste("y ~ aGroup*(",paste(terms,collapse="+"),") + ",newTerms)
        }else{
            newForm <- paste("y ~ aGroup*(",paste(terms[-gfInd],collapse="+"),") + ",newTerms)
        }
    }else{
        stop("Something went wrong!")
    }


    if(inter == 0){
        newForm <- paste(newForm,"-1")
    }
    newForm <- as.formula(newForm,env=parent.frame())

    if(scale==TRUE){
        if(pred == FALSE){
            stdobj = standardize::standardize(newForm,newMf,family="binomial")
            newMf = stdobj$data
            ret$stdobj = stdobj
        }else{
            newMf = predict(stdobj,newMf)
            newMf$y = u$y
        }
    } 
    
    ret$crFormula <- newForm
    ret$y <- newMf$y
    ret$X <- model.matrix(newForm,data=newMf)
       
    if(is.null(penalty.matrix)){
        penalty.matrix = pen.matrix(newForm,newMf,gfList$rangeKey,gfList$sdUkey)
    }
    
    
    if(length(lambda) == 1){
        if(penrand == TRUE){
            lambda = rep(lambda,ncol(penalty.matrix))
        }else{
            rkey = unique(gfList$rangeKey[!is.na(gfList$rangeKey)])
            skey = unique(gfList$sdUkey[!is.na(gfList$sdUkey)])
            tic = length(rkey) + length(skey)
            lambda = c(rep(lambda,(ncol(ret$X)-1)),rep(0,tic))
        }
    }else{
        if(length(lambda) != ncol(penalty.matrix)){
            stop("Length of lambda wrong")
        }
    }

    
    
    ##I realllly could have planned everything better
    tmbL = list()
    tmbL$data = list()
    tmbL$parm = list()
    tmbL$other = list()

    tmbL$other$repped <- u
    tmbL$other$response <- response
    tmbL$other$age <- age
    tmbL$other$newMf <- newMf
    tmbL$other$nullDev <- nullDev
    if(scale == TRUE){
        tmbL$other$stdobj = stdobj
    }
    
    if(!is.null(atilde)){
        tmbL$other$atilde <- subcohort
    }
    tmbL$data$y = ret$y
    tmbL$data$X = ret$X
    tmbL$data$ages = length(unique(age))
    tmbL$data$size = size
    tmbL$data$lambda = lambda
    tmbL$data$P = penalty.matrix
    tmbL$data$model = ret$model

    tmbL$parm$beta = if(!is.null(starting.param)){starting.param$beta} else{rep(0,ncol(ret$X))}
    if(ret$model == 1 || ret$model == 2){
        tmbL$data$A = ret$A
        tmbL$data$cohort = ret$cohort
        tmbL$data$ranXind = ret$ranXind
        tmbL$data$Qlinks = ret$Qlinks
        tmbL$data$fem = ret$fem
        tmbL$parm$ranX = matrix(0,nrow=ncol(ret$A),ncol=length(ret$Qlinks[,1]))
        if(!is.null(starting.param$ranX)){
            tmbL$parm$ranX = matrix(starting.param$ranX,ncol=length(ret$Qlinks[,1]))
        }
        
        if(ret$model == 1){
            tmbL$data$range_fraction = ret$range_fraction
            tmbL$parm$log_ranges = if(!is.null(starting.param)){starting.param$log_ranges} else {rep(0,length(unique(ret$Qlinks[,1])))}
            tmbL$parm$log_sigma_us = if(!is.null(starting.param)){starting.param$log_sigma_us} else{rep(0,length(unique(ret$Qlinks[,2])))}
        }
        if(ret$model == 2){
            tmbL$parm$log_kappas = if(!is.null(starting.param)){starting.param$log_kappas} else{rep(0,length(unique(ret$Qlinks[,1])))}
            tmbL$parm$log_taus = if(!is.null(starting.param)){starting.param$log_taus}else{rep(0,length(unique(ret$Qlinks[,2])))}
        }
    }
        
    if(ret$model == 3){
        tmbL$data$A = ret$A
        tmbL$data$cohort = ret$cohort
        if(!is.null(gfList$year)){
            tmbL$parm$rhoT = if(!is.null(starting.param)){starting.param$rhoT}else{c(0.1,0.1)}
            Y = length(unique(newMf$year))
            tmbL$data$year = newMf$year
        }else{
            tmbL$data$cohort = ret$cohort
            Y = 1
            tmbL$parm$rhoT = if(!is.null(starting.param)){starting.param$rhoT}else{c(0.1,0)}
            tmbL$map = list(rhoT=as.factor(c(1,NA)))
            tmbL$data$year = as.factor(rep(1,length(newMf$y)))
        }
        aa = length(minAge:maxAge)-1
        aay = aa*Y
        tmbL$data$fem = ret$fem
        tmbL$data$range_fraction = ret$range_fraction
        if(!is.null(gfList$ageKey)){
            ak = gfList$ageKey
            if(nrow(ak) != Y){
                stop("ageKey rows must be equal to number of years")
            }
            if(ncol(ak) != aa){
                stop("ageKey columns must be equal to number of ages-1")
            }

            tmbL$data$ageTimeKey = ak
            aay = length(unique(ak))
        }else{
            tmbL$data$ageTimeKey = matrix(0:(aay-1),nrow=Y,ncol=aa,byrow=TRUE)
        }
        tmbL$parm$log_range = if(!is.null(starting.param)){starting.param$log_range}else{0}
        tmbL$parm$log_sigma_u = if(!is.null(starting.param)){starting.param$log_sigma_u}else{0}
        tmbL$parm$ranX = tmbL$parm$ranX = matrix(0,nrow=ncol(ret$A),ncol=aay)
        if(!is.null(starting.param$ranX)){
            tmbL$parm$ranX = matrix(starting.param$ranX,ncol=aay)
        }
        
        
        


    }

        if(ret$model == 4){
        tmbL$data$A = ret$A
        tmbL$data$cohort = ret$cohort
        if(!is.null(gfList$year)){
            tmbL$parm$rhoT = if(!is.null(starting.param)){starting.param$rhoT}else{c(0.1,0.1)}
            Y = length(unique(newMf$year))
            tmbL$data$year = newMf$year
        }else{
            tmbL$data$cohort = ret$cohort
            Y = 1
            tmbL$parm$rhoT = if(!is.null(starting.param)){starting.param$rhoT}else{c(0.1,0)}
            tmbL$map = list(rhoT=as.factor(c(1,NA)))
            tmbL$data$year = as.factor(rep(1,length(newMf$y)))
        }
        aa = length(minAge:maxAge)-1
        aay = aa*Y
        tmbL$data$fem = ret$fem
        tmbL$data$range_fraction = ret$range_fraction
        if(!is.null(gfList$ageKey)){
            ak = gfList$ageKey
            if(nrow(ak) != Y){
                stop("ageKey rows must be equal to number of years")
            }
            if(ncol(ak) != aa){
                stop("ageKey columns must be equal to number of ages-1")
            }

            tmbL$data$ageTimeKey = ak
            aay = length(unique(ak))
        }else{
            tmbL$data$ageTimeKey = matrix(0:(aay-1),nrow=Y,ncol=aa,byrow=TRUE)
        }
        tmbL$parm$log_kappa = if(!is.null(starting.param)){starting.param$log_kappa}else{0}
        tmbL$parm$log_tau = if(!is.null(starting.param)){starting.param$log_tau}else{0}
        tmbL$parm$ranX = tmbL$parm$ranX = matrix(0,nrow=ncol(ret$A),ncol=aay)
        if(!is.null(starting.param$ranX)){
            tmbL$parm$ranX = matrix(starting.param$ranX,ncol=aay)
        }
        
        
        


    }
        
                        
            
    
        
    attr(tmbL,"ageV") <- unique(age)
    tmbL
}

#' Function to handle the gaussian field construction
#'
#' This function specifies the gaussian field used within a spatialALK formula
#'
#' @param loc_x The x coordinates of the spatial data
#' @param loc_y The y coordinates of the spatial data
#' @param model specify whether the model is basic spde or barrier
#' @param mesh The mesh used to build the precision matrix
#' @param barrier.triangles the set of triangles in the barrier of the mesh, only for barrier model
#' @param rangeKey optionally specify which age groups to share ranges in the precision matrices, default is all ages unique
#' @param sdUkey optionally specify which age groups to share sigma_u's in the precision matrices, default is all ages unique
#' @param ageKey for spatAR, specify age groups share same settings
#' @param year factor variable indicating the year to use
#' 
#' 
#' @export
#'
gf <- function(loc_x,loc_y,model="barrier",mesh=NULL,barrier.triangles=NULL,rangeKey=NULL,sdUkey=NULL,
               ageKey=NULL,year=NULL){

    
    
    gf <- list()
    ##gf$loc_x = loc_x
    ##gf$loc_y = loc_y
    gf$locxn = deparse(substitute(loc_x))
    gf$locyn = deparse(substitute(loc_y))
    if(model == "barrier"){
        gf$model = 1L
        if(is.null(barrier.triangles)){
            stop("Set of barrier triangles needed to use barrier model!")
        }
        ##Yes, three colons are probably bad, but AFAIK I can't get this anywhere else!
        gf$fem <- INLA:::inla.barrier.fem(mesh,barrier.triangles)
    }else if(model == "spde"){
        gf$model = 2L
        ##R-INLA in TMB only supports alpha = 2!
        spde <- INLA::inla.spde2.matern(mesh=mesh,alpha=2)
        gf$fem <- spde$param.inla[c("M0","M1","M2")]
    }else if(model == "spatAR"){
        gf$model = 3L
        if(is.null(barrier.triangles)){
            stop("Set of barrier triangles needed to use barrier model!")
        }
        
        
        gf$fem <- INLA:::inla.barrier.fem(mesh,barrier.triangles)
        if(!missing(year)){
            gf$year = deparse(substitute(year))
        }
    }else if(model == "spatARNoB"){
        gf$model = 4L
        spde <- INLA::inla.spde2.matern(mesh=mesh,alpha=2)
        gf$fem <- spde$param.inla[c("M0","M1","M2")]
    }else{
        stop("Only 'barrier', 'spatAR', 'spatARNoB' and 'spde' supported.")
    }

    if(!is.null(rangeKey)){
        gf$rangeKey = rangeKey
    }
    if(!is.null(sdUkey)){
        gf$sdUkey = sdUkey
    }
    if(!is.null(ageKey)){
        gf$ageKey = ageKey
    }
    


    gf$mesh <- mesh

    gf
}

contr.ave <- function(n,contrasts=TRUE){
    if(length(n) <= 1L){
        if(is.numeric(n) && length(n)  == 1L && n > 1L)
            levels <- seq_len(n)
        else stop("not enough degrees of freedom to define contrasts")
    }
    else levels <- n
    levels <- as.character(levels)
    if(contrasts){
        n <- length(levels)
        contr <- diag(n)-1/(n+1)
        colnames(contr) <- NULL
        contr
    }
    else diag(levels)
}



pen.matrix <- function(object,data,rangeKey=NULL,sdUkey=NULL){
    X <- model.matrix(object,data)
    tt <- terms(object)
    ttl <- attr(tt,"term.labels")
    assm <- attr(X,"assign")
    ##Throw away the intercept terms
    keep <- which(assm != 0)
    X <- X[,keep] 
    assm = assm[keep]
    runs = rle(assm)
    pen = matrix(0,nrow=ncol(X),ncol=ncol(X))
    ind = 1
    for(i in 1:length(runs$values)){
        lab = ttl[runs$values[i]]
        sects = strsplit(lab,":")[[1]]
        dd = data[,sects,drop=FALSE]
        idu = ind+runs$lengths[[i]]-1
        if(all(sapply(dd,is.factor))){
            pen[ind:idu,ind:idu] = contr.ave(runs$lengths[i])
            ind = ind + runs$lengths[i]
        }else{
            vars = apply(X[,ind:idu,drop=FALSE],2,var)
            pen[ind:idu,ind:idu] = diag(vars,nrow=length(vars),ncol=length(vars))
            ind = ind + runs$lengths[i]
        }
    }

    if(!is.null(rangeKey) & !is.null(sdUkey))
    {
        rkey = unique(rangeKey[!is.na(rangeKey)])
        skey = unique(sdUkey[!is.na(sdUkey)])
        tot = ncol(X)+length(rkey)+length(skey)
        dig = diag(tot)
        dig[1:ncol(X),1:ncol(X)] = pen
        pen = dig
    }
    
    
    pen
        
}

