#'Generate an emperical ALK from the observed proportions in each age group
#'
#' Note the shitty fix for DFO data. 
#'
#' @param age a vector containing the ages
#' @param length a vector containing the lengths
#' @param minAge pool ages below this together
#' @param maxAge pool ages above this together
#' @param lengthBins the desired vector of length bins for the ALK
#' @param proportion if true gets the ALK in proportion form
#' @export
getALKraw <- function(age,length,minAge,maxAge,lengthBins,proportion=FALSE)
{
    ##Shitty fix for lengthbin problem for rv survey data -1 is cause 3cm blah blah I hate length bins
    lengthPlus <- c(lengthBins-1,Inf)
    sizeGroup <- cut(length,breaks=lengthPlus,right=FALSE)
    sizeGroup <- as.factor(sizeGroup)
    modAge <- age
    modAge[modAge < minAge] = minAge
    modAge[modAge > maxAge] = maxAge
    modAge <- as.factor(modAge)
    ages <- length(levels(modAge))
    levels(modAge)[ages] <- paste(levels(modAge)[ages],"+",sep="")
    ALK = xtabs(~sizeGroup + modAge,data=cbind(age,length))
    if(proportion==TRUE){
        ALK = ALK/rowSums(ALK,na.rm=TRUE)
    }
    rownames(ALK) <- lengthBins
    class(ALK) <- "rawALK"
    attr(ALK,"prop") <- proportion
    ALK
}


#'Plot a spatial field
#'
#' This function returns a ggplot2 plot of a spatial field
#'
#' @param field field to plot
#' @param mesh mesh to project on to
#' @param xlim x limits of plotting area
#' @param ylim y limits of plotting area
#' @param poly.barrier optional barrier.polygon to
#' @param dims dimension of the projection
#' 
#' @export
#'
plotSpField <- function(field,mesh,xlim,ylim,poly.barrier=NULL,dims=c(300,300)){
    stopifnot(nrow(field) == mesh$n)

    loc = data.frame(x=mesh$loc[,1],y=mesh$loc[,2])
    loc = sf::st_as_sf(loc,coords=c("x","y"))
    b = sf::st_bbox(loc)    
    
    if(missing(xlim)){
        xlim = c(b$xmin,b$xmax)
    }
    if(missing(ylim)){
        ylim = c(b$ymin,b$ymax)
    }

    proj = INLA::inla.mesh.projector(mesh,xlim=xlim,ylim=ylim,dims=dims)
    field.proj = INLA::inla.mesh.project(proj,field)

    f = expand.grid(proj$x,proj$y)
    f$z = as.vector(field.proj)
    names(f) = c("x","y","z")

    rast = ggplot2::geom_raster(ggplot2::aes(x,y,fill=z),data=f)

    if(is.null(poly.barrier)){
        plot = ggplot2::ggplot() + rast  + ggplot2::scale_fill_viridis_c(option='plasma')
    }else{
        sfbarrier = sf::st_as_sf(poly.barrier)
        sfbarrier = st_crop(sfbarrier,xmin=xlim[1],xmax=xlim[2],ymin=ylim[1],ymax=ylim[2])
        
        plot = ggplot2::ggplot() + ggplot2::geom_sf(data=sfbarrier)
        plot$layers = c(rast,plot$layers)
        plot = plot  + ggplot2::scale_fill_viridis_c(option='plasma')
    }
    

    plot
}

#'The correlation of a precision Q at a single location
#'
#' The spatial correlation of a single point
#' 
#' @param Q the precision matrix Q to plot
#' @param location the location to plot
#' @param mesh the mesh used
#'
#' @export
#'
PointCorr <- function(Q,location,mesh){
    Sigma = INLA::inla.qinv(Q)
    var = Matrix::diag(Sigma)
    sd = sqrt(var)

    A = INLA::inla.spde.make.A(mesh=mesh,loc=matrix(c(location[1],location[2]),1,2))
    id.node = which.max(A[1,])

    Inode = rep(0,nrow(Q))
    Inode[id.node] = 1
    covar = solve(Q,Inode)

    corr = drop(matrix(covar))/(sd*sd[id.node])
    corr
}



#'Plot a predicted ALK
#'
#' @param x the predALK object to plot
#' @param ... the options to pass onto matplot
#' @param rawALK the rawALK proportions
#' object to plot with the predicted one
#' @export
plot.predALK <- function(x,lwd=2,type="l",ylab="Proportion",
                        xlab="Length (CM)",col=NULL,lty=seq_len(length(colors)),...,rawALK=NULL,plotFreq=FALSE,autoLen=TRUE){
    if(is.null(col)){
        colors <- seq_len(ncol(x))
    }else{
        colors <- col
    }

    if(autoLen == TRUE){
            layout(matrix(c(1,2),nrow=1),width=c(4,1))
            par(mar=c(5,4,4,0))
    }

    matplot(rownames(x),x,lwd=lwd,type=type,ylab=ylab,xlab=xlab,
            lty=lty,col=colors,...)

        if(!is.null(rawALK)){
        prop <- attr(rawALK,"prop")
        realRaw <- rawALK
        if(prop == FALSE){
            rawALK <- rawALK/rowSums(rawALK,na.rm=TRUE)
        }
        ##Don't want to plot zeroes
        rawALK[rawALK == 0] <- NA
        labs <- data.frame(label=rep(colnames(rawALK),nrow(rawALK)),
                           y=as.vector(t(rawALK)),rawN=as.vector(t(realRaw)))
        labs$x <- sort(as.numeric(rep(rownames(rawALK),ncol(rawALK))))
        labs$cCols <- sub("*\\+","",labs$label)
        labs$cCols <- as.numeric(labs$cCols)
        labs$trueCols <- colors[labs$cCols]
        if(plotFreq == FALSE){
            text(x=labs$x,y=labs$y,labels=labs$label,col=labs$trueCols)
        }else if(plotFreq == TRUE & prop == FALSE){
            text(x=labs$x,y=labs$y,labels=labs$label,col=labs$trueCols)
            text(x=labs$x,y=labs$y,labels=labs$rawN,col="orange",pos=1,cex=0.65,adj=0.5)
        }else{
            stop("Can't plot frequencies of proportions")
        }


        }
    if(autoLen == TRUE){

        par(mar=c(5,0,4,2))
        plot(c(0,1),type="n",axes=F,xlab="",ylab="")
        if(plotFreq == FALSE){
            legend("center",legend=colnames(x),col=colors,lty=lty,text.col=colors,lwd=lwd)
        }else{
            legend("center",legend=c(colnames(x),"Num. Aged"),col=colors,"orange",lty=lty,text.col=c(colors,"orange"),lwd=lwd)
        }
    }
}
