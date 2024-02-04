SosDiscRobust <- function (x, ...) UseMethod("SosDiscRobust")

SosDiscRobust.formula <- function(formula, data=NULL, ..., subset, na.action)
{
    m <- match.call(expand.dots = FALSE)
    m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    grouping <- model.response(m)
    x <- model.matrix(Terms, m)
    xint <- match("(Intercept)", colnames(x), nomatch=0)
    if(xint > 0)
    x <- x[, -xint, drop=FALSE] # always without intercept
    #if(length(colnames(x))==0){ colnames(x) <- paste0("V", 1:ncol(x))}
    res <- SosDiscRobust.default(x, grouping, ...)

    ##    res$terms <- Terms

    ## fix up call to refer to the generic, but leave arg name as 'formula'
    cl <- match.call()
    cl[[1]] <- as.name("rsos")
    res@call <- cl

    ##    res$contrasts <- attr(x, "contrasts")
    ##    res$xlevels <- .getXlevels(Terms, m)
    ##    res$na.action <- attr(m, "na.action")

    res
}

SosDiscRobust.default <- function(x,
                         grouping,
                         prior = proportions,
                         lambda,                             # Sparsity parameter for L1 norm penalty
                         Q=length(unique(grouping))-1,       # number of optimal scoring coefficient vectors, Q has to be smaller than the number of groups
                         alpha=0.5,                          # robustness parameter in sparseLTS (for initial estimation)
                         maxit=100,                          # number of iterations for the estimation of optimal scoring coefficients and case weights
                         tol=1.0e-4,
                         trace=FALSE, ...)
{

###### functions
    ## this is stright from nnet:::formula
    class.ind <- function(cl)
    {
        Ik <- diag(length(levels(cl)))
        x <- Ik[as.numeric(cl),]
        dimnames(x) <- list(names(cl), levels(cl))

        x
    }

    orth.Q <- function(dp,Qj,theta)
    {
        DtQ <- t(dp)%*%Qj # t(C)
        theta <- theta - Qj%*%as.matrix(solve(t(Qj)%*%DtQ))%*%t(DtQ)%*%theta
        thetan <- sqrt(as.numeric(t(theta)%*%dp%*%theta))
        theta*1/thetan
    }

    rtheta <- function(K,dp)
    {
        jj <- rnorm(K)
        jj/sqrt(sum(jj^2)*dp)
    }

######### end functions


    if(is.null(dim(x)))
        stop("x is not a matrix")

    xcall <- match.call()
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)

    if(length(colnames(x))==0)
    {
      stop("data matrix x has no column names.")
      #  colnames(x) <- paste0("V", 1:ncol(x))
    }

    grouping <- as.factor(grouping)
    xmed <- apply(x, 2, median)
    xscale <- apply(x, 2, mad)
    xs <- scale(x[,xscale!=0], center=xmed[xscale!=0], scale=xscale[xscale!=0])
    predNames <- colnames(x)

    grpx <- .getGrouping(grouping, n)
    if(!missing(prior))
    {
        if(any(prior < 0) || round(sum(prior), 5) != 1)
            stop("invalid prior")
        if(length(prior) != grpx$ng)
            stop("prior is of incorrect length")
        prior <- prior[grpx$counts > 0]
    } else
        proportions <- grpx$proportions

    K <- grpx$ng
    Y <- class.ind(grpx$grouping)

    if(missing(lambda))
        lambda <- 0

    stopifnot(lambda >= 0)

    if(missing(Q))
        Q <- length(unique(grouping))-1
    if(Q > K-1)
        stop("at most Q=K-1 variates allowed")

    dpi <- matrix(0, nrow=grpx$ng, ncol=grpx$ng)
    diag(dpi) <- grpx$proportions
    ydp <- scale(Y,FALSE,grpx$proportions) # scale y with prior probabilities, D^{-1}

    b <- matrix(0,sum(xscale!=0),Q)
    theta <- matrix(0,K,Q)
    Qj <- matrix(1,K,1)
    rss <- rep(0,Q)
    wslts <- matrix(ncol=Q, nrow=n)
    wrew <- matrix(ncol=Q, nrow=n)

    for(j in 1:Q)
    {
        #### Robust Initial coefficient estimators with LTS and L1:
        thetaj <- rtheta(K,grpx$proportions)   # random theta, fulfills 1st side constraint
        thetaj <- orth.Q(dpi,Qj,thetaj)        # fulfills 2nd side constrinat, orth to previous thetas


        for(srep in 1:2)    # repeat initial coefficient estimation two times for better results
        {
            Yc <- Y %*% thetaj		
            sLTSmod <- sparseLTS(xs, Yc, lambda=lambda, mode="lambda", intercept=FALSE, alpha=alpha, normalize=FALSE) # no internal cv
            betas <- sLTSmod$coefficients # reweighted sLTS

            if(all(betas==0))
            {
                b[,j] <- matrix(0,p,1)
                break
            }
            yhatj=as.matrix(xs)%*%betas

            # suppressWarnings( # issues warning because solution is not unique
            #     qr1 <- rq(yhatj ~ -1+ydp, tau = 0.5)
            # )

            # thetaj <- orth.Q(dpi,Qj,qr1$coefficients)

            qr1 <- .l1reg(A=ydp, b=yhatj, p=1)
            thetaj <- orth.Q(dpi,Qj,qr1$x)

            # rmRSS2 <- mean((as.matrix(xs)%*%betas-Y%*%thetaj)^2)+lambda*sum(abs(betas))
            rmRSS <- sum((as.matrix(xs) %*% betas - Y %*% thetaj)^2 * sLTSmod$w) / sum(sLTSmod$w) + lambda * sum(abs(betas))
            if (trace)
            {
                cat('\nstarting rep: ', srep, ' weighted lasso cost: ',rmRSS,' |b|_1: ', sum(abs(betas)),' sum of weights: ', sum(sLTSmod$w),'\n')
            }
        }

        #### Reweighting Step
        if(!all(betas==0))
        {
            if (trace)
                cat('\n re-weighting:\n')

            ## residuals
            ri <- Y %*% thetaj - as.matrix(xs) %*% betas
            wi <- .hampel(ri, q1=qnorm(0.95), q2=qnorm(0.975), q3=qnorm(0.999), as.factor(Yc))
            rmRSS <- 1e6
            RSSold <- Inf
            ite <- 0
            while (abs(RSSold-rmRSS)/rmRSS > tol & ite < maxit)
            {
                RSSold <- rmRSS
                ite <- ite + 1

                ## 1. Estimate beta:
                Yc <- Y%*%thetaj
                betas <- .solvebeta(sqrt(wi)*xs, sqrt(wi)*Yc, paras=c(0, sum(wi)*lambda), sparse="penalty")

                if(all(betas==0))
                {
                    b[,j] <- matrix(0,p,1)
                    break
                }

                ## 2. Estimate theta:
                yhatj <- as.matrix(xs)%*%betas
                dpi <- 1/sum(wi)*t(sqrt(wi)*Y)%*%(sqrt(wi)*Y)

                thetaj <- orth.Q(dpi,Qj,solve(dpi)%*%t(sqrt(wi)*Y)%*%(sqrt(wi)*yhatj))

                ## 3. update residuals:
                ri <- Y%*%thetaj - yhatj
                wi <- .hampel(ri,q1=qnorm(0.95),q2=qnorm(0.975),q3=qnorm(0.999),as.factor(Yc))


                # rmRSS2 <- mean((as.matrix(xs)%*%betas-y%*%thetaj)^2)+lambda*sum(abs(betas))
                rmRSS <- sum((as.matrix(xs)%*%betas-Y%*%thetaj)^2*wi)/sum(wi)+lambda*sum(abs(betas))
                if (trace)
                    cat('\nite: ', ite, ' weighted lasso cost: ', rmRSS, ' |b|_1: ', sum(abs(betas)), ' sum of weights: ', sum(wi), '\n')
            }

            if(!all(betas==0))
            {
                rss[j] <- rmRSS
                Qj <- cbind(Qj,thetaj)
                theta[,j] <- thetaj
                b[,j] <- betas
                wslts[,j] <- sLTSmod$wt	
                wrew[,j] <- wi
            }
        }
    }

    if (trace)
    {
        cat('final update, total weighted lasso cost: ', sum(rss), ' |b|_1: ', sum(abs(b)),'\n')
    }

    #### calcualte classification rule with Linda
    if (all(b==0))
    { # all values in b and sl are zero
        warning("no non-zero coefficients - try other regularization parameter")
        b <- matrix(0,p,Q)
        sl <- matrix(0,n,1)
        colnames(sl) <- paste("score", 1:ncol(sl), sep = "")
        notZero <- matrix(TRUE,p,1) # fake NotZero for alignment purpose
        varnames <- colnames(x)

        lobj<-list(
            prior=prior,
            means = NULL,
            x = x,
            covw = matrix(0,1,1)
        )
        #origP <- ncol(x)
        wm2 <- rep(1,n)
        mahadist2 <- rep(NA,n)
    } else
    {
        # remove predictors which are not included (do not have non-zero parameter estimates)
        notZero <- apply(b, 1, function(x) any(x != 0))
        b <- b[notZero,]
        #origP <- ncol(x)
        xs <- xs[, notZero, drop = FALSE]
        varnames <- colnames(xs)

        ### remove directions with only zero elements (this can be caused by a too high weight on L1-penalty)
        if (is.vector(b))
        {
            notZeroC <- (b!=0)
            b <- as.matrix(b[notZeroC])
        } else
        {
            notZeroC <- apply(b,2,function(x) any(x!=0))
            b <- as.matrix(b[,notZeroC])		
        }

        sl <- xs %*% b
        colnames(sl) <- paste("score", 1:ncol(sl), sep = "")
        lobj <- Linda(as.matrix(sl), grpx$grouping, prior=prior, method="mcd")
        if(ncol(sl) == 1)
        {
            mahadist2 <- vector(length=n)
            for(i in 1:grpx$ng)
            {
                yi <- grpx$lev[i]
                mahadist2[grpx$grouping==yi] <- scale(sl[grpx$grouping==yi,],center=lobj@center[rownames(lobj@center)==yi,], scale=sqrt(as.numeric(lobj@cov)))^2
            }

            delta <- 0.975
            wm2 <- as.numeric(mahadist2<=qchisq(delta,df=1))		#!!!!
            #wm2 <- as.numeric(abs(mahadist2)<=qnorm(delta))		
        } else
        {
            mahadist2 <- vector(length=n)
            for(i in 1:grpx$ng)
            {
                yi <- grpx$lev[i]
                mahadist2[grpx$grouping==yi] <- mahalanobis(sl[grpx$grouping==yi,], center=lobj@center[rownames(lobj@center)==yi,], cov=lobj@cov)
            }
            delta <- 0.975
            wm2 <- as.numeric(mahadist2<=qchisq(delta,df=ncol(sl)))		
        }
    }

    ret <- new("SosDiscRobust",
        call = xcall,
        prior=prior,
        counts=grpx$counts,

        beta = b,
        theta = theta,
        lambda = lambda,
        varnames = varnames,

        # varIndex = which(colnames(x)%in%varnames),
        # origP = ncol(x),
        # rss = rss,

        fit = lobj,
        center = xmed,
        scale = xscale,

        mahadist2=mahadist2,
        wlinda = wm2,

        X=x,
        grp=grpx$grouping)

    return (ret)
}

########## internal functions
.hampel <- function(x, q1=qnorm(0.95), q2=qnorm(0.975), q3=qnorm(0.999), y)
{
  # caluculates hampel weights for each group coded in y separately
  # x .......... vector of residuals
  # q1,q2,q3 ... cut off values of hampel's weight function
  # y .......... group conding (factor)

    w <- vector(length=length(x))

    for(i in levels(y))
    {
        wtmp <- vector(length=sum(y==i))
        resi <- x[y==i]
        resis <- robStandardize(resi)
        wtmp[which(abs(resis) <= q1)] <- 1
        wtmp[which(abs(resis) > q1 & abs(resis) <= q2)] <- q1/abs(resis[which(abs(resis) > q1 & abs(resis) <= q2)])
        wtmp[which(abs(resis) > q2 & abs(resis) <= q3)] <- q1 *
            (q3 - abs(resis[which(abs(resis) > q2 & abs(resis) <= q3)]))/
            (q3 - q2) * 1/abs(resis[which(abs(resis) > q2 & abs(resis) <= q3)])

        wtmp[which(abs(resis) > q3)] <- 0
        w[y==i] <- wtmp
    }

    return(w)
}

#################################################################################################################
setMethod("predict", "SosDisc", function(object, newdata, ...)
{
    if(all(object@beta == 0))
    {
        return(new("PredictSosDisc",
            classification <- rep(NA, nrow(newdata)),
            w <- rep(NA, nrow(newdata)),
            mahadist2 <-  matrix(NA, nrow=nrow(newdata), ncol= length(object@fit@prior)) ))
    } else
    {
        if(missing(newdata))
            newdata <- object@X

      if(length(colnames(newdata))==0){
        stop("Column names for newdata are missing")
      }
        if(!is.matrix(newdata))
            newdata <- as.matrix(newdata)
        newdata <- scale(newdata,object@center, object@scale)

        #if(!is.null(object@varnames))
        #{
            newdata <- newdata[, object@varnames, drop = FALSE]
        #} else
        #{
        #    if(ncol(newdata) != object@origP) stop("dimensions of training and testing X different")
        #    newdata <- newdata[, object@varIndex, drop = FALSE]
        #}

        Yfacs <- object@fit@grp
        xnew <- newdata %*% object@beta
        pred <- predict(object@fit,xnew)
        mahadist2 <- matrix(nrow=nrow(newdata), ncol=length(unique(Yfacs)))
        for(i in 1:length(unique(Yfacs)))
        {
            yi <- unique(Yfacs)[i]
            mahadist2[,i] <- mahalanobis(xnew, center=object@fit@center[rownames(object@fit@center)==yi,], cov=object@fit@cov)
        }

        delta <- 0.975
        w <- apply(mahadist2, 1, min) <= qchisq(delta, df=ncol(object@beta))
        return(new("PredictSosDisc",
            classification=pred@classification,
            w=w,
            mahadist2=mahadist2)
        )
    }
})

setMethod("show", "SosDisc", function(object)
{
    digits <- max(3, getOption("digits") - 3)

    if(!is.null(cl <- object@call))
    {
        names(cl)[2] <- ""
        cat("Call:\n")
        dput(cl)
    }

    spar <- paste(format(object@lambda, digits = digits), "L1 bound")

    classInfo <- paste(unique(object@fit@grp), collapse = ", ")

    cat("lambda =", spar,
      "\nclasses =", classInfo,
      "\nPrior Probabilities of Groups: ", object@fit@prior,
      "\n\n")

    # top <- if(!is.null(object@varnames)) object@varnames else paste("Predictor", object@varIndex, sep = "")
    top <- object@varnames
    varOrder <- if(is.matrix(object@beta)) order(apply(abs(object@beta), 1, sum)) else order(abs(object@beta))
    top <- top[varOrder]
    top <- top[1:min(5, length(top))]
    top <- paste(top, collapse = ", ")

    if(nrow(object@beta) > 5)
    {
        cat("Top 5 predictors (out of ",
            length(object@varnames),
            "):\n\t",
            top,
            sep = "")
    } else {
        cat("Predictors:\t",
            top,
            "\n",
            sep = "")
    }
    cat("\n")
    invisible(object)
})

setMethod("summary", "SosDisc", function(object, ...){
  digits <- max(3, getOption("digits") - 3)

  if(!is.null(cl <- object@call))
  {
    names(cl)[2] <- ""
    cat("Call:\n")
    dput(cl)
  }

  spar <- paste(format(object@lambda, digits = digits), "L1 bound")

  classInfo <- paste(unique(object@fit@grp), collapse = ", ")

  cat("lambda =", spar,
      "\nclasses =", classInfo,
      "\nPrior Probabilities of Groups: ", object@fit@prior,
      "\n\n")

  # top <- if(!is.null(object@varnames)) object@varnames else paste("Predictor", object@varIndex, sep = "")
  top <- object@varnames
  varOrder <- if(is.matrix(object@beta)) order(apply(abs(object@beta), 1, sum)) else order(abs(object@beta))
  top <- top[varOrder]
  top <- top[1:min(5, length(top))]
  top <- paste(top, collapse = ", ")

  if(nrow(object@beta) > 5)
  {
    cat("Top 5 predictors (out of ",
        length(object@varnames),
        "):\n\t",
        top,
        sep = "")
  } else {
    cat("Predictors:\t",
        top,
        "\n",
        sep = "")
  }
  cat("\n\nConfusion table\n")
  tab <- table(object@fit@grp,predict(object@fit)@classification)
  print(tab)
  cat("\nMisclassification rate: ", 1-mean(diag(tab)/apply(tab,1,sum)), "\n")
  invisible(object)

})

setMethod(f="plot", signature="SosDisc", definition=function(x, ind=c(1,2),...){
  if(ncol(x@fit@X)==1){
    ind <- 1
    plot(x@fit@X[,ind], col=x@fit@grp, ylab=expression(paste(X,beta)))
  }
  else {
    plot(as.data.frame(x@fit@X[,ind]), col=x@fit@grp)
  }

})
