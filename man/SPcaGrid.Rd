\name{SPcaGrid}
\alias{SPcaGrid}
\alias{SPcaGrid.formula}
\alias{SPcaGrid.default}
\title{Sparse Robust Principal Components based on Projection Pursuit (PP): GRID search Algorithm }
\description{
    Computes an approximation of the PP-estimators for sparse and robust PCA using the grid search algorithm in the plane.
}
\usage{
    SPcaGrid(x, ...)
    \method{SPcaGrid}{default}(x, k = 0, kmax = ncol(x), method = c ("mad", "sd", "qn", "Qn"), 
    lambda = 1, scale=FALSE, na.action = na.fail, trace=FALSE, \dots)
    \method{SPcaGrid}{formula}(formula, data = NULL, subset, na.action, \dots)
}
\arguments{
  \item{formula}{a formula with no response variable, referring only to
    numeric variables.}
  \item{data}{an optional data frame (or similar: see
    \code{\link{model.frame}}) containing the variables in the
    formula \code{formula}.}
  \item{subset}{an optional vector used to select rows (observations) of the
    data matrix \code{x}.}
  \item{na.action}{a function which indicates what should happen
    when the data contain \code{NA}s.  The default is set by
    the \code{na.action} setting of \code{\link{options}}, and is
    \code{\link{na.fail}} if that is unset. The default is \code{\link{na.omit}}.}
  \item{\dots}{arguments passed to or from other methods.}
  \item{x}{a numeric matrix (or data frame) which provides
    the data for the principal components analysis.}
  \item{k}{number of principal components to compute. If \code{k} is missing, 
    or \code{k = 0}, the algorithm itself will determine the number of 
    components by finding such \code{k} that \eqn{l_k/l_1 >= 10.E-3} and 
    \eqn{\Sigma_{j=1}^k l_j/\Sigma_{j=1}^r l_j >= 0.8}. 
    It is preferable to investigate the scree plot in order to choose the number 
    of components and then run again. Default is \code{k=0}. }
  \item{kmax}{maximal number of principal components to compute.
    Default is \code{kmax=10}. If \code{k} is provided, \code{kmax} 
    does not need to be specified, unless \code{k} is larger than 10.}
  \item{method}{ the scale estimator used to detect the direction with the 
  largest variance. Possible values are \code{"sd"}, \code{"mad"} and 
  \code{"Qn"}. \code{"mad"} is the 
  default value.}
  \item{lambda}{ the sparseness constraint's strength(\code{sPCAgrid} only).
    A single value for all components, or a vector of length \code{k} with 
    different values for each component can be specified.
    See \code{\link[pcaPP]{opt.TPO}} for the choice of this argument. }
  \item{scale}{a value indicating whether and how the variables should be 
  scaled. If \code{scale = FALSE} (default) or \code{scale = NULL} no scaling is 
  performed (a vector of 1s is returned in the \code{scale} slot). 
  If \code{scale = TRUE} the data are scaled to have unit variance. Alternatively it can 
  be a function like \code{sd} or \code{mad} or a vector of length equal 
  the number of columns of \code{x}. The value is passed to the underlying function
  and the result returned is stored in the \code{scale} slot. 
  Default is \code{scale = FALSE}}
  \item{trace}{whether to print intermediate results. Default is \code{trace = FALSE}}
}



\details{
\code{SPcaGrid}, serving as a constructor for objects of class \code{\link{SPcaGrid-class}} 
is a generic function with "formula" and "default" methods. For details see 
\code{\link[pcaPP]{sPCAgrid}} and the relevant references.
}
\value{
  An S4 object of class \code{\link{SPcaGrid-class}} which is a subclass of \code{\link[rrcov]{PcaGrid-class}} which in turn is a subclass of the 
  virtual class \code{\link[rrcov]{PcaRobust-class}}. 
}


\references{
    C. Croux, P. Filzmoser, M. Oliveira, (2007).
    Algorithms for Projection-Pursuit Robust Principal Component Analysis,
    \emph{Chemometrics and Intelligent Laboratory Systems}, Vol. 87, pp. 218-225.
    
    C. Croux, P. Filzmoser, H. Fritz (2013).
    Robust Sparse Principal Component Analysis,
    \emph{Technometrics} \bold{55}(2), pp. 202--2014,
    \doi{10.1080/00401706.2012.727746}.   
  
    V. Todorov, P. Filzmoser (2013). Comparing classical and robust sparse PCA.
    In R Kruse, M Berthold, C Moewes, M Gil, P Grzegorzewski, O Hryniewicz (eds.), 
    \emph{Synergies of Soft Computing and Statistics for Intelligent Data Analysis}, 
    volume 190 of \emph{Advances in Intelligent Systems and Computing}, pp. 283--291. 
    Springer, Berlin; New York. ISBN 978-3-642-33041-4,
    \doi{10.1007/978-3-642-33042-1_31}.
}

\author{ Valentin Todorov \email{valentin.todorov@chello.at} 
}

\examples{

data(bus)
bus <- as.matrix(bus)

## calculate MADN for each variable
xmad <- apply(bus, 2, mad)
cat("\nMin, Max of MADN: ", min(xmad), max(xmad), "\n")

## calculate MADN for each variable
xqn <- apply(bus, 2, Qn)
cat("\nMin, Max of Qn: ", min(xqn), max(xqn), "\n")


## MADN vary between 0 (for variable 9) and 34. Therefore exclude
##  variable 9 and divide the remaining variables by their MADNs.
bus1 <- bus[, -c(9)]
p <- ncol(bus1)

madbus <- apply(bus1, 2, mad)
bus2 <- sweep(bus1, 2, madbus, "/", check.margin = FALSE)

xsd <- apply(bus1, 2, sd)
bus.sd <- sweep(bus1, 2, xsd, "/", check.margin = FALSE)

xqn <- apply(bus1, 2, Qn)
bus.qn <- sweep(bus1, 2, xqn, "/", check.margin = FALSE)

\dontrun{
spc <- SPcaGrid(bus2, lambda=0, method="sd", k=p, kmax=p)
rspc <- SPcaGrid(bus2, lambda=0, method="Qn", k=p, kmax=p)
summary(spc)
summary(rspc)
screeplot(spc, type="line", main="Classical PCA", sub="PC", cex.main=2)
screeplot(rspc, type="line", main="Robust PCA", sub="PC", cex.main=2)

##  find lambda

K <- 4
lambda.sd <- 1.64
    to.sd <- .tradeoff(bus2, k=K, lambda.max=2.5, lambda.n=100, method="sd")
    plot(to.sd, type="b", xlab="lambda", ylab="Explained Variance (percent)")
    abline(v=lambda.sd, lty="dotted")
 
spc.sd.p <- SPcaGrid(bus2, lambda=lambda.sd, method="sd", k=p)
.CPEV(spc.sd.p, k=K)
spc.sd <- SPcaGrid(bus2, lambda=lambda.sd, method="sd", k=K)
getLoadings(spc.sd)[,1:K]
plot(spc.sd)

lambda.qn <- 2.06
    to.qn <- .tradeoff(bus2, k=K, lambda.max=2.5, lambda.n=100, method="Qn")
    plot(to.qn, type="b", xlab="lambda", ylab="Explained Variance (percent)")
    abline(v=lambda.qn, lty="dotted")

spc.qn.p <- SPcaGrid(bus2, lambda=lambda.qn, method="Qn", k=p)
.CPEV(spc.qn.p, k=K)
spc.qn <- SPcaGrid(bus2, lambda=lambda.qn, method="Qn", k=K)
getLoadings(spc.qn)[,1:K]
plot(spc.qn)
}

## DD-plots
##
## Not run:
\dontrun{
usr <- par(mfrow=c(2,2))
plot(SPcaGrid(bus2, lambda=0, method="sd", k=4), id.n.sd=0, main="Standard PCA")
plot(SPcaGrid(bus2, lambda=0, method="Qn", k=4), id.n.sd=0, ylim=c(0,20))

plot(SPcaGrid(bus2, lambda=1.64, method="sd", k=4), id.n.sd=0, main="Stdandard sparse PCA")
plot(SPcaGrid(bus2, lambda=3.07, method="Qn", k=4), id.n.sd=0, main="Robust sparse PCA")

par(usr)
## End (Not run)
}
}

\keyword{robust}
\keyword{multivariate}
