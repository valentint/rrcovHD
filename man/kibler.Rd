\name{kibler}
\alias{kibler}
\alias{kibler.orig}
\docType{data}
\title{
    1985 Auto Imports Database
}

\description{
 The original data set \code{kibler.orig} consists of three types of entities: 
 (a) the specification of an auto in terms of various characteristics, 
 (b) its assigned insurance risk rating and
 (c) its normalized losses in use as compared to other cars.  
 
The second rating corresponds to the degree to which the auto is more 
    risky than its price indicates. Cars are initially assigned a risk 
    factor symbol associated with its price.   Then, if it is more 
    risky (or less), this symbol is adjusted by moving it up (or down) 
    the scale.  Actuarians call this process "symboling".  A value 
    of +3 indicates that the auto is risky, -3 that it is probably pretty safe.
    
The third factor is the relative average loss payment per insured 
    vehicle year.  This value is normalized for all autos within a particular 
    size classification (two-door small, station wagons, sports/speciality, etc...), 
    and represents the average loss per car per year.
}
\usage{data(kibler)}
\format{
  A data frame with 195 observations on the following 14 variables. The original 
  data set (also available as \code{kibler.orig}) contains 205 cases and 
  26 variables of which 15 continuous, 1 integer and 10 nominal. The 
  non-numeric variables and variables with many missing values were removed. 
  Cases with missing values were removed too. 
  \describe{
    \item{\code{symboling}}{a numeric vector}
    \item{\code{wheel-base}}{a numeric vector}
    \item{\code{length}}{a numeric vector}
    \item{\code{width}}{a numeric vector}
    \item{\code{height}}{a numeric vector}
    \item{\code{curb-weight}}{a numeric vector}
    \item{\code{bore}}{a numeric vector}
    \item{\code{stroke}}{a numeric vector}
    \item{\code{compression-ratio}}{a numeric vector}
    \item{\code{horsepower}}{a numeric vector}
    \item{\code{peak-rpm}}{a numeric vector}
    \item{\code{city-mpg}}{a numeric vector}
    \item{\code{highway-mpg}}{a numeric vector}
    \item{\code{price}}{a numeric vector}
  }
}
\details{

The original data set contains 205 cases and 26 variables of which 15 continuous, 1 integer and 10 nominal. The non-numeric variables and
variables with many missing values were removed. Cases with missing values were removed too. Thus the data set remains with 195 cases and 14 variables.

}
\source{
    www.cs.umb.edu/~rickb/files/UCI/
}
\references{
    Kibler, D., Aha, D.W. and Albert, M. (1989).  Instance-based prediction
       of real-valued attributes.  \emph{Computational Intelligence}, Vo.l 5,
       51-57. 
}
\examples{
data(kibler)
x.sd <- apply(kibler,2,sd)
xsd <- sweep(kibler, 2, x.sd, "/", check.margin = FALSE)
apply(xsd, 2, sd)

x.mad <- apply(kibler, 2, mad)
xmad <- sweep(kibler, 2, x.mad, "/", check.margin = FALSE)
apply(xmad, 2, mad)

x.qn <- apply(kibler, 2, Qn)
xqn <- sweep(kibler, 2, x.qn, "/", check.margin = FALSE)
apply(xqn, 2, Qn)


## Display the scree plot of the classical and robust PCA
screeplot(PcaClassic(xsd))
screeplot(PcaGrid(xqn))

#########################################
##
## DD-plots
##
\dontrun{
usr <- par(mfrow=c(2,2))
plot(SPcaGrid(xsd, lambda=0, method="sd", k=4), main="Standard PCA")    # standard
plot(SPcaGrid(xqn, lambda=0, method="Qn", k=4))                         # robust, non-sparse

plot(SPcaGrid(xqn, lambda=1,43, method="sd", k=4), main="Stdandard sparse PCA")  # sparse
plot(SPcaGrid(xqn, lambda=2.36, method="Qn", k=4), main="Robust sparse PCA")     # robust sparse
par(usr)

#########################################
##  Table 2 in Croux et al
##  - to compute EV=Explained variance and Cumulative EV we
##      need to get all 14 eigenvalues
##
rpca <- SPcaGrid(xqn, lambda=0, k=14)
srpca <- SPcaGrid(xqn, lambda=2.36, k=14)
tab <- cbind(round(getLoadings(rpca)[,1:4], 2), round(getLoadings(srpca)[,1:4], 2))

vars1 <- getEigenvalues(rpca);  vars1 <- vars1/sum(vars1)
vars2 <- getEigenvalues(srpca); vars2 <- vars2/sum(vars2)
cvars1 <- cumsum(vars1)
cvars2 <- cumsum(vars2)
ev <- round(c(vars1[1:4], vars2[1:4]),2)
cev <- round(c(cvars1[1:4], cvars2[1:4]),2)
rbind(tab, ev, cev)
}

}
\keyword{datasets}
