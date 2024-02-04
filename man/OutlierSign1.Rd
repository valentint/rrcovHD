\name{OutlierSign1}
\alias{OutlierSign1}
\alias{OutlierSign1.formula}
\alias{OutlierSign1.default}

\title{
Outlier identification in high dimensions using the SIGN1 algorithm 
}

\description{
    Fast algorithm for identifying multivariate outliers in high-dimensional 
    and/or large datasets, using spatial signs, see Filzmoser, Maronna, and Werner (CSDA, 2007). 
    The computation of the distances is based on Mahalanobis distances. 
}

\usage{
    OutlierSign1(x, ...)
    \method{OutlierSign1}{default}(x, grouping, qcrit = 0.975, trace=FALSE, \dots)
    \method{OutlierSign1}{formula}(formula, data, \dots, subset, na.action)
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
  \item{x}{a matrix or data frame. }
  \item{grouping}{grouping variable:  a factor specifying the class for each observation.}
  \item{qcrit}{a numeric value between 0 and 1 indicating the quantile to be used as critical value for outlier detection (default to 0.975).}
  \item{trace}{whether to print intermediate results. Default is \code{trace = FALSE}}
}
\details{
    Based on the robustly sphered and normed data, robust principal components are 
    computed. These are used for computing the covariance matrix which is the basis 
    for Mahalanobis distances. A critical value from the chi-square distribution 
    is then used as outlier cutoff. }
\value{
  An S4 object of class \code{\linkS4class{OutlierSign1}} which 
  is a subclass of the virtual class \code{\linkS4class{Outlier}}.
}
\references{
    P. Filzmoser, R. Maronna and M. Werner (2008).
    Outlier identification in high dimensions,
    \emph{Computational Statistics & Data Analysis}, Vol. 52 1694--1711.
    
    Filzmoser P & Todorov V (2013).
    Robust tools for the imperfect world,
    \emph{Information Sciences} \bold{245}, 4--20.
    \doi{10.1016/j.ins.2012.10.017}.
}

\author{ Valentin Todorov \email{valentin.todorov@chello.at} 
}

\seealso{
     \code{\linkS4class{OutlierSign1}}, \code{\linkS4class{OutlierSign2}}, \code{\linkS4class{Outlier}} 
}

\examples{

data(hemophilia)
obj <- OutlierSign1(gr~.,data=hemophilia)
obj

getDistance(obj)            # returns an array of distances
getClassLabels(obj, 1)      # returns an array of indices for a given class
getCutoff(obj)              # returns an array of cutoff values (for each class, usually equal)
getFlag(obj)                #  returns an 0/1 array of flags
plot(obj, class=2)          # standard plot function
}
\keyword{robust}
\keyword{multivariate}
