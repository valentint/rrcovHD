\name{OutlierMahdist}
\alias{OutlierMahdist}
\alias{OutlierMahdist.formula}
\alias{OutlierMahdist.default}
\title{
Outlier identification using robust (mahalanobis) distances based on robust multivariate location and covariance matrix 
}
\description{
    This function uses the Mahalanobis distance as a basis for 
    multivariate outlier detection. The standard method for multivariate 
    outlier detection is robust estimation of the parameters in the 
    Mahalanobis distance and the comparison with a critical value 
    of the Chi2 distribution (Rousseeuw and Van Zomeren, 1990). 
}
\usage{
    OutlierMahdist(x, ...)
    \method{OutlierMahdist}{default}(x, grouping, control, trace=FALSE, \dots)
    \method{OutlierMahdist}{formula}(formula, data, \dots, subset, na.action)
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
  \item{control}{ a control object (S4) for one of the available control classes,
    e.g. \code{\link{CovControlMcd-class}}, \code{\link{CovControlOgk-class}},
    \code{\link{CovControlSest-class}}, etc.,
    containing estimation options. The class of this object defines
    which estimator will be used. Alternatively a character string can be specified
    which names the estimator - one of auto, sde, mcd, ogk, m, mve, sfast, surreal,
    bisquare, rocke. If 'auto' is specified or the argument is missing, the
    function will select the estimator (see below for details)}

  \item{trace}{whether to print intermediate results. Default is \code{trace = FALSE}}
}

\details{
    If the data set consists of two or more classes 
    (specified by the grouping variable \code{grouping}) the proposed method iterates
    through the classes present in the data, separates each class from the rest and
    identifies the outliers relative to this class, thus treating both types of outliers,
    the mislabeled and the abnormal samples in a homogenous way.
  
    The estimation method is selected by the control object \code{control}. 
    If a character string naming an estimator is specified, a
    new control object will be created and used (with default estimation options).
    If this argument is missing or a character string
    'auto' is specified, the function will select the robust estimator
    according to the size of the dataset - for details see \code{\link{CovRobust}}.
}
\value{
  An S4 object of class \code{\linkS4class{OutlierMahdist}} which 
  is a subclass of the virtual class \code{\linkS4class{Outlier}}.
}
\references{
    P. J. Rousseeuw and B. C. Van Zomeren (1990). 
    Unmasking multivariate outliers and leverage points. 
    \emph{Journal of the American Statistical Association}. 
    Vol. 85(411), pp. 633-651.

    P. J. Rousseeuw and A. M. Leroy (1987).
    \emph{Robust Regression and Outlier Detection.} Wiley.

    P. J. Rousseeuw and K. van Driessen (1999)
    A fast algorithm for the minimum covariance determinant estimator.
    \emph{Technometrics} \bold{41}, 212--223.

    Todorov V & Filzmoser P (2009).
    An Object Oriented Framework for Robust Multivariate Analysis.
    \emph{Journal of Statistical Software}, \bold{32}(3), 1--47,
    \doi{10.18637/jss.v032.i03}.   
    
    Filzmoser P & Todorov V (2013).
    Robust tools for the imperfect world,
    \emph{Information Sciences} \bold{245}, 4--20.
    \doi{10.1016/j.ins.2012.10.017}.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} 
}
\examples{

data(hemophilia)
obj <- OutlierMahdist(gr~.,data=hemophilia)
obj

getDistance(obj)            # returns an array of distances
getClassLabels(obj, 1)      # returns an array of indices for a given class
getCutoff(obj)              # returns an array of cutoff values (for each class, usually equal)
getFlag(obj)                #  returns an 0/1 array of flags
plot(obj, class=2)          # standard plot function
}
\keyword{robust}
\keyword{multivariate}

