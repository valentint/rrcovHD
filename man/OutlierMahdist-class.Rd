\name{OutlierMahdist-class}
\Rdversion{1.1}
\docType{class}
\alias{OutlierMahdist-class}
\alias{getCutoff,OutlierMahdist-method}
\alias{getDistance,OutlierMahdist-method}

\title{Class \code{OutlierMahdist} - Outlier identification using 
    robust (mahalanobis) distances based on robust multivariate 
    location and covariance matrix}
\description{
    Holds the results of outlier identification using robust mahalanobis 
    distances computed by robust multivarite location and covarince matrix.   
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("OutlierMahdist", ...)} but the
    usual way of creating \code{OutlierMahdist} objects is a call to the function
    \code{OutlierMahdist()} which serves as a constructor.
}
\section{Slots}{
  \describe{
    \item{\code{covobj}:}{A list containing the robust estimates 
        of multivariate location and covariance matrix for each class}
    \item{\code{call}:}{Object of class \code{"language"}}
    \item{\code{counts}:}{Number of observations in each class}
    \item{\code{grp}:}{Grouping variable}
    \item{\code{wt}:}{Weights}
    \item{\code{flag}:}{0/1 flags identifying the outliers}
    \item{\code{method}:}{Method used to compute the robust estimates 
        of multivariate location and covariance matrix}
    \item{\code{singularity}:}{a list with singularity 
    information for the covariance matrix (or \code{NULL} 
    of not singular)}
  }
}
\section{Extends}{
Class \code{"\linkS4class{Outlier}"}, directly.
}
\section{Methods}{
  \describe{
    \item{getCutoff}{Return the cutoff value used to identify outliers }
    \item{getDistance}{Return a vector containing the computed distances }
	 }
}
\references{
    Todorov V & Filzmoser P (2009).
    An Object Oriented Framework for Robust Multivariate Analysis.
    \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
    \doi{10.18637/jss.v032.i03}.   
  
    Filzmoser P & Todorov V (2013).
    Robust tools for the imperfect world,
    \emph{Information Sciences} \bold{245}, 4--20.
    \doi{10.1016/j.ins.2012.10.017}.
  
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} }
\seealso{
     \code{\link{OutlierMahdist}}, \code{\link{Outlier-class}}
}
\examples{
showClass("OutlierMahdist")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}

