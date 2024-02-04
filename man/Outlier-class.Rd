\name{Outlier-class}
\Rdversion{1.1}
\docType{class}
\alias{Outlier-class}
\alias{getClassLabels,Outlier-method}
\alias{getDistance,Outlier-method}
\alias{getFlag,Outlier-method}
\alias{getOutliers,Outlier-method}
\alias{getWeight,Outlier-method}
\alias{plot,Outlier,missing-method}
\alias{show,Outlier-method}

\title{Class \code{"Outlier"} -- a base class for outlier identification}
\description{
    The class \code{Outlier} represents the results of outlier identification. 
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{call}:}{Object of class \code{"language"}}
    \item{\code{counts}:}{Number of observations in each class}
    \item{\code{grp}:}{Grouping variable}
    \item{\code{wt}:}{Vector of weights}
    \item{\code{flag}:}{0/1 flags identifying the outliers}
    \item{\code{method}:}{A character string specifying the method used to 
    identify the outliers. In case of \code{\linkS4class{OutlierMahdist}} 
    class this is the name of the robust estimator
        of multivariate location and covariance matrix used}
    \item{\code{singularity}:}{a list with singularity 
    information for the covariance matrix (or \code{NULL} 
    if not singular)}
  }
}
\section{Methods}{
  \describe{
    \item{getClassLabels}{Returns a vector with indices for a given class}
    \item{getDistance}{Returns a vector containing the computed distances }
    \item{getFlag}{Returns the flags identifying the outliers}
    \item{getOutliers}{Returns a vector with the indices of the identified outliers}
    \item{getWeight}{Returns a vector of weights}
    \item{plot}{}
    \item{show}{}
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
\examples{
showClass("Outlier")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}

