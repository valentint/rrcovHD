\name{OutlierPCDist-class}
\Rdversion{1.1}
\docType{class}
\alias{OutlierPCDist-class}
\alias{getCutoff,OutlierPCDist-method}
\alias{getDistance,OutlierPCDist-method}

\title{Class \code{"OutlierPCDist"} - Outlier identification in high dimensions using using the PCDIST algorithm}
\description{
    The function implements a simple, automatic outlier detection method suitable
    for high dimensional data that treats each class independently and uses 
    a statistically principled threshold for outliers. The algorithm can 
    detect both mislabeled and abnormal samples without reference to other classes.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("OutlierPCDist", ...)} but the
    usual way of creating \code{OutlierPCDist} objects is a call to the function
    \code{OutlierPCDist()} which serves as a constructor.
}
\section{Slots}{
  \describe{
    \item{\code{covobj}:}{A list containing intermediate results of the PCDIST algorithm for each class}
    \item{\code{k}:}{Number of selected PC}
    \item{\code{call}, \code{counts}, \code{grp}, \code{wt},
        \code{flag}, \code{method}, \code{singularity}:}{from the \code{"\linkS4class{Outlier}"} class.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{Outlier}"}, directly.
}

\section{Methods}{
  \describe{
    \item{getCutoff}{Return the cutoff value used to identify outliers }
	 }
}
\references{
    A.D. Shieh and Y.S. Hung  (2009).
    Detecting Outlier Samples in Microarray Data,
    \emph{Statistical Applications in Genetics and Molecular Biology} Vol. 8.

    Filzmoser P & Todorov V (2013).
    Robust tools for the imperfect world,
    \emph{Information Sciences} \bold{245}, 4--20.
    \doi{10.1016/j.ins.2012.10.017}.
}

\author{ Valentin Todorov \email{valentin.todorov@chello.at} 
}

\seealso{
     \code{\link{OutlierPCDist}}, \code{\linkS4class{Outlier}} 
}
\examples{
showClass("OutlierPCDist")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}

