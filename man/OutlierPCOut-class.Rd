\name{OutlierPCOut-class}
\Rdversion{1.1}
\docType{class}
\alias{OutlierPCOut-class}
\alias{getCutoff,OutlierPCOut-method}
\alias{getDistance,OutlierPCOut-method}
\alias{plot,OutlierPCOut,missing-method}

\title{Class \code{"OutlierPCOut"} - Outlier identification in high dimensions using using the PCOUT algorithm}

\description{
    Holds the results of outlier identification using the PCOUT algorithm.   
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("OutlierPCOut", ...)} but the
    usual way of creating \code{OutlierPCOut} objects is a call to the function
    \code{OutlierPCOut()} which serves as a constructor.
}
\section{Slots}{
  \describe{
    \item{\code{covobj}:}{A list containing intermediate results of the PCOUT algorithm for each class}
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
    \item{getDistance}{Return a vector containing the computed distances }
    \item{plot}{Plot the results of the outlier detection process}
	 }
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
\author{ Valentin Todorov \email{valentin.todorov@chello.at} }
\seealso{
     \code{\link{OutlierPCOut}}, \code{"\linkS4class{Outlier}"}
}
\examples{
showClass("OutlierMahdist")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}


