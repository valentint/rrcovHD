\name{OutlierSign2-class}
\Rdversion{1.1}
\docType{class}
\alias{OutlierSign2-class}
\alias{getCutoff,OutlierSign2-method}
\alias{getDistance,OutlierSign2-method}

\title{Class \code{"OutlierSign2"} - Outlier identification in high dimensions using the SIGN2 algorithm}
\description{
    Fast algorithm for identifying multivariate outliers in high-dimensional 
    and/or large datasets, using spatial signs, see Filzmoser, Maronna, and Werner (CSDA, 2007). 
    The computation of the distances is based on principal components. 
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("OutlierSign2", ...)} but the
    usual way of creating \code{OutlierSign2} objects is a call to the function
    \code{OutlierSign2()} which serves as a constructor.
}
\section{Slots}{
  \describe{
    \item{\code{covobj}:}{A list containing intermediate results of the SIGN2 algorithm for each class}
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

\author{ Valentin Todorov \email{valentin.todorov@chello.at} 
}

\seealso{
     \code{\link{OutlierSign2}}, \code{\linkS4class{OutlierSign1}}, \code{\linkS4class{Outlier}} 
}

\examples{
showClass("OutlierSign2")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
