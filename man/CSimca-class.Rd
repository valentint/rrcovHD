\name{CSimca-class}
\Rdversion{1.1}
\docType{class}
\alias{CSimca-class}

\title{Class \code{"CSimca"} - classification in high dimensions based on the (classical) SIMCA method}
\description{
    The class \code{CSimca} represents the SIMCA algorithm for classification in high dimensions. 
    The objects of class \code{CSImca} contain the results
    of the SIMCA method.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CSImca", ...)} but the 
    usual way of creating \code{CSimca} objects is a call to the function
    \code{CSimca()} which serves as a constructor.
}
\section{Slots}{
  \describe{
    \item{\code{call}:}{the (matched) function call.}
    \item{\code{prior}:}{prior probabilities used, default to group proportions}
    \item{\code{counts}:}{number of observations in each class}
    \item{\code{pcaobj}:}{A list of Pca objects - one for each group}
    \item{\code{k}:}{Object of class \code{"numeric"} number of (choosen) principal components }
    \item{\code{flag}:}{Object of class \code{"Uvector"} The observations whose score distance is larger 
    than cutoff.sd or whose orthogonal distance is larger than cutoff.od can be considered 
    as outliers and receive a flag equal to zero.
    The regular observations receive a flag 1 }
    \item{\code{X}:}{the training data set (same as the input parameter x of the constructor function)}
    \item{\code{grp}:}{grouping variable:  a factor specifying the class for each observation.}   
  }
}
\section{Extends}{
Class \code{"\linkS4class{Simca}"}, directly.
}
\section{Methods}{
No methods defined with class "CSimca" in the signature.
}
\references{
    Vanden Branden K, Hubert M (2005) Robust classification in high 
    dimensions based on the SIMCA method. Chemometrics and 
    Intellegent Laboratory Systems 79:10--21

    Todorov V & Filzmoser P (2009),
    An Object Oriented Framework for Robust Multivariate Analysis.
    \emph{Journal of Statistical Software}, \bold{32}(3), 1--47,
    \doi{10.18637/jss.v032.i03}.   

    Todorov V & Filzmoser P (2014),
    Software Tools for Robust Analysis of High-Dimensional Data.
    \emph{Austrian Journal of Statistics}, \bold{43}(4),  255--266,
    \doi{10.17713/ajs.v43i4.44}.     
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} 
}

\examples{
showClass("CSimca")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}

