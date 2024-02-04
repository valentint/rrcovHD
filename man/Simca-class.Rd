\name{Simca-class}
\Rdversion{1.1}
\docType{class}
\alias{Simca-class}
\alias{predict,Simca-method}
\alias{show,Simca-method}
\alias{summary,Simca-method}

\title{Class \code{"Simca"} - virtual base class for all classic and robust SIMCA 
classes representing classification in high dimensions based on the SIMCA method}
\description{
The class \code{Simca} searves as a base class for deriving all other 
classes representing the results of the classical and robust SIMCA methods}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
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
\section{Methods}{
  \describe{
    \item{predict}{\code{signature(object = "Simca")}: calculates prediction using the results in 
        \code{object}. An optional data frame or matrix in which to look for variables with which 
        to predict. If omitted, the training data set is used. If the original fit used a formula or 
        a data frame or a matrix with column names, newdata must contain columns with the 
        same names. Otherwise it must contain the same number of columns, 
        to be used in the same order. }
    \item{show}{\code{signature(object = "Simca")}: prints the results }
    \item{summary}{\code{signature(object = "Simca")}: prints summary information }
	 }
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
showClass("Simca")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}


