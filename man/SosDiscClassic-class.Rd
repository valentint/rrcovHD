\name{SosDiscClassic-class}
\Rdversion{1.1}
\docType{class}
\alias{SosDiscClassic-class}

\title{Class \code{SosDiscClassic} - sparse multigroup classification by the optimal scoring approach}

\description{
    Sparse multigroup classification by the optimal scoring approach.
}

\section{Objects from the Class}{
    Objects can be created by calls of the form \code{new("SosDiscClassic", ...)} but the 
    usual way of creating \code{SosDiscClassic} objects is a call to the function
    \code{SosDiscRobust()} which serves as a constructor.
}
\section{Slots}{
  \describe{
    \item{\code{call}:}{The (matched) function call.}
    \item{\code{prior}:}{Prior probabilities; same as input parameter.}
    \item{\code{counts}:}{Number of observations in each class.}
    \item{\code{beta}:}{Object of class \code{"matrix"}: Q coefficient vectors of the predictor matrix from optimal scoring (see Details); 
      rows corespond to variables listed in \code{varnames}.}
    \item{\code{theta}:}{Object of class \code{"matrix"}: Q coefficient vectors of the dummy matrix for class coding from optimal scoring (see Details).}
    \item{\code{lambda}:}{Non-negative tuning paramer from L1 norm penaly; same as input parameter}
    \item{\code{varnames}:}{Character vector: Names of included predictor variables 
      (variables where at least one beta coefficient is non-zero).}
    \item{\code{center}:}{Centering vector of the input predictors (coordinate wise median).}
    \item{\code{scale}:}{Scaling vector of the input predictors (mad).}
    \item{\code{fit}:}{Object of class \code{"Linda"}: Linda model (robust LDA model) estimated in the low dimensional subspace \eqn{X[\beta_1,...,\beta_Q]} (see Details)}
    \item{\code{mahadist2}:}{These will go later to Linda object: squared robust Mahalanobis distance 
      (calculated with estimates from Linda, with common covariance structure of all groups) 
      of each observation to its group center in the low dimensional subspace \eqn{X[\beta_1,...,\beta_Q]} (see Details).}
    \item{\code{wlinda}:}{These will go later to Linda object: 0-1 weights derived from \code{mahadist2};
      observations where the squred robust Mahalanobis distance is larger than the 0.975 quantile 
      of the chi-square distribution with Q degrees of freedom resive weight zero.}
    \item{\code{X}:}{The training data set (same as the input parameter \code{x} of the constructor function)}
    \item{\code{grp}:}{Grouping variable: a factor specifying the class for each observation (same as the input parameter \code{grouping})}   
  }
}
\section{Extends}{
    Class \code{"\linkS4class{SosDisc}"}, directly.
}
\section{Methods}{
    No methods defined with class "SosDiscClassic" in the signature.
}
\references{
    Clemmensen L, Hastie T, Witten D & Ersboll B (2012),
    Sparse discriminant analysis.
    \emph{Technometrics}, \bold{53}(4), 406--413. 

    Ortner I, Filzmoser P & Croux C (2020),
    Robust and sparse multigroup classification by the optimal scoring approach.
    Data Mining and Knowledge Discovery \bold{34}, 723--741.
    \doi{10.1007/s10618-019-00666-8}.     
}
\author{ 
    Irene Ortner \email{irene.ortner@applied-statistics.at} and Valentin Todorov \email{valentin.todorov@chello.at} 
}

\examples{
    showClass("SosDiscClassic")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
\keyword{classification}
\keyword{sparse}

