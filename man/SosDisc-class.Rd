\name{SosDisc-class}
\Rdversion{1.1}
\docType{class}
\alias{SosDisc-class}
\alias{predict, SosDisc-method}
\alias{show, SosDisc-method}
\alias{summary, SosDisc-method}

\title{
    Class \code{"SosDisc"} - virtual base class for all classic and robust SosDisc 
    classes representing the results of the robust and sparse multigroup classification 
    by the optimal scoring approach
}
\description{
    Robust and sparse multigroup classification by the optimal scoring approach.
    The class \code{SosDisc} searves as a base class for deriving all other 
    classes representing the results of the robust and sparse multigroup classification 
    by the optimal scoring approach.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
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
\details{
The sparse optimal scoring problem (Clemmensen et al, 2011):
for \eqn{h=1,....,Q}
\deqn{
\min_{\beta_h,\theta_h} \frac{1}{n} \|Y \theta_h - X \beta_h \|_2^2   + \lambda \|\beta_h\|_1
}{
min{\beta_h,\theta_h} 1/n ||Y \theta_h - X \beta_h ||_2^2   + \lambda ||\beta_h||_1
}
subject to
\deqn{
\frac{1}{n} \theta_h^T Y^T Y\theta_h=1, \quad \theta_h^T Y^T Y \theta_l=0 \quad \forall l<h.
}{
1/n \theta_h^T Y^T Y \theta_h = 1,  \theta_h^T Y^T Y \theta_l = 0 for all l<h,
}

where \eqn{X} deontes the robustly centered and scaled input matrix \code{x} (or alternativly the predictors from \code{formular}) and \eqn{Y} is an dummy matrix coding die classmemberships from \code{grouping}.

For each \eqn{h} this problem can be solved interatively for \eqn{\beta_h} and \eqn{\theta_h}. In order to obtain robust estimates, \eqn{\beta_h} is estimated with reweighted sparse least trimmed squares regression (Alfons et al, 2013) and \eqn{\theta_h} with least absolut deviation regression in the first two iterations. To speed up the following repetitions an iterative down-weighting of observations with large residuals is combined with the iterative estimation of the optimal scoring coefficients with their classical estimates.

The classification model is estimated on the low dimensional sparse subspace \eqn{X[\beta_1,...,\beta_Q]} with robust LDA (\code{\link[rrcov]{Linda}}).
}
\section{Methods}{
  \describe{
    \item{predict}{\code{signature(object = "SosDisc")}: calculates prediction using the results in 
        \code{object}. An optional data frame or matrix in which to look for variables with which 
        to predict. If omitted, the training data set is used. If the original fit used a formula or 
        a data frame or a matrix with column names, newdata must contain columns with the 
        same names. }
    \item{show}{\code{signature(object = "SosDisc")}: prints the results }
    \item{summary}{\code{signature(object = "SosDisc")}: prints summary information }
    }
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
    showClass("SosDisc")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
\keyword{classification}
\keyword{sparse}

