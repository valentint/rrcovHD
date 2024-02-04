\name{SosDiscRobust}
\alias{SosDiscRobust}
\alias{SosDiscRobust.default}
\alias{SosDiscRobust.formula}
\title{
    Robust and sparse multigroup classification by the optimal scoring approach
}
\description{
    Robust and sparse multigroup classification by the optimal scoring approach is robust against outliers, provides a low-dimensional and sparse representation of the predictors and is also applicable if the number of variables exeeds the number of observations.
}
\usage{
SosDiscRobust(x, ...)
\method{SosDiscRobust}{default}(x, grouping, prior=proportions, 
    lambda, Q=length(unique(grouping))-1, alpha=0.5, maxit=100, 
    tol = 1.0e-4, trace=FALSE, \dots)
\method{SosDiscRobust}{formula}(formula, data = NULL, \dots, subset, na.action)
}
\arguments{
  \item{formula}{A formula of the form \code{y~x}, it describes the response
  and the predictors. The formula can be more complicated, such as
  \code{y~log(x)+z} etc (see \code{\link{formula}} for more details).
  The response should
  be a factor representing the response variable, or any vector
  that can be coerced to such (such as a logical variable).}
  \item{data}{An optional data frame (or similar: see
    \code{\link{model.frame}}) containing the variables in the
    formula \code{formula}.}
  \item{subset}{An optional vector used to select rows (observations) of the
    data matrix \code{x}.}
  \item{na.action}{A function which indicates what should happen
    when the data contain \code{NA}s.  The default is set by
    the \code{na.action} setting of \code{\link{options}}, and is
    \code{\link{na.fail}} if that is unset. The default is \code{\link{na.omit}}.}
  \item{x}{A matrix or data frame containing the explanatory variables (training set); 
    colnames of x have to be provided. }
  \item{grouping}{Grouping variable:  a factor specifying the class for each observation.}
  \item{prior}{Prior probabilities, a vector of positive numbers that sum up to 1;
    default to the class proportions for the training set.}
  \item{lambda}{A non-negative tuning parameter for L1 norm penalty introducing sparsity on the 
    optimal scoring coefficients \eqn{\boldsymbol{\beta}_h} (see Details). 
    If the number of variables exceeds the number of observations \code{lambda} has to be positive.}
  \item{Q}{Number of optimal scoring coefficient vectors; \code{Q} has to be smaller than the number of groups. 
    Defaults to number of groups - 1.}
  \item{alpha}{Robustness parameter used in sparseLTS (for initial estimation, see Details). Default \code{alpha=0.5}.}
  \item{maxit}{Number of iterations for the estimation of optimal scoring coefficients and case weights. Default \code{maxit=100}.}
  \item{tol}{Tolerance for convergence of the normed weighted change in the residual sum of squares
    for the estiamtion of optimal scoring coefficeints. Default is \code{tol=1.0e-4}.}
  \item{trace}{Whether to print intermediate results. Default is \code{trace = FALSE}.}
  \item{\dots}{Arguments passed to or from other methods.}
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
\frac{1}{n} \theta_h^T Y^T Y\theta_h=1, \quad \theta_h^T Y^T Y \theta_l=0 \quad \forall l<h,
}{
1/n \theta_h^T Y^T Y \theta_h = 1,  \theta_h^T Y^T Y \theta_l = 0 for all l<h.
}

where \eqn{X} deontes the robustly centered and scaled input matrix \code{x} (or alternativly the predictors from \code{formular}) and \eqn{Y} is an dummy matrix coding die classmemberships from \code{grouping}.

For each \eqn{h} this problem can be solved interatively for \eqn{\beta_h} and \eqn{\theta_h}. In order to obtain robust estimates, \eqn{\beta_h} is estimated with reweighted sparse least trimmed squares regression (Alfons et al, 2013) and \eqn{\theta_h} with least absolut deviation regression in the first two iterations. To speed up the following repetitions an iterative down-weighting of observations with large residuals is combined with the iterative estimation of the optimal scoring coefficients with their classical estimates.

The classification model is estimated on the low dimensional sparse subspace \eqn{X[\beta_1,...,\beta_Q]} with robust LDA (\code{\link[rrcov]{Linda}}).
}
\value{
  An S4 object of class \code{\link{SosDiscRobust-class}} which is a subclass of of the 
  virtual class \code{\link{SosDisc-class}}. 
}
\references{
    Clemmensen L, Hastie T, Witten D & Ersboll B (2011),
    Sparse discriminant analysis.
    \emph{Technometrics}, \bold{53}(4), 406--413. 
    
    Alfons A, Croux C & Gelper S (2013),
    Sparse least trimmed squares regression for analysing high-dimensional large data sets.
    \emph{The Annals of Applied Statistics}, \bold{7}(1), 226--248.

    Ortner I, Filzmoser P & Croux C (2020),
    Robust and sparse multigroup classification by the optimal scoring approach.
    Data Mining and Knowledge Discovery \bold{34}, 723--741.
    \doi{10.1007/s10618-019-00666-8}.     
}


\author{ 
Irene Ortner \email{irene.ortner@applied-statistics.at} and Valentin Todorov \email{valentin.todorov@chello.at} 
}

\examples{

## EXAMPLE 1 ######################################
data(olitos)
grind <- which(colnames(olitos)=="grp")

set.seed(5008642)
mod <- SosDiscRobust(grp~., data=olitos, lambda=0.3, maxIte=30, Q=3, tol=1e-2)

pred <- predict(mod, newdata=olitos[,-grind])

summary(mod)
plot(mod, ind=c(1:3))


## EXAMPLE 2 ######################################
##

\dontrun{
library(sparseLDA)
data(penicilliumYES)

## for demonstration only:
set.seed(5008642)
X <- penicilliumYES$X[, sample(1:ncol(penicilliumYES$X), 100)]

## takes a subsample of the variables
## to have quicker computation time

colnames(X) <- paste0("V",1:ncol(X))
y <- as.factor(c(rep(1,12), rep(2,12), rep(3,12)))

set.seed(5008642)
mod <- SosDiscRobust(X, y, lambda=1, maxit=5, Q=2, tol=1e-2)

summary(mod)
plot(mod)
}

}
\keyword{robust}
\keyword{multivariate}
\keyword{classification}
\keyword{sparse}

