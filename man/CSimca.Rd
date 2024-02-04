\name{CSimca}
\alias{CSimca}
\alias{CSimca.formula}
\alias{CSimca.default}
\title{
    Classification in high dimensions based on the (classical) SIMCA method
}
\description{
    CSimca performs the (classical) SIMCA method. This method classifies
    a data matrix x with a known group structure. To reduce the dimension on 
    each group a PCA analysis is performed. Afterwards a classification
    rule is developped to determine the assignment of new observations. 
}
\usage{
CSimca(x, ...)
\method{CSimca}{default}(x, grouping, prior=proportions, k, kmax = ncol(x), 
    tol = 1.0e-4, trace=FALSE, \dots)
\method{CSimca}{formula}(formula, data = NULL, \dots, subset, na.action)
}
\arguments{
  \item{formula}{a formula of the form \code{y~x}, it describes the response
  and the predictors. The formula can be more complicated, such as
  \code{y~log(x)+z} etc (see \code{\link{formula}} for more details).
  The response should
  be a factor representing the response variable, or any vector
  that can be coerced to such (such as a logical variable).}
  \item{data}{an optional data frame (or similar: see
    \code{\link{model.frame}}) containing the variables in the
    formula \code{formula}.}
  \item{subset}{an optional vector used to select rows (observations) of the
    data matrix \code{x}.}
  \item{na.action}{a function which indicates what should happen
    when the data contain \code{NA}s.  The default is set by
    the \code{na.action} setting of \code{\link{options}}, and is
    \code{\link{na.fail}} if that is unset. The default is \code{\link{na.omit}}.}
  \item{x}{a matrix or data frame containing the explanatory variables (training set). }
  \item{grouping}{grouping variable:  a factor specifying the class for each observation.}
  \item{prior}{prior probabilities, default to the class proportions for the training set.}
  \item{tol}{tolerance}
  \item{k}{number of principal components to compute. If \code{k} is missing, 
    or \code{k = 0}, the algorithm itself will determine the number of 
    components by finding such \code{k} that \eqn{l_k/l_1 >= 10.E-3} and 
    \eqn{\Sigma_{j=1}^k l_j/\Sigma_{j=1}^r l_j >= 0.8}. 
    It is preferable to investigate the scree plot in order to choose the number 
    of components and then run again. Default is \code{k=0}. }
  \item{kmax}{maximal number of principal components to compute.
    Default is \code{kmax=10}. If \code{k} is provided, \code{kmax} 
    does not need to be specified, unless \code{k} is larger than 10.}
  \item{trace}{whether to print intermediate results. Default is \code{trace = FALSE}}
  \item{\dots}{arguments passed to or from other methods.}
}

\details{
\code{CSimca}, serving as a constructor for objects of class \code{\link{CSimca-class}} 
is a generic function with "formula" and "default" methods.

SIMCA is a two phase procedure consisting of PCA performed on each group 
separately for dimension reduction followed by classification rules built 
in the lower dimensional space (note that the dimension in 
each group can be different). In original SIMCA new observations are 
classified by means of their deviations from the different PCA models.
Here (and also in the robust versions implemented in this package) the classification
rules will be obtained using two popular distances arising from PCA - 
orthogonal distances (OD) and score distances (SD). For the definition of these distances,
the definition of the cutoff values and the standartization of the distances see 
Vanden Branden K, Hubert M (2005) and Todorov and Filzmoser (2009).
}
\value{
  An S4 object of class \code{\link{CSimca-class}} which is a subclass of of the 
  virtual class \code{\link{Simca-class}}. 
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
data(pottery)
dim(pottery)        # 27 observations in 2 classes, 6 variables
head(pottery)

## Build the SIMCA model. Use RSimca for a robust version
cs <- CSimca(origin~., data=pottery)
cs
summary(cs)


## generate a sample from the pottery data set -
##  this will be the "new" data to be predicted
smpl <- sample(1:nrow(pottery), 5)
test <- pottery[smpl, -7]          # extract the test sample. Remove the last (grouping) variable
print(test)


## predict new data
pr <- predict(cs, newdata=test)

pr@classification 
}
\keyword{robust}
\keyword{multivariate}
