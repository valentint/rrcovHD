\name{OutlierPCOut}
\alias{OutlierPCOut}
\alias{OutlierPCOut.formula}
\alias{OutlierPCOut.default}
\title{
Outlier identification in high dimensions using the PCOUT algorithm 
}
\description{
    The function implements a computationally fast procedure for identifying 
    outliers that is particularly effective in high dimensions. This algorithm 
    utilizes simple properties of principal components to identify outliers in 
    the transformed space, leading to significant computational advantages 
    for high-dimensional data. This approach requires considerably less 
    computational time than existing methods for outlier detection, and is 
    suitable for use on very large data sets. It is also capable of analyzing 
    the data situation commonly found in certain biological applications in which 
    the number of dimensions is several orders of magnitude larger than the number of observations. 
}


\usage{
    OutlierPCOut(x, ...)
    \method{OutlierPCOut}{default}(x, grouping, explvar=0.99, trace=FALSE, \dots)
    \method{OutlierPCOut}{formula}(formula, data, \dots, subset, na.action)
}

\arguments{
  \item{formula}{a formula with no response variable, referring only to
    numeric variables.}
  \item{data}{an optional data frame (or similar: see
    \code{\link{model.frame}}) containing the variables in the
    formula \code{formula}.}
  \item{subset}{an optional vector used to select rows (observations) of the
    data matrix \code{x}.}
  \item{na.action}{a function which indicates what should happen
    when the data contain \code{NA}s.  The default is set by
    the \code{na.action} setting of \code{\link{options}}, and is
    \code{\link{na.fail}} if that is unset. The default is \code{\link{na.omit}}.}
  \item{\dots}{arguments passed to or from other methods.}
  \item{x}{a matrix or data frame. }
  \item{grouping}{grouping variable:  a factor specifying the class for each observation.}
  \item{explvar}{ a numeric value between 0 and 1 indicating how much variance
  should be covered by the robust PCs (default to 0.99) }    
  \item{trace}{whether to print intermediate results. Default is \code{trace = FALSE}}
}

\details{
    If the data set consists of two or more classes 
    (specified by the grouping variable \code{grouping}) the proposed method iterates
    through the classes present in the data, separates each class from the rest and
    identifies the outliers relative to this class, thus treating both types of outliers,
    the mislabeled and the abnormal samples in a homogenous way.
}
\value{
  An S4 object of class \code{\linkS4class{OutlierPCOut}} which 
  is a subclass of the virtual class \code{\linkS4class{Outlier}}.
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
     \code{\linkS4class{OutlierPCOut}}, \code{\linkS4class{Outlier}} 
}

\examples{

data(hemophilia)
obj <- OutlierPCOut(gr~.,data=hemophilia)
obj

getDistance(obj)            # returns an array of distances
getClassLabels(obj, 1)      # returns an array of indices for a given class
getCutoff(obj)              # returns an array of cutoff values (for each class, usually equal)
getFlag(obj)                #  returns an 0/1 array of flags
plot(obj, class=2)          # standard plot function
}
\keyword{robust}
\keyword{multivariate}
