\name{OutlierPCDist}
\alias{OutlierPCDist}
\alias{OutlierPCDist.formula}
\alias{OutlierPCDist.default}

\title{
Outlier identification in high dimensions using the PCDIST algorithm 
}
\description{
    The function implements a simple, automatic outlier detection method suitable
    for high dimensional data that treats each class independently and uses 
    a statistically principled threshold for outliers. The algorithm can 
    detect both mislabeled and abnormal samples without reference to other classes.
}
\usage{
    OutlierPCDist(x, ...)
    \method{OutlierPCDist}{default}(x, grouping, control, k, explvar, trace=FALSE, \dots)
    \method{OutlierPCDist}{formula}(formula, data, \dots, subset, na.action)
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

    \item{control}{ a control object (S4) for one of the available control classes,
    e.g. \code{\link{CovControlMcd-class}}, \code{\link{CovControlOgk-class}},
    \code{\link{CovControlSest-class}}, etc.,
    containing estimation options. The class of this object defines
    which estimator will be used. Alternatively a character string can be specified
    which names the estimator - one of auto, sde, mcd, ogk, m, mve, sfast, surreal,
    bisquare, rocke. If 'auto' is specified or the argument is missing, the
    function will select the estimator (see below for details)}
   
    \item{k}{Number of components to select for PCA. If missing, the number 
    of components will be calculated automatically}
    \item{explvar}{Minimal explained variance to be used for calculation of 
    the number of components in PCA. If \code{explvar} is not provided, 
    automatic dimensionality selection using profile likelihood, as 
    proposed by Zhu and Ghodsi will be used.}
    \item{trace}{whether to print intermediate results. Default is \code{trace = FALSE}}
}

\details{
    If the data set consists of two or more classes 
    (specified by the grouping variable \code{grouping}) the proposed method iterates
    through the classes present in the data, separates each class from the rest and
    identifies the outliers relative to this class, thus treating both types of outliers,
    the mislabeled and the abnormal samples in a homogenous way.
    
    The first step of the algorithm is dimensionality reduction using (classical) PCA. 
    The number of components to select can be provided by the user but if missing, 
    the number of components will be calculated either using the provided minimal 
    explained variance or by the automatic dimensionality selection using 
    profile likelihood, as proposed by Zhu and Ghodsi.
}

\value{
  An S4 object of class \code{\linkS4class{OutlierPCDist}} which 
  is a subclass of the virtual class \code{\linkS4class{Outlier}}.
}

\references{
    A.D. Shieh and Y.S. Hung  (2009).
    Detecting Outlier Samples in Microarray Data,
    \emph{Statistical Applications in Genetics and Molecular Biology} \bold{8}.

    M. Zhu, and A. Ghodsi (2006). Automatic dimensionality selection from the scree
    plot via the use of profile likelihood. \emph{Computational Statistics & Data Analysis},
    \bold{51}, pp. 918--930.

    Filzmoser P & Todorov V (2013).
    Robust tools for the imperfect world,
    \emph{Information Sciences} \bold{245}, 4--20.
    \doi{10.1016/j.ins.2012.10.017}.
}

\author{ Valentin Todorov \email{valentin.todorov@chello.at} 
}

\seealso{
     \code{\linkS4class{OutlierPCDist}}, \code{\linkS4class{Outlier}} 
}

\examples{

data(hemophilia)
obj <- OutlierPCDist(gr~.,data=hemophilia)
obj

getDistance(obj)            # returns an array of distances
getClassLabels(obj, 1)      # returns an array of indices for a given class
getCutoff(obj)              # returns an array of cutoff values (for each class, usually equal)
getFlag(obj)                #  returns an 0/1 array of flags
plot(obj, class=2)          # standard plot function
}
\keyword{robust}
\keyword{multivariate}

