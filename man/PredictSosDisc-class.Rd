\name{PredictSosDisc-class}
\Rdversion{1.1}
\docType{class}
\alias{PredictSosDisc-class}
\alias{show,PredictSosDisc-method}

\title{Class \code{"PredictSosDisc"} - prediction of \code{"SosDisc"} objects}
\description{ The prediction of a "SosDisc" object }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("PredictSosDisc", ...)}
but most often by invoking \code{predict()} on a \code{"SosDisc"} object. 
They contain values meant for printing by \code{show()}}
\section{Slots}{
  \describe{
    \item{\code{classification}:}{Object of class \code{"factor"} representing the predicted classification }
    \item{\code{mahadist2}:}{A \code{"matrix"} containing the squared robust Mahalanobis distances to each group center in the subspace (see Details).}
    \item{\code{w}:}{A \code{"vector"} containing the weights derived from robust Mahalanobis distances to the closest group center (see Details).}
  }
}
\section{Details}{
For the prediction of the class membership a two step approach is taken. 
First, the \code{newdata} are scaled and centered (by obj@scale and obj@center) 
and  multiplied by \code{obj@beta} for dimension reduction. Then the 
classification of the transformed data is obtained by prediction with 
the \code{Linda} object obj@fit. The Mahalanobis distances to the closest 
group center in this subspace is used to derive case weights \code{w}. 
Observations where the squared robust mahalanobis distance is larger 
than the 0.975 quantile of the chi-square distribution with Q degrees 
of freedom receive weight zero, all others weight one.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "PredictSosDisc")}: Prints the results.}
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
\seealso{
  \code{\link{SosDisc-class}}
}
\examples{
    showClass("PredictSosDisc")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
\keyword{classification}
\keyword{sparse}

