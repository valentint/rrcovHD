\name{SummarySosDisc-class}
\Rdversion{1.1}
\docType{class}
\alias{SummarySosDisc-class}
\alias{show,SummarySosDisc-method}

\title{Class \code{"SummarySosDisc"} - summary of \code{"SosDisc"} objects}
\description{ 
    Contains summary information about a \code{SosDisc} object representing 
    the results of the robust and sparse multigroup classification 
    by the optimal scoring approach.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SummarySosDisc", ...)},
but most often by invoking \code{summary()} on an \code{"SosDisc"} object. 
They contain values meant for printing by \code{show()}.
}
\section{Slots}{
    \describe{
        \item{\code{obj}:}{Object of class \code{"SosDisc"}}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "SummarySosDisc")}: display the object }
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
showClass("SummarySosDisc")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}
\keyword{classification}
\keyword{sparse}

