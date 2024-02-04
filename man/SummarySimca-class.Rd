\name{SummarySimca-class}
\Rdversion{1.1}
\docType{class}
\alias{SummarySimca-class}
\alias{show,SummarySimca-method}

\title{Class \code{"SummarySimca"} - summary of \code{"Simca"} objects}
\description{ 
    Contains summary information about a \code{Simca} object - 
    classification in high dimensions based on the SIMCA method 
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SummarySimca", ...)},
but most often by invoking \code{summary()} on an \code{"Simca"} object. 
They contain values meant for printing by \code{show()}.
}
\section{Slots}{
    \describe{
        \item{\code{simcaobj}:}{Object of class \code{"Simca"}}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "SummarySimca")}: display the object }
     }
}
\references{
    Todorov V & Filzmoser P (2009),
    An Object Oriented Framework for Robust Multivariate Analysis.
    \emph{Journal of Statistical Software}, \bold{32}(3), 1--47,
    \doi{10.18637/jss.v032.i03}.   

    Todorov V & Filzmoser P (2014),
    Software Tools for Robust Analysis of High-Dimensional Data.
    \emph{Austrian Journal of Statistics}, \bold{43}(4),  255--266,
    \doi{10.17713/ajs.v43i4.44}.       
}
\author{Valentin Todorov \email{valentin.todorov@chello.at}}
\seealso{
  \code{\link{Simca-class}}
}
\examples{
showClass("SummarySimca")
}
\keyword{classes}
\keyword{robust}
\keyword{multivariate}

