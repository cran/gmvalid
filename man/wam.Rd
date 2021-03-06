\name{wam}
\alias{wam}
\docType{data}
\title{ Women and Mathematics }
\description{
  The data set describes the attitude of high school students toward
  mathematics achievement and related topics.
}
\usage{data(wam)}
\format{
  A data frame with 1190 observations on the following 6 variables.
  \describe{
    \item{\code{lecture}}{a factor with levels \code{1} (yes) \code{2} (no) }
    \item{\code{gender}}{a factor with levels \code{1} (female) \code{2} (male) }
    \item{\code{school}}{a factor with levels \code{1} (suburban) \code{2} (urban) }
    \item{\code{work}}{a factor with levels \code{1} (will need math) \code{2} (won't) }
    \item{\code{subject}}{a factor with levels \code{1} (science) \code{2} (arts) }
    \item{\code{plans}}{a factor with levels \code{1} (college) \code{2} (job) }
  }
}
\details{
  1190 high school students of eight schools in New Jersey responded to 
  the questionnaire. The survey was performed by Lacampagne (1979).
}
\source{ 
  Fowlkes, EB, Freeny, AE and Landwehr, JM (1988) \emph{Evaluating logistic models for large contingency tables.}
  Journal of the American Statistical Association, 83:611-622.
}
\references{
  Lacampagne, CB (1979) \emph{An Evaluation of the Women and Mathematics (WAM) Program and Associated Sex-Related Differences in the Teaching, 
  Learning and Counseling of Mathematics.}
  Ed.D. Thesis, Columbia University Teachers College, USA.
}
\examples{
  data(wam)
\dontrun{
## until CoCo is available again ...
  gm.analysis(wam, program="coco")
}
  sel = gm.gms(wam)
}
\keyword{datasets}
