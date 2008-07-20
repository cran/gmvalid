\name{gm.chi}
\alias{gm.chi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Conditional chi-squared test }
\description{
  Conducts a conditional chi-squared test for the null hypothesis
  X independent Y given Z.
}
\usage{
gm.chi(data, X = 1, Y = 2, Z = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ Data frame or a table (array). Variables should have names, \code{data} has to be discrete. }
  \item{X}{ Number of column of the first variable in data. }
  \item{Y}{ Number of column of the second variable in data. }
  \item{Z}{ Number of column of the variable in condition in data. }
}
\details{
  The conditional chi-squared test is conducted like a marginal chi-squared test where the
  expected cell counts are calculated as:
  E(i,j,k) = E(i,+,k)*E(+,j,k)/E(+,+,k)
}
\value{
  \item{ chi.squared }{ Value of the chi-squared statistic. }
  \item{ DF }{ Degrees of freedom. }
  \item{ p.value }{ P-value of the test. }
}
\references{ 
  Agresti A (1990)
  \emph{Categorial Data Analysis}
  Wiley, New York.
}
\author{   
  Ronja Foraita, Fabian Sobotka \cr
  Bremen Institute for Prevention Research and Social Medicine \cr
  (BIPS)  \url{http://www.bips.uni-bremen.de}
   }

\seealso{ \code{\link{gm.gamma}}, \code{\link[stats]{chisq.test}} }
\examples{
  data(dp)
  gm.chi(dp,1,2,3)
  gm.chi(dp,2,3,1)
  gm.chi(dp,1,3,2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ attribute }
