\name{gm.gamma}
\alias{gm.gamma}
\title{ Conditional Gamma coefficient estimation and confidence intervals }
\description{
  Calculates conditional gamma coefficients of two binary or ordinal variables X and Y given a 
  set of variables \{A,B,...\}.
}
\usage{
  gm.gamma(X = 0, Y = 0, data = 0, conditions = 0,
        type = c("conditional", "single", "marginal"), conf.level = 0.95)
}
\arguments{
  \item{X}{ Index of the variable's position in \code{data} or a vector. }
  \item{Y}{ Index of the variable's position in \code{data} or a vector. }
  \item{data}{ Data frame or table }
  \item{conditions}{ Variable indices in \code{data} or a data frame of the conditioning variables. }
  \item{type}{ Character string specifying the method,
                must be one of "conditional" (default), "single" or "marginal". 
               May be abbreviated. }
  \item{conf.level}{ Confidence level of the interval (default is 0.95). }
}
\details{
    The gamma coefficient is a correlation measure between two
    binary or ordinal variables. It ranges between -1 and 1, where -1 or respectively
    +1 stands for a purely negative or positive monotone relation. The relation has
    not to be of linear nature!

    If \code{type} = "conditional", the conditional gamma coefficients are calculated and
    if \code{type} = "marginal", the marginal gamma coefficients are computed. \cr
    If \code{X} and \code{Y} are given, the "single" gamma coefficient between both
    variables are computed.\cr
    If \code{X} or \code{Y} are zero, the function computes all possible conditional gamma coefficients.
           
    Confidence intervals are calculated using the asymptotic variance given in Olszak and Ritschard (1995).
}
\value{
  A matrix containing the gamma estimate(s), standard deviation(s), confidence interval(s) and p-value(s).
}
\references{ 
  Davis JA (1967) 
  \emph{A partial coefficient for Goodman and Kruskal's gamma.}
  Journal of the American Statistical Association, 62:189-193.

  Olszak M, Ritschard G (1995) 
  \emph{The behaviour of nominal and ordinal partial association measures.}
  The Statistician, 44(2):195-212.
 }
\author{ 
  Ronja Foraita, Fabian Sobotka \cr
  Bremen Institute for Prevention Research and Social Medicine \cr
  (BIPS)  \url{http://www.bips.uni-bremen.de}
   }
\seealso{ \code{\link{gm.or}}, \code{\link{gm.rr}}}
\examples{
  data(dp)

  ### Conditional Gamma by victime
  gm.gamma(1,3,conditions=2,data=dp)
  ### the same
  gm.gamma(dp$Defendants.Race,dp$Death.Penalty,data=dp,conditions=dp$Victims.Race)
  
  ### Stratified Gamma
  dp.black <- data.frame(victime=dp$Victims.Race[dp$Victims.Race=="black"],
                        killer=dp$Defendants.Race[dp$Victims.Race=="black"],
                        death.penalty=dp$Death.Penalty[dp$Victims.Race=="black"])
  dp.white <- data.frame(victime=dp$Victims.Race[dp$Victims.Race=="white"],
                        killer=dp$Defendants.Race[dp$Victims.Race=="white"],
                        death.penalty=dp$Death.Penalty[dp$Victims.Race=="white"])  
  table(dp.black[,c(2,3,1)])
  table(dp.white[,c(2,3,1)])  

  gm.gamma(2,3,data=dp.black)  
  gm.gamma(2,3,data=dp.white)  
  
  ### Marginal Gamma
  gm.gamma(1,3,data=dp)

  ### Analyse complete data set
  gm.gamma(data=dp,type="m")
  
  ### Plot model
  gamma <- gm.gamma(data=dp)
   #> all edges
  mat <- matrix(NA,nrow=3,ncol=3)
  mat[upper.tri(mat)] <- gamma[,1]
  gm.plot(model="abc",data.analysis=mat)
   #> only significant edges
  mat <- matrix(NA,nrow=3,ncol=3)   
  tmp <- vector()
  for( i in 1:dim(gamma)[1] ) ifelse(gamma[i,5]<0.05, tmp[i] <- gamma[i,1], tmp[i] <-NA)
  mat[upper.tri(mat)] <- tmp
  gm.plot(model="ab,bc",data.analysis=mat)
}
\keyword{ htest }
