\name{gm.or, gm.rr}
\alias{gm.or}
\alias{gm.rr}
\title{ Stratified odds ratios and risk ratios }
\description{
  Computes (stratified) odds ratios or risk ratios and their asymptotic confidence interval. 
}
\usage{
gm.or(X, Y, data = 0, conditions = 0,
        reference = c("last", "biggest", "first"), conf.level = 0.95)

gm.rr(X, Y, data = 0, conditions = 0,
        reference = c("last", "biggest", "first"), conf.level = 0.95)
}
\arguments{
  \item{X}{ Index of the variable's position in \code{data} or a vector. }
  \item{Y}{ Index of the variable's position in \code{data} or a vector. }
  \item{data}{ Data frame or table }
  \item{conditions}{ Variable indices in \code{data} or a data frame of the conditional variables. }
  \item{reference}{ Character string specifying the reference category.
                    Must be one of "last" (default), "biggest" (largest) or "first",
                    may be abbreviated. }
  \item{conf.level}{ Confidence level of the interval (default is 0.95). }
}
\details{
  Calculates odds ratios by conditional maximum likelihood estimation (Fisher) for stratified odds ratios and
  odds ratios by unconditional maximum likelihood estimation (Wald) for marginal odds ratios. 
  Confidence intervals are calculated using exact methods.
  
  Calculates risk ratios by unconditional maximum likelihood estimation (Wald). 
  Confidence intervals are calculated using normal approximation.
  
  Is based on the functions \code{oddsratio.fisher}, \code{oddsratio.wald} and \code{riskratio.wald} (package: epitools).
}
\value{

  A matrix containing estimate(s), confidence interval(s) and p-value(s).
}
\references{ 
  \url{http://www.epitools.net}
  
  Rothman KJ, Greenland S (1998) 
  \emph{Modern Epidemiology.} (2nd)
  Lippincott-Raven Publisher
 }
\author{ 
  Ronja Foraita, Fabian Sobotka \cr
  Bremen Institute for Prevention Research and Social Medicine \cr
  (BIPS)  \url{http://www.bips.uni-bremen.de}
 }
\seealso{ \code{\link[epitools]{oddsratio}}, \code{\link[epitools]{riskratio}} }
\examples{
  group  <- c("treatment","placebo1","placebo2")
  target <- c("low","medium","high")
  mat    <- matrix(c(78,35,53,77,10,89,16,119,32),nrow=3,ncol=3,byrow=TRUE, 
                    dimnames=list("group"=group,"target"=target))
  treat  <- data.frame(expand.table(mat))
  table(treat)
  
  ### Marginal OR
  gm.or(1,2,treat,reference="f")
  gm.or(treat$target,treat$group)
  
  ### Stratified OR
  data <- gm.modelsim(1000,"ab,bcd",list(c(1,1),c(1,1),c(1,1),c(1,1)))  
  gm.or(1,2,conditions=c(3,4),data=data)
  
  ### Marginal RR
  gm.rr(1,2,treat,reference="f")
  gm.rr(treat$target,treat$group)
  
  ### Stratified RR
  data <- gm.modelsim(1000,"ab,bcd",list(c(1,1),c(1,1),c(1,1),c(1,1)))  
  gm.rr(1,2,conditions=c(3,4),data=data)
  
  ### ALSO
  gm.or(X=data$a,Y=data$b,conditions=data$d)

}
\keyword{ htest }
