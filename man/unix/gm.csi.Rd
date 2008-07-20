\name{gm.csi}
\alias{gm.csi}
\title{ Conditional Synergy Index }
\description{
  Computes the conditional synergy index of two dichotomous variables in relation to 
a binary group variable on the penetrance scale.
}
\usage{
gm.csi(X, Y, group, data=0, reference=c(1,1,2), pen = NULL, conf.level = 0.95)
}
\arguments{
  \item{X}{ Index of the variable's position in \code{data} or a vector. }
  \item{Y}{ Index of the variable's position in \code{data} or a vector. }
  \item{group}{ Binary group or outcome variable addressed as variable index in 
                \code{data} or as vector. By default, 2 is the reference category.}
  \item{data}{ Data frame or table. }
  \item{reference}{ Vector to define the \code{reference} categories of \code{X}, \code{Y} and \code{group}. 
                    By default, 2 is the reference category for group.}
  \item{pen}{ If FALSE the penetrance P( group | X,Y ) is estimated from the data. }
  \item{conf.level}{ Confidence level of the interval (default is 0.95). }
}
\details{
    Foraita's conditional synergy index (CSI) is an interaction measure on the penetrance scale
    between two dichotomous variables on a  binary outcome. 
    The index equals 1 under additivity, CSI > 1 in the case of antagonism and 
    CSI < 1 in the presence of synergy.
           
    The confidence intervals are calculated using an asymptotic variance given in Foraita (2007). 
}
\value{
A list containing:
  \item{penetrance}{ a matrix containing the joint penetrance. If pen = NULL, the penetrance will be
                     estimated from the data.  }
        
  \item{measure}{a matrix containing the estimate, variance, standard error,
                  confidence interval and p-value. Figures in brackets show the
                  reference category respectively the category under consideration.}
}
\references{ 
Foraita R (2007)
\emph{A conditional synergy index to assess biological interaction.}
Working Paper. Please send an e-mail to \email{foraita@bips.uni-bremen.de}.
 }
\author{ 
  Ronja Foraita, Fabian Sobotka \cr
  Bremen Institute for Prevention Research and Social Medicine \cr
  (BIPS)  \url{http://www.bips.uni-bremen.de}
   }
\seealso{ \code{\link{gm.si}} }
\examples{
  data(idd35)
  gm.csi(1,2,3,data=idd35)

  ### >> constructing an additive and multiplicative penetrance
  x <- c(0.1,0.4)
  y <- c(0.05,0.5)
  add.pen <- outer(x,y,FUN="+")
  mult.pen <- outer(x,y)
  het.pen <- outer(x,y,FUN="+") - outer(x,y)

  ### >> Function that samples data using the penetrance 
  make.data <- function(R,pen,category) 
    {
      s.vec <- sample(c(1,2,3,4),R,replace=TRUE,prob=as.vector(pen))
      fact.1 <- fact.2 <- vector()
      for( i in 1:R ) {
        ifelse( s.vec[i] == 1 || s.vec[i] == 3 , fact.1[i] <- 1, fact.1[i] <- 2 ) 
        ifelse( s.vec[i] == 1 || s.vec[i] == 2 , fact.2[i] <- 1, fact.2[i] <- 2 ) 
      }
      cbind(X=fact.1,Y=fact.2,group=rep(category,R))  
    }

  ### >>> Building datasets with affected and unaffected subjects   
  add.aff <- make.data(200,add.pen,2)
  add.uaf <- make.data(200,1-add.pen,1)  
  add.df <- as.data.frame(rbind(add.uaf,add.aff))
  
  mult.aff <- make.data(200,mult.pen,2)
  mult.uaf <- make.data(200,1-mult.pen,1)  
  mult.df <- as.data.frame(rbind(mult.uaf,mult.aff))
  
  het.aff <- make.data(200,het.pen,2)
  het.uaf <- make.data(200,1-het.pen,1)  
  het.df <- as.data.frame(rbind(het.uaf,het.aff))
   
  gm.csi(1,2,3,add.df,pen=add.pen)   # Additivity
  gm.csi(1,2,3,mult.df,pen=mult.pen) # Synergy
  gm.csi(1,2,3,het.df,pen=het.pen)   # Antagonism
}
\keyword{ htest }
