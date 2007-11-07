### Name: gm.csi
### Title: Conditional Synergy Index
### Aliases: gm.csi
### Keywords: htest

### ** Examples

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



