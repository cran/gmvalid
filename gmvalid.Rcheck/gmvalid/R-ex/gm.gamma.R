### Name: gm.gamma
### Title: Conditional Gamma coefficient estimation and confidence
###   intervals
### Aliases: gm.gamma
### Keywords: htest

### ** Examples

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



