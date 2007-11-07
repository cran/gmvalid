### Name: gm.or, gm.rr
### Title: Stratified odds ratios or risk ratios
### Aliases: gm.or gm.rr
### Keywords: htest

### ** Examples

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




