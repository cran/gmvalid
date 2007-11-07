### Name: gm.boot.mim
### Title: Graphical model validation using the bootstrap (MIM).
### Aliases: gm.boot.mim
### Keywords: nonparametric multivariate models

### ** Examples

  ### should provide good results because of simulated data
  gm.a <- gm.modelsim(2000,"ABC,CDE")
  gm.boot.mim(50,gm.a)
  
  ### on real data sets a forward bootstrap seems to have better results
  data(wynder)
  gm.boot.mim(100,wynder,strategy="f",calculations=c("s","e"),options="u")
  
  ### with model given
  data(wam)
  gm.boot.mim(10,wam,model="a,bcde,cdef")
  
  ### EH-strategy
  gm.boot.mim(50,wam,strategy="eh",model="a,bc,de,f - abcde,bcdef")



