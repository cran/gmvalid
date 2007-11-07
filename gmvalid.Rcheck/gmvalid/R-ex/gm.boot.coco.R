### Name: gm.boot.coco
### Title: Graphical model validation using the bootstrap (CoCo).
### Aliases: gm.boot.coco
### Keywords: nonparametric multivariate models

### ** Examples

  ### should provide good results because of simulated data
  gm <- gm.modelsim(2000,"ABC,CDE")
  gm.boot.coco(50,gm,recursive=TRUE)
  
  ### on real data sets a forward bootstrap seems to have better results
  data(wynder)
  gm.boot.coco(100,wynder,strategy="f",calculations=c("s","e"),decomposable.mode=TRUE)
  
  ### with a given model
  data(wam)
  gm.boot.coco(10,wam,model="ab,bcde,cdef")



