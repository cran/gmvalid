### Name: gm.coco
### Title: Graphical model selection (CoCo)
### Aliases: gm.coco
### Keywords: models multivariate

### ** Examples

  data(wam)
  gm.coco(wam)
  ### giving many options to the strategy
  gm.coco(wam,recursive=TRUE,follow=TRUE,decomposable.mode=TRUE,
            coherent=TRUE,IC=TRUE,BIC=TRUE)
  
  ### giving base
  gm.coco(wam,strategy="e",model=c("ab,cd","ae,be"))
  gm.coco(wam,strategy="f",model="abc,cd,de,f")



