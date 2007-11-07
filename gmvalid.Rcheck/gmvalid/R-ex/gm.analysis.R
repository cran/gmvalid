### Name: gm.analysis
### Title: Analyze a data set
### Aliases: gm.analysis
### Keywords: multivariate nonparametric graphs hplot models

### ** Examples

  data(wam)
  gm.analysis(wam)
  
  ### showing various options in action
  gm.analysis(wam,program="c",strategy="f",edge.measure="b",
    boot.N=50,plot.significant=FALSE,recursive=TRUE,follow=TRUE,decomposable.mode=TRUE)
              
  gm.analysis(wam,edge.measure="p",options="u")



