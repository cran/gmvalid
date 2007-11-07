### Name: gm.cv
### Title: Cross Validation for Graphical (Chain) Models
### Aliases: gm.cv
### Keywords: models regression

### ** Examples

  
  ABC <- gm.modelsim(500,"ABC,CD")
  out <- gm.cv(5,data=ABC, strategy="f")
  out
  
  ### DAG using a stepwise selection
  out.dag <- gm.cv(3,data=ABC,option="j",chain="d|b|c|a")  
  
  ### Chain graph using BIC as selection criteria and allowing for 
  ### non-decomposable models
  cg <- gm.modelsim(1000,"ABD,BCE")  
  out.cg <- gm.cv(3,data=cg,option="bu",chain="cb|de|a")
  
  ## Not run: 
##D gm.cv(3,data=ABC,chain="DBD|A") # you have to use lowercase letters
##D             gm.cv(3,data=ABC,chain="dca|b") # a is supposed to be outcome variable 
##D                                   # and thus have to be in the very right block    
##D             
## End(Not run)



