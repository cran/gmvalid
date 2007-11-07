### Name: wynder
### Title: Alcohol, Smoking and Oral Cancer
### Aliases: wynder
### Keywords: datasets

### ** Examples

  data(wynder)
  tab <- table(wynder)
  dimnames(tab) <- list(c("<1","1-6",">6"), c("<15","16-34",">34"), c("controls","cases"))  
  names(dimnames(tab)) <- c("Alcohol (unit/day)","Smoking (cigarettes/day)","Group")  
  tab



