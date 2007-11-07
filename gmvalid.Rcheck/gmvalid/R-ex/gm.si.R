### Name: gm.si
### Title: Synergy Index
### Aliases: gm.si
### Keywords: htest

### ** Examples

  data(wynder)
  gm.si(1,2,3,wynder)

  # Smoking and alcohol in relation to oral cancer among male veterans under age 60.
  # (from "Modern Epidemiology")
  oral <- array(c(20,3,18,8,12,6,166,225),dim=c(2,2,2), 
            dimnames=list(Group=c("control","cases"),
            Smoker=c("no","yes"),Alcohol=c("no","yes")))
  oral.df <- expand.table(oral)
  # grouping variable is first in data frame
  gm.si(2,3,1,oral.df)
  
  # Effects must be ascending in respect to the reference category
  show.effect <- array(c(1,7,2,7,7,12,106,48),dim=c(2,2,2),
                        dimnames=list(A=1:2,B=1:2,C=1:2))
  # produces NaN
  gm.si(1,2,3,expand.table(show.effect))
  # > re-ordering variable B helps
  gm.si(1,2,3,expand.table(show.effect),reference=c(1,2,2))




