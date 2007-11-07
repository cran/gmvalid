### Name: gm.mim
### Title: Graphical model selection (MIM)
### Aliases: gm.mim
### Keywords: models multivariate

### ** Examples

  data(wam)
  gm.mim(wam)
  ### giving strategy
  gm.mim(wam,strategy="e")
  
  ### giving minimal and maximal model
  gm.mim(wam,strategy="e",model="a,bc,de,f - abcd,cdef")
  ### giving block structure
  gm.mim(wam,strategy="f",model="a,b,c|abc,de|abcd,ef",chain="abc|de|f",options="BNU")



