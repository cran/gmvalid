### * <HEADER>
###
attach(NULL, name = "CheckExEnv")
assign("nameEx", 
       local({
	   s <- "__{must remake R-ex/*.R}__"
           function(new) {
               if(!missing(new)) s <<- new else s
           }
       }),
       pos = "CheckExEnv")
## Add some hooks to label plot pages for base and grid graphics
assign("base_plot_hook",
       function() {
           pp <- par(c("mfg","mfcol","oma","mar"))
           if(all(pp$mfg[1:2] == c(1, pp$mfcol[2]))) {
               outer <- (oma4 <- pp$oma[4]) > 0; mar4 <- pp$mar[4]
               mtext(sprintf("help(\"%s\")", nameEx()), side = 4,
                     line = if(outer)max(1, oma4 - 1) else min(1, mar4 - 1),
              outer = outer, adj = 1, cex = .8, col = "orchid", las=3)
           }
       },
       pos = "CheckExEnv")
assign("grid_plot_hook",
       function() {
           pushViewport(viewport(width=unit(1, "npc") - unit(1, "lines"),
                                 x=0, just="left"))
           grid.text(sprintf("help(\"%s\")", nameEx()),
                     x=unit(1, "npc") + unit(0.5, "lines"),
                     y=unit(0.8, "npc"), rot=90,
                     gp=gpar(col="orchid"))
       },
       pos = "CheckExEnv")
setHook("plot.new",     get("base_plot_hook", pos = "CheckExEnv"))
setHook("persp",        get("base_plot_hook", pos = "CheckExEnv"))
setHook("grid.newpage", get("grid_plot_hook", pos = "CheckExEnv"))
assign("cleanEx",
       function(env = .GlobalEnv) {
	   rm(list = ls(envir = env, all.names = TRUE), envir = env)
           RNGkind("default", "default")
	   set.seed(1)
   	   options(warn = 1)
	   .CheckExEnv <- as.environment("CheckExEnv")
	   delayedAssign("T", stop("T used instead of TRUE"),
		  assign.env = .CheckExEnv)
	   delayedAssign("F", stop("F used instead of FALSE"),
		  assign.env = .CheckExEnv)
	   sch <- search()
	   newitems <- sch[! sch %in% .oldSearch]
	   for(item in rev(newitems))
               eval(substitute(detach(item), list(item=item)))
	   missitems <- .oldSearch[! .oldSearch %in% sch]
	   if(length(missitems))
	       warning("items ", paste(missitems, collapse=", "),
		       " have been removed from the search path")
       },
       pos = "CheckExEnv")
assign("ptime", proc.time(), pos = "CheckExEnv")
grDevices::postscript("gmvalid-Ex.ps")
assign("par.postscript", graphics::par(no.readonly = TRUE), pos = "CheckExEnv")
options(contrasts = c(unordered = "contr.treatment", ordered = "contr.poly"), pager="console")
options(warn = 1)    
library('gmvalid')

assign(".oldSearch", search(), pos = 'CheckExEnv')
assign(".oldNS", loadedNamespaces(), pos = 'CheckExEnv')
cleanEx(); nameEx("dp")
### * dp

flush(stderr()); flush(stdout())

### Name: dp
### Title: Death penalty example of Simpson's paradox
### Aliases: dp
### Keywords: datasets

### ** Examples

data(dp)
## Graphical model analysis shows that 'defendants' race' is 
## independent from 'death penalty' given 'victims' race'.
gm.analysis(dp,program="coco",recursive=TRUE)



cleanEx(); nameEx("gm.analysis")
### * gm.analysis

flush(stderr()); flush(stdout())

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



cleanEx(); nameEx("gm.boot.coco")
### * gm.boot.coco

flush(stderr()); flush(stdout())

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



cleanEx(); nameEx("gm.boot.mim")
### * gm.boot.mim

flush(stderr()); flush(stdout())

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



cleanEx(); nameEx("gm.coco")
### * gm.coco

flush(stderr()); flush(stdout())

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



cleanEx(); nameEx("gm.csi")
### * gm.csi

flush(stderr()); flush(stdout())

### Name: gm.csi
### Title: Conditional Synergy Index
### Aliases: gm.csi
### Keywords: htest

### ** Examples

  data(idd35)
  gm.csi(1,2,3,data=idd35)

  ### >> constructing an additive and multiplicative penetrance
  x <- c(0.1,0.4)
  y <- c(0.05,0.5)
  add.pen <- outer(x,y,FUN="+")
  mult.pen <- outer(x,y)
  het.pen <- outer(x,y,FUN="+") - outer(x,y)

  ### >> Function that samples data using the penetrance 
  make.data <- function(R,pen,category) 
    {
      s.vec <- sample(c(1,2,3,4),R,replace=TRUE,prob=as.vector(pen))
      fact.1 <- fact.2 <- vector()
      for( i in 1:R ) {
        ifelse( s.vec[i] == 1 || s.vec[i] == 3 , fact.1[i] <- 1, fact.1[i] <- 2 ) 
        ifelse( s.vec[i] == 1 || s.vec[i] == 2 , fact.2[i] <- 1, fact.2[i] <- 2 ) 
      }
      cbind(X=fact.1,Y=fact.2,group=rep(category,R))  
    }

  ### >>> Building datasets with affected and unaffected subjects   
  add.aff <- make.data(200,add.pen,2)
  add.uaf <- make.data(200,1-add.pen,1)  
  add.df <- as.data.frame(rbind(add.uaf,add.aff))
  
  mult.aff <- make.data(200,mult.pen,2)
  mult.uaf <- make.data(200,1-mult.pen,1)  
  mult.df <- as.data.frame(rbind(mult.uaf,mult.aff))
  
  het.aff <- make.data(200,het.pen,2)
  het.uaf <- make.data(200,1-het.pen,1)  
  het.df <- as.data.frame(rbind(het.uaf,het.aff))
   
  gm.csi(1,2,3,add.df,pen=add.pen)   # Additivity
  gm.csi(1,2,3,mult.df,pen=mult.pen) # Synergy
  gm.csi(1,2,3,het.df,pen=het.pen)   # Antagonism



cleanEx(); nameEx("gm.cv")
### * gm.cv

flush(stderr()); flush(stdout())

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



cleanEx(); nameEx("gm.gamma")
### * gm.gamma

flush(stderr()); flush(stdout())

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



cleanEx(); nameEx("gm.generate")
### * gm.generate

flush(stderr()); flush(stdout())

### Name: gm.generate
### Title: Random data frames of binary variables with given marginals
### Aliases: gm.generate
### Keywords: datagen distribution

### ** Examples

gm.generate(10,c(.5,.2,.2))
gm.generate(15,c(.5,.5,.5,.5,.5,.5))



cleanEx(); nameEx("gm.mim")
### * gm.mim

flush(stderr()); flush(stdout())

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



cleanEx(); nameEx("gm.modelsim")
### * gm.modelsim

flush(stderr()); flush(stdout())

### Name: gm.modelsim
### Title: Random data frames with given dependence model and marginals
### Aliases: gm.modelsim
### Keywords: datagen graphs

### ** Examples

    gm.modelsim(100,"AB,AC")
    table( gm.modelsim(100,"a,b,c") )
    
    tmp.df <- gm.modelsim(10000,"abf,cd,cf,bdeg,bfg")
    
    # with given number of categories
    tmp.df <- gm.modelsim(1000,"AB,C",list(c(1,1,1),c(1,1),c(1,1,1)))

    # with given number of categories and marginals
    tmp.df <- gm.modelsim(1000,"ABC",list(c(0.3,0.3,0.4),c(0.6,.4),c(0.25,0.25,0.5)))
    table(tmp.df)

    ## Not run: 
##D tmp.df <- gm.modelsim(100,"ABC",list(3,2,3))# (number of categories will be 2 x 2 x 2 )
##D             gm.modelsim(100,"123")
##D             
## End(Not run)



cleanEx(); nameEx("gm.or")
### * gm.or

flush(stderr()); flush(stdout())

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




cleanEx(); nameEx("gm.plot")
### * gm.plot

flush(stderr()); flush(stdout())

### Name: gm.plot
### Title: Plot graphical models
### Aliases: gm.plot
### Keywords: hplot graphs

### ** Examples

  gm.plot("ABC,CDE")

  gm.plot("VBA,EVC")
  
  gm.plot(c("ABC,CDE","AB,BC,CD,DE","ABC,DEF,GHI"))
  
  gm.plot("AB,AC",FALSE,matrix(0.5,nrow=3,ncol=3))



cleanEx(); nameEx("gm.si")
### * gm.si

flush(stderr()); flush(stdout())

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




cleanEx(); nameEx("gm.sim.ixj")
### * gm.sim.ixj

flush(stderr()); flush(stdout())

### Name: gm.sim.ixj
### Title: Random (i x j)-way dependency tables with given marginals
### Aliases: gm.sim.ixj
### Keywords: datagen

### ** Examples

    gm.sim.ixj(1000,c(1,1,1),c(1,1,1,1,1))
    gm.sim.ixj(1000,c(.2,.3,.4,.1),c(.5,.2,.3)) 
    
    round(gm.sim.ixj(30,c(1,1),c(1,1)))       
    
    tab <- round(gm.sim.ixj(500,c(.5,.5),c(.5,.5)))
    chisq.test(tab)   



cleanEx(); nameEx("gmvalid-package")
### * gmvalid-package

flush(stderr()); flush(stdout())

### Encoding: latin1

### Name: gmvalid-package
### Title: Validation of graphical models
### Aliases: gmvalid-package gmvalid
### Keywords: datagen models multivariate nonparametric graphs htest

### ** Examples

    ### Generates a data frame given a dependence model
    gm.a <- gm.modelsim(1000,"ABC,CDE")
    
    ### Modelselection with graphical output
    gm.analysis(gm.a)   
    
    ### Model validation using the bootstrap 
    gm.boot.coco(100,gm.a,recursive=TRUE,follow=TRUE)

    ### Model prediction using cross validation
    gm.cv(3,data=gm.a,strategy="f",options="b")
    
    ### Testing interaction on the penetrance scale
    ### using the conditional synergy index (CSI)
    gm.csi(1,2,3,data=gm.a)

    ### Testing interaction on a additivity scale
    ### using the synergy index (S)
    gm.si(1,2,3,data=gm.a)   

    ### Gamma Coefficient B indpendent D given C
    gm.gamma(2,4,data=gm.a,conditions=3)




cleanEx(); nameEx("idd35")
### * idd35

flush(stderr()); flush(stdout())

### Name: idd35
### Title: Type 1 Diabetes susceptibility loci Idd3 and Idd5
### Aliases: idd35
### Keywords: datasets

### ** Examples

  data(idd35)
  table(idd35)



cleanEx(); nameEx("wam")
### * wam

flush(stderr()); flush(stdout())

### Name: wam
### Title: Women and Mathematics
### Aliases: wam
### Keywords: datasets

### ** Examples

  data(wam)
  gm.analysis(wam, program="coco")



cleanEx(); nameEx("wynder")
### * wynder

flush(stderr()); flush(stdout())

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



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
