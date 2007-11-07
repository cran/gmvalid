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



