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



