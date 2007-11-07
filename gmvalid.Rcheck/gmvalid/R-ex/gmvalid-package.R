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




