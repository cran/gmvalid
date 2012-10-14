gm.analysis <-
function (data, edge.measure = c("gamma.cond", "marg.gamma", 
    "cor", "boot", "cov", "p.value"), program = c("gms", "coco", 
    "mim"), strategy = c("backwards", "forwards", "eh"), plot.significant = TRUE, 
    boot.N = 100, ...) 
{
    if (is.array(data)) {
        require(epitools, quietly = TRUE)
        if (!length(dimnames(data))) 
            for (i in length(dim(data)):1) dimnames(data)[[i]] = as.character(c(1:dim(data)[i]))
        data = data.frame(expand.table(data))
    }
    var.names = cbind(letters[1:dim(data)[2]], dimnames(data)[2][[1]])
    dimnames(data)[2][[1]] = letters[1:dim(data)[2]]
    program = match.arg(program)
    edge.measure = match.arg(edge.measure)
    strategy = match.arg(strategy)
    if (program == "gms") 
        model = gm.gms(data, strategy, ...)$model
    else if (program == "coco") {
        require(CoCo)
        selection.output = gm.coco(data, strategy, ...)
    }
    else if (program == "mim") {
        require(mimR)
        selection.output = gm.mim(data, strategy, ...)
    }
    model = selection.output$accepted
    analysis = matrix(0, nrow = dim(data)[2], ncol = dim(data)[2])
    if (edge.measure == "gamma.cond") {
        elements = dimnames(data)[[2]]
        dimnames(analysis) = list(elements, elements)
        index = 1
        gamma.ana = gm.gamma(data = data)
        for (i in 1:(dim(data)[2] - 1)) for (j in (i + 1):dim(data)[2]) {
            analysis[i, j] = gamma.ana[index, 1]
            index = index + 1
        }
    }
    else if (edge.measure == "marg.gamma") {
        elements = dimnames(data)[[2]]
        dimnames(analysis) = list(elements, elements)
        index = 1
        gamma.ana = gm.gamma(data = data, type = "marginal")
        for (i in 1:(dim(data)[2] - 1)) for (j in (i + 1):dim(data)[2]) {
            analysis[i, j] = gamma.ana[index, 1]
            index = index + 1
        }
    }
    else if (edge.measure == "cor") 
        analysis = cor(data)
    else if (edge.measure == "cov") 
        analysis = cov(data)
    else if (edge.measure == "boot") {
        if (program == "coco") 
            analysis = gm.boot.coco(boot.N, data, strategy, "edge", 
                ...)$"bootstrapped edges"
        else if (program == "mim") 
            analysis = gm.boot.mim(boot.N, data, strategy, "edge", 
                ...)$"bootstrapped edges"
    }
    else if (edge.measure == "p.value" && program == "mim") {
        analysis = selection.output$"p values"
        if (plot.significant == FALSE) {
            plot.significant = TRUE
            warning("p.values only available for significant edges.")
        }
    }
    else if (edge.measure == "p.value" && program == "coco") 
        warning("CoCo does not support the output of p.values.")
    gm.plot(model, data.analysis = analysis, significant = plot.significant)
    list(strategy = strategy, model = model, edge.measure = edge.measure, 
        analysis = analysis, `variable names` = var.names)
}
