gm.activate.stop <-
function (start.model, data = data) 
{
    MODEL <- gm.make.unique.cliques(start.model)
    statout <- list()
    for (i in 1:length(start.model)) {
        if (length(start.model[[i]]) > 1) 
            statout[[i]] <- gm.gamma(data = data[, start.model[[i]]])
    }
    STAT <- statout
    return(list(measure = STAT, model = MODEL))
}
