gm.addEdge <-
function (data, start.model, clique, conf.level = 0.95, output = TRUE) 
{
    model <- start.model[[clique]]
    if (length(model) > 1 && clique < length(start.model)) {
        bvars = NULL
        for (i in (clique + 1):length(start.model)) bvars <- c(bvars, 
            start.model[[i]])
        for (i in 1:length(model)) if (length(which(bvars == 
            model[i])) > 0) 
            bvars = bvars[-which(bvars == model[i])]
        bvars = unique(bvars)
        result <- matrix(nrow = length(bvars) * length(model), 
            ncol = 7)
        index = 1
        for (i in 1:length(bvars)) for (j in 1:length(model)) {
            r <- gm.gamma(X = j, Y = i, data = data, conditions = model[-j], 
                conf.level = conf.level)
            r <- cbind(r, model[j], bvars[i])
            result[index, ] <- r
            index = index + 1
        }
        nrows <- seq(1, dim(result)[1])
        kante.hin <- nrows[result[, 5] <= 1 - conf.level]
        result.out <- matrix(0, nrow = nrow(result), ncol = (ncol(result) + 
            1))
        result.out[1:nrow(result), 1:ncol(result)] <- result
        result.out[kante.hin, (ncol(result) + 1)] <- 1
        if (output) 
            print(result.out)
        ifelse(length(kante.hin) > 0, zeile <- kante.hin, zeile <- NA)
        ifelse(length(kante.hin) > 0, edge.name <- result.out[kante.hin, 
            6:7], edge.name <- NA)
        ifelse(length(kante.hin) > 0, p.value <- result.out[kante.hin, 
            5], p.value <- NA)
        OUT <- list(zeile = zeile, edge.name = edge.name, p.value = p.value, 
            result = result.out)
    }
    else if (length(model) == 1 && clique < length(start.model)) {
        bvars = NULL
        for (i in (clique + 1):length(start.model)) bvars <- c(bvars, 
            start.model[[i]])
        bvars = unique(bvars)
        result <- matrix(nrow = length(bvars), ncol = 5)
        for (i in 1:length(bvars)) {
            result[i, ] <- gm.gamma(X = data[, model], Y = data[, 
                bvars[i]], conf.level = conf.level)
        }
        result <- cbind(result, rep(model, dim(result)[1]), bvars)
        dimnames(result) <- list(1:dim(result)[1], c("estimate", 
            "SE", "lower", "upper", "p.value", "knot1", "knot2"))
        nrows <- seq(1, dim(result)[1])
        kante.hin <- nrows[result[, 5] <= 1 - conf.level]
        result.out <- matrix(0, nrow = nrow(result), ncol = (ncol(result) + 
            1))
        result.out[1:nrow(result), 1:ncol(result)] <- result
        result.out[kante.hin, (ncol(result) + 1)] <- 1
        if (output) 
            print(result.out)
        ifelse(length(kante.hin) > 0, zeile <- kante.hin, zeile <- NA)
        ifelse(length(kante.hin) > 0, edge.name <- result.out[kante.hin, 
            6:7], edge.name <- NA)
        ifelse(length(kante.hin) > 0, p.value <- result.out[kante.hin, 
            5], p.value <- NA)
        OUT <- list(zeile = zeile, edge.name = edge.name, p.value = p.value, 
            result = result.out)
    }
    else {
        OUT <- list(zeile = NA, edge.name = NA, p.value = NA, 
            result = NA)
    }
    OUT
}
