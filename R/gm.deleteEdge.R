gm.deleteEdge <-
function (data, model, conf.level = 0.95, output = TRUE) 
{
    dimName <- combinations(length(model), 2, v = model)
    result <- gm.gamma(data = data[, model], conf.level = conf.level)
    nrows <- seq(1, dim(result)[1])
    kante.weg <- nrows[result[, 5] > 1 - conf.level]
    if (max(result[, 5], na.rm = TRUE) <= 1 - conf.level) 
        kante.weg <- vector()
    result.out <- matrix(0, nrow = nrow(result), ncol = (ncol(result) + 
        1))
    dimnames(result.out) <- list(dimnames(result)[[1]], c(dimnames(result)[[2]], 
        "edge.out"))
    result.out[1:nrow(result), 1:ncol(result)] <- result
    result.out[kante.weg, (ncol(result) + 1)] <- 1
    if (output) 
        print(result.out)
    ifelse(length(kante.weg) > 0, edge.out <- kante.weg, edge.out <- NA)
    ifelse(length(kante.weg) > 0, edge.name <- dimName[kante.weg, 
        ], edge.name <- NA)
    ifelse(length(kante.weg) > 0, p.value <- result[kante.weg, 
        5], p.value <- NA)
    OUT <- list(edge.out = edge.out, edge.name = edge.name, p.value = p.value, 
        result = result.out)
    return(OUT)
}
