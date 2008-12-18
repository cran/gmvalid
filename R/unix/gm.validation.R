`gm.validation` <-
function (data, N = 0, Umax = 0.5, 
    conf.level = 0.95, ...) 
{
    if (is.data.frame(data) || is.array(data)) {
        boot.out <- gm.boot.coco(N, data, ...)
    }
    else boot.out <- data
    result <- list()
    result$"original model" = paste(sort(strsplit(boot.out$"original model", 
        ",")[[1]]), collapse = ",")
    result$"mode model" = boot.out$"bootstrapped models"[1, ]
    edge.square = 0
    edge.stab = 0
    edges = boot.out$"edge frequencies"
    modelmatrix = edges
    for (i in 1:dim(edges)[1]) for (j in 1:dim(edges)[2]) modelmatrix[i, 
        j] = ifelse(edges[i, j] >= Umax, 1, 0)
    m.name = .gm.modelparse(modelmatrix)
    result$"mean model" = m.name
    for (i in 1:(dim(edges)[1] - 1)) for (j in (i + 1):dim(edges)[2]) {
        edge.square = edge.square + min((1 - Umax)^2/0.25 * edges[i, 
            j]^2, Umax^2/0.25 * (1 - edges[i, j])^2)/((1 - Umax)^2/0.25 * 
            Umax^2)
        edge.stab = edge.stab + min((1 - Umax)/0.5 * edges[i, 
            j], Umax/0.5 * (1 - edges[i, j]))/2/(Umax/0.5 * (1 - 
            Umax))
    }
    edge.square = edge.square/sum(1:(dim(edges)[1] - 1))
    edge.stab = edge.stab/sum(1:(dim(edges)[1] - 1))
    result$MEU = edge.stab
    result$MSEU = edge.square
    differ.result = c(NULL)
    elements = boot.out$"variable names"[, 1]
    original = modelmatrix
    N = sum(boot.out$"edge differences"$abs)
    for (k in 1:dim(boot.out$"bootstrapped models")[1]) {
        diff = as.character(sum(abs(original - .gm.matrixparse(boot.out$"bootstrapped models"[k, 
            1]))))
        if (length(which(names(differ.result) == diff)) == 0) {
            tmp = names(differ.result)
            differ.result = c(differ.result, as.numeric(boot.out$"bootstrapped models"[k, 
                2]))
            names(differ.result) = c(tmp, diff)
        }
        else differ.result[which(names(differ.result) == diff)] = differ.result[which(names(differ.result) == 
            diff)] + as.numeric(boot.out$"bootstrapped models"[k, 
            2])
    }
    differ.result = differ.result[order(names(differ.result))]
    result$"edge differences" = differ.result
    var.sum = 0
    expect = 0
    for (i in 1:length(differ.result)) {
        var.sum = var.sum + differ.result[i] * N * as.numeric(names(differ.result[i]))^2
        expect = expect + differ.result[i] * as.numeric(names(differ.result[i]))
    }
    var.sum = as.numeric(var.sum/(N - 1))
    result$"total possible edges" = sum(1:(dim(edges)[1] - 1))
    result$"model std" = sqrt(var.sum)
    result$"MED" = as.numeric(expect)
    differ.result = differ.result * N
    ii = 1
    c.interval = 0
    while (c.interval <= conf.level * sum(differ.result)) {
        c.interval = c.interval + differ.result[ii]
        ii = ii + 1
    }
    result$"bootstrap percentile 95" = as.numeric(names(differ.result[ii - 
        1]))
    result$"variable names" = boot.out$"variable names"
    result
}
