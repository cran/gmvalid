gm.chi <-
function (data, X = 1, Y = 2, Z = 3) 
{
    if (is.array(data)) {
        if (!length(dimnames(data))) 
            for (i in length(dim(data)):1) dimnames(data)[[i]] = as.character(c(1:dim(data)[i]))
        data = data.frame(expand.table(data))
        if (dim(data)[2] < 2) 
            stop("Dimension missmatch!")
    }
    data = data[, c(X, Y, Z)]
    data = table(data)
    p <- data
    chi <- 0
    for (k in 1:dim(p)[3]) {
        s_k = sum(data[, , k])
        for (j in 1:dim(p)[2]) {
            s_j = sum(data[, j, k])
            for (i in 1:dim(p)[1]) {
                p[i, j, k] = sum(data[i, , k]) * s_j/s_k
                if (p[i, j, k] > 0) 
                  chi = chi + (p[i, j, k] - data[i, j, k])^2/p[i, 
                    j, k]
            }
        }
    }
    p.value = 1 - pchisq(chi, (dim(p)[1] - 1) * (dim(p)[2] - 
        1) * dim(p)[3])
    df = (dim(p)[1] - 1) * (dim(p)[2] - 1) * dim(p)[3]
    list(chi.squared = chi, DF = df, p.value = p.value)
}
