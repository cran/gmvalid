gm.sim.ixj <-
function (N, pa, pb) 
{
    ifelse(length(pa) >= length(pb), base <- pa, base <- pb)
    ifelse(length(pa) >= length(pb), pr <- pb, pr <- pa)
    m <- matrix(nrow = length(base), ncol = length(pr), byrow = TRUE)
    for (i in length(dim(m)):1) dimnames(m)[[i]] = as.character(c(1:dim(m)[i]))
    flag <- 1
    while (flag == 1) {
        for (i in 1:length(base)) m[i, ] <- gm.sim.row.i(base[i], 
            pr)
        ifelse(any(is.na(m)) == TRUE, flag <- 1, flag <- 0)
        if (flag == 0) {
            test = expand.table(N * m)
            ifelse(abs(cor(as.numeric(test[, 1]), as.numeric(test[, 
                2]))) < 0.5, flag <- 1, flag <- 0)
        }
    }
    if (length(pa) >= length(pb)) 
        result = N * m
    else {
        result = N * t(m)
    }
    result
}
