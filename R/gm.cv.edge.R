gm.cv.edge <-
function (pvalue, k, data) 
{
    knodes <- dim(data)[2]
    diagMat <- function(x) x[row(x) < col(x)]
    together <- function(x) paste(x[1], x[2], sep = "")
    pv.v <- lapply(pvalue, diagMat)
    pv.m <- matrix(0, nrow = length(pv.v), ncol = length(pv.v[[1]]))
    for (i in 1:length(pv.v)) {
        pv.m[i, ] <- pv.v[[i]]
    }
    tmp <- matrix(1:knodes^2, byrow = TRUE, ncol = knodes)
    tmp.v <- diagMat(tmp)
    OUT <- matrix(NA, ncol = length(tmp.v), nrow = length(pv.v))
    for (i in 1:length(pv.v)) {
        elliott.smith <- matrix(c(tmp.v, pv.v[[i]]), byrow = TRUE, 
            nrow = 2)
        OUT[i, ] <- elliott.smith[2, order(elliott.smith[1, ])]
    }
    dimName <- combinations(knodes, 2, letters)
    dimnames(OUT) <- list(1:k, apply(dimName, 1, together))
    MEAN <- apply(OUT, 2, mean, na.rm = TRUE)
    SE <- apply(OUT, 2, gm.cv.sd, na.rm = TRUE)
    STATISTIC <- list(mean = MEAN, se = SE, pvalue = OUT)
}
