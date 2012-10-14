gm.cliquerecursive <-
function (x) 
{
    result = NULL
    if (all(as.logical(x) == upper.tri(x))) 
        result = c(result, ",", dimnames(x)[[1]][1:dim(x)[1]])
    else if (dim(x)[1] > 2) 
        for (i in 1:dim(x)[2]) {
            y = x[, -i]
            y = y[-i, ]
            result = c(result, Recall(y))
        }
    result
}
