gm.cv.divide <-
function (k, data) 
{
    N <- dim(data)[1]
    ds.names <- paste("ds", 1:k, sep = "")
    DS <- list(length = k)
    OUT <- list()
    mod <- N%%k
    inv <- N%/%k
    sVec <- 1:N
    for (j in 1:(k)) {
        if (j <= mod) {
            assign(paste("ds", j, sep = ""), sample(sVec, inv + 
                1))
            OUT[[j]] <- eval(parse(text = (paste("ds", j, sep = ""))))
            sVec <- setdiff(sVec, OUT[[j]])
        }
        else {
            assign(paste("ds", j, sep = ""), sample(sVec, inv))
            OUT[[j]] <- eval(parse(text = (paste("ds", j, sep = ""))))
            sVec <- setdiff(sVec, OUT[[j]])
        }
    }
    return(lapply(OUT, sort))
}
