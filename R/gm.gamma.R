`gm.gamma` <-
function (X = 0, Y = 0, data = 0, conditions = 0, type = c("conditional", 
    "single", "marginal"), conf.level = 0.95) 
{
    if (all(X == 0 && Y == 0 && data == 0)) 
        stop("Data, please!")
    result = NULL
    dname = NULL
    type = match.arg(type)
    if (X != 0 && Y != 0) 
        type = "single"
    if (type == "single" && all(X == 0 || Y == 0)) 
        stop("For a single calculation you have to give both X and Y!")
    require(epitools, quietly = TRUE)
    if (any(data != 0) && length(X) == 1 && length(Y) == 1) {
        if (is.array(data)) {
            if (!length(dimnames(data))) 
                for (i in length(dim(data)):1) dimnames(data)[[i]] = as.character(c(1:dim(data)[i]))
            data = data.frame(expand.table(data))
            if (dim(data)[2] < 2) 
                stop("Dimension missmatch!")
        }
    }
    if (length(X) > 1 && length(Y) > 1 && length(X) == length(Y)) {
        if (all(conditions == 0)) {
            data = data.frame(X, Y)
        }
        else {
            data = data.frame(X, Y, conditions)
            conditions = 3:dim(data)[2]
        }
        X = 1
        Y = 2
    }
    if (type == "conditional") 
        for (i in 1:(dim(data)[2] - 1)) for (j in (i + 1):dim(data)[2]) {
            partial.g = .gm.gamma.single(data, i, j, c(1:dim(data)[2])[-c(i, 
                j)], conf.level)
            result = c(result, partial.g)
            dim(result) = c(5, length(result)/5)
            dname = c(dname, dimnames(partial.g)[[1]])
        }
    else if (type == "marginal") 
        for (i in 1:(dim(data)[2] - 1)) for (j in (i + 1):dim(data)[2]) {
            partial.g = .gm.gamma.single(data, i, j, 0, conf.level)
            result = c(result, partial.g)
            dim(result) = c(5, length(result)/5)
            dname = c(dname, dimnames(partial.g)[[1]])
        }
    else if (type == "single") {
        result = t(.gm.gamma.single(data, X, Y, conditions, conf.level))
        dname = dimnames(result)[[2]]
    }
    result = t(result)
    dimnames(result) = list(dname, c("estimate", "SE", "lower", 
        "upper", "p.value"))
    if (type == "conditional" || any(conditions != 0)) 
        names(dimnames(result)) = c("Hypothesis", paste(c("Conditional Gamma with ", 
            conf.level * 100, "% C.I."), collapse = ""))
    else names(dimnames(result)) = c("Hypothesis", paste(c("Marginal Gamma with ", 
        conf.level * 100, "% C.I."), collapse = ""))
    result
}
