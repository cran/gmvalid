gm.or <-
function (X, Y, data = 0, conditions = 0, reference = c("last", 
    "biggest", "first"), conf.level = 0.95) 
{
    require(epitools, quietly = TRUE)
    if (any(data != 0) && length(X) == 1 && length(Y) == 1) {
        if (is.array(data)) {
            if (!length(dimnames(data))) 
                for (i in length(dim(data)):1) dimnames(data)[[i]] = as.character(c(1:dim(data)[i]))
            data = data.frame(expand.table(data))
            if (dim(data)[2] <= 1) 
                stop("Dimension missmatch!")
        }
        data = data[, c(X, Y, conditions)]
    }
    else {
        if (length(X) < 2 && length(Y) < 2 || length(X) != length(Y)) 
            stop("The data has false dimensions.")
        if (all(conditions == 0)) {
            data = data.frame(X, Y)
        }
        else {
            data = data.frame(X, Y, conditions)
            conditions = 3:dim(data)[2]
        }
    }
    reference = match.arg(reference)
    if (reference == "last") {
        reference = NULL
        for (i in 1:dim(data)[2]) {
            reference = c(reference, length(levels(as.factor(data[, 
                i]))))
        }
    }
    else if (reference == "biggest") {
        reference = NULL
        for (i in 1:dim(data)[2]) {
            reference = c(reference, which(levels(as.factor(data[, 
                i])) == names(table(data[, i])[which(table(data[, 
                i]) == max(table(data[, i])))])[1]))
        }
    }
    else if (reference == "first") {
        reference = matrix(1, nrow = 1, ncol = dim(data)[2])
    }
    dname = dimnames(data)[[2]]
    data = table(data)
    res_table = array(0, dim = c(2, 2))
    result = matrix(0, nrow = (prod(dim(data)[2:length(dim(data))]) * 
        (dim(data)[1] - 1)), ncol = 4)
    res.index = 1
    dim.names = NULL
    if (all(conditions == 0)) 
        for (i in (1:dim(data)[1])[-reference[1]]) for (j in 1:dim(data)[2]) {
            res_table[1, 1] = data[reference[1], reference[2]]
            res_table[2, 1] = data[i, reference[2]]
            res_table[1, 2] = data[reference[1], j]
            res_table[2, 2] = data[i, j]
            if (any(res_table > 0)) {
                or = oddsratio.wald(res_table, conf.level = conf.level)
                result[res.index, ] = c(or$measure[2, ], or$p.value[2, 
                  1])
                if (j == reference[2]) 
                  result[res.index, 2:4] = NA
            }
            else result[res.index, ] = c(NA, NA, NA, NA)
            res.index = res.index + 1
            if (dim(data)[1] == 2) 
                dim.names = c(dim.names, paste(dname[2], "(", 
                  dimnames(data)[[2]][j], ")", sep = ""))
            else dim.names = c(dim.names, paste(dname[1], "(", 
                dimnames(data)[[1]][i], ")", " ~ ", dname[2], 
                "(", dimnames(data)[[2]][j], ")", sep = ""))
        }
    else {
        index = matrix(1, nrow = 1, ncol = length(conditions))
        looping = length(conditions)
        flag = 0
        next.dim = FALSE
        while (any(index < dim(data)[-c(1, 2)]) || flag == 0) {
            coord = paste(c("", index), collapse = ",")
            cname = NULL
            for (h in 1:length(index)) cname = c(cname, dname[h + 
                2], "(", dimnames(data)[[h + 2]][index[h]], ") ")
            cname = paste(cname, collapse = "")
            data.tmp = eval(parse(text = paste("data[,", coord, 
                "]", sep = "")))
            for (i in (1:dim(data)[1])[-reference[1]]) for (j in 1:dim(data.tmp)[2]) {
                res_table[1, 1] = data.tmp[reference[1], reference[2]]
                res_table[2, 1] = data.tmp[i, reference[2]]
                res_table[1, 2] = data.tmp[reference[1], j]
                res_table[2, 2] = data.tmp[i, j]
                if (any(res_table > 0)) {
                  or = oddsratio.fisher(res_table, conf.level = conf.level)
                  result[res.index, ] = c(or$measure[2, ], or$p.value[2, 
                    1])
                  if (j == reference[2]) 
                    result[res.index, 2:4] = NA
                }
                else result[res.index, ] = c(NA, NA, NA, NA)
                if (dim(data)[1] == 2) 
                  dim.names = c(dim.names, paste(dname[2], "(", 
                    dimnames(data)[[2]][j], ")", " stratified by ", 
                    cname, sep = ""))
                else dim.names = c(dim.names, paste(dname[1], 
                  "(", dimnames(data)[[1]][i], ") ~ ", dname[2], 
                  "(", dimnames(data)[[2]][j], ") stratified by ", 
                  cname, sep = ""))
                res.index = res.index + 1
            }
            if (all(index == dim(data)[-c(1, 2)])) {
                flag = 1
                break
            }
            while (looping > 0 && index[looping] == dim(data)[looping + 
                2]) {
                looping = looping - 1
                next.dim = TRUE
            }
            if (next.dim) {
                index[looping] = index[looping] + 1
                index[(looping + 1):length(index)] = 1
                looping = length(index)
                next.dim = FALSE
            }
            else if (looping > 0) 
                index[looping] = index[looping] + 1
        }
    }
    dimnames(result) = list(dim.names, c("estimate", "lower", 
        "upper", "p.value"))
    if (dim(data)[1] == 2) 
        nnd.name.1 = paste(c("between ", dname[1], "(", dimnames(data)[[1]][reference[1]], 
            ") and"), collapse = "")
    else nnd.name.1 = ""
    if (all(conditions == 0)) 
        nnd.name.2 = paste(c("Marginal Odds Ratio with ", conf.level * 
            100, "% C.I."), collapse = "")
    else nnd.name.2 = paste(c("Stratified Odds Ratio with ", 
        conf.level * 100, "% C.I."), collapse = "")
    names(dimnames(result)) = c(nnd.name.1, nnd.name.2)
    result
}
