gm.csi <-
function (X, Y, group, data = 0, reference = c(1, 1, 2), pen = NULL, 
    conf.level = 0.95) 
{
    require(gtools, quietly = TRUE)
    if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || 
        conf.level < 0 || conf.level > 1)) 
        stop("'conf.level' must be a single number between 0 and 1")
    if (length(reference) < 3) 
        stop("Please give reference categories for all 3 variables.")
    if (any(data != 0) && length(X) == 1 && length(Y) == 1) {
        if (is.array(data)) {
            if (!length(dimnames(data))) 
                for (i in length(dim(data)):1) dimnames(data)[[i]] = as.character(c(1:dim(data)[i]))
            data = data.frame(expand.table(data))
            if (dim(data)[2] < 3) 
                stop("Dimension missmatch!")
        }
        data = data[, c(X, Y, group)]
    }
    else {
        if (length(X) < 2 && length(Y) < 2 && length(group) < 
            2 || length(X) != length(Y) || length(group) != length(X)) 
            stop("The data has false dimensions.")
        data = data.frame(X, Y, group)
    }
    dname = dimnames(data)[[2]]
    csi <- vector()
    var.csi <- vector()
    tab <- table(data)
    if (dim(tab)[3] > 2) 
        stop("Group variable has more than 2 dimensions.")
    A <- as.vector(table(data[, 1]))
    B <- as.vector(table(data[, 2]))
    AB <- as.vector(table(data[, c(1, 2)]))
    ifelse(is.null(pen), penetrance <- FALSE, penetrance <- TRUE)
    if (penetrance == FALSE) {
        pen <- tab[, , reference[3]]/(tab[, , 1] + tab[, , 2])
        if (any(is.na(pen) == TRUE)) {
            penetrance[pen == "NaN"] <- 0
        }
    }
    dim.i <- dim(pen)[1]
    dim.j <- dim(pen)[2]
    z <- qnorm(0.5 + 0.5 * conf.level)
    ifelse(dim.i > 2, combi.i <- combinations(dim.i - 1, dim.i - 
        2) + 1, combi.i <- FALSE)
    ifelse(dim.j > 2, combi.j <- combinations(dim.i - 1, dim.i - 
        2) + 1, combi.j <- FALSE)
    index <- 0
    for (k in 1:(dim.i - 1)) {
        ifelse(combi.i == FALSE, i <- "", i <- combi.i[k, ])
        for (l in 1:(dim.j - 1)) {
            index <- index + 1
            ifelse(combi.j == FALSE, j <- "", j <- combi.i[l, 
                ])
            ifelse(dim.i > 2, fA <- A[-i], fA <- A)
            pA <- fA/sum(A)
            ifelse(dim.j > 2, fB <- B[-j], fB <- B)
            pB <- fB/sum(B)
            ifelse(dim.j > 2, fAB <- AB[-j], fAB <- AB)
            pAB <- fAB/sum(AB)
            if (is.numeric(i) & is.numeric(j)) {
                coord <- paste(c(i, j), collapse = ",")
                coord.pr <- paste(c(-i, -j), ",", sep = "", collapse = "")
            }
            else if (is.numeric(i) & !is.numeric(j)) {
                coord <- paste(c(i, ""), collapse = ",")
                coord.pr <- paste(c(-i, ""), ",", sep = "", collapse = "")
            }
            else if (!is.numeric(i) & is.numeric(j)) {
                coord <- paste(c("", j), collapse = ",")
                coord.pr <- paste(c("", -j), ",", sep = "", collapse = "")
            }
            else if (!is.numeric(i) & !is.numeric(j)) {
                coord <- paste(c("", ""), collapse = ",")
                coord.pr <- paste(c("", ""), ",", sep = "", collapse = "")
            }
            D11 <- eval(parse(text = paste("gm.delta2(1,1,", 
                coord, ",pen*sum(tab),A,B,AB)", sep = "")))/eval(parse(text = paste("gm.delta2(1,1,", 
                coord, ",(1-pen)*sum(tab),A,B,AB)", sep = "")))
            D12 <- eval(parse(text = paste("gm.delta2(1,2,", 
                coord, ",pen*sum(tab),A,B,AB)", sep = "")))/eval(parse(text = paste("gm.delta2(1,2,", 
                coord, ",(1-pen)*sum(tab),A,B,AB)", sep = "")))
            D21 <- eval(parse(text = paste("gm.delta2(2,1,", 
                coord, ",pen*sum(tab),A,B,AB)", sep = "")))/eval(parse(text = paste("gm.delta2(2,1,", 
                coord, ",(1-pen)*sum(tab),A,B,AB)", sep = "")))
            D22 <- eval(parse(text = paste("gm.delta2(2,2,", 
                coord, ",pen*sum(tab),A,B,AB)", sep = "")))/eval(parse(text = paste("gm.delta2(2,2,", 
                coord, ",(1-pen)*sum(tab),A,B,AB)", sep = "")))
            Delta <- matrix(c(D11, D21, D12, D22), nrow = 2)
            d <- matrix(c(eval(parse(text = paste("gm.delta2(1,1,", 
                coord, ",1-pen,pA,pB,pAB)", sep = ""))), eval(parse(text = paste("gm.delta2(2,1,", 
                coord, ",1-pen,pA,pB,pAB)", sep = ""))), eval(parse(text = paste("gm.delta2(1,2,", 
                coord, ",1-pen,pA,pB,pAB)", sep = ""))), eval(parse(text = paste("gm.delta2(2,2,", 
                coord, ",1-pen,pA,pB,pAB)", sep = "")))), nrow = 2)
            csi[index] <- (D11 + D22)/(D12 + D21)
            pr <- eval(parse(text = paste("tab[", coord.pr, "]/sum(tab[", 
                coord.pr, "])", sep = "")))
            c1.ar <- array(0, dim = c(2, 2, 2))
            c1.ar[1, 1, ] <- gm.var.csi.part1(1, 1, Delta, d, 
                pr)
            c1.ar[2, 2, ] <- gm.var.csi.part1(2, 2, Delta, d, 
                pr)
            c2.ar <- array(0, dim = c(2, 2, 2))
            c2.ar[1, 1, 1] <- gm.var.csi.part2(1, 1, 1, Delta, 
                d, pr)
            c2.ar[2, 2, 1] <- gm.var.csi.part2(2, 2, 1, Delta, 
                d, pr)
            c2.ar[1, 1, 2] <- gm.var.csi.part2(1, 1, 2, Delta, 
                d, pr)
            c2.ar[2, 2, 2] <- gm.var.csi.part2(2, 2, 2, Delta, 
                d, pr)
            c3.ar <- array(0, dim = c(2, 2, 2))
            c3.ar[1, 1, 1] <- gm.var.csi.part3(1, 1, 1, Delta, 
                d, pr)
            c3.ar[2, 2, 1] <- gm.var.csi.part3(2, 2, 1, Delta, 
                d, pr)
            c3.ar[1, 1, 2] <- gm.var.csi.part3(1, 1, 2, Delta, 
                d, pr)
            c3.ar[2, 2, 2] <- gm.var.csi.part3(2, 2, 2, Delta, 
                d, pr)
            concor <- c1.ar + c2.ar - csi[index] * c3.ar
            d1.ar <- array(0, dim = c(2, 2, 2))
            d1.ar[1, 2, 1] <- gm.var.csi.part3(1, 2, 1, Delta, 
                d, pr)
            d1.ar[2, 1, 1] <- gm.var.csi.part3(2, 1, 1, Delta, 
                d, pr)
            d1.ar[1, 2, 2] <- gm.var.csi.part3(1, 2, 2, Delta, 
                d, pr)
            d1.ar[2, 1, 2] <- gm.var.csi.part3(2, 1, 2, Delta, 
                d, pr)
            d2.ar <- array(0, dim = c(2, 2, 2))
            d2.ar[1, 2, 1] <- gm.var.csi.part4(1, 2, 1, Delta, 
                d, pr)
            d2.ar[2, 1, 1] <- gm.var.csi.part4(2, 1, 1, Delta, 
                d, pr)
            d2.ar[1, 2, 2] <- gm.var.csi.part4(1, 2, 2, Delta, 
                d, pr)
            d2.ar[2, 1, 2] <- gm.var.csi.part4(2, 1, 2, Delta, 
                d, pr)
            d3.ar <- array(0, dim = c(2, 2, 2))
            d3.ar[1, 2, 1] <- gm.var.csi.part2(1, 2, 1, Delta, 
                d, pr)
            d3.ar[2, 1, 1] <- gm.var.csi.part2(2, 1, 1, Delta, 
                d, pr)
            d3.ar[1, 2, 2] <- gm.var.csi.part2(1, 2, 2, Delta, 
                d, pr)
            d3.ar[2, 1, 2] <- gm.var.csi.part2(2, 1, 2, Delta, 
                d, pr)
            discor <- d1.ar - csi[index] * (d2.ar + d3.ar)
            var.csi[index] <- sum(pr * (1 - pr) * (concor + discor)^2)/(sum(tab) * 
                (Delta[1, 2] + Delta[2, 1])^2)
        }
    }
    nv <- (dim.i - 1) * (dim.j - 1)
    ci.mat <- matrix(nrow = nv, ncol = 6)
    for (lauf in 1:nv) {
        CI.lower <- csi[lauf] - z * sqrt(var.csi)[lauf]
        CI.upper <- csi[lauf] + z * sqrt(var.csi)[lauf]
        pvalue <- 2 * pnorm(-abs(csi[lauf] - 1), sd = sqrt(var.csi))
        ci.mat[lauf, ] <- c(csi[lauf], var.csi[lauf], sqrt(var.csi)[lauf], 
            CI.lower, CI.upper, pvalue)
    }
    row.name <- vector()
    for (i in seq(1, (dim.i - 1))) {
        for (j in seq(1, (dim.j - 1))) {
            row.name <- rbind(row.name, c(i, j) + 1)
        }
    }
    dim.names = NULL
    for (i in 1:dim(row.name)[1]) dim.names = c(dim.names, paste(dname[1], 
        "(", row.name[i, 1], ") ~ ", dname[2], "(", row.name[i, 
            2], ")", sep = ""))
    dimnames(ci.mat) = list(dim.names, c("estimate", "VAR", "SE", 
        "lower", "upper", "p.value"))
    names(dimnames(ci.mat)) = c(paste(dname[1], "(", reference[1], 
        ") ~ ", dname[2], "(", reference[2], ") grouped by ", 
        dname[3], "(", reference[3], ")", sep = ""), paste(c("Conditional Synergy Index with ", 
        conf.level * 100, "% C.I."), collapse = ""))
    list(penetrance = pen, measure = ci.mat)
}
