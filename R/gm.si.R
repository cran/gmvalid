`gm.si` <-
function (X, Y, group, data = 0, reference = c(1, 1, 2), conf.level = 0.95) 
{
    require(epitools, quietly = TRUE)
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
    S.ij <- function(i, j) (OR[i, j] - 1)/(OR[i, 1] + OR[1, j] - 
        2)
    S.co <- function(i, j) (RR[i, j] - 1)/(RR[i, 1] + RR[1, j] - 
        2)
    variance.OR <- function(i, j, tab) {
        bd <- 1/tab[1, 1, 1] + 1/tab[1, 1, 2]
        var.OR[i, j] <- OR[i, j]^2 * (bd + 1/tab[i, j, 1] + 1/tab[i, 
            j, 2])
    }
    cov.1 <- function(i, j, k, l) if (i == k & j == l) 
        var.OR[i, j]
    else OR[i, j] * OR[k, l] * bd
    cov.2 <- function(i, j, k, l) OR[i, j] * (OR[k, 1] + OR[1, 
        l]) * bd
    cov.3 <- function(i, j, k, l) OR[k, l] * (OR[i, 1] + OR[1, 
        j]) * bd
    cov.4 <- function(i, j, k, l) cov.1(1, j, 1, l) + cov.1(1, 
        j, k, 1) + cov.1(i, 1, 1, l) + cov.1(i, 1, k, 1)
    se.S <- function(i, j) {
        sqrt((var.RR[i, j] + var.RR[1, 1])/(RR[i, j] - RR[1, 
            1])^2 + (var.RR[1, j] + var.RR[i, 1] + 4 * var.RR[1, 
            1])/(RR[1, j] + RR[i, 1] - 2 * RR[1, 1])^2 - 4 * 
            var.RR[1, 1]/((RR[i, j] - RR[1, 1]) * (RR[i, 1] + 
            RR[1, j] - 2 * RR[1, 1])))
    }
    se.cc.S <- function(i, j) {
        sqrt(var.OR[i, j]/(OR[i, j] - 1)^2 + (var.OR[1, j] + 
            var.OR[i, 1] + 2 * cov.1(1, 2, 2, 1))/(OR[1, j] + 
            OR[i, 1] - 2)^2 - 2 * cov.2(2, 2, 2, 2)/((OR[i, j] - 
            1) * (OR[1, j] + OR[i, 1] - 2)))
    }
    cov.cc <- function(i, j, k, l) {
        cov.1(i, j, k, l)/((OR[i, 1] + OR[1, j] - 2) * (OR[k, 
            1] + OR[1, l] - 2)) - ((OR[i, j] - 1) * cov.3(i, 
            j, k, l))/((OR[i, 1] + OR[1, j] - 2)^2 * (OR[k, 1] + 
            OR[1, l] - 2)) - ((OR[k, l] - 1) * cov.2(i, j, k, 
            l))/((OR[i, 1] + OR[1, j] - 2) * (OR[k, 1] + OR[1, 
            l] - 2)^2) + ((OR[i, j] - 1) * (OR[k, l] - 1) * cov.4(i, 
            j, k, l))/((OR[i, 1] + OR[1, j] - 2)^2 * (OR[k, 1] + 
            OR[1, l] - 2)^2)
    }
    cov.cc.2 <- function(i, j, k, l) {
        var.OR[i, j]/(OR[i, 1] + OR[1, j] - 2)^2 - 2 * ((OR[i, 
            j] - 1) * cov.2(i, j, k, l))/(OR[i, 1] + OR[1, j] - 
            2)^3 + ((OR[i, j] - 1)^2 * (var.OR[1, j] + var.OR[i, 
            1] + 2 * cov.1(1, j, i, 1)))/(OR[i, 1] + OR[1, j] - 
            2)^4
    }
    tab <- table(data)
    if (dim(tab)[3] > 2) 
        stop("Group variable has more than 2 dimensions.")
    bd <- 1/tab[reference[1], reference[2], 1] + 1/tab[reference[1], 
        reference[2], 2]
    R <- tab[, , reference[3]]/table(data[-group])
    RR <- R/R[reference[1], reference[2]]
    M <- tab[, , 1] + tab[, , 2]
    O <- tab[, , reference[3]]/tab[, , c(1, 2)[-reference[3]]]
    OR <- O/O[reference[1], reference[2]]
    dim.i <- dim(tab)[1]
    dim.j <- dim(tab)[2]
    var.RR <- RR/M
    var.OR <- matrix(nrow = dim.i, ncol = dim.j)
    for (i in 1:dim.i) {
        for (j in 1:dim.j) {
            var.OR[i, j] <- variance.OR(i, j, tab)
        }
    }
    cov.vec <- vector()
    index <- 0
    for (i in 2:dim.i) {
        for (j in 2:dim.j) {
            for (k in 2:dim.i) {
                for (l in 2:dim.j) {
                  index <- index + 1
                  if (i == k && j == l) {
                    cov.vec[index] <- cov.cc.2(i, j, k, l)
                  }
                  else {
                    cov.vec[index] <- cov.cc(i, j, k, l)
                  }
                }
            }
        }
    }
    cov.mat <- matrix(cov.vec, byrow = TRUE, nrow = (dim.i - 
        1) * (dim.j - 1), ncol = (dim.i - 1) * (dim.j - 1))
    dimnames(cov.mat) <- list(seq(1, ((dim.i - 1) * (dim.j - 
        1))), seq(1:((dim.i - 1) * (dim.j - 1))))
    S <- vector(length = (dim.i - 1) * (dim.j - 1))
    S.cohort <- vector(length = (dim.i - 1) * (dim.j - 1))
    k <- 0
    for (i in 2:dim.i) {
        for (j in 2:dim.j) {
            k <- k + 1
            S[k] <- S.ij(i, j)
            S.cohort[k] <- S.co(i, j)
        }
    }
    overall.S <- sum(S^3/diag(cov.mat), na.rm = TRUE)/sum(S^2/diag(cov.mat), 
        na.rm = TRUE)
    overall.S.cohort <- sum(S.cohort^3/diag(cov.mat), na.rm = TRUE)/sum(S.cohort^2/diag(cov.mat), 
        na.rm = TRUE)
    w <- S^2/diag(cov.mat)
    var.overall.S <- sum(outer(w, w) * cov.mat)/sum(outer(w, 
        w))
    se.overall.S <- sqrt(var.overall.S/overall.S^2)
    alpha <- 1 - conf.level
    z <- qnorm(1 - alpha/2)
    pvalue <- ci.low.vec <- ci.upp.vec <- vector()
    combinations <- merge(seq(2, dim.i), seq(2, dim.j))
    for (elliott.smith in 1:((dim.i - 1) * (dim.j - 1))) {
        ci.low.vec[elliott.smith] <- exp(log(S[elliott.smith]) - 
            z * sqrt(diag(cov.mat)[elliott.smith])/S[elliott.smith])
        ci.upp.vec[elliott.smith] <- exp(log(S[elliott.smith]) + 
            z * sqrt(diag(cov.mat)[elliott.smith])/S[elliott.smith])
        pvalue[elliott.smith] <- plnorm(exp(abs(-log(S[elliott.smith]))), 
            sd = sqrt(diag(cov.mat)[elliott.smith])/S[elliott.smith], 
            lower.tail = FALSE)
    }
    out.matrix.kxl <- matrix(nrow = ((dim.i - 1) * (dim.j - 1)), 
        ncol = 5)
    out.matrix.kxl[, 1] <- S
    out.matrix.kxl[, 2] <- sqrt(diag(cov.mat))/S
    out.matrix.kxl[, 3] <- ci.low.vec
    out.matrix.kxl[, 4] <- ci.upp.vec
    out.matrix.kxl[, 5] <- pvalue
    if (dim.i == 2 && dim.j == 2) {
        CI.lower <- exp(log(S) - z * se.S(2, 2))
        CI.upper <- exp(log(S) + z * se.S(2, 2))
        out.matrix.kxl <- rbind(out.matrix.kxl, c(S.cohort, se.S(2, 
            2), CI.lower, CI.upper, plnorm(exp(abs(-log(S))), 
            sd = se.S(2, 2), lower.tail = FALSE)))
        dimnames(out.matrix.kxl) <- list(c("Case/Control", "Cohort"), 
            c("estimates", "SE", "lower", "upper", "p.value"))
        names(dimnames(out.matrix.kxl)) = c(paste(dname[1], "(", 
            reference[1], ") ~ ", dname[2], "(", reference[2], 
            ") grouped by ", dname[3], "(", reference[3], ")", 
            sep = ""), paste(c("Synergy Index with ", conf.level * 
            100, "% C.I."), collapse = ""))
        OUT <- list(OddsRatio = OR, measure = out.matrix.kxl)
    }
    else {
        CI.cc.lower <- exp(log(overall.S) - z * se.overall.S)
        CI.cc.upper <- exp(log(overall.S) + z * se.overall.S)
        pvalue[length(pvalue) + 1] <- plnorm(overall.S, sd = se.overall.S)
        row.name <- vector()
        for (i in seq(1, (dim.i - 1))) {
            for (j in seq(1, (dim.j - 1))) {
                row.name <- rbind(row.name, c(i, j) + 1)
            }
        }
        dim.names = NULL
        for (i in 1:dim(row.name)[1]) dim.names = c(dim.names, 
            paste(dname[1], "(", row.name[i, 1], ") ~ ", dname[2], 
                "(", row.name[i, 2], ")", sep = ""))
        out.matrix.kxl <- rbind(out.matrix.kxl, c(overall.S, 
            se.overall.S, CI.cc.lower, CI.cc.upper, pvalue[length(pvalue)]))
        dim.names = c(dim.names, "overall")
        dimnames(out.matrix.kxl) <- list(dim.names, c("estimates", 
            "SE", "lower", "upper", "p.value"))
        names(dimnames(out.matrix.kxl)) = c(paste(dname[1], "(", 
            reference[1], ") ~ ", dname[2], "(", reference[2], 
            ") grouped by ", dname[3], "(", reference[3], ")", 
            sep = ""), paste(c("Synergy Index with ", conf.level * 
            100, "% C.I. (Case-Control design)"), collapse = ""))
        OUT <- list(OddsRatio = OR, covariance = cov.mat, measure = out.matrix.kxl)
    }
    OUT
}
