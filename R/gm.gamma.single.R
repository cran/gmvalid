gm.gamma.single <-
function (data, X, Y, conditions, conf.level) 
{
    data <- data[, c(X, Y, conditions)]
    data <- table(data)
    N <- sum(data)
    pc.all <- NULL
    pd.all <- NULL
    data.all <- NULL
    if (all(conditions == 0)) {
        conc <- gm.concordant(data)
        disc <- gm.discordant(data)
        all.pairs <- N * (N - 1)
        for (i in 1:dim(data)[1]) for (j in 1:dim(data)[2]) {
            pc <- 0
            pd <- 0
            if (i < dim(data)[1] && j < dim(data)[2]) 
                pc <- sum(data[(i + 1):dim(data)[1], (j + 1):dim(data)[2]])
            if (i > 1 && j > 1) 
                pc <- pc + sum(data[1:(i - 1), 1:(j - 1)])
            if (i < dim(data)[1] && j > 1) 
                pd <- sum(data[(i + 1):dim(data)[1], 1:(j - 1)])
            if (i > 1 && j < dim(data)[2]) 
                pd <- pd + sum(data[1:(i - 1), (j + 1):dim(data)[2]])
            pd.all <- c(pd.all, pd)
            pc.all <- c(pc.all, pc)
            data.all <- c(data.all, data[i, j])
        }
    }
    else {
        index <- matrix(1, nrow = 1, ncol = length(conditions))
        looping <- length(conditions)
        flag <- 0
        next.dim <- FALSE
        conc <- 0
        disc <- 0
        all.pairs <- 0
        while (any(index < dim(data)[-c(1, 2)]) || flag == 0) {
            coord <- paste(c("", index), collapse = ",")
            coord2 <- paste(matrix("", nrow = 1, ncol = length(index) + 
                1), collapse = ",")
            for (i in 1:dim(data)[1]) for (j in 1:dim(data)[2]) {
                pc <- 0
                pd <- 0
                if (i < dim(data)[1] && j < dim(data)[2]) 
                  pc <- sum(eval(parse(text = paste("data[(i+1):dim(data)[1],(j+1):dim(data)[2]", 
                    coord2, "]", sep = ""))))
                if (i > 1 && j > 1) 
                  pc <- pc + sum(eval(parse(text = paste("data[1:(i-1),1:(j-1)", 
                    coord2, "]", sep = ""))))
                if (i < dim(data)[1] && j > 1) 
                  pd <- sum(eval(parse(text = paste("data[(i+1):dim(data)[1],1:(j-1)", 
                    coord2, "]", sep = ""))))
                if (i > 1 && j < dim(data)[2]) 
                  pd <- pd + sum(eval(parse(text = paste("data[1:(i-1),(j+1):dim(data)[2]", 
                    coord2, "]", sep = ""))))
                pd.all <- c(pd.all, pd)
                pc.all <- c(pc.all, pc)
                data.all <- c(data.all, eval(parse(text = paste("data[i,j", 
                  coord, "]", sep = ""))))
            }
            coord <- paste(c("", coord), collapse = ",")
            part.data <- eval(parse(text = paste("data[", coord, 
                "]", sep = "")))
            part.data.length <- sum(part.data)
            all.pairs <- all.pairs + part.data.length * (part.data.length - 
                1)
            conc <- conc + gm.concordant(part.data)
            disc <- disc + gm.discordant(part.data)
            if (all(index == dim(data)[-c(1, 2)])) {
                flag <- 1
                break
            }
            while (looping > 0 && index[looping] == dim(data)[looping + 
                2]) {
                looping <- looping - 1
                next.dim <- TRUE
            }
            if (next.dim) {
                index[looping] <- index[looping] + 1
                index[(looping + 1):length(index)] <- 1
                looping <- length(index)
                next.dim <- FALSE
            }
            else if (looping > 0) 
                index[looping] <- index[looping] + 1
        }
    }
    gamma <- (conc - disc)/(conc + disc)
    conc.p <- 2 * conc/all.pairs
    disc.p <- 2 * disc/all.pairs
    var.sum <- sum(data.all * (pd.all * conc.p - pc.all * disc.p)^2)
    var <- (16/(N * (conc.p + disc.p))^4) * var.sum
    se <- sqrt(var)
    result <- matrix(c(gamma, se, gamma - qnorm(mean(c(1, conf.level))) * 
        se, gamma + qnorm(mean(c(1, conf.level))) * se, pnorm(-abs(gamma), 
        sd = se)), nrow = 1, ncol = 5)
    dname <- names(dimnames(data))
    if (all(conditions == 0)) 
        dimnames(result) <- list(paste(dname[1], "~", dname[2], 
            sep = " "), c("estimate", "SE", "lower", "upper", 
            "p.value"))
    else {
        cname <- paste(dname[-c(1, 2)], collapse = "")
        dimnames(result) <- list(paste(dname[1], "~", dname[2], 
            "|", cname, sep = " "), c("estimate", "SE", "lower", 
            "upper", "p.value"))
    }
    result
}
