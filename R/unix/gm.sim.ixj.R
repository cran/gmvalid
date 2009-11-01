`gm.sim.ixj` <-
function (N, pa, pb) 
{
    ifelse(length(pa) >= length(pb), base <- pa, base <- pb)
    ifelse(length(pa) >= length(pb), pr <- pb, pr <- pa)
    m <- matrix(nrow = length(base), ncol = length(pr), byrow = TRUE)
    for (i in length(dim(m)):1) dimnames(m)[[i]] = as.character(c(1:dim(m)[i]))
    flag <- 1
    while (flag == 1) {
        for (i in 1:length(base)) m[i, ] <- .gm.sim.row.i(base[i], 
            pr)
        ifelse(any(is.na(m)) == TRUE, flag <- 1, flag <- 0)
        if (flag == 0) {
            test = data.frame(expand.table(N * m))
            ifelse(abs(cor(test)[2, 1]) < 0.5, flag <- 1, flag <- 0)
        }
    }
    if (length(pa) >= length(pb)) 
        result = N * m
    else {
        result = N * t(m)
    }
    result
}




`.gm.activate.stop` <-
function (start.model, data = data) 
{
    MODEL <- .gm.make.unique.cliques(start.model)
    statout <- list()
    for (i in 1:length(start.model)) {
        if (length(start.model[[i]]) > 1) 
            statout[[i]] <- gm.gamma(data = data[, start.model[[i]]])
    }
    STAT <- statout
    return(list(measure = STAT, model = MODEL))
}
`.gm.addEdge` <-
function (data, start.model, clique, conf.level = 0.95, output = TRUE) 
{
    model <- start.model[[clique]]
    if (length(model) > 1 && clique < length(start.model)) {
        bvars = NULL
        for (i in (clique + 1):length(start.model)) bvars <- c(bvars, 
            start.model[[i]])
        for (i in 1:length(model)) if (length(which(bvars == 
            model[i])) > 0) 
            bvars = bvars[-which(bvars == model[i])]
        bvars = unique(bvars)
        result <- matrix(nrow = length(bvars) * length(model), 
            ncol = 7)
        index = 1
        for (i in 1:length(bvars)) for (j in 1:length(model)) {
            r <- gm.gamma(X = j, Y = i, data = data, conditions = model[-j], 
                conf.level = conf.level)
            r <- cbind(r, model[j], bvars[i])
            result[index, ] <- r
            index = index + 1
        }
        nrows <- seq(1, dim(result)[1])
        kante.hin <- nrows[result[, 5] <= 1 - conf.level]
        result.out <- matrix(0, nrow = nrow(result), ncol = (ncol(result) + 
            1))
        result.out[1:nrow(result), 1:ncol(result)] <- result
        result.out[kante.hin, (ncol(result) + 1)] <- 1
        if (output) 
            print(result.out)
        ifelse(length(kante.hin) > 0, zeile <- kante.hin, zeile <- NA)
        ifelse(length(kante.hin) > 0, edge.name <- result.out[kante.hin, 
            6:7], edge.name <- NA)
        ifelse(length(kante.hin) > 0, p.value <- result.out[kante.hin, 
            5], p.value <- NA)
        OUT <- list(zeile = zeile, edge.name = edge.name, p.value = p.value, 
            result = result.out)
    }
    else if (length(model) == 1 && clique < length(start.model)) {
        bvars = NULL
        for (i in (clique + 1):length(start.model)) bvars <- c(bvars, 
            start.model[[i]])
        bvars = unique(bvars)
        result <- matrix(nrow = length(bvars), ncol = 5)
        for (i in 1:length(bvars)) {
            result[i, ] <- gm.gamma(X = data[, model], Y = data[, 
                bvars[i]], conf.level = conf.level)
        }
        result <- cbind(result, rep(model, dim(result)[1]), bvars)
        dimnames(result) <- list(1:dim(result)[1], c("estimate", 
            "SE", "lower", "upper", "p.value", "knot1", "knot2"))
        nrows <- seq(1, dim(result)[1])
        kante.hin <- nrows[result[, 5] <= 1 - conf.level]
        result.out <- matrix(0, nrow = nrow(result), ncol = (ncol(result) + 
            1))
        result.out[1:nrow(result), 1:ncol(result)] <- result
        result.out[kante.hin, (ncol(result) + 1)] <- 1
        if (output) 
            print(result.out)
        ifelse(length(kante.hin) > 0, zeile <- kante.hin, zeile <- NA)
        ifelse(length(kante.hin) > 0, edge.name <- result.out[kante.hin, 
            6:7], edge.name <- NA)
        ifelse(length(kante.hin) > 0, p.value <- result.out[kante.hin, 
            5], p.value <- NA)
        OUT <- list(zeile = zeile, edge.name = edge.name, p.value = p.value, 
            result = result.out)
    }
    else {
        OUT <- list(zeile = NA, edge.name = NA, p.value = NA, 
            result = NA)
    }
    OUT
}
`.gm.calculate.delta` <-
function (x, y, r) 
{
    delta_x <- x[2] - x[1]
    delta_y <- y[2] - y[1]
    x_null <- r/sqrt(delta_x^2 + delta_y^2) * delta_x
    if (y[1] < y[2]) 
        y_null <- -sqrt(abs(r^2 - x_null^2))
    else if (y[1] > y[2]) 
        y_null <- sqrt(abs(r^2 - x_null^2))
    else y_null <- 0
    return(c(x_null, y_null))
}
`.gm.concordant` <-
function (x) 
{
    mat.lr <- function(r, c) {
        lr <- x[(r.x > r) & (c.x > c)]
        sum(lr)
    }
    r.x <- row(x)
    c.x <- col(x)
    sum(x * mapply(mat.lr, r = r.x, c = c.x))
}
`.gm.deleteEdge` <-
function (data, model, conf.level = 0.95, output = TRUE) 
{
    dimName <- combinations(length(model), 2, v = model)
    result <- gm.gamma(data = data[, model], conf.level = conf.level)
    nrows <- seq(1, dim(result)[1])
    kante.weg <- nrows[result[, 5] > 1 - conf.level]
    if (max(result[, 5], na.rm = TRUE) <= 1 - conf.level) 
        kante.weg <- vector()
    result.out <- matrix(0, nrow = nrow(result), ncol = (ncol(result) + 
        1))
    dimnames(result.out) <- list(dimnames(result)[[1]], c(dimnames(result)[[2]], 
        "edge.out"))
    result.out[1:nrow(result), 1:ncol(result)] <- result
    result.out[kante.weg, (ncol(result) + 1)] <- 1
    if (output) 
        print(result.out)
    ifelse(length(kante.weg) > 0, edge.out <- kante.weg, edge.out <- NA)
    ifelse(length(kante.weg) > 0, edge.name <- dimName[kante.weg, 
        ], edge.name <- NA)
    ifelse(length(kante.weg) > 0, p.value <- result[kante.weg, 
        5], p.value <- NA)
    OUT <- list(edge.out = edge.out, edge.name = edge.name, p.value = p.value, 
        result = result.out)
    return(OUT)
}
`.gm.delta` <-
function (k, l, i, j, pen, A, B) 
{
    if (!missing(i)) {
        coord.A <- paste(c(-i, l), collapse = ",")
    }
    else if (missing(i)) {
        coord.A <- paste(c("", l), collapse = ",")
    }
    if (!missing(j)) {
        coord.B <- paste(c(k, -j), collapse = ",")
    }
    else if (missing(i)) {
        coord.B <- paste(c(k, ""), collapse = ",")
    }
    sum(A * eval(parse(text = paste("pen[", coord.A, "]", sep = "")))) + 
        sum(B * eval(parse(text = paste("pen[", coord.B, "]", 
            sep = "")))) - pen[k, l]
}
`.gm.delta2` <-
function (k, l, i, j, pen, A, B, AB) 
{
    if (!missing(i)) {
        coord.A <- paste(c(-i, l), collapse = ",")
    }
    else if (missing(i)) {
        coord.A <- paste(c("", l), collapse = ",")
    }
    if (!missing(j)) {
        coord.B <- paste(c(k, -j), collapse = ",")
    }
    else if (missing(i)) {
        coord.B <- paste(c(k, ""), collapse = ",")
    }
    if (!missing(j)) {
        coord.AB <- paste(c(k, -j), collapse = ",")
    }
    else if (missing(i)) {
        coord.AB <- paste(c(k, ""), collapse = ",")
    }
    m <- k * l
    if (k == 2 && l == 1) 
        m <- 3
    mm <- c(4, 3, 2, 1)
    n <- mm[m]
    (AB[m] - AB[n]) * pen[k, l] + B[-k] * pen[-k, l] + A[-l] * 
        pen[k, -l]
}
`.gm.discordant` <-
function (x) 
{
    mat.ll <- function(r, c) {
        ll <- x[(r.x > r) & (c.x < c)]
        sum(ll)
    }
    r.x <- row(x)
    c.x <- col(x)
    sum(x * mapply(mat.ll, r = r.x, c = c.x))
}
`.gm.draw.text` <-
function (label, x, y, xy_null = c(0, 0), color, alignment = c("left", 
    "bottom"), fontsize) 
{
    x_label <- x[1] + xy_null[1]
    y_label <- y[1] + xy_null[2]
    x_delta <- x[2] - x[1]
    y_delta <- y[2] - y[1]
    angle = atan(y_delta/x_delta) * (180/pi)
    if (angle == -90) 
        angle = 90
    if (is.numeric(label)) 
        label <- as.numeric(as.integer(label * 10000)/10000)
    pushViewport(viewport(x = x_label, y = y_label, width = 0, 
        height = , angle = angle, name = "vp1", just = alignment))
    grid.text(label, x = 0, y = unit(0.75, "mm"), just = alignment, 
        gp = gpar(fontsize = fontsize, col = color))
    popViewport()
}
`.gm.gamma.backward` <-
function (data, start.model, onestep, headlong, conf.level) 
{
    durchgang <- 1
    flag <- 0
    while (flag == 0) {
        cat("Durchgang : ", durchgang, "\n")
        edge.in.more.than.one.clique <- vector()
        result <- list()
        pv.mat <- matrix(ncol = 5)
        no.clique <- length(start.model)
        edge.out <- edge.name <- pv <- vert1 <- vert2 <- vector()
        for (i in 1:no.clique) {
            if (length(start.model[[i]]) > 1) {
                result[[i]] <- .gm.deleteEdge(data = data, model = start.model[[i]], 
                  conf.level = conf.level, output = FALSE)
                edge.out <- result[[i]]$edge.out
                edge.name <- result[[i]]$edge.name
                for (k in 1:length(edge.out)) if (!is.na(edge.out[k])) {
                  pv <- result[[i]]$result[edge.out[k], 5]
                  vert1 <- split(t(edge.name), 1:2)[[1]][k]
                  vert2 <- split(t(edge.name), 1:2)[[2]][k]
                  pv.mat <- rbind(pv.mat, cbind(i, edge.out[k], 
                    pv, vert1, vert2))
                }
                else {
                  pv.mat <- rbind(pv.mat, cbind(i, edge.out[k], 
                    NA, NA, NA))
                }
            }
            else {
                pv.mat <- rbind(pv.mat, cbind(i, NA, NA, NA, 
                  NA))
            }
        }
        pv.mat <- pv.mat[order(pv.mat[, 3], na.last = TRUE, decreasing = TRUE), 
            ]
        if (onestep == TRUE) {
            model.matrix = matrix(1, nrow = dim(data)[2], ncol = dim(data)[2])
            model.matrix = upper.tri(model.matrix)
            i = 1
            while (!is.na(pv.mat[i, 4])) {
                selected.edge <- c(pv.mat[i, 4], pv.mat[i, 5])
                selected.edge <- sort(selected.edge)
                model.matrix[selected.edge[1], selected.edge[2]] = 0
                i = i + 1
            }
            model.string = .gm.modelparse(model.matrix)
            model.string <- strsplit(model.string, ",")[[1]]
            start.model = list()
            for (i in 1:length(model.string)) start.model[[i]] <- match(strsplit(model.string[i], 
                "")[[1]], letters)
            OUT <- .gm.activate.stop(start.model, data)
            flag <- 1
        }
        else if (length(pv.mat[, 3][!is.na(pv.mat[, 3])]) > 0) {
            if (headlong == TRUE) {
                selected <- sample(length(pv.mat[, 3][!is.na(pv.mat[, 
                  3])]), 1)
            }
            else selected <- 1
            selected.clique <- pv.mat[selected, 1]
            selected.edge <- c(pv.mat[selected, 4], pv.mat[selected, 
                5])
            if (length(start.model) > 1 && flag == 0) {
                for (lauf.model in 1:length(start.model)) {
                  edge.in.more.than.one.clique <- c(edge.in.more.than.one.clique, 
                    lauf.model[all(selected.edge %in% start.model[[lauf.model]]) == 
                      TRUE])
                }
            }
            else {
                edge.in.more.than.one.clique <- selected.clique
            }
        }
        else {
            selected.clique <- NA
            selected.edge <- NA
            OUT <- .gm.activate.stop(start.model, data = data)
            flag <- 1
        }
        if (flag == 0) {
            for (blur in 1:length(edge.in.more.than.one.clique)) {
                edge.match <- match(selected.edge, start.model[[edge.in.more.than.one.clique[blur]]])
                start.model[[length(start.model) + 1]] <- start.model[[edge.in.more.than.one.clique[blur]]][-edge.match[1]]
                start.model[[length(start.model) + 1]] <- start.model[[edge.in.more.than.one.clique[blur]]][-edge.match[2]]
            }
            start.model <- unique(start.model[-edge.in.more.than.one.clique])
            merke <- vector()
            lauf <- 0
            l <- order(sapply(start.model, length))
            ml <- start.model[l]
            for (i in 1:(length(ml) - 1)) {
                for (j in (i + 1):length(ml)) {
                  lauf <- lauf + 1
                  ifelse(all(ml[[i]] %in% ml[[j]]), merke[lauf] <- i, 
                    merke[lauf] <- NA)
                }
            }
            if (any(!is.na(merke))) {
                merke <- merke[!is.na(merke)]
                start.model <- ml[-merke]
            }
        }
        if (all(sapply(start.model, length) == 1) && flag == 
            0) {
            OUT <- .gm.activate.stop(start.model, data = data)
            flag <- 1
        }
        else {
            durchgang <- durchgang + 1
        }
    }
    OUT
}
`.gm.gamma.forward` <-
function (data, start.model, onestep, headlong, conf.level = 0.95) 
{
    durchgang <- 1
    flag <- 0
    model.matrix = matrix(0, nrow = dim(data)[2], ncol = dim(data)[2])
    while (flag == 0) {
        cat("Durchgang : ", durchgang, "\n")
        vert.in.clique <- vector()
        result <- list()
        pv.mat <- matrix(ncol = 5)
        no.clique <- length(start.model)
        edge.name <- pv <- vert1 <- vert2 <- vector()
        for (i in 1:no.clique) {
            if (length(start.model[[i]]) > 0) {
                result[[i]] <- .gm.addEdge(data = data, start.model = start.model, 
                  clique = i, conf.level = conf.level, output = FALSE)
                zeile <- result[[i]]$zeile
                edge.name <- result[[i]]$edge.name
                edge.name <- matrix(edge.name, ncol = 2)
                if (all(!is.na(zeile))) {
                  pv <- result[[i]]$p.value
                  pv.mat <- rbind(pv.mat, cbind(i, zeile, pv, 
                    edge.name))
                }
                else {
                  pv.mat <- rbind(pv.mat, cbind(i, zeile, NA, 
                    NA, NA))
                }
            }
            else {
                pv.mat <- rbind(pv.mat, cbind(i, NA, NA, NA, 
                  NA))
            }
        }
        pv.mat <- pv.mat[order(pv.mat[, 3], na.last = TRUE, decreasing = FALSE), 
            ]
        if (onestep == TRUE) {
            i = 1
            while (!is.na(pv.mat[i, 4])) {
                selected.edge <- c(pv.mat[i, 4], pv.mat[i, 5])
                selected.edge <- sort(selected.edge)
                model.matrix[selected.edge[1], selected.edge[2]] = 1
                i = i + 1
            }
            model.string = .gm.modelparse(model.matrix)
            model.string <- strsplit(model.string, ",")[[1]]
            start.model = list()
            for (i in 1:length(model.string)) start.model[[i]] <- match(strsplit(model.string[i], 
                "")[[1]], letters)
            OUT <- .gm.activate.stop(start.model, data)
            flag <- 1
        }
        else if (length(pv.mat[, 3][!is.na(pv.mat[, 3])]) > 0) {
            if (headlong == TRUE) {
                selected <- sample(length(pv.mat[, 3][!is.na(pv.mat[, 
                  3])]), 1)
            }
            else selected <- 1
            selected.clique <- pv.mat[selected, 1]
            selected.edge <- c(pv.mat[selected, 4], pv.mat[selected, 
                5])
            selected.edge <- sort(selected.edge)
            kk = 2
            while (model.matrix[selected.edge[1], selected.edge[2]] == 
                1) {
                selected.edge <- c(pv.mat[kk, 4], pv.mat[kk, 
                  5])
                if (any(is.na(selected.edge))) {
                  flag = 1
                  break
                }
                selected.edge <- sort(selected.edge)
                kk = kk + 1
            }
            if (length(start.model) > 1 && flag == 0) {
                model.matrix[selected.edge[1], selected.edge[2]] = 1
                model.string = .gm.modelparse(model.matrix)
                model.string <- strsplit(model.string, ",")[[1]]
                start.model = list()
                for (i in 1:length(model.string)) start.model[[i]] = match(strsplit(model.string[i], 
                  "")[[1]], letters)
            }
            else {
                selected.clique <- NA
                selected.edge <- NA
                OUT <- .gm.activate.stop(start.model, data = data)
                flag <- 1
            }
        }
        else {
            selected.clique <- NA
            selected.edge <- NA
            OUT <- .gm.activate.stop(start.model, data = data)
            flag <- 1
        }
        if (flag == 0 && sapply(start.model, length) == dim(data)[2]) {
            OUT <- .gm.activate.stop(start.model, data = data)
            flag <- 1
        }
        else {
            durchgang <- durchgang + 1
        }
    }
    OUT
}
`.gm.gamma.single` <-
function (data, X, Y, conditions, conf.level) 
{
    data <- data[, c(X, Y, conditions)]
    data <- table(data)
    N <- sum(data)
    pc.all <- NULL
    pd.all <- NULL
    data.all <- NULL
    if (all(conditions == 0)) {
        conc <- .gm.concordant(data)
        disc <- .gm.discordant(data)
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
            conc <- conc + .gm.concordant(part.data)
            disc <- disc + .gm.discordant(part.data)
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
`.gm.make.unique.cliques` <-
function (start.model) 
{
    letter.func <- function(x) letters[x]
    letter.list <- lapply(start.model, letter.func)
    letter.list <- lapply(letter.list, paste, collapse = "")
    letter.list <- unique(letter.list)
    model <- paste(matrix(sapply(letter.list, paste, collapse = ""), 
        nrow = 1, ncol = length(letter.list)), collapse = ",")
}
`.gm.matrixparse` <-
function (result) 
{
    m = strsplit(result, "")[[1]]
    elements = unique(m)
    if (any(elements == ",")) 
        elements = elements[-which(elements == ",")]
    elements = sort(elements)
    dep.table = matrix(0, nrow = length(elements), ncol = length(elements))
    dimnames(dep.table) = list(elements, elements)
    for (i in 1:length(m)) {
        if (length(which(elements == m[i])) == 0) {
        }
        else {
            j = 1
            while (i + j <= length(m) && length(which(elements == 
                m[i + j])) > 0) {
                pos1 = which(elements == m[i])
                pos2 = which(elements == m[i + j])
                dep.table[pos1, pos2] = 1
                j = j + 1
            }
        }
    }
    dep.table
}
`.gm.modelparse` <-
function (data) 
{
    v <- dim(data)[1]
    for (i in 1:(v - 1)) {
        for (j in (i + 1):v) {
            data[j, i] <- data[i, j]
        }
    }
    maxvar <- v
    SetSize <- function(X) {
        i <- 0
        while ((X[(i + 1)] > 0) && (i < maxvar)) {
            i <- i + 1
        }
        out <- i
    }
    Neighbours <- function(u, G) {
        j <- 0
        Y <- rep(0, maxvar)
        for (i in 1:maxvar) {
            if (G[i, u] == 1) {
                j <- j + 1
                Y[j] <- i
            }
        }
        out <- Y
    }
    AndSet <- function(X, Y) {
        nx <- SetSize(X)
        ny <- SetSize(Y)
        k <- 0
        Z <- rep(0, maxvar)
        if ((nx > 0) && (ny > 0)) {
            for (i in 1:nx) {
                for (j in 1:ny) {
                  if (X[i] == Y[j]) {
                    k <- k + 1
                    Z[k] <- X[i]
                  }
                }
            }
        }
        out <- Z
    }
    BK <- function(R, P, X, G) {
        result <- NULL
        if ((SetSize(P) == 0) && (SetSize(X) == 0)) {
            result <- rbind(result, R)
        }
        else {
            k <- SetSize(P)
            if (k == 0) {
            }
            else {
                for (i in k:1) {
                  u <- P[i]
                  P <- c(P[-i], 0)
                  nr <- SetSize(R)
                  if (nr == 0) {
                    Rnew <- c(u, rep(0, (maxvar - 1)))
                  }
                  else {
                    Rnew <- c(R[1:nr], u, rep(0, (maxvar - nr - 
                      1)))
                  }
                  N <- Neighbours(u, G)
                  Pnew <- AndSet(P, N)
                  Xnew <- AndSet(X, N)
                  result = rbind(result, BK(Rnew, Pnew, Xnew, 
                    G))
                  nx <- SetSize(X)
                  if (nx == 0) {
                    X <- c(u, rep(0, (maxvar - 1)))
                  }
                  else {
                    X <- c(X[1:nx], u, rep(0, (maxvar - nx - 
                      1)))
                  }
                }
            }
        }
        result
    }
    X <- BK(rep(0, maxvar), 1:maxvar, rep(0, maxvar), data)
    result <- NULL
    for (i in 1:dim(X)[1]) {
        result <- c(result, ",", letters[sort(X[i, ])])
    }
    result <- paste(result[-1], collapse = "")
    result <- paste(sort(strsplit(result, ",")[[1]]), collapse = ",")
    result
}
`.gm.power.set` <-
function (v, first) 
{
    if (length(v) < 2) 
        return(c(paste(c(first, v), collapse = ""), paste(c(first), 
            collapse = "")))
    else {
        return(c(Recall(v[2:length(v)], first), Recall(v[2:length(v)], 
            c(first, v[1]))))
    }
}
`.gm.sim.row.i` <-
function (base, pr) 
{
    c00 <- runif(1, min = 0.001, max = min(base, pr[1]))
    rest <- base - c00
    pr.i <- length(pr) - 1
    sim <- vector()
    i <- 1
    while (pr.i > 0) {
        if (pr.i > 1) {
            sim[i] <- runif(1, min = 0.001, max = min(rest, pr[i + 
                1]))
        }
        else {
            sim[i] <- rest
        }
        rest <- rest - sim[i]
        pr.i <- pr.i - 1
        i <- i + 1
    }
    row <- c(c00, sim)
    return(row)
}
`.gm.var.csi.part1` <-
function (k, l, Delta, d, pr) 
{
    if (k == 1) 
        m <- 2
    else if (k == 2) 
        m <- 1
    if (l == 1) 
        n <- 2
    else if (l == 2) 
        n <- 1
    out <- ((-pr[m, n, 2] + Delta[m, n] * pr[m, n, 1])/sum(pr[m, 
        n, ]) + (pr[k, n, 2] + Delta[k, n] * pr[k, n, 1])/sum(pr[k, 
        n, ]) + (pr[m, l, 2] + Delta[m, l] * pr[m, l, 1])/sum(pr[m, 
        l, ]))/d[m, n]
}
`.gm.var.csi.part2` <-
function (k, l, j, Delta, d, pr) 
{
    ind <- vector()
    if (k == 1) 
        m <- 2
    else if (k == 2) 
        m <- 1
    if (l == 1) 
        n <- 2
    else if (l == 2) 
        n <- 1
    if (j == 1) 
        y <- 2
    else if (j == 2) 
        y <- 1
    ifelse(y == 2, ind <- 1, ind <- 0)
    out <- ((-1)^(1 - ind) * (ind + (1 - ind) * Delta[k, l] - 
        pr[k, l, y] * (sum(pr[m, n, ])/sum(pr[k, l, ])^2) * (1 + 
            Delta[k, l])))/d[k, l]
}
`.gm.var.csi.part3` <-
function (k, l, j, Delta, d, pr) 
{
    ind <- vector()
    if (k == 1) 
        m <- 2
    else if (k == 2) 
        m <- 1
    if (l == 1) 
        n <- 2
    else if (l == 2) 
        n <- 1
    if (j == 1) 
        y <- 2
    else if (j == 2) 
        y <- 1
    ifelse(y == 2, ind <- 1, ind <- 0)
    out <- (((-1)^(1 - ind) * (ind + (1 - ind) * Delta[k, n] - 
        pr[k, l, y] * (sum(pr[m, l, ])/sum(pr[k, l, ])^2) * (1 + 
            Delta[k, n])))/d[k, n] + ((-1)^(1 - ind) * (ind + 
        (1 - ind) * Delta[m, l] - pr[k, l, y] * (sum(pr[k, n, 
        ])/sum(pr[k, l, ])^2) * (1 + Delta[m, l])))/d[m, l])
}
`.gm.var.csi.part4` <-
function (k, l, j, Delta, d, pr) 
{
    ind <- vector()
    ind.g <- vector()
    if (k == 1) 
        m <- 2
    else if (k == 2) 
        m <- 1
    if (l == 1) 
        n <- 2
    else if (l == 2) 
        n <- 1
    if (j == 1) 
        y <- 2
    else if (j == 2) 
        y <- 1
    ifelse(y == 2, ind <- 1, ind <- 0)
    ifelse(k > l, ind.g <- 1, ind.g <- 0)
    out <- ((ind.g * (pr[m, n, 2] - Delta[m, n] * pr[m, n, 1])/sum(pr[m, 
        n, ]))/d[m, n] + ((1 - ind.g) * (-pr[m, n, 2] + Delta[m, 
        n] * pr[m, n, 1])/sum(pr[m, n, ]))/d[m, n] + ((pr[k, 
        n, 2] - Delta[m, n] * pr[k, n, 1])/sum(pr[k, n, ]) + 
        (pr[m, l, 2] - Delta[m, n] * pr[m, l, 1])/sum(pr[m, l, 
            ]))/d[m, n])
}
