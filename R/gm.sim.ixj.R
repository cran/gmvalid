`gm.sim.ixj` <-
function (N, pa, pb) 
{
    ifelse(length(pa) >= length(pb), base <- pa, base <- pb)
    ifelse(length(pa) >= length(pb), pr <- pb, pr <- pa)
    m <- matrix(nrow = length(base), ncol = length(pr), byrow = TRUE)
    flag <- 1
    while (flag == 1) {
        for (i in 1:length(base)) m[i, ] <- .gm.sim.row.i(base[i], 
            pr)
        ifelse(any(is.na(m)) == TRUE, flag <- 1, flag <- 0)
        ifelse(chisq.test(m * N)$p.value > 0.01, flag <- 1, flag <- 0)
    }
    if (length(pa) >= length(pb)) 
        return(N * m)
    else return(N * t(m))
}





`.gm.calculate.delta` <-
function (x, y, r) 
{
    delta_x <- x[2] - x[1]
    delta_y <- y[2] - y[1]
    x_null <- r/sqrt(delta_x^2 + delta_y^2) * delta_x
    if (y[1] < y[2]) 
        y_null <- -sqrt(r^2 - x_null^2)
    else if (y[1] > y[2]) 
        y_null <- sqrt(r^2 - x_null^2)
    else y_null <- 0
    return(c(x_null, y_null))
}
`.gm.chain.selectCV` <-
function (k, ds, data, conf.level = conf.level, strategy = strategy, 
    option.vector = option.vector, chain = chain) 
{
    EOUT <- list()
    POUT <- list()
    for (var.i in 1:dim(data)[2]) {
        data[[names(data)[var.i]]] = as.factor(data[[names(data)[var.i]]])
    }
    chain.elements <- strsplit(chain, "|")[[1]]
    chain.stop <- seq(1:length(chain.elements))[chain.elements %in% 
        setdiff(chain.elements, "|") == FALSE]
    chain.lauf <- seq(1:length(chain.elements))[chain.elements %in% 
        setdiff(chain.elements, "|") == TRUE]
    if (length(chain.lauf) != dim(data)[2]) 
        stop("number of knodes in the chain doesn't correspond to data")
    cckk <- list()
    chain.stop.2 = c(0, chain.stop, length(chain.elements) + 
        1)
    for (i in 1:(length(chain.stop.2) - 1)) cckk[[i]] = chain.elements[(chain.stop.2[i] + 
        1):(chain.stop.2[i + 1] - 1)]
    cckk <- lapply(cckk, sort)
    chain.knodes <- sapply(cckk, paste, collapse = "")[-length(cckk)]
    if (length(cckk) > 1) {
        c.k.o <- NULL
        for (lauf.4 in (length(chain.knodes) - 1):1) {
            c.k.o <- c(c.k.o, paste(chain.knodes[order(chain.knodes[1:(lauf.4 + 
                1)])], collapse = ""))
            c.k.o <- paste(strsplit(c.k.o, "")[[1]][order(strsplit(c.k.o, 
                "")[[1]])], collapse = "")
        }
    }
    else {
        c.k.o <- chain.knodes
    }
    for (i in 1:k) {
        nrows <- dim(data)[2]
        Edges <- matrix(0, nrow = nrows, ncol = nrows, dimnames = list(letters[1:nrows], 
            letters[1:nrows]))
        pval <- matrix(NA, nrow = nrows, ncol = nrows, dimnames = list(letters[1:nrows], 
            letters[1:nrows]))
        df_train <- data[unlist(ds[-i]), ]
        crit.level <- paste("CritLevel ", 1 - conf.level)
        modelselect <- paste("satmodel; stepwise ", strategy, 
            option.vector)
        if (strategy == "f") 
            modelselect <- paste("maine; stepwise ", strategy, 
                option.vector)
        set.chain <- paste("setblocks ", chain, "; blockmode +")
        toMIM(as.gmData.data.frame(df_train))
        mim.cmd(eval(crit.level), look.nice = FALSE)
        mim.cmd(eval(set.chain), look.nice = FALSE)
        mim.cmd(eval(modelselect), look.nice = FALSE)
        model <- mim.cmd("print b", look.nice = FALSE)
        model <- model[-(1:grep("is", model))]
        m1 <- model[model != "|"]
        m2 <- strsplit(m1, ",")
        if (m2[[1]] == "1") 
            m2 <- m2[-seq(from = 1, to = length(m2), by = 2)]
        m2 <- m2[seq(from = length(m2), to = 1)]
        for (lauf.5 in 1:(length(m2) - 1)) {
            if (c.k.o[lauf.5] %in% m2[[lauf.5]]) 
                m2[[lauf.5]] <- m2[[lauf.5]][-match(c.k.o[lauf.5], 
                  m2[[lauf.5]])]
        }
        m2 <- m2[seq(from = length(m2), to = 1)]
        m <- strsplit(unlist(m2), NULL)
        m[[length(m)]] <- setdiff(m[[length(m)]], ".")
        for (elliott.smith in 1:length(m)) {
            platzhalter <- vector()
            tmp <- vector()
            if (length(m[[elliott.smith]]) == 1) 
                next
            if (length(m[[elliott.smith]]) < 3) {
                tmp <- match(m[[elliott.smith]], letters)
                Edges[tmp[1], tmp[2]] <- 1
                platzhalter <- paste("testdel ", eval(unlist(m2)[elliott.smith]))
                p1 <- mim.cmd(eval(platzhalter), look.nice = FALSE)
                if (p1[1] != "Edge") {
                  p <- as.numeric(p1[length(p1)])
                  pval[tmp[1], tmp[2]] <- p
                }
            }
            if (length(m[[elliott.smith]]) >= 3) {
                nr <- choose(length(m[[elliott.smith]]), 2)
                tmpM <- combinations(length(m[[elliott.smith]]), 
                  2, m[[elliott.smith]])
                for (k in 1:nr) {
                  tmp <- match(tmpM[k, ], letters)
                  Edges[tmp[1], tmp[2]] <- 1
                  platzhalter <- paste("testdel ", paste(tmpM[k, 
                    1], tmpM[k, 2], sep = ""))
                  p1 <- mim.cmd(eval(platzhalter), look.nice = FALSE)
                  if (p1[1] != "Edge") {
                    p <- as.numeric(p1[length(p1)])
                    pval[tmp[1], tmp[2]] <- p
                  }
                }
                rm(tmpM)
            }
        }
        EOUT[[i]] <- Edges
        POUT[[i]] <- round(pval, digits = 5)
        mim.cmd("clear")
    }
    OUT <- list(edges = EOUT, pvalue = POUT)
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
`.gm.divideCV` <-
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
`.gm.edgeCV` <-
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
    SE <- apply(OUT, 2, sd, na.rm = TRUE)
    STATISTIC <- list(mean = MEAN, se = SE, pvalue = OUT)
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
`.gm.jointCV` <-
function (pvalue, k, ds, data, conf.level = 0.95, chain = FALSE) 
{
    pvalue.bak <- pvalue
    for (var.i in 1:dim(data)[2]) {
        data[[names(data)[var.i]]] = as.factor(data[[names(data)[var.i]]])
    }
    risk <- list()
    risk.dec <- list()
    test.knodes <- list()
    success <- vector()
    gm <- matrix(nrow = k, ncol = 3, dimnames = list(1:k, c("selected model", 
        "number of edges", "success")))
    knodes <- dim(data)[2]
    dimName <- combinations(knodes, 2, letters)
    edgeName <- apply(dimName, 1, paste, collapse = "")
    if (chain == FALSE) {
        influences <- paste(letters[2:knodes], collapse = "")
        setblock <- paste("setblock ", influences, "|a; blockmode +", 
            sep = "")
    }
    else {
        setblock <- paste("setblocks ", chain, "; blockmode +")
    }
    eo <- grep("a", edgeName)
    edgeOut <- edgeName[1:length(eo)]
    for (i in 1:k) {
        df_train <- data[unlist(ds[-i]), ]
        df_test <- data[ds[[i]], ]
        toMIM(as.gmData.data.frame(df_train))
        pvalue[is.na(pvalue)] <- 1
        e <- ifelse(pvalue[i, ] > 1 - conf.level, TRUE, FALSE)
        edgeDEL <- names(pvalue[i, ])[e]
        mim.cmd(eval(setblock), look.nice = FALSE)
        mim.cmd("satmod", look.nice = FALSE)
        platzhalter <- paste("del ", paste(edgeDEL, col = "", 
            sep = ""), sep = " ")
        if (length(edgeDEL) != 0) {
            mim.cmd(eval(platzhalter), look.nice = FALSE)
        }
        model <- mim.cmd("print", look.nice = FALSE)
        model <- model[-(1:grep("is", model))]
        if (length(model) == 3) {
            model <- model[seq(1, length(model), 2)]
        }
        else if (length(model) > 3) {
            model <- model[seq(2, length(model), 2)]
        }
        model <- c(model[1], strsplit(model[2], ",")[[1]][-1])
        gm[i, 1] <- paste(model, collapse = " | ")
        m <- strsplit(model[length(model)], ",")[[1]]
        m <- m[charmatch("a", m)]
        test.edges <- strsplit(m, NULL)[[1]]
        knode <- setdiff(test.edges, c(".", "a"))
        test.knodes[[i]] <- match(knode, letters)
        if (length(test.knodes[[i]]) <= 0) {
            risk[[i]] <- NA
            risk.dec[[i]] <- NA
            gm[i, 3] <- NA
            gm[i, 2] <- length(test.knodes[[i]])
        }
        else if (length(test.knodes[[i]]) == 1) {
            tab.xy <- table(df_train[, test.knodes[[i]]])
            tab <- table(df_train[, c(test.knodes[[i]], 1)])
            unaff <- tab[, 1]/tab.xy
            aff <- tab[, 2]/tab.xy
            risk[[i]] <- aff/unaff
            risk.dec[[i]] <- ifelse(risk[[i]] >= 1, 2, 1)
            risk.dec[[i]][is.na(risk.dec[[i]])] <- 1.5
            pred.mat <- data.frame(Var1 = 1:length(risk.dec[[i]]), 
                pred = c(risk.dec[[i]]))
            df_test$Var1 <- df_test[, test.knodes[[i]][1]]
            pred.merge <- merge(df_test, pred.mat, by = c("Var1"))
            diff <- abs(as.numeric(pred.merge[[names(df_test)[1]]]) - 
                pred.merge$pred)
            gm[i, 3] <- 1 - sum(diff)/length(diff)
            gm[i, 2] <- length(test.knodes[[i]])
        }
        else {
            tab.xy <- table(df_train[, test.knodes[[i]]])
            tab <- table(df_train[, c(test.knodes[[i]], 1)])
            choord <- paste(matrix("", nrow = 1, ncol = length(dim(tab))), 
                collapse = ",")
            unaff <- eval(parse(text = paste("tab[", choord, 
                "1]", sep = "")))/tab.xy
            aff <- eval(parse(text = paste("tab[", choord, "2]", 
                sep = "")))/tab.xy
            risk[[i]] <- aff/unaff
            risk.dec[[i]] <- ifelse(risk[[i]] >= 1, 2, 1)
            risk.dec[[i]][is.na(risk.dec[[i]])] <- 1.5
            for (lauf.dim in 1:length(dim(risk.dec[[i]]))) {
                dimnames(risk.dec[[i]])[[lauf.dim]] <- 1:dim(risk.dec[[i]])[lauf.dim]
            }
            if (length(dim(risk.dec[[i]])) > 1) {
                tmp <- as.table(array(1, dim = dim(risk.dec[[i]])))
                dimnames(tmp) <- dimnames(risk.dec[[i]])
                p <- expand.table(tmp)
            }
            pred <- as.vector(aperm(risk.dec[[i]], length(dim(risk.dec[[i]])):1))
            pred.mat <- cbind(p, pred)
            for (lauf.dim in 1:length(dim(risk.dec[[i]]))) {
                levels(df_test[, test.knodes[[i]][lauf.dim]]) <- seq(from = 1, 
                  to = length(levels(df_test[, test.knodes[[i]][lauf.dim]])))
            }
            pred.merge <- merge(df_test, pred.mat, by = names(df_test)[test.knodes[[i]]])
            diff <- abs(as.numeric(pred.merge[[names(df_test)[1]]]) - 
                pred.merge$pred)
            gm[i, 3] <- 1 - sum(diff)/length(diff)
            gm[i, 2] <- length(test.knodes[[i]])
        }
    }
    success.mat <- matrix(nrow = k, ncol = 2, dimnames = list(gm[, 
        1], c("No. edges -> a", "success")))
    names(dimnames(success.mat)) <- c("selected models", paste(c("Success probability with ", 
        conf.level * 100, "% C.I."), collapse = ""))
    success.mat[, 1] <- as.numeric(gm[, 2])
    success.mat[, 2] <- as.numeric(gm[, 3])
    OUT <- list(pvalue = pvalue.bak, ratio = risk, risk = risk.dec, 
        success = success.mat[order(success.mat[, 2], decreasing = TRUE), 
            ])
}
`.gm.power.set` <-
function (v, first, dot = FALSE) 
{
    if (dot == TRUE && length(v) < 2) 
        return(c(paste(c("", first, v), collapse = ":"), paste(c("", 
            first), collapse = ":")))
    else if (length(v) < 2) 
        return(c(paste(c(first, v), collapse = ""), paste(c(first), 
            collapse = "")))
    else {
        return(c(Recall(v[2:length(v)], first, dot), Recall(v[2:length(v)], 
            c(first, v[1]), dot)))
    }
}
`.gm.selectCV` <-
function (k, ds, data, conf.level = conf.level, strategy = strategy, 
    option.vector = option.vector) 
{
    EOUT <- list()
    POUT <- list()
    for (var.i in 1:dim(data)[2]) {
        data[[names(data)[var.i]]] = as.factor(data[[names(data)[var.i]]])
    }
    for (i in 1:k) {
        nrows <- dim(data)[2]
        Edges <- matrix(0, nrow = nrows, ncol = nrows, dimnames = list(letters[1:nrows], 
            letters[1:nrows]))
        pval <- matrix(NA, nrow = nrows, ncol = nrows, dimnames = list(letters[1:nrows], 
            letters[1:nrows]))
        df_train <- data[unlist(ds[-i]), ]
        crit.level <- paste("CritLevel ", 1 - conf.level)
        modelselect <- paste("satmodel; step ", strategy, option.vector)
        if (strategy == "f") 
            modelselect <- paste("maine; step ", strategy, option.vector)
        toMIM(as.gmData.data.frame(df_train))
        mim.cmd(eval(crit.level), look.nice = FALSE)
        mim.cmd(eval(modelselect), look.nice = FALSE)
        model <- mim.cmd("print m", look.nice = FALSE)
        m1 <- strsplit(model[length(model)], ",")[[1]]
        m <- strsplit(m1, NULL)
        m[[length(m)]] <- setdiff(m[[length(m)]], ".")
        for (elliott.smith in 1:length(m)) {
            platzhalter <- vector()
            tmp <- vector()
            if (length(m[[elliott.smith]]) == 1) 
                next
            if (length(m[[elliott.smith]]) < 3) {
                tmp <- match(m[[elliott.smith]], letters)
                Edges[tmp[1], tmp[2]] <- 1
                platzhalter <- paste("testdel ", eval(m1[elliott.smith]))
                p1 <- mim.cmd(eval(platzhalter), look.nice = FALSE)
                p <- as.numeric(p1[length(p1)])
                pval[tmp[1], tmp[2]] <- p
            }
            if (length(m[[elliott.smith]]) >= 3) {
                nr <- choose(length(m[[elliott.smith]]), 2)
                tmpM <- combinations(length(m[[elliott.smith]]), 
                  2, m[[elliott.smith]])
                for (k in 1:nr) {
                  tmp <- match(tmpM[k, ], letters)
                  Edges[tmp[1], tmp[2]] <- 1
                  platzhalter <- paste("testdel ", paste(tmpM[k, 
                    1], tmpM[k, 2], sep = ""))
                  p1 <- mim.cmd(eval(platzhalter), look.nice = FALSE)
                  p <- as.numeric(p1[length(p1)])
                  pval[tmp[1], tmp[2]] <- p
                }
                rm(tmpM)
            }
        }
        EOUT[[i]] <- Edges
        POUT[[i]] <- round(pval, digits = 5)
        mim.cmd("clear")
    }
    OUT <- list(edges = EOUT, pvalue = POUT)
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
function (k, l,Delta,d,pr) 
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
function (k, l, j,Delta,d,pr) 
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
function (k, l, j,Delta,d,pr) 
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
function (k, l, j,Delta,d,pr) 
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

