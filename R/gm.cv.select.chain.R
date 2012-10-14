gm.cv.select.chain <-
function (k, ds, data, conf.level = conf.level, strategy = strategy, 
    option.vector = option.vector, chain = chain, show.output = FALSE) 
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
        mim.cmd(eval(set.chain))
        mim.cmd(eval(modelselect), look.nice = show.output)
        model <- mim.cmd("print b", look.nice = FALSE)
        model <- model[-(1:grep("is", model))]
        m1 <- model[model != "|"]
        m2 <- strsplit(m1, ",")
        if (any(m2[[1]] == "1")) 
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
                if (is.na(tmp[1])) 
                  tmp <- match(m[[elliott.smith]], LETTERS)
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
                  if (is.na(tmp[1])) 
                    tmp <- match(tmpM[k, ], LETTERS)
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
