gm.cv.select <-
function (k, ds, data, conf.level = conf.level, strategy = strategy, 
    option.vector = option.vector, show.output = show.output) 
{
    EOUT <- list()
    POUT <- list()
    for (var.i in 1:dim(data)[2]) {
        data[[names(data)[var.i]]] = as.factor(data[[names(data)[var.i]]])
    }
    crit.level <- paste("CritLevel ", 1 - conf.level)
    for (i in 1:k) {
        nrows <- dim(data)[2]
        Edges <- matrix(0, nrow = nrows, ncol = nrows, dimnames = list(letters[1:nrows], 
            letters[1:nrows]))
        pval <- matrix(NA, nrow = nrows, ncol = nrows, dimnames = list(letters[1:nrows], 
            letters[1:nrows]))
        df_train <- data[unlist(ds[-i]), ]
        if (strategy == "f") {
            modelselect <- paste("maine; step ", strategy, option.vector)
        }
        else {
            modelselect <- paste("satmodel; step ", strategy, 
                option.vector)
        }
        toMIM(as.gmData.data.frame(df_train))
        mim.cmd(eval(crit.level), look.nice = FALSE)
        mim.cmd(eval(modelselect), look.nice = show.output)
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
                  pval[tmp[1], tmp[2]] <- as.numeric(p1[length(p1)])
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
