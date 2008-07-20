`gm.cv` <-
function (k, data, outcome = 1, strategy = c("backwards", "forwards"), 
    chain = FALSE, options = "", conf.level = 0.95) 
{
    require(mimR, quietly = TRUE)
    if ((!missing(k) && (length(k) != 1 || is.na(k))) || k < 
        1) 
        stop("k must be a single number > 0")
    if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || 
        conf.level < 0 || conf.level > 1)) 
        stop("'conf.level' must be a single number between 0 and 1")
    if (chain != FALSE && any(is.na(match(setdiff(strsplit(chain, 
        "|")[[1]], "|"), letters)))) 
        stop("Use lowercase letters in 'chain'.")
    if (chain != FALSE && is.na(match(strsplit(chain, "")[[1]], 
        "a")[length(match(strsplit(chain, "")[[1]], "a"))])) 
        stop("outcome variable 'a' has to be in the very right in 'chain'.")
    if (is.array(data)) {
        require(epitools, quietly = TRUE)
        if (!length(dimnames(data))) 
            for (i in length(dim(data)):1) dimnames(data)[[i]] = as.character(c(1:dim(data)[i]))
        data = data.frame(expand.table(data))
    }
    if (k >= dim(data)[1]) 
        stop("k must be smaller than number of observations")
    data = data[, c(outcome, c(1:dim(data)[2])[-outcome])]
    var.names = cbind(letters[1:dim(data)[2]], dimnames(data)[2][[1]])
    if (dim(table(data))[1] != 2) 
        stop("outcome variable 'a' has to be binary.")
    names(data) <- letters[1:length(data)]
    strategy <- match.arg(strategy)
    if (strategy == "backwards") 
        strategy <- ""
    else if (strategy == "forwards") 
        strategy <- "f"
    RANDOM <- .gm.cv.divide(k, data)
    if (chain == FALSE) {
        SELECTION <- .gm.cv.select(k, RANDOM, data, conf.level = conf.level, 
            strategy = strategy, option.vector = options)
    }
    else {
        ch = strsplit(chain, "")[[1]]
        s = unique(ch)
        p1 = which(s == "|")
        if (length(p1)) 
            s = s[-p1]
        p2 = which(s == ",")
        if (length(p2)) 
            s = s[-p2]
        s = sort(s)
        for (i in 1:length(ch)) {
            if (ch[i] != "|" && ch[i] != ",") 
                ch[i] = letters[which(s == ch[i])]
        }
        chain = paste(ch, collapse = "")
        SELECTION <- .gm.cv.select.chain(k, RANDOM, data, conf.level = conf.level, 
            strategy = strategy, option.vector = options, chain = chain)
    }
    JUDGE <- .gm.cv.edge(SELECTION$pvalue, k, data)
    JOINT <- .gm.cv.joint(JUDGE$pvalue, k, RANDOM, data, conf.level = conf.level, 
        chain = chain)
    names(dimnames(JOINT$pvalue)) <- list("k-folds", "possible edges")
    mim.cmd(paste(c("brmodel", dimnames(JOINT$success)[[1]][1]), 
        collapse = " "), look.nice = FALSE)
    mim.cmd("gr")
    JOINT$"variable names" = var.names
    print(JOINT$success)
    print(JOINT$statistics)
    print(JOINT$"variable names")
    return(JOINT)
}
