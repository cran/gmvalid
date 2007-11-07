`gm.cv` <-
function (k, data, outcome = 1, strategy = c("backwards", "forwards"), 
    chain = FALSE, options = "", conf.level = 0.95) 
{
    if (k >= dim(data)[1]) 
        stop("k must be smaller than number of observations")
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
    data = data[, c(outcome, c(1:dim(data)[2])[-outcome])]
    strategy <- match.arg(strategy)
    if (strategy == "backwards") 
        strategy <- ""
    else if (strategy == "forwards") 
        strategy <- "f"
    RANDOM <- .gm.divideCV(k, data)
    if (chain == FALSE) {
        SELECTION <- .gm.selectCV(k, RANDOM, data, conf.level = conf.level, 
            strategy = strategy, option.vector = options)
    }
    else {
        SELECTION <- .gm.chain.selectCV(k, RANDOM, data, conf.level = conf.level, 
            strategy = strategy, option.vector = options, chain = chain)
    }
    JUDGE <- .gm.edgeCV(SELECTION$pvalue, k, data)
    JOINT <- .gm.jointCV(JUDGE$pvalue, k, RANDOM, data, conf.level = conf.level)
    names(dimnames(JOINT$pvalue)) <- list("k-folds", "possible edges")
    print(JOINT$success)
    return(JOINT)
}
