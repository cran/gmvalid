gm.mim <-
function (data, strategy = c("backwards", "forwards", "eh", "combined"), 
    model = FALSE, chain = FALSE, options = "", tests = FALSE) 
{
    require(mimR, quietly = TRUE)
    strategy = match.arg(strategy)
    rejected = NULL
    accepted = NULL
    sat.tests = NULL
    if (is.array(data)) {
        require(epitools, quietly = TRUE)
        if (!length(dimnames(data))) 
            for (i in length(dim(data)):1) dimnames(data)[[i]] = as.character(c(1:dim(data)[i]))
        data = data.frame(expand.table(data))
    }
    var.names = cbind(letters[1:dim(data)[2]], dimnames(data)[2][[1]])
    dimnames(data)[2][[1]] = letters[1:dim(data)[2]]
    for (i in 1:dim(data)[2]) data[[names(data)[i]]] = as.factor(data[[names(data)[i]]])
    toMIM(as.gmData(data))
    if (strategy == "backwards") {
        if (chain != FALSE) {
            mim.cmd(paste(c("setblocks ", chain), collapse = ""))
            mim.cmd("blockmode +")
            if (model != FALSE) 
                mim.cmd(paste(c("brmodel ", model), collapse = ""))
            else {
                mim.cmd("satmod")
                model = chain
            }
        }
        else if (model != FALSE) 
            mim.cmd(paste(c("model ", model), collapse = ""))
        else {
            mim.cmd("satmod")
            model = paste(letters[1:dim(data)[2]], collapse = "")
        }
        mim.text = mim.cmd(paste(c("stepwise ", options), collapse = ""))
        if (chain != FALSE) {
            accepted = mim.cmd("pr")
            accepted = accepted[7:length(accepted)]
            accepted[1:length(accepted)%%2 == 0] = "|"
            accepted = paste(accepted, collapse = "")
        }
        else accepted = toString(mim.text[length(mim.text)])
        if (accepted != paste(letters[1:dim(data)[2]], collapse = "")) {
            sat.tests = mim.cmd("test")
            sat.tests = sat.tests[8:length(sat.tests)]
        }
        if (chain != FALSE) 
            mim.cmd("blockmode -")
    }
    else if (strategy == "forwards") {
        if (chain != FALSE) {
            mim.cmd(paste(c("setblocks ", chain), collapse = ""))
            mim.cmd("blockmode +")
            if (model != FALSE) 
                mim.cmd(paste(c("brmodel ", model), collapse = ""))
            else {
                mim.cmd("maine")
                model = paste(letters[1:dim(data)[2]], collapse = ",")
            }
        }
        else if (model != FALSE) 
            mim.cmd(paste(c("Model ", model), collapse = ""))
        else {
            mim.cmd("maine")
            model = paste(letters[1:dim(data)[2]], collapse = ",")
        }
        mim.text = mim.cmd(paste(c("stepwise f", options), collapse = ""))
        if (chain != FALSE) {
            accepted = mim.cmd("pr")
            accepted = accepted[7:length(accepted)]
            accepted[1:length(accepted)%%2 == 0] = "|"
            accepted = paste(accepted, collapse = "")
        }
        else accepted = toString(mim.text[length(mim.text)])
        if (accepted != paste(letters[1:dim(data)[2]], collapse = "")) {
            sat.tests = mim.cmd("test")
            sat.tests = sat.tests[-(1:(grep("LR:", sat.tests) - 
                1))]
        }
        if (chain != FALSE) 
            mim.cmd("blockmode -")
    }
    else if (strategy == "eh") {
        if (chain != FALSE) {
            mim.cmd(paste(c("setblocks ", chain), collapse = ""))
            mim.cmd("blockmode +")
        }
        if (model == FALSE) 
            model = paste(c(paste(letters[1:dim(data)[2]], collapse = ","), 
                " - ", paste(letters[1:dim(data)[2]], collapse = "")), 
                collapse = "")
        mim.text = mim.cmd(paste(c("initsearch ", model, ";startsearch ", 
            options), collapse = ""))
        accepted = mim.cmd("ehshow a")
        rejected = mim.cmd("ehshow r")
        sat.tests = c(sat.tests, accepted[5:length(accepted)][1:(length(accepted) - 
            4)%%7 != 1])
        sat.tests = c(sat.tests, rejected[5:length(rejected)][1:(length(rejected) - 
            4)%%7 != 1])
        accepted = accepted[5:length(accepted)][1:(length(accepted) - 
            4)%%7 == 1]
        rejected = rejected[5:length(rejected)][1:(length(rejected) - 
            4)%%7 == 1]
        if (chain != FALSE) {
            mim.cmd("blockmode -")
        }
    }
    else if (strategy == "combined" && chain == FALSE) {
        mim.cmd(paste(c("model ", gm.screening(data)$model), 
            collapse = ""))
        mim.text = mim.cmd(paste(c("stepwise j", options), collapse = ""))
        accepted = toString(mim.text[length(mim.text)])
        if (accepted != paste(letters[1:dim(data)[2]], collapse = "")) {
            sat.tests = mim.cmd("test")
            sat.tests = sat.tests[8:length(sat.tests)]
        }
    }
    else stop("No valid strategy!")
    if (length(sat.tests) > 0 && length(sat.tests)%%2 == 0) {
        sat.tests = sat.tests[1:length(sat.tests)%%2 == 0]
        sat.tests = as.numeric(sat.tests)
        dim(sat.tests) = c(length(sat.tests)/3, 3)
        if (length(rejected) > 0) 
            dimnames(sat.tests) = list(c(accepted, rejected), 
                c("LR", "df", "p.value"))
        else dimnames(sat.tests) = list(accepted, c("LR", "df", 
            "p.value"))
    }
    p.list = list()
    if (chain == FALSE && tests == TRUE) 
        for (ll in 1:length(accepted)) {
            p.matrix = matrix(NA, nrow = dim(data)[2], ncol = dim(data)[2])
            dimnames(p.matrix) = list(letters[1:dim(data)[2]], 
                letters[1:dim(data)[2]])
            p.test = mim.cmd(paste(c("model ", accepted[ll], 
                "; step ou"), collapse = ""), look.nice = FALSE)
            if (any(p.test == "+")) 
                p.test = p.test[-which(p.test == "+")]
            if (length(which(p.test == "P")) > 0) {
                p.test = p.test[-(1:grep("P", p.test)[2])]
                edges = grep("[", p.test, fixed = TRUE, value = TRUE)
                p.values = NULL
                if (length(grep("[", p.test, fixed = TRUE)) > 
                  1) 
                  for (i in 1:(length(grep("[", p.test, fixed = TRUE)) - 
                    1)) {
                    p.values = c(p.values, as.numeric(p.test[grep("[", 
                      p.test, fixed = TRUE)[i + 1] - 1]))
                  }
                p.values = c(p.values, as.numeric(p.test[length(p.test)]))
                for (i in 1:length(edges)) {
                  e = strsplit(edges[i], "")[[1]]
                  p.matrix[which(letters == e[2]), which(letters == 
                    e[3])] = p.values[i]
                }
            }
            p.test = mim.cmd(paste(c("model ", accepted[ll], 
                "; step ufo"), collapse = ""), look.nice = FALSE)
            if (any(p.test == "+")) 
                p.test = p.test[-which(p.test == "+")]
            if (length(which(p.test == "P")) > 0) {
                p.test = p.test[-(1:grep("P", p.test)[2])]
                edges = grep("[", p.test, fixed = TRUE, value = TRUE)
                p.values = NULL
                if (length(grep("[", p.test, fixed = TRUE)) > 
                  1) 
                  for (i in 1:(length(grep("[", p.test, fixed = TRUE)) - 
                    1)) {
                    p.values = c(p.values, as.numeric(p.test[grep("[", 
                      p.test, fixed = TRUE)[i + 1] - 1]))
                  }
                p.values = c(p.values, as.numeric(p.test[length(p.test)]))
                for (i in 1:length(edges)) {
                  e = strsplit(edges[i], "")[[1]]
                  p.matrix[which(letters == e[2]), which(letters == 
                    e[3])] = p.values[i]
                }
            }
            p.list[[accepted[ll]]] = p.matrix
        }
    if (tests == TRUE) 
        res = list(accepted = accepted, rejected = rejected, 
            base = model, strategy = strategy, `tests against saturated` = sat.tests, 
            `p values` = p.list, `variable names` = var.names)
    else res = list(accepted = accepted, rejected = rejected, 
        base = model, strategy = strategy, `variable names` = var.names)
    res
}
