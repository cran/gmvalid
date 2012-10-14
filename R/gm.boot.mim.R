gm.boot.mim <-
function (N, data, strategy = c("backwards", "forwards", "eh", 
    "combined"), calculations = c("diff", "edge", "clique"), 
    model = FALSE, options = "") 
{
    require(mimR, quietly = TRUE)
    result.model = NULL
    result.freq = NULL
    strategy = match.arg(strategy)
    calculations = match.arg(calculations, several.ok = TRUE)
    if (is.array(data)) {
        require(epitools, quietly = TRUE)
        if (!length(dimnames(data))) 
            for (i in length(dim(data)):1) dimnames(data)[[i]] = as.character(c(1:dim(data)[i]))
        data = data.frame(expand.table(data))
    }
    var.names = cbind(letters[1:dim(data)[2]], dimnames(data)[2][[1]])
    dimnames(data)[2][[1]] = letters[1:dim(data)[2]]
    length.data = dim(data)[1]
    elements = dimnames(data)[[2]]
    for (i in 1:length(elements)) if (any(LETTERS == elements[i])) 
        elements[i] = letters[which(LETTERS == elements[i])]
    for (i in 1:dim(data)[2]) data[[names(data)[i]]] = as.factor(data[[names(data)[i]]])
    if (strategy == "backwards") {
        for (i in 1:N) {
            if (N > 3 && i%%(round(N/4)) == 1) {
                cat(paste(c("Run number ", i, "\n"), collapse = ""))
                flush.console()
            }
            toMIM(as.gmData(data[sample(length.data, length.data, 
                replace = TRUE), ]))
            if (model != FALSE) 
                mim.cmd(paste(c("Model ", model), collapse = ""))
            else {
                mim.cmd("satmod")
            }
            mim.text = mim.cmd(paste(c("stepwise ", options), 
                collapse = ""), look.nice = FALSE)
            m.name = toString(mim.text[length(mim.text)])
            m.name = strsplit(m.name, ",")[[1]]
            m.name = sort(m.name)
            m.name = paste(m.name, collapse = ",")
            if (!any(result.model == m.name)) {
                result.model = c(result.model, m.name)
                result.freq = c(result.freq, 1)
            }
            else {
                index = which(result.model == m.name)
                result.freq[index] = result.freq[index] + 1
            }
        }
    }
    else if (strategy == "forwards") {
        for (i in 1:N) {
            if (N > 3 && i%%(round(N/4)) == 1) {
                print(paste(c("Run number ", i), collapse = ""))
                flush.console()
            }
            toMIM(as.gmData(data[sample(length.data, length.data, 
                replace = TRUE), ]))
            if (model != FALSE) 
                mim.cmd(paste(c("Model ", model), collapse = ""))
            else {
                mim.cmd("maine")
            }
            mim.text = mim.cmd(paste(c("stepwise f", options), 
                collapse = ""), look.nice = FALSE)
            m.name = toString(mim.text[length(mim.text)])
            m.name = strsplit(m.name, ",")[[1]]
            m.name = sort(m.name)
            m.name = paste(m.name, collapse = ",")
            if (!any(result.model == m.name)) {
                result.model = c(result.model, m.name)
                result.freq = c(result.freq, 1)
            }
            else {
                index = which(result.model == m.name)
                result.freq[index] = result.freq[index] + 1
            }
        }
    }
    else if (strategy == "eh") {
        for (i in 1:N) {
            if (N > 3 && i%%(round(N/4)) == 1) {
                print(paste(c("Run number ", i), collapse = ""))
                flush.console()
            }
            toMIM(as.gmData(data[sample(length.data, length.data, 
                replace = TRUE), ]))
            if (model == FALSE) 
                model = paste(c(paste(letters[1:dim(data)[2]], 
                  collapse = ","), " - ", paste(letters[1:dim(data)[2]], 
                  collapse = "")), collapse = "")
            mim.text = mim.cmd(paste(c("initsearch ", model, 
                ";startsearch ", options), collapse = ""), look.nice = FALSE)
            m.name = mim.cmd("ehshow a", look.nice = FALSE, return.look.nice = TRUE)
            m.name = m.name[1:length(m.name)%%2 == 0]
            for (j in 1:length(m.name)) if (!any(result.model == 
                m.name)) {
                result.model = c(result.model, m.name)
                result.freq = c(result.freq, 1)
            }
            else {
                index = which(result.model == m.name)
                result.freq[index] = result.freq[index] + 1
            }
        }
    }
    else if (strategy == "combined") {
        for (i in 1:N) {
            if (N > 3 && i%%(round(N/4)) == 1) {
                print(paste(c("Run number ", i), collapse = ""))
                flush.console()
            }
            sam = data[sample(length.data, length.data, replace = TRUE), 
                ]
            toMIM(as.gmData(sam))
            mim.cmd(paste(c("Model ", gm.screening(sam)$model), 
                collapse = ""))
            mim.text = mim.cmd(paste(c("stepwise j", options), 
                collapse = ""), look.nice = FALSE)
            m.name = toString(mim.text[length(mim.text)])
            m.name = strsplit(m.name, ",")[[1]]
            m.name = sort(m.name)
            m.name = paste(m.name, collapse = ",")
            if (!any(result.model == m.name)) {
                result.model = c(result.model, m.name)
                result.freq = c(result.freq, 1)
            }
            else {
                index = which(result.model == m.name)
                result.freq[index] = result.freq[index] + 1
            }
        }
    }
    result = cbind(result.model, result.freq)
    if (length(which(calculations == "clique"))) {
        result.clique = NULL
        for (h in 1:length(result)) {
            m = strsplit(result[h, 1], ",")[[1]]
            for (hh in 1:length(m)) if (!any(names(result.clique) == 
                m[hh])) {
                eval(parse(text = paste("result.clique = c(result.clique,\"", 
                  m[hh], "\" = ", result[h, 2], ")", sep = "")))
            }
            else {
                index = which(names(result.clique) == m[hh])
                result.clique[index] = result.clique[index] + 
                  as.numeric(result[h, 2])
            }
        }
    }
    if (length(which(calculations == "edge"))) {
        dep.table = matrix(0, nrow = dim(data)[2], ncol = dim(data)[2])
        for (k in 1:length(result)) {
            dep.table = dep.table + gm.matrixparse(result[k, 
                1]) * as.numeric(result[k, 2])
        }
        dimnames(dep.table) = list(elements, elements)
    }
    if (strategy != "eh" && length(which(calculations == "diff"))) {
        toMIM(as.gmData(data))
        if (strategy == "backwards") {
            mim.text = mim.cmd(paste(c("satmod; stepwise ", options), 
                collapse = ""), look.nice = FALSE)
            m.name = toString(mim.text[length(mim.text)])
        }
        else if (strategy == "forwards") {
            mim.text = mim.cmd(paste(c("maine; stepwise f", options), 
                collapse = ""), look.nice = FALSE)
            m.name = toString(mim.text[length(mim.text)])
        }
        else if (strategy == "eh") {
            mim.text = mim.cmd(paste(c("initsearch; startsearch ", 
                options), collapse = ""), look.nice = FALSE)
            m.name = mim.cmd("ehshow a", look.nice = FALSE, return.look.nice = TRUE)
            m.name = m.name[1:length(m.name)%%2 == 0]
        }
        result.original = m.name
        original = gm.matrixparse(m.name)
        list.result = list(more = c(NULL), less = c(NULL), abs = c(NULL))
        for (k in 1:length(result)) {
            booted = gm.matrixparse(result[k, 1])
            more = as.character(sum((booted - original)[(booted - 
                original) > 0]))
            less = as.character(sum((original - booted)[(original - 
                booted) > 0]))
            diff = as.character(sum(abs(original - booted)))
            if (length(which(names(list.result[[1]]) == more)) == 
                0) {
                tmp = names(list.result[[1]])
                list.result[[1]] = c(list.result[[1]], as.numeric(result[k, 
                  2]))
                names(list.result[[1]]) = c(tmp, more)
            }
            else list.result[[1]][which(names(list.result[[1]]) == 
                more)] = list.result[[1]][which(names(list.result[[1]]) == 
                more)] + as.numeric(result[k, 2])
            if (length(which(names(list.result[[2]]) == less)) == 
                0) {
                tmp = names(list.result[[2]])
                list.result[[2]] = c(list.result[[2]], as.numeric(result[k, 
                  2]))
                names(list.result[[2]]) = c(tmp, less)
            }
            else list.result[[2]][which(names(list.result[[2]]) == 
                less)] = list.result[[2]][which(names(list.result[[2]]) == 
                less)] + as.numeric(result[k, 2])
            if (length(which(names(list.result[[3]]) == diff)) == 
                0) {
                tmp = names(list.result[[3]])
                list.result[[3]] = c(list.result[[3]], as.numeric(result[k, 
                  2]))
                names(list.result[[3]]) = c(tmp, diff)
            }
            else list.result[[3]][which(names(list.result[[3]]) == 
                diff)] = list.result[[3]][which(names(list.result[[3]]) == 
                diff)] + as.numeric(result[k, 2])
        }
    }
    summie = sum(as.numeric(result[, 2]))
    result[, 2] = as.numeric(result[, 2])/summie
    if (length(result.model) > 1) 
        result = list(`bootstrapped models` = result[order(result[, 
            2], decreasing = TRUE), ])
    else result = list(`bootstrapped models` = result)
    if (length(which(calculations == "clique"))) 
        result$"bootstrapped cliques" = sort(result.clique/summie, 
            decreasing = TRUE)
    if (length(which(calculations == "edge"))) 
        result$"edge frequencies" = dep.table/summie
    if (strategy != "eh" && length(which(calculations == "diff"))) {
        result$"original model" = result.original
        result$"edge differences" = lapply(list.result, sort, 
            decreasing = TRUE)
    }
    result$replications = N
    result$"variable names" = var.names
    result
}
