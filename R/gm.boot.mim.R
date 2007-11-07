`gm.boot.mim` <-
function (N, data, strategy = c("backwards", "forwards", "eh"), 
    calculations = c("subgraph", "diff", "edge", "clique"), model = FALSE, 
    options = "") 
{
    require(mimR, quietly = TRUE)
    result = NULL
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
            if (!any(names(result) == m.name)) {
                eval(parse(text = paste("result = c(result,\"", 
                  m.name, "\" = 1)", sep = "")))
            }
            else {
                index = which(names(result) == m.name)
                result[index] = result[index] + 1
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
            if (!any(names(result) == m.name)) {
                eval(parse(text = paste("result = c(result,\"", 
                  m.name, "\" = 1)", sep = "")))
            }
            else {
                index = which(names(result) == m.name)
                result[index] = result[index] + 1
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
            for (j in 1:length(m.name)) if (!any(names(result) == 
                m.name[j])) {
                eval(parse(text = paste("result = c(result,\"", 
                  m.name[j], "\" = 1)", sep = "")))
            }
            else {
                index = which(names(result) == m.name[j])
                result[index] = result[index] + 1
            }
        }
    }
    if (length(which(calculations == "clique"))) {
        result.clique = NULL
        for (h in 1:length(result)) {
            m = strsplit(names(result)[h], ",")[[1]]
            for (hh in 1:length(m)) if (!any(names(result.clique) == 
                m[hh])) {
                eval(parse(text = paste("result.clique = c(result.clique,\"", 
                  m[hh], "\" = ", result[h], ")", sep = "")))
            }
            else {
                index = which(names(result.clique) == m[hh])
                result.clique[index] = result.clique[index] + 
                  result[h]
            }
        }
    }
    if (length(which(calculations == "subgraph"))) {
        result.subgraph = NULL
        for (h in 1:length(result)) {
            m = strsplit(names(result)[h], "")[[1]]
            dep_table = NULL
            i = 1
            while (i <= length(m)) {
                if (length(which(elements == m[i])) == 0) {
                  i = i + 1
                }
                else {
                  j = 0
                  clique = NULL
                  while (i + j <= length(m) && length(which(elements == 
                    m[i + j])) > 0) {
                    pos = which(elements == m[i + j])%%27
                    clique = c(clique, pos)
                    j = j + 1
                  }
                  i = i + j
                  dep_table = c(dep_table, 0, clique)
                }
            }
            dep_table = c(dep_table, 0)
            zeros = which(dep_table == 0)
            zeros = zeros - 1
            model.result = NULL
            for (k in 1:(length(zeros) - 1)) {
                check.result = model.result
                ps = .gm.power.set(m[(zeros[k] + 1):(zeros[k + 
                  1] - 1)], NULL)
                vv = array(result[h], dim = length(ps))
                names(vv) = ps
                vv = vv[-which(names(vv) == "")]
                for (j in 1:length(vv)) {
                  if (length(which(names(model.result) == names(vv)[j])) == 
                    0) {
                    tmp = names(model.result)
                    model.result = c(model.result, vv[j])
                    names(model.result) = c(tmp, names(vv)[j])
                  }
                }
            }
            check.result.3 = result.subgraph
            for (j in 1:length(model.result)) {
                if (length(which(names(result.subgraph) == names(model.result)[j])) == 
                  0) {
                  tmp = names(result.subgraph)
                  result.subgraph = c(result.subgraph, model.result[j])
                  names(result.subgraph) = c(tmp, names(model.result)[j])
                }
                else if (check.result.3[which(names(result.subgraph) == 
                  names(model.result)[j])] == result.subgraph[which(names(result.subgraph) == 
                  names(model.result)[j])]) 
                  result.subgraph[which(names(result.subgraph) == 
                    names(model.result)[j])] = result.subgraph[which(names(result.subgraph) == 
                    names(model.result)[j])] + model.result[j]
            }
        }
        i = 1
        while (i <= length(result.subgraph)) {
            j = 1
            length.2 = length(result.subgraph)
            while (j <= length(result.subgraph)) {
                w = which(.gm.power.set(strsplit(names(result.subgraph)[j], 
                  "")[[1]], NULL) == names(result.subgraph)[i])
                if (length(w) > 0 && result.subgraph[i] == result.subgraph[j] && 
                  i != j) {
                  result.subgraph = result.subgraph[-i]
                }
                j = j + 1
            }
            if (length.2 == length(result.subgraph)) 
                i = i + 1
        }
    }
    if (length(which(calculations == "edge"))) {
        dep.table = matrix(0, nrow = dim(data)[2], ncol = dim(data)[2])
        for (k in 1:length(result)) {
            m = strsplit(names(result)[k], "")[[1]]
            dep.table.i = matrix(0, nrow = dim(data)[2], ncol = dim(data)[2])
            for (i in 1:length(m)) {
                if (length(which(elements == m[i])) == 0) {
                }
                else {
                  j = 1
                  while (i + j <= length(m) && length(which(elements == 
                    m[i + j])) > 0) {
                    pos1 = which(elements == m[i])
                    pos2 = which(elements == m[i + j])
                    dep.table.i[pos1, pos2] = result[k]
                    j = j + 1
                  }
                }
            }
            dep.table = dep.table + dep.table.i
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
        m = strsplit(m.name, "")[[1]]
        original = matrix(0, nrow = dim(data)[2], ncol = dim(data)[2])
        for (i in 1:length(m)) {
            if (length(which(elements == m[i])) == 0) {
            }
            else {
                j = 1
                while (i + j <= length(m) && length(which(elements == 
                  m[i + j])) > 0) {
                  pos1 = which(elements == m[i])
                  pos2 = which(elements == m[i + j])
                  original[pos1, pos2] = 1
                  j = j + 1
                }
            }
        }
        list.result = list("more edges" = c(NULL), "less edges" = c(NULL), 
            differences = c(NULL))
        for (k in 1:length(result)) {
            m = strsplit(names(result)[k], "")[[1]]
            booted = matrix(0, nrow = dim(data)[2], ncol = dim(data)[2])
            for (i in 1:length(m)) {
                if (length(which(elements == m[i])) == 0) {
                }
                else {
                  j = 1
                  while (i + j <= length(m) && length(which(elements == 
                    m[i + j])) > 0) {
                    pos1 = which(elements == m[i])
                    pos2 = which(elements == m[i + j])
                    booted[pos1, pos2] = 1
                    j = j + 1
                  }
                }
            }
            more = as.character(sum((booted - original)[(booted - 
                original) > 0]))
            less = as.character(sum((original - booted)[(original - 
                booted) > 0]))
            diff = as.character(sum(abs(original - booted)))
            if (length(which(names(list.result[[1]]) == more)) == 
                0) {
                tmp = names(list.result[[1]])
                list.result[[1]] = c(list.result[[1]], result[k])
                names(list.result[[1]]) = c(tmp, more)
            }
            else list.result[[1]][which(names(list.result[[1]]) == 
                more)] = list.result[[1]][which(names(list.result[[1]]) == 
                more)] + result[k]
            if (length(which(names(list.result[[2]]) == less)) == 
                0) {
                tmp = names(list.result[[2]])
                list.result[[2]] = c(list.result[[2]], result[k])
                names(list.result[[2]]) = c(tmp, less)
            }
            else list.result[[2]][which(names(list.result[[2]]) == 
                less)] = list.result[[2]][which(names(list.result[[2]]) == 
                less)] + result[k]
            if (length(which(names(list.result[[3]]) == diff)) == 
                0) {
                tmp = names(list.result[[3]])
                list.result[[3]] = c(list.result[[3]], result[k])
                names(list.result[[3]]) = c(tmp, diff)
            }
            else list.result[[3]][which(names(list.result[[3]]) == 
                diff)] = list.result[[3]][which(names(list.result[[3]]) == 
                diff)] + result[k]
        }
    }
    summie = sum(result)
    result = list("bootstrapped models" = sort(result/summie, 
        decreasing = TRUE))
    if (length(which(calculations == "subgraph"))) 
        result$"bootstrapped subgraphs" = sort(result.subgraph/summie, 
            decreasing = TRUE)
    if (length(which(calculations == "clique"))) 
        result$"bootstrapped cliques" = sort(result.clique/summie, 
            decreasing = TRUE)
    if (length(which(calculations == "edge"))) 
        result$"bootstrapped edges" = dep.table/summie
    if (strategy != "eh" && length(which(calculations == "diff"))) {
        result$"original model" = result.original
        result$"differences from original data set" = lapply(list.result, 
            sort, decreasing = TRUE)
    }
    result$"variable names" = var.names
    result
}
