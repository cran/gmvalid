`gm.boot.coco` <-
function (N, data, strategy = c("backwards", "forwards", "eh"), 
    calculations = c("subgraph", "diff", "edge", "clique"), model = FALSE, 
    ...) 
{
    require(CoCo, quietly = TRUE)
    result = NULL
    strategy = match.arg(strategy)
    calculations = match.arg(calculations, several.ok = TRUE)
    if (!is.array(data)) 
        for (i in 1:dim(data)[2]) data[[names(data)[i]]] = as.factor(data[[names(data)[i]]])
    else {
        require(epitools, quietly = TRUE)
        if (!length(dimnames(data))) 
            for (i in length(dim(data)):1) dimnames(data)[[i]] = as.character(c(1:dim(data)[i]))
        data = data.frame(expand.table(data))
    }
    var.names = cbind(letters[1:dim(data)[2]], dimnames(data)[2][[1]])
    dimnames(data)[2][[1]] = letters[1:dim(data)[2]]
    length.data = dim(data)[1]
    elements = dimnames(data)[[2]]
    CoCoObject = makeCoCo(silent = TRUE)
    if (strategy == "backwards") {
        for (i in 1:N) {
            if (N > 3 && i%%(round(N/4)) == 1) {
                cat(paste(c("Run number ", i, "\n"), collapse = ""))
                flush.console()
            }
            capture.output(enterTable(table(data[sample(length.data, 
                length.data, replace = TRUE), ]), object = CoCoObject))
            if (model != FALSE) 
                enterModel(model)
            else enterModel("*")
            capture.output(backward(object = CoCoObject, ...))
            m.name = returnModel("last", split.string = TRUE)
            m.name = paste(sort(m.name), collapse = ",")
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
            capture.output(enterTable(table(data[sample(length.data, 
                length.data, replace = TRUE), ]), object = CoCoObject))
            if (model != FALSE) 
                enterModel(model)
            else enterModel(".")
            capture.output(forward(object = CoCoObject, ...))
            m.name = returnModel("last", split.string = TRUE)
            m.name = paste(sort(m.name), collapse = ",")
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
            capture.output(enterTable(table(data[sample(length.data, 
                length.data, replace = TRUE), ]), object = CoCoObject))
            if (length(model) > 1) {
                ehForceFix(model[1], fix = "in")
                ehForceFix("*", fix = "out")
                ehSetBase(model[2])
            }
            else if (model != FALSE) 
                ehSetBase(model)
            else
                ehSetBase("*")
            capture.output(eh(object = CoCoObject, ...))
            ehExtract(object = CoCoObject, "accepted")
            makeCurrent(1)
            m.name = returnModel("current", split.string = TRUE)
            m.name = paste(sort(m.name), collapse = ",")
            makeCurrent("next")
            tmp.check = returnModel("last", split.string = TRUE)
            tmp.check = paste(sort(tmp.check), collapse = ",")
            while (m.name[length(m.name)] != tmp.check) {
                mod = returnModel("current", split.string = TRUE)
                m.name = c(m.name, paste(sort(mod), collapse = ","))
                makeCurrent("next")
            }
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
    endCoCo(silent = TRUE)
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
        dot.flag = 0
        for (h in 1:length(result)) {
            m = strsplit(names(result)[h], "")[[1]]
            if (any(m == ":")) {
                dot.flag = 1
                m = strsplit(names(result)[h], ",")[[1]]
                m = strsplit(m, ":")
                for (i in 1:length(m)) for (j in 1:length(m[[i]])) m[[i]][j] = paste(":", 
                  m[[i]][j], sep = "")
                mm = m[[1]][-which(m[[1]] == ":")]
                if (length(m) > 1) 
                  for (i in 2:length(m)) mm = c(mm, ",", m[[i]][-which(m[[i]] == 
                    ":")])
                m = mm
                elements = unique(m)
                elements = elements[-which(elements == ",")]
            }
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
                    pos = which(elements == m[i + j])
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
                if (dot.flag == 1) 
                  w = which(.gm.power.set(strsplit(names(result.subgraph)[j], 
                    ":")[[1]], NULL, dot = TRUE) == names(result.subgraph)[i])
                else w = which(.gm.power.set(strsplit(names(result.subgraph)[j], 
                  "")[[1]], NULL, dot = FALSE) == names(result.subgraph)[i])
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
    elements = dimnames(data)[[2]]
    if (length(which(calculations == "edge"))) {
        dep.table = matrix(0, nrow = dim(data)[2], ncol = dim(data)[2])
        for (k in 1:length(result)) {
            m = strsplit(names(result)[k], "")[[1]]
            if (any(m == ":")) {
                m = strsplit(names(result)[k], ",")[[1]]
                m = strsplit(m, ":")
                mm = m[[1]][-(m[[1]] == "")]
                if (length(m) > 1) 
                  for (i in 2:length(m)) mm = c(mm, ",", m[[i]][-(m[[i]] == 
                    "")])
                m = mm
            }
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
        CoCoObject = makeCoCo(silent = TRUE)
        capture.output(enterDataFrame(data))
        if (strategy == "backwards") {
            if (model != FALSE) 
                enterModel(model)
            else enterModel("*")
            capture.output(backward(object = CoCoObject, ...))
        }
        else if (strategy == "forwards") {
            if (model != FALSE) 
                enterModel(model)
            else enterModel(".")
            capture.output(forward(object = CoCoObject, ...))
        }
        else if (strategy == "eh") {
            capture.output(eh(object = CoCoObject, ...))
            ehExtract("accepted")
        }
        m.name = returnModel("last", split.string = TRUE)
        m.name = paste(m.name, collapse = ",")
        result.original = m.name
        endCoCo(silent = TRUE)
        m = strsplit(m.name, "")[[1]]
        if (any(m == ":")) {
            m = strsplit(m.name, ",")[[1]]
            m = strsplit(m, ":")
            mm = m[[1]][-(m[[1]] == "")]
            if (length(m) > 1) 
                for (i in 2:length(m)) mm = c(mm, ",", m[[i]][-(m[[i]] == 
                  "")])
            m = mm
        }
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
            if (any(m == ":")) {
                m = strsplit(names(result)[k], ",")[[1]]
                m = strsplit(m, ":")
                mm = m[[1]][-(m[[1]] == "")]
                if (length(m) > 1) 
                  for (i in 2:length(m)) mm = c(mm, ",", m[[i]][-(m[[i]] == 
                    "")])
                m = mm
            }
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