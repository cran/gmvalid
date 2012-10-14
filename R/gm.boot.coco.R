gm.boot.coco <-
function (N, data, strategy = c("backwards", "forwards", "combined"), 
    calculations = c("diff", "edge", "clique"), model = FALSE, 
    criterion = c("lr", "aic", "bic"), ...) 
{
    require(CoCo, quietly = TRUE)
    result.model = NULL
    result.freq = NULL
    strategy = match.arg(strategy)
    calculations = match.arg(calculations, several.ok = TRUE)
    criterion = match.arg(criterion)
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
    CoCoObject = makeCoCo()
    if (strategy == "backwards") {
        for (i in 1:N) {
            if (N > 3 && i%%(round(N/4)) == 1) {
                cat(paste(c("Run number ", i, "\n"), collapse = ""))
                flush.console()
            }
            capture.output(enterTable(table(data[sample(length.data, 
                length.data, replace = TRUE), ]), object = CoCoObject))
            if (criterion == "bic") 
                capture.output(optionsCoCo(ic = TRUE, bic = TRUE, 
                  object = CoCoObject))
            else if (criterion == "aic") 
                capture.output(optionsCoCo(ic = TRUE, object = CoCoObject))
            if (model != FALSE) 
                enterModel(model)
            else enterModel("*")
            capture.output(backward(object = CoCoObject, ...))
            m.name = returnModel("last", split.string = TRUE)
            m.name = paste(sort(m.name), collapse = ",")
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
            capture.output(enterTable(table(data[sample(length.data, 
                length.data, replace = TRUE), ]), object = CoCoObject))
            if (criterion == "bic") 
                capture.output(optionsCoCo(ic = TRUE, bic = TRUE, 
                  object = CoCoObject))
            else if (criterion == "aic") 
                capture.output(optionsCoCo(ic = TRUE, object = CoCoObject))
            if (model != FALSE) 
                enterModel(model)
            else enterModel(".")
            capture.output(forward(object = CoCoObject, ...))
            m.name = returnModel("last", split.string = TRUE)
            m.name = paste(sort(m.name), collapse = ",")
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
    else if (strategy == "combined") {
        for (i in 1:N) {
            if (N > 3 && i%%(round(N/4)) == 1) {
                print(paste(c("Run number ", i), collapse = ""))
                flush.console()
            }
            sam = data[sample(length.data, length.data, replace = TRUE), 
                ]
            capture.output(enterTable(table(sam), object = CoCoObject))
            if (criterion == "bic") 
                capture.output(optionsCoCo(ic = TRUE, bic = TRUE, 
                  object = CoCoObject))
            else if (criterion == "aic") 
                capture.output(optionsCoCo(ic = TRUE, object = CoCoObject))
            enterModel(gm.screening(sam)$model)
            capture.output(backward(object = CoCoObject, ...))
            makeBase("last")
            capture.output(forward(object = CoCoObject, ...))
            m.name = returnModel("last", split.string = TRUE)
            m.name = paste(sort(m.name), collapse = ",")
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
    endCoCo()
    result = cbind(result.model, result.freq)
    if (length(which(calculations == "clique"))) {
        result.clique = NULL
        for (h in 1:length(result.model)) {
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
        for (k in 1:length(result.model)) {
            dep.table = dep.table + gm.matrixparse(result[k, 
                1]) * as.numeric(result[k, 2])
        }
        dimnames(dep.table) = list(elements, elements)
    }
    if (length(which(calculations == "diff"))) {
        CoCoObject = makeCoCo()
        capture.output(enterDataFrame(data))
        if (criterion == "bic") 
            capture.output(optionsCoCo(ic = TRUE, bic = TRUE, 
                object = CoCoObject))
        else if (criterion == "aic") 
            capture.output(optionsCoCo(ic = TRUE, object = CoCoObject))
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
        else if (strategy == "combined") {
            enterModel(gm.screening(data)$model)
            capture.output(backward(object = CoCoObject, ...))
        }
        m.name = returnModel("last", split.string = TRUE)
        m.name = paste(m.name, collapse = ",")
        result.original = m.name
        endCoCo()
        original = gm.matrixparse(m.name)
        list.result = list(more = c(NULL), less = c(NULL), abs = c(NULL))
        for (k in 1:length(result.model)) {
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
    if (length(which(calculations == "diff"))) {
        result$"original model" = result.original
        result$"edge differences" = lapply(list.result, sort, 
            decreasing = TRUE)
    }
    result$replications = N
    result$"variable names" = var.names
    result
}
