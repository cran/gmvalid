`gm.coco` <-
function (data, strategy = c("backwards", "forwards", "eh"), 
    model = FALSE, ...) 
{
    require(CoCo, quietly = TRUE)
    strategy = match.arg(strategy)
    result.base = NA
    result.rejec = NULL
    result.accep = NA
    if (is.array(data)) {
        require(epitools, quietly = TRUE)
        if (!length(dimnames(data))) 
            for (i in length(dim(data)):1) dimnames(data)[[i]] = as.character(c(1:dim(data)[i]))
        data = data.frame(expand.table(data))
    }
    var.names = cbind(letters[1:dim(data)[2]], dimnames(data)[2][[1]])
    dimnames(data)[2][[1]] = letters[1:dim(data)[2]]
    for (i in 1:dim(data)[2]) data[[names(data)[i]]] = as.factor(data[[names(data)[i]]])
    CoCoObject = makeCoCo()
    enterDataFrame(data)
    if (strategy == "backwards") {
        if (model != FALSE) 
            enterModel(model)
        else enterModel("*")
        backward(object = CoCoObject, ...)
        result.base = returnModel("base", split.string = TRUE)
        result.base = sort(result.base)
        result.base = paste(result.base, collapse = ",")
        makeCurrent("next")
        if (returnModel("current") != returnModel("last")) {
            result.rejec = returnModel("current", split.string = TRUE)
            result.rejec = sort(result.rejec)
            result.rejec = paste(result.rejec, collapse = ",")
            makeCurrent("next")
        }
        result.accep = returnModel("last", split.string = TRUE)
        result.accep = sort(result.accep)
        result.accep = paste(result.accep, collapse = ",")
        while (TRUE) {
            mod = returnModel("current", split.string = TRUE)
            mod = sort(mod)
            mod = paste(mod, collapse = ",")
            if (mod == result.accep) 
                break
            result.rejec = c(result.rejec, mod)
            makeCurrent("next")
        }
    }
    else if (strategy == "forwards") {
        if (model != FALSE) 
            enterModel(model)
        else enterModel(".")
        forward(object = CoCoObject, ...)
        result.base = returnModel("base", split.string = TRUE)
        result.base = sort(result.base)
        result.base = paste(result.base, collapse = ",")
        makeCurrent("next")
        if (returnModel("current") != returnModel("last")) {
            result.rejec = returnModel("current", split.string = TRUE)
            result.rejec = sort(result.rejec)
            result.rejec = paste(result.rejec, collapse = ",")
            makeCurrent("next")
        }
        result.accep = returnModel("last", split.string = TRUE)
        result.accep = sort(result.accep)
        result.accep = paste(result.accep, collapse = ",")
        while (TRUE) {
            mod = returnModel("current", split.string = TRUE)
            mod = sort(mod)
            mod = paste(mod, collapse = ",")
            if (mod == result.accep) 
                break
            result.rejec = c(result.rejec, mod)
            makeCurrent("next")
        }
    }
    else if (strategy == "eh") {
        if (length(model) > 1) {
            ehForceFix(model[1], fix = "in")
            ehForceFix("*", fix = "out")
            ehSetBase(model[2])
        }
        else if (model != FALSE) 
            ehSetBase(model)
        else
            ehSetBase("*")
        result.base = returnModel("base",split.string=TRUE)
        result.base = sort(result.base)
        result.base = paste(result.base,collapse=",")
        eh(object = CoCoObject, ...)
        ehExtract("accepted")
        makeCurrent(1)
        result.accep = returnModel("current", split.string = TRUE)
        result.accep = sort(result.accep)
        result.accep = paste(result.accep, collapse = ",")
        makeCurrent("next")
        tmp.check = returnModel("last", split.string = TRUE)
        tmp.check = sort(tmp.check)
        tmp.check = paste(tmp.check, collapse = ",")
        while (result.accep[length(result.accep)] != tmp.check) {
            mod = returnModel("current", split.string = TRUE)
            mod = sort(mod)
            mod = paste(mod, collapse = ",")
            result.accep = c(result.accep, mod)
            makeCurrent("next")
        }
        ehExtract("rejected")
        makeCurrent("next")
        result.rejec = returnModel("current", split.string = TRUE)
        result.rejec = sort(result.rejec)
        result.rejec = paste(result.rejec, collapse = ",")
        makeCurrent("next")
        tmp.check = returnModel("last", split.string = TRUE)
        tmp.check = sort(tmp.check)
        tmp.check = paste(tmp.check, collapse = ",")
        while (result.rejec[length(result.rejec)] != tmp.check) {
            mod = returnModel("current", split.string = TRUE)
            mod = sort(mod)
            mod = paste(mod, collapse = ",")
            result.rejec = c(result.rejec, mod)
            makeCurrent("next")
        }
    }
    else stop("No valid strategy!")
    tests = 0
    if (length(result.rejec)) {
        tests = matrix(0, nrow = length(result.accep) * length(result.rejec), 
            ncol = 3)
        dimnames(tests) = list(c(1:(length(result.accep) * length(result.rejec))), 
            c("df", "deviance", "p.value"))
        index = 1
        for (i in 1:length(result.accep)) for (j in 1:length(result.rejec)) if (isSubmodel(model.1 = result.accep[i], 
            model.2 = result.rejec[j])) {
            tests[index, ] = returnDeviance(model.1 = result.accep[i], 
                model.2 = result.rejec[j])[c(3, 11, 12)]
            dimnames(tests)[[1]][index] = paste(c(result.accep[i], 
                " ~ ", result.rejec[j]), collapse = "")
            index = index + 1
        }
        else if (isSubmodel(model.1 = result.rejec[j], model.2 = result.accep[i])) {
            tests[index, ] = returnDeviance(model.1 = result.rejec[j], 
                model.2 = result.accep[i])[c(3, 11, 12)]
            dimnames(tests)[[1]][index] = paste(c(result.accep[i], 
                " ~ ", result.rejec[j]), collapse = "")
            index = index + 1
        }
    }
    if (all(tests == 0) && all(result.accep != returnModel("*", 
        split.string = TRUE))) {
        tests = matrix(0, nrow = length(result.accep) + length(result.rejec), 
            ncol = 3)
        dimnames(tests) = list(c(1:(length(result.accep) + length(result.rejec))), 
            c("df", "deviance", "p.value"))
        index = 1
        for (i in 1:length(result.accep)) {
            tests[index, ] = returnDeviance(model.1 = result.accep[i], 
                model.2 = "*")[c(3, 11, 12)]
            dimnames(tests)[[1]][index] = paste(c(result.accep[i], 
                " ~ SATURATED"), collapse = "")
            index = index + 1
        }
        if (length(result.rejec)) 
            for (i in 1:length(result.rejec)) {
                tests[index, ] = returnDeviance(model.1 = result.rejec[i], 
                  model.2 = "*")[c(3, 11, 12)]
                dimnames(tests)[[1]][index] = paste(c(result.rejec[i], 
                  " ~ SATURATED"), collapse = "")
                index = index + 1
            }
    }
    endCoCo()
    list(accepted = result.accep, rejected = result.rejec, base = result.base, 
        strategy = strategy, tests = tests, "variable names" = var.names)
}
