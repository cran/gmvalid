gm.coco <-
function (data, strategy = c("backwards", "forwards", "eh", "combined"), 
    model = FALSE, eh.strategy = c("smallest", "alternating", 
        "rough"), criterion = c("lr", "aic", "bic"), tests = FALSE, 
    ...) 
{
    require(CoCo, quietly = TRUE)
    strategy = match.arg(strategy)
    eh.strategy = match.arg(eh.strategy)
    criterion = match.arg(criterion)
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
    if (criterion == "bic") 
        optionsCoCo(ic = TRUE, bic = TRUE, object = CoCoObject)
    else if (criterion == "aic") 
        optionsCoCo(ic = TRUE, object = CoCoObject)
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
        else ehSetBase("*")
        result.base = returnModel("base", split.string = TRUE)
        result.base = sort(result.base)
        result.base = paste(result.base, collapse = ",")
        eh(object = CoCoObject, strategy = eh.strategy, ...)
        ehExtract("accepted")
        makeCurrent("next")
        result.accep = returnModel("current", split.string = TRUE)
        result.accep = sort(result.accep)
        result.accep = paste(result.accep, collapse = ",")
        tmp.check = returnModel("last", split.string = TRUE)
        tmp.check = sort(tmp.check)
        tmp.check = paste(tmp.check, collapse = ",")
        while (result.accep[length(result.accep)] != tmp.check) {
            makeCurrent("next")
            mod = returnModel("current", split.string = TRUE)
            mod = sort(mod)
            mod = paste(mod, collapse = ",")
            result.accep = c(result.accep, mod)
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
    else if (strategy == "combined") {
        result.base = gm.screening(data)$model
        enterModel(result.base)
        backward(object = CoCoObject, ...)
        makeBase("last")
        forward(object = CoCoObject, ...)
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
    else stop("No valid strategy!")
    model.tests = 0
    if (length(result.rejec) && tests == TRUE) {
        model.tests = matrix(0, nrow = length(result.accep) * 
            length(result.rejec), ncol = 3)
        dimnames(model.tests) = list(c(1:(length(result.accep) * 
            length(result.rejec))), c("df", "deviance", "p.value"))
        index = 1
        for (i in 1:length(result.accep)) for (j in 1:length(result.rejec)) if (isSubmodel(model.1 = result.accep[i], 
            model.2 = result.rejec[j])) {
            model.tests[index, ] = returnDeviance(model.1 = result.accep[i], 
                model.2 = result.rejec[j])[c(3, 11, 12)]
            dimnames(model.tests)[[1]][index] = paste(c(result.accep[i], 
                " ~ ", result.rejec[j]), collapse = "")
            index = index + 1
        }
        else if (isSubmodel(model.1 = result.rejec[j], model.2 = result.accep[i])) {
            model.tests[index, ] = returnDeviance(model.1 = result.rejec[j], 
                model.2 = result.accep[i])[c(3, 11, 12)]
            dimnames(model.tests)[[1]][index] = paste(c(result.accep[i], 
                " ~ ", result.rejec[j]), collapse = "")
            index = index + 1
        }
    }
    if (all(model.tests == 0) && all(result.accep != returnModel("*", 
        split.string = TRUE)) && tests == TRUE) {
        model.tests = matrix(0, nrow = length(result.accep) + 
            length(result.rejec), ncol = 3)
        dimnames(model.tests) = list(c(1:(length(result.accep) + 
            length(result.rejec))), c("df", "deviance", "p.value"))
        index = 1
        for (i in 1:length(result.accep)) {
            model.tests[index, ] = returnDeviance(model.1 = result.accep[i], 
                model.2 = "*")[c(3, 11, 12)]
            dimnames(model.tests)[[1]][index] = paste(c(result.accep[i], 
                " ~ SATURATED"), collapse = "")
            index = index + 1
        }
        if (length(result.rejec)) 
            for (i in 1:length(result.rejec)) {
                model.tests[index, ] = returnDeviance(model.1 = result.rejec[i], 
                  model.2 = "*")[c(3, 11, 12)]
                dimnames(model.tests)[[1]][index] = paste(c(result.rejec[i], 
                  " ~ SATURATED"), collapse = "")
                index = index + 1
            }
    }
    endCoCo()
    if (tests == TRUE) 
        list(accepted = result.accep, rejected = result.rejec, 
            base = result.base, strategy = strategy, tests = model.tests, 
            `variable names` = var.names)
    else list(accepted = result.accep, rejected = result.rejec, 
        base = result.base, strategy = strategy, `variable names` = var.names)
}
