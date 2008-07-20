`gm.gms` <-
function (data, strategy = c("backwards", "forwards", "combined"), 
    model = FALSE, onestep = FALSE, headlong = FALSE, conf.level = 0.95) 
{
    if (dim(data)[2] <= 1) 
        stop("Model selection needs at least two variables")
    if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || 
        conf.level < 0 || conf.level > 1)) 
        stop("'conf.level' must be a single number between 0 and 1")
    require(gtools, quietly = TRUE)
    start.model <- NA
    strategy = match.arg(strategy)
    model.list = list()
    if (model != FALSE) {
        model <- strsplit(model, ",")[[1]]
        start.model = list()
        for (i in 1:length(model)) start.model[[i]] = match(strsplit(model[i], 
            "")[[1]], letters)
    }
    if (strategy == "backwards") {
        model.list <- list(seq(1, dim(data)[2]))
        if (is.na(start.model)) 
            start.model <- model.list
        OUT <- .gm.gamma.backward(data = data, start.model = start.model, 
            onestep = onestep, headlong = headlong, conf.level = conf.level)
    }
    else if (strategy == "forwards") {
        model.list <- eval(parse(text = paste("list(", paste(seq(1, 
            dim(data)[2]), collapse = ",", sep = ""), ")")))
        if (is.na(start.model)) 
            start.model <- model.list
        OUT <- .gm.gamma.forward(data = data, start.model = start.model, 
            onestep = onestep, headlong = headlong, conf.level = conf.level)
    }
    else if (strategy == "combined") {
        model = gm.screening(data, conf.level)
        OUT <- .gm.gamma.backward(data = data, start.model = model, 
            onestep = onestep, headlong = headlong, conf.level = conf.level)
        OUT <- .gm.gamma.forward(data = data, start.model = OUT$accepted, 
            onestep = onestep, headlong = headlong, conf.level = conf.level)
    }
    OUT
}
