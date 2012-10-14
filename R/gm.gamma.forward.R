gm.gamma.forward <-
function (data, start.model, onestep, headlong, conf.level = 0.95) 
{
    durchgang <- 1
    flag <- 0
    model.matrix = matrix(0, nrow = dim(data)[2], ncol = dim(data)[2])
    while (flag == 0) {
        cat("Durchgang : ", durchgang, "\n")
        vert.in.clique <- vector()
        result <- list()
        pv.mat <- matrix(ncol = 5)
        no.clique <- length(start.model)
        edge.name <- pv <- vert1 <- vert2 <- vector()
        for (i in 1:no.clique) {
            if (length(start.model[[i]]) > 0) {
                result[[i]] <- gm.addEdge(data = data, start.model = start.model, 
                  clique = i, conf.level = conf.level, output = FALSE)
                zeile <- result[[i]]$zeile
                edge.name <- result[[i]]$edge.name
                edge.name <- matrix(edge.name, ncol = 2)
                if (all(!is.na(zeile))) {
                  pv <- result[[i]]$p.value
                  pv.mat <- rbind(pv.mat, cbind(i, zeile, pv, 
                    edge.name))
                }
                else {
                  pv.mat <- rbind(pv.mat, cbind(i, zeile, NA, 
                    NA, NA))
                }
            }
            else {
                pv.mat <- rbind(pv.mat, cbind(i, NA, NA, NA, 
                  NA))
            }
        }
        pv.mat <- pv.mat[order(pv.mat[, 3], na.last = TRUE, decreasing = FALSE), 
            ]
        if (onestep == TRUE) {
            i = 1
            while (!is.na(pv.mat[i, 4])) {
                selected.edge <- c(pv.mat[i, 4], pv.mat[i, 5])
                selected.edge <- sort(selected.edge)
                model.matrix[selected.edge[1], selected.edge[2]] = 1
                i = i + 1
            }
            model.string = gm.modelparse(model.matrix)
            model.string <- strsplit(model.string, ",")[[1]]
            start.model = list()
            for (i in 1:length(model.string)) start.model[[i]] <- match(strsplit(model.string[i], 
                "")[[1]], letters)
            OUT <- gm.activate.stop(start.model, data)
            flag <- 1
        }
        else if (length(pv.mat[, 3][!is.na(pv.mat[, 3])]) > 0) {
            if (headlong == TRUE) {
                selected <- sample(length(pv.mat[, 3][!is.na(pv.mat[, 
                  3])]), 1)
            }
            else selected <- 1
            selected.clique <- pv.mat[selected, 1]
            selected.edge <- c(pv.mat[selected, 4], pv.mat[selected, 
                5])
            selected.edge <- sort(selected.edge)
            kk = 2
            while (model.matrix[selected.edge[1], selected.edge[2]] == 
                1) {
                selected.edge <- c(pv.mat[kk, 4], pv.mat[kk, 
                  5])
                if (any(is.na(selected.edge))) {
                  flag = 1
                  break
                }
                selected.edge <- sort(selected.edge)
                kk = kk + 1
            }
            if (length(start.model) > 1 && flag == 0) {
                model.matrix[selected.edge[1], selected.edge[2]] = 1
                model.string = gm.modelparse(model.matrix)
                model.string <- strsplit(model.string, ",")[[1]]
                start.model = list()
                for (i in 1:length(model.string)) start.model[[i]] = match(strsplit(model.string[i], 
                  "")[[1]], letters)
            }
            else {
                selected.clique <- NA
                selected.edge <- NA
                OUT <- gm.activate.stop(start.model, data = data)
                flag <- 1
            }
        }
        else {
            selected.clique <- NA
            selected.edge <- NA
            OUT <- gm.activate.stop(start.model, data = data)
            flag <- 1
        }
        if (flag == 0 && sapply(start.model, length) == dim(data)[2]) {
            OUT <- gm.activate.stop(start.model, data = data)
            flag <- 1
        }
        else {
            durchgang <- durchgang + 1
        }
    }
    OUT
}
