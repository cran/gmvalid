gm.gamma.backward <-
function (data, start.model, onestep, headlong, conf.level) 
{
    durchgang <- 1
    flag <- 0
    while (flag == 0) {
        cat("Durchgang : ", durchgang, "\n")
        edge.in.more.than.one.clique <- vector()
        result <- list()
        pv.mat <- matrix(ncol = 5)
        no.clique <- length(start.model)
        edge.out <- edge.name <- pv <- vert1 <- vert2 <- vector()
        for (i in 1:no.clique) {
            if (length(start.model[[i]]) > 1) {
                result[[i]] <- gm.deleteEdge(data = data, model = start.model[[i]], 
                  conf.level = conf.level, output = FALSE)
                edge.out <- result[[i]]$edge.out
                edge.name <- result[[i]]$edge.name
                for (k in 1:length(edge.out)) if (!is.na(edge.out[k])) {
                  pv <- result[[i]]$result[edge.out[k], 5]
                  vert1 <- split(t(edge.name), 1:2)[[1]][k]
                  vert2 <- split(t(edge.name), 1:2)[[2]][k]
                  pv.mat <- rbind(pv.mat, cbind(i, edge.out[k], 
                    pv, vert1, vert2))
                }
                else {
                  pv.mat <- rbind(pv.mat, cbind(i, edge.out[k], 
                    NA, NA, NA))
                }
            }
            else {
                pv.mat <- rbind(pv.mat, cbind(i, NA, NA, NA, 
                  NA))
            }
        }
        pv.mat <- pv.mat[order(pv.mat[, 3], na.last = TRUE, decreasing = TRUE), 
            ]
        if (onestep == TRUE) {
            model.matrix = matrix(1, nrow = dim(data)[2], ncol = dim(data)[2])
            model.matrix = upper.tri(model.matrix)
            i = 1
            while (!is.na(pv.mat[i, 4])) {
                selected.edge <- c(pv.mat[i, 4], pv.mat[i, 5])
                selected.edge <- sort(selected.edge)
                model.matrix[selected.edge[1], selected.edge[2]] = 0
                i = i + 1
            }
            print(model.matrix)
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
            if (length(start.model) > 1 && flag == 0) {
                for (lauf.model in 1:length(start.model)) {
                  edge.in.more.than.one.clique <- c(edge.in.more.than.one.clique, 
                    lauf.model[all(selected.edge %in% start.model[[lauf.model]]) == 
                      TRUE])
                }
            }
            else {
                edge.in.more.than.one.clique <- selected.clique
            }
        }
        else {
            selected.clique <- NA
            selected.edge <- NA
            OUT <- gm.activate.stop(start.model, data = data)
            flag <- 1
        }
        if (flag == 0) {
            for (blur in 1:length(edge.in.more.than.one.clique)) {
                edge.match <- match(selected.edge, start.model[[edge.in.more.than.one.clique[blur]]])
                start.model[[length(start.model) + 1]] <- start.model[[edge.in.more.than.one.clique[blur]]][-edge.match[1]]
                start.model[[length(start.model) + 1]] <- start.model[[edge.in.more.than.one.clique[blur]]][-edge.match[2]]
            }
            start.model <- unique(start.model[-edge.in.more.than.one.clique])
            merke <- vector()
            lauf <- 0
            l <- order(sapply(start.model, length))
            ml <- start.model[l]
            for (i in 1:(length(ml) - 1)) {
                for (j in (i + 1):length(ml)) {
                  lauf <- lauf + 1
                  ifelse(all(ml[[i]] %in% ml[[j]]), merke[lauf] <- i, 
                    merke[lauf] <- NA)
                }
            }
            if (any(!is.na(merke))) {
                merke <- merke[!is.na(merke)]
                start.model <- ml[-merke]
            }
        }
        if (all(sapply(start.model, length) == 1) && flag == 
            0) {
            OUT <- gm.activate.stop(start.model, data = data)
            flag <- 1
        }
        else {
            durchgang <- durchgang + 1
        }
    }
    OUT
}
