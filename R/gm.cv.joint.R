gm.cv.joint <-
function (pvalue, k, ds, data, conf.level = 0.95, chain = FALSE) 
{
    pvalue.bak <- pvalue
    for (var.i in 1:dim(data)[2]) {
        data[[names(data)[var.i]]] = as.factor(data[[names(data)[var.i]]])
    }
    risk <- list()
    risk.dec <- list()
    stats <- list()
    test.knodes <- list()
    success <- vector()
    gm <- matrix(nrow = k, ncol = 4, dimnames = list(1:k, c("selected model", 
        "number of edges", "success", "block")))
    knodes <- dim(data)[2]
    dimName <- combinations(knodes, 2, letters)
    edgeName <- apply(dimName, 1, paste, collapse = "")
    if (chain == FALSE) {
        influences <- paste(letters[2:knodes], collapse = "")
        setblock <- paste("setblock ", influences, "|a; blockmode +", 
            sep = "")
    }
    else {
        setblock <- paste("setblocks ", chain, "; blockmode +")
    }
    eo <- grep("a", edgeName)
    edgeOut <- edgeName[1:length(eo)]
    for (i in 1:k) {
        df_train <- data[unlist(ds[-i]), ]
        df_test <- data[ds[[i]], ]
        toMIM(as.gmData.data.frame(df_train))
        pvalue[is.na(pvalue)] <- 1
        e <- ifelse(pvalue[i, ] > 1 - conf.level, TRUE, FALSE)
        edgeDEL <- names(pvalue[i, ])[e]
        platzhalter <- paste("del ", paste(edgeDEL, col = "", 
            sep = ""), sep = " ")
        mim.cmd("satmod", look.nice = FALSE)
        if (length(edgeDEL) != 0) {
            mim.cmd(eval(platzhalter), look.nice = FALSE)
        }
        out.model <- mim.cmd("print", look.nice = FALSE)
        out.model <- out.model[-(1:grep("is:", out.model))]
        mim.cmd(eval(setblock), look.nice = FALSE)
        mim.cmd("satmod", look.nice = FALSE)
        if (length(edgeDEL) != 0) {
            mim.cmd(eval(platzhalter), look.nice = FALSE)
        }
        model <- mim.cmd("print", look.nice = FALSE)
        model <- model[-(1:grep("is", model))]
        if (length(model) == 3) {
            model <- model[seq(1, length(model), 2)]
        }
        else if (length(model) > 3) {
            model <- model[seq(2, length(model), 2)]
        }
        model.split <- strsplit(model[length(model)], ",")[[1]]
        model[length(model)] <- model.split[charmatch("a", model.split)]
        model[length(model)] <- strsplit(model[length(model)], 
            ".", fixed = TRUE)[[1]]
        gm[i, 1] <- out.model
        gm[i, 4] <- paste(model, collapse = " | ")
        test.edges <- strsplit(model[length(model)], NULL)[[1]]
        knode <- setdiff(test.edges, c(".", "a"))
        test.knodes[[i]] <- match(knode, letters)
        if (length(test.knodes[[i]]) <= 0) {
            risk[[i]] <- NA
            risk.dec[[i]] <- NA
            gm[i, 3] <- NA
            gm[i, 2] <- length(test.knodes[[i]])
            stats[[i]] <- NA
        }
        else if (length(test.knodes[[i]]) == 1) {
            tab.xy <- table(df_train[, test.knodes[[i]]])
            tab <- table(df_train[, c(test.knodes[[i]], 1)])
            unaff <- tab[, 1]/tab.xy
            aff <- tab[, 2]/tab.xy
            risk[[i]] <- aff/unaff
            risk[[i]] = matrix(risk[[i]], ncol = dim(risk[[i]]))
            for (lauf.dim in length(dim(risk[[i]])):1) {
                dimnames(risk[[i]])[[lauf.dim]] <- 1:dim(risk[[i]])[lauf.dim]
            }
            risk.dec[[i]] <- ifelse(risk[[i]] >= 1, 2, 1)
            for (lauf.dim in length(dim(risk.dec[[i]])):1) {
                dimnames(risk.dec[[i]])[[lauf.dim]] <- 1:dim(risk.dec[[i]])[lauf.dim]
            }
            pred.mat <- data.frame(Var1 = 1:length(risk.dec[[i]]), 
                pred = c(risk.dec[[i]]))
            df_test$Var1 <- df_test[, test.knodes[[i]][1]]
            pred.merge <- merge(df_test, pred.mat, by = c("Var1"))
            diff <- abs(as.numeric(pred.merge[[names(df_test)[1]]]) - 
                pred.merge$pred)
            gm[i, 3] <- 1 - sum(diff)/length(diff)
            gm[i, 2] <- length(test.knodes[[i]])
            TAB.T <- table(pred.merge[[names(df_test)[1]]], pred.merge$pred)
            TN <- TAB.T[1, 1]
            FP <- TAB.T[1, 2]
            FN <- TAB.T[2, 1]
            TP <- TAB.T[2, 2]
            sens <- TP/(TP + FN)
            spec <- TN/(TN + FP)
            prec <- TP/(TP + FP)
            kapp <- 2 * (TP * TN + FP * FN)/((TP + FN) * (FN + 
                TN) + (TP + FP) * (FP + TN))
            stats[[i]] <- matrix(c(sens, spec, prec, kapp), nrow = 4, 
                ncol = 1, dimnames = list(c("Sensitivity", "Specificity", 
                  "Precision", "Cohen's Kappa"), ""))
        }
        else {
            tab.xy <- table(df_train[, test.knodes[[i]]])
            tab <- table(df_train[, c(test.knodes[[i]], 1)])
            choord <- paste(matrix("", nrow = 1, ncol = length(dim(tab))), 
                collapse = ",")
            unaff <- eval(parse(text = paste("tab[", choord, 
                "1]", sep = "")))/tab.xy
            aff <- eval(parse(text = paste("tab[", choord, "2]", 
                sep = "")))/tab.xy
            risk[[i]] <- aff/unaff
            risk[[i]] = array(as.vector(risk[[i]]), dim = dim(risk[[i]]), 
                dimnames = dimnames(risk[[i]]))
            risk.dec[[i]] <- ifelse(risk[[i]] >= 1, 2, 1)
            for (lauf.dim in length(dim(risk.dec[[i]])):1) {
                dimnames(risk.dec[[i]])[[lauf.dim]] <- 1:dim(risk.dec[[i]])[lauf.dim]
            }
            if (length(dim(risk.dec[[i]])) > 1) {
                tmp <- array(1, dim = dim(risk.dec[[i]]))
                dimnames(tmp) <- dimnames(risk.dec[[i]])
                p <- expand.table(tmp)
            }
            pred <- as.vector(aperm(risk.dec[[i]], length(dim(risk.dec[[i]])):1))
            pred.mat <- cbind(p, pred)
            for (lauf.dim in 1:length(test.knodes[[i]])) {
                levels(df_test[, test.knodes[[i]][lauf.dim]]) <- seq(from = 1, 
                  to = length(levels(df_test[, test.knodes[[i]][lauf.dim]])))
            }
            pred.merge <- merge(df_test, pred.mat, by = names(df_test)[test.knodes[[i]]])
            diff <- abs(as.numeric(pred.merge[[names(df_test)[1]]]) - 
                pred.merge$pred)
            gm[i, 3] <- 1 - sum(diff, na.rm = TRUE)/length(diff)
            gm[i, 2] <- length(test.knodes[[i]])
            TAB.T <- table(pred.merge[[names(df_test)[1]]], pred.merge$pred)
            TN <- TAB.T[1, 1]
            FP <- TAB.T[1, 2]
            FN <- TAB.T[2, 1]
            TP <- TAB.T[2, 2]
            sens <- TP/(TP + FN)
            spec <- TN/(TN + FP)
            prec <- TP/(TP + FP)
            kapp <- 2 * (TP * TN + FP * FN)/((TP + FN) * (FN + 
                TN) + (TP + FP) * (FP + TN))
            stats[[i]] <- matrix(c(sens, spec, prec, kapp), nrow = 4, 
                ncol = 1, dimnames = list(c("Sensitivity", "Specificity", 
                  "Precision", "Cohen's Kappa"), ""))
        }
    }
    success.mat <- matrix(nrow = k, ncol = 2, dimnames = list(gm[, 
        1], c("# edges -> a", "success")))
    block.mat <- matrix(nrow = k, ncol = 2, dimnames = list(gm[, 
        4], c("# edges -> a", "success")))
    names(dimnames(success.mat)) <- c("selected graphs", "Success probability")
    success.mat[, 1] <- block.mat[, 2] <- as.numeric(gm[, 2])
    success.mat[, 2] <- block.mat[, 2] <- as.numeric(gm[, 3])
    if (any(!is.na(success.mat[, 2]))) {
        risk = risk[[which.max(success.mat[, 2])]]
        risk.dec = risk.dec[[which.max(success.mat[, 2])]]
        stats <- stats[[which.max(success.mat[, 2])]]
    }
    else {
        risk <- NA
        risk.dec <- NA
        stats <- NA
    }
    success.mat <- success.mat[order(success.mat[, 2], decreasing = TRUE), 
        ]
    block.mat <- block.mat[order(block.mat[, 2], decreasing = TRUE), 
        ]
    OUT <- list(pvalue = pvalue.bak, ratio = risk, risk = risk.dec, 
        success = success.mat, statistics = stats, chain = chain, 
        graph = dimnames(block.mat)[[1]][[1]])
}
