gm.modelparse <-
function (data) 
{
    v <- dim(data)[1]
    for (i in 1:(v - 1)) {
        for (j in (i + 1):v) {
            data[j, i] <- data[i, j]
        }
    }
    maxvar <- v
    SetSize <- function(X) {
        i <- 0
        while ((X[(i + 1)] > 0) && (i < maxvar)) {
            i <- i + 1
        }
        out <- i
    }
    Neighbours <- function(u, G) {
        j <- 0
        Y <- rep(0, maxvar)
        for (i in 1:maxvar) {
            if (G[i, u] == 1) {
                j <- j + 1
                Y[j] <- i
            }
        }
        out <- Y
    }
    AndSet <- function(X, Y) {
        nx <- SetSize(X)
        ny <- SetSize(Y)
        k <- 0
        Z <- rep(0, maxvar)
        if ((nx > 0) && (ny > 0)) {
            for (i in 1:nx) {
                for (j in 1:ny) {
                  if (X[i] == Y[j]) {
                    k <- k + 1
                    Z[k] <- X[i]
                  }
                }
            }
        }
        out <- Z
    }
    BK <- function(R, P, X, G) {
        result <- NULL
        if ((SetSize(P) == 0) && (SetSize(X) == 0)) {
            result <- rbind(result, R)
        }
        else {
            k <- SetSize(P)
            if (k == 0) {
            }
            else {
                for (i in k:1) {
                  u <- P[i]
                  P <- c(P[-i], 0)
                  nr <- SetSize(R)
                  if (nr == 0) {
                    Rnew <- c(u, rep(0, (maxvar - 1)))
                  }
                  else {
                    Rnew <- c(R[1:nr], u, rep(0, (maxvar - nr - 
                      1)))
                  }
                  N <- Neighbours(u, G)
                  Pnew <- AndSet(P, N)
                  Xnew <- AndSet(X, N)
                  result = rbind(result, BK(Rnew, Pnew, Xnew, 
                    G))
                  nx <- SetSize(X)
                  if (nx == 0) {
                    X <- c(u, rep(0, (maxvar - 1)))
                  }
                  else {
                    X <- c(X[1:nx], u, rep(0, (maxvar - nx - 
                      1)))
                  }
                }
            }
        }
        result
    }
    X <- BK(rep(0, maxvar), 1:maxvar, rep(0, maxvar), data)
    result <- NULL
    for (i in 1:dim(X)[1]) {
        result <- c(result, ",", letters[sort(X[i, ])])
    }
    result <- paste(result[-1], collapse = "")
    result <- paste(sort(strsplit(result, ",")[[1]]), collapse = ",")
    result
}
