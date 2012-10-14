gm.delta2 <-
function (k, l, i, j, pen, A, B, AB) 
{
    if (!missing(i)) {
        coord.A <- paste(c(-i, l), collapse = ",")
    }
    else if (missing(i)) {
        coord.A <- paste(c("", l), collapse = ",")
    }
    if (!missing(j)) {
        coord.B <- paste(c(k, -j), collapse = ",")
    }
    else if (missing(i)) {
        coord.B <- paste(c(k, ""), collapse = ",")
    }
    if (!missing(j)) {
        coord.AB <- paste(c(k, -j), collapse = ",")
    }
    else if (missing(i)) {
        coord.AB <- paste(c(k, ""), collapse = ",")
    }
    m <- k * l
    if (k == 2 && l == 1) 
        m <- 3
    mm <- c(4, 3, 2, 1)
    n <- mm[m]
    (AB[m] - AB[n]) * pen[k, l] + B[-k] * pen[-k, l] + A[-l] * 
        pen[k, -l]
}
