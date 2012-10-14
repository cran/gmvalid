gm.generate <-
function (N, p = c(0.5, 0.5, 0.5)) 
{
    if (any(p < 0) || any(p > 1)) 
        stop("P shall be a vector of probabilities!")
    y = matrix(0, nrow = N, ncol = length(p))
    for (i in 1:N) for (j in 1:length(p)) {
        y[i, j] = rbinom(1, 1, p[j]) + 1
    }
    data.frame(y)
}
