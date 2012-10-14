gm.var.csi.part4 <-
function (k, l, j, Delta, d, pr) 
{
    ind <- vector()
    ind.g <- vector()
    if (k == 1) 
        m <- 2
    else if (k == 2) 
        m <- 1
    if (l == 1) 
        n <- 2
    else if (l == 2) 
        n <- 1
    if (j == 1) 
        y <- 2
    else if (j == 2) 
        y <- 1
    ifelse(y == 2, ind <- 1, ind <- 0)
    ifelse(k > l, ind.g <- 1, ind.g <- 0)
    out <- ((ind.g * (pr[m, n, 2] - Delta[m, n] * pr[m, n, 1])/sum(pr[m, 
        n, ]))/d[m, n] + ((1 - ind.g) * (-pr[m, n, 2] + Delta[m, 
        n] * pr[m, n, 1])/sum(pr[m, n, ]))/d[m, n] + ((pr[k, 
        n, 2] - Delta[m, n] * pr[k, n, 1])/sum(pr[k, n, ]) + 
        (pr[m, l, 2] - Delta[m, n] * pr[m, l, 1])/sum(pr[m, l, 
            ]))/d[m, n])
}
