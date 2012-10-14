gm.var.csi.part1 <-
function (k, l, Delta, d, pr) 
{
    if (k == 1) 
        m <- 2
    else if (k == 2) 
        m <- 1
    if (l == 1) 
        n <- 2
    else if (l == 2) 
        n <- 1
    out <- ((-pr[m, n, 2] + Delta[m, n] * pr[m, n, 1])/sum(pr[m, 
        n, ]) + (pr[k, n, 2] + Delta[k, n] * pr[k, n, 1])/sum(pr[k, 
        n, ]) + (pr[m, l, 2] + Delta[m, l] * pr[m, l, 1])/sum(pr[m, 
        l, ]))/d[m, n]
}
