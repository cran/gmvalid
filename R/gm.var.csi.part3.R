gm.var.csi.part3 <-
function (k, l, j, Delta, d, pr) 
{
    ind <- vector()
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
    out <- (((-1)^(1 - ind) * (ind + (1 - ind) * Delta[k, n] - 
        pr[k, l, y] * (sum(pr[m, l, ])/sum(pr[k, l, ])^2) * (1 + 
            Delta[k, n])))/d[k, n] + ((-1)^(1 - ind) * (ind + 
        (1 - ind) * Delta[m, l] - pr[k, l, y] * (sum(pr[k, n, 
        ])/sum(pr[k, l, ])^2) * (1 + Delta[m, l])))/d[m, l])
}
