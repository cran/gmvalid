gm.sim.row.i <-
function (base, pr) 
{
    c00 <- runif(1, min = 0.001, max = min(base, pr[1]))
    rest <- base - c00
    pr.i <- length(pr) - 1
    sim <- vector()
    i <- 1
    while (pr.i > 0) {
        if (pr.i > 1) {
            sim[i] <- runif(1, min = 0.001, max = min(rest, pr[i + 
                1]))
        }
        else {
            sim[i] <- rest
        }
        rest <- rest - sim[i]
        pr.i <- pr.i - 1
        i <- i + 1
    }
    row <- c(c00, sim)
    return(row)
}
