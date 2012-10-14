gm.calculate.delta <-
function (x, y, r) 
{
    delta_x <- x[2] - x[1]
    delta_y <- y[2] - y[1]
    x_null <- r/sqrt(delta_x^2 + delta_y^2) * delta_x
    if (y[1] < y[2]) 
        y_null <- -sqrt(abs(r^2 - x_null^2))
    else if (y[1] > y[2]) 
        y_null <- sqrt(abs(r^2 - x_null^2))
    else y_null <- 0
    return(c(x_null, y_null))
}
