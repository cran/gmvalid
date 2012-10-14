gm.delta <-
function (k, l, i, j, pen, A, B) 
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
    sum(A * eval(parse(text = paste("pen[", coord.A, "]", sep = "")))) + 
        sum(B * eval(parse(text = paste("pen[", coord.B, "]", 
            sep = "")))) - pen[k, l]
}
