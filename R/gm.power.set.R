gm.power.set <-
function (v, first) 
{
    if (length(v) < 2) 
        return(c(paste(c(first, v), collapse = ""), paste(c(first), 
            collapse = "")))
    else {
        return(c(Recall(v[2:length(v)], first), Recall(v[2:length(v)], 
            c(first, v[1]))))
    }
}
