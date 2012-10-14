gm.cv.try_default <-
function (expr, default = NA, quiet = FALSE) 
{
    result <- default
    if (quiet) {
        tryCatch(result <- expr, error = function(e) {
        })
    }
    else {
        try(result <- expr)
    }
    result
}
