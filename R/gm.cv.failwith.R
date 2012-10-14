gm.cv.failwith <-
function (default = NULL, f, quiet = FALSE) 
{
    function(...) gm.cv.try_default(f(...), default, quiet = quiet)
}
