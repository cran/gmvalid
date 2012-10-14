gm.make.unique.cliques <-
function (start.model) 
{
    letter.func <- function(x) letters[x]
    letter.list <- lapply(start.model, letter.func)
    letter.list <- lapply(letter.list, paste, collapse = "")
    letter.list <- unique(letter.list)
    model <- paste(matrix(sapply(letter.list, paste, collapse = ""), 
        nrow = 1, ncol = length(letter.list)), collapse = ",")
}
