gm.modelsim <-
function (N, model, categories = 0) 
{
    require(epitools, quietly = TRUE)
    m = strsplit(model, "")[[1]]
    elements = c(m[1])
    for (i in 2:length(m)) if (length(which(elements == m[i])) == 
        0 && length(which(c(letters, LETTERS) == m[i])) > 0) 
        elements = c(elements, m[i])
    elements = sort(elements)
    if (length(elements) != length(categories) && categories != 
        0) 
        stop("Not enough category definitions!")
    dep.table.1 = NULL
    i = 1
    while (i <= length(m)) {
        if (length(which(c(letters, LETTERS) == m[i])) == 0) {
            i = i + 1
        }
        else {
            j = 0
            clique = NULL
            while (i + j <= length(m) && length(which(c(letters, 
                LETTERS) == m[i + j])) > 0) {
                pos = which(c(LETTERS, 0, letters) == m[i + j])%%27
                clique = c(clique, pos)
                j = j + 1
            }
            i = i + j
            dep.table.1 = c(dep.table.1, 0, clique)
        }
    }
    dep.table.1 = c(dep.table.1, 0)
    zeros = which(dep.table.1 == 0)
    for (i in 1:(length(zeros) - 1)) dep.table.1[(zeros[i] + 
        1):(zeros[i + 1] - 1)] = sort(dep.table.1[(zeros[i] + 
        1):(zeros[i + 1] - 1)])
    dep.table = NULL
    for (i in 1:(length(zeros) - 1)) {
        if (zeros[i + 1] - zeros[i] == 2) 
            dep.table = c(dep.table, 0, dep.table.1[zeros[i] + 
                1])
        else {
            for (j in (zeros[i] + 1):(zeros[i + 1] - 2)) for (k in (j + 
                1):(zeros[i + 1] - 1)) {
                dep.table = c(dep.table, 0, dep.table.1[j], dep.table.1[k])
            }
        }
    }
    dep.table = c(dep.table, 0)
    zeros = which(dep.table == 0)
    if (length(categories[[1]]) < 2) {
        categories = as.list(categories)
        for (i in 1:length(elements)) categories[[i]] = c(0.5, 
            0.5)
    }
    dimension = NULL
    for (i in 1:length(categories)) dimension = c(dimension, 
        length(categories[[i]]))
    y = array(1, dim = dimension)
    coord = NULL
    in.braces <- function(i, flag = 1) {
        result = which(elements == letters[dep.table[zeros[i] + 
            flag]])
        if (length(result) == 0) 
            result = which(elements == LETTERS[dep.table[zeros[i] + 
                flag]])
        result
    }
    for (i in 1:(length(zeros) - 1)) {
        if (zeros[i + 1] - zeros[i] == 2) {
            coord2 = array(c(","), dim = length(elements) - 1)
            coord = NULL
            if (dep.table[zeros[i] + 1] > 1) 
                coord = coord2[1:(in.braces(i) - 1)]
            coord = c(coord, 1)
            if (length(elements) > (in.braces(i))) 
                coord = c(coord, coord2[(in.braces(i)):(length(elements) - 
                  1)])
            for (k in 1:length(categories[[in.braces(i)]])) {
                coord2 = paste(coord, collapse = "")
                eval(parse(text = paste("y[", coord2, "] = categories[[in.braces(i)]][k] * y[", 
                  coord2, "]", sep = "")))
                coord[in.braces(i)] = as.character(k + 1)
            }
        }
        else if (zeros[i + 1] - zeros[i] == 3) {
            table = gm.sim.ixj(N, categories[[in.braces(i)]], 
                categories[[in.braces(i, 2)]])
            coord2 = array(c(","), dim = length(elements) - 1)
            coord = NULL
            if (in.braces(i) > 1) 
                coord = coord2[1:(min(in.braces(i), in.braces(i, 
                  2)) - 1)]
            coord = c(coord, 1)
            coord = c(coord, coord2[min(in.braces(i), in.braces(i, 
                2)):(max(in.braces(i), in.braces(i, 2)) - 1)], 
                1)
            if (length(elements) > max(in.braces(i), in.braces(i, 
                2))) 
                coord = c(coord, coord2[max(in.braces(i), in.braces(i, 
                  2)):(length(elements) - 1)])
            for (k in 1:length(categories[[in.braces(i)]])) {
                coord[min(in.braces(i), in.braces(i, 2))] = as.character(k)
                for (kk in 1:length(categories[[in.braces(i, 
                  2)]])) {
                  coord[max(in.braces(i), in.braces(i, 2)) + 
                    1] = as.character(kk)
                  coord2 = paste(coord, collapse = "")
                  eval(parse(text = paste("y[", coord2, "] = table[k,kk] * y[", 
                    coord2, "]", sep = "")))
                }
            }
        }
        y = y/sum(y)
    }
    y = round(N * y)
    for (i in length(dim(y)):1) dimnames(y)[[i]] = as.character(c(1:dim(y)[i]))
    b = expand.table(y)
    bb = data.frame(b)
    dimnames(bb) = list(c(1:dim(bb)[1]), elements)
    return(bb)
}
