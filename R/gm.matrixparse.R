gm.matrixparse <-
function (result) 
{
    m = strsplit(result, "")[[1]]
    elements = unique(m)
    if (any(elements == ",")) 
        elements = elements[-which(elements == ",")]
    elements = sort(elements)
    dep.table = matrix(0, nrow = length(elements), ncol = length(elements))
    dimnames(dep.table) = list(elements, elements)
    for (i in 1:length(m)) {
        if (length(which(elements == m[i])) == 0) {
        }
        else {
            j = 1
            while (i + j <= length(m) && length(which(elements == 
                m[i + j])) > 0) {
                pos1 = which(elements == m[i])
                pos2 = which(elements == m[i + j])
                dep.table[pos1, pos2] = 1
                j = j + 1
            }
        }
    }
    dep.table
}
