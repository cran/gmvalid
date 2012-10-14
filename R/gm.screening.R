gm.screening <-
function (data, conf.level = 0.95) 
{
    if (is.array(data)) {
        require(epitools, quietly = TRUE)
        if (!length(dimnames(data))) 
            for (i in length(dim(data)):1) dimnames(data)[[i]] = as.character(c(1:dim(data)[i]))
        data = data.frame(expand.table(data))
    }
    var.names = cbind(letters[1:dim(data)[2]], dimnames(data)[2][[1]])
    dimnames(data)[2][[1]] = letters[1:dim(data)[2]]
    for (i in 1:dim(data)[2]) data[[names(data)[i]]] = as.factor(data[[names(data)[i]]])
    g = gm.gamma(data = data, type = "m")
    m = matrix(0, nrow = dim(data)[2], ncol = dim(data)[2])
    index = 1
    for (i in 1:(dim(data)[2] - 1)) for (j in (i + 1):dim(data)[2]) {
        p = min(g[index, 5], chisq.test(data[, i], data[, j])$p.value, 
            na.rm = TRUE)
        if (!is.finite(p) || p <= 1 - conf.level) {
            m[i, j] = 1
        }
        else {
            for (k in 1:dim(data)[2]) if (i != k && j != k) {
                p = min(gm.gamma(X = i, Y = j, data = data, conditions = k, 
                  type = "s")[5], gm.chi(table(data[, c(i, j, 
                  k)]))$p.value, na.rm = TRUE)
                if (!is.finite(p) || p <= 1 - conf.level) {
                  m[i, j] = 1
                  break
                }
            }
        }
        index = index + 1
    }
    for (i in 1:(dim(data)[2] - 2)) for (j in (i + 1):(dim(data)[2] - 
        1)) {
        if (m[i, j] == 1) {
            for (k in (j + 1):dim(data)[2]) {
                if (m[i, k] == 1 && m[j, k] == 1) {
                  p1 = min(gm.gamma(X = i, Y = j, data = data, 
                    conditions = k, type = "s")[5], gm.chi(X = i, 
                    Y = j, data = data, Z = k)$p.value, na.rm = TRUE)
                  p2 = min(gm.gamma(X = j, Y = k, data = data, 
                    conditions = i, type = "s")[5], gm.chi(X = j, 
                    Y = k, data = data, Z = i)$p.value, na.rm = TRUE)
                  p3 = min(gm.gamma(X = k, Y = i, data = data, 
                    conditions = j, type = "s")[5], gm.chi(X = k, 
                    Y = i, data = data, Z = j)$p.value, na.rm = TRUE)
                  if (all(is.finite(c(p1, p2, p3)))) 
                    if (p1 > 1 - conf.level) 
                      if (p2 > 1 - conf.level) 
                        if (p3 > 1 - conf.level) {
                          m[i, j] = 0
                          m[j, k] = 0
                          m[i, k] = 0
                        }
                        else if (p1 > p2) 
                          m[i, j] = 0
                        else m[j, k] = 0
                      else if (p3 > 1 - conf.level) 
                        if (p1 > p3) 
                          m[i, j] = 0
                        else m[i, k] = 0
                      else m[i, j] = 0
                    else if (p2 > 1 - conf.level) 
                      if (p3 > 1 - conf.level) 
                        if (p2 > p3) 
                          m[j, k] = 0
                        else m[i, k] = 0
                      else m[j, k] = 0
                    else if (p3 > 1 - conf.level) 
                      m[i, k] = 0
                }
            }
        }
    }
    list(mat = m, model = gm.modelparse(m))
}
