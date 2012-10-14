gm.draw.text <-
function (label, x, y, xy_null = c(0, 0), color, alignment = c("left", 
    "bottom"), fontsize) 
{
    x_label <- x[1] + xy_null[1]
    y_label <- y[1] + xy_null[2]
    x_delta <- x[2] - x[1]
    y_delta <- y[2] - y[1]
    angle = atan(y_delta/x_delta) * (180/pi)
    if (angle == -90) 
        angle = 90
    if (is.numeric(label)) 
        label <- as.numeric(as.integer(label * 10000)/10000)
    pushViewport(viewport(x = x_label, y = y_label, width = 0, 
        height = , angle = angle, name = "vp1", just = alignment))
    grid.text(label, x = 0, y = unit(0.75, "mm"), just = alignment, 
        gp = gpar(fontsize = fontsize, col = color))
    popViewport()
}
