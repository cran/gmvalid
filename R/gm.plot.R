gm.plot <- function(model, significant=TRUE, data.analysis=0)
 
#plotting a graphical model

{
                fontsize = 8
                lwd = 10

  if(is.list(data.analysis) && length(model) != length(data.analysis))
    stop("There has to be a data matrix in the list for every model.")

  #load the GRID-library if not already loaded
  require( grid, quietly=TRUE )

  windows()

  row.grid = floor(sqrt(length(model)))
  col.grid = ceiling(sqrt(length(model)))
  
  if(length(model) > row.grid*col.grid)
    row.grid = row.grid + 1
  index = 1

  # Open viewport
  grid.newpage()
  grid.rect(gp=gpar(fill="grey"))
  pushViewport(viewport(layout=grid.layout(row.grid,col.grid)))

  for(rrr in 1:row.grid)
  for(ccc in 1:col.grid)
  {
  m = strsplit(model[index],"")[[1]]
  if(any(m == ":"))
  {
    m = strsplit(model[index],",")[[1]]
    m = strsplit(m,":")
    mm = m[[1]][-(m[[1]] == "")]

    if(length(m) > 1)
      for(i in 2:length(m))
        mm = c(mm,",",m[[i]][-(m[[i]] == "")])

    m = mm
    elements = unique(m)

    if(any(elements == ","))
      elements = elements[-which(elements == ",")]
  }
  else
  {
    elements = c(m[1])
    for(i in 2:length(m))
      if(length(which(elements == m[i])) == 0 && length(which(c(letters,LETTERS) == m[i])) > 0)
        elements = c(elements,m[i])
  }

  elements = sort(elements)

  circle_factor = length(elements)/2/pi
  radius = 0.3 / length(elements)

  dep_table = matrix(0,nrow=length(elements),ncol=length(elements))

  if(is.list(data.analysis) || is.matrix(data.analysis))
  {
    if(is.list(data.analysis))
      data.ana.tmp = data.analysis[[model[index]]]
    else
      data.ana.tmp = data.analysis
  }

  for(i in 1:length(m))
  {
    if(length(which(elements == m[i])) == 0) 
    {
    }
    else
    {
      j=1

      while(i+j <= length(m) && length(which(elements == m[i + j])) > 0)
      {
        pos1 = which(elements == m[i])
        pos2 = which(elements == m[i+j])
        
        if(pos1 > pos2)
        {
          tmp = pos1
          pos1 = pos2
          pos2 = tmp
        }

        dep_table[pos1,pos2] = 1

        j = j + 1
      }
    }
  }

  pushViewport(viewport(layout.pos.row=rrr,layout.pos.col=ccc))
  pushViewport( viewport( name="Model", width=1, height=row.grid/col.grid ) )
  grid.rect(gp=gpar(fill="white"))

  for(i in 1:dim(dep_table)[1])
  {
    x1 = cos(i/circle_factor)/20*dim(dep_table)[1] + 0.5
    y1 = sin(i/circle_factor)/20*dim(dep_table)[1] + 0.5

    grid.circle( x=x1, y=y1, r=radius, gp=gpar( fill="white", col="black" ) )
    .gm.draw.text ( label = elements[i], x = c(x1-0.01,x1+0.01), y = c(y1,y1), xy_null = c( -0.005, -0.015 ),
                color = "black", fontsize = fontsize + 2 )
  }

  for(i in 1:(dim(dep_table)[1]))
    for(j in (1:dim(dep_table)[2])[-i])
    {
        if(dep_table[i,j] == 0 && !significant)
        {
          l.type = 3
        }
        else if(dep_table[i,j] == 0 && significant)
        {
          l.type = 0
        }
        else
        {
          l.type = 1
        }

        result <- .gm.calculate.delta( c( cos(i/circle_factor)/4+0.5, cos(j/circle_factor)/4+0.5 ), c( sin(i/circle_factor)/4+0.5, sin(j/circle_factor)/4+0.5 ), radius )

        x <- c( cos(i/circle_factor)/20*dim(dep_table)[1]+0.5 + result[1], cos(j/circle_factor)/20*dim(dep_table)[1]+0.5 - result[1] )
        y <- c( sin(i/circle_factor)/20*dim(dep_table)[1]+0.5 - result[2], sin(j/circle_factor)/20*dim(dep_table)[1]+0.5 + result[2] )

        if(is.list(data.analysis) || is.matrix(data.analysis))
        {
          grid.lines ( x = x, y = y, arrow = arrow( length = unit( .005, "cm" ),type="open"),
                gp = gpar( fill="white", col="black", lty = l.type, lwd=abs(lwd*data.ana.tmp[i,j]) ) )

          if(dep_table[i,j] != 0  || (!significant && j > i))
          if(x[1] < x[2])
            .gm.draw.text ( label = data.ana.tmp[i,j], x = x, y = y, xy_null = c( result[2]*abs(lwd*data.ana.tmp[i,j])/50, result[1]*abs(lwd*data.ana.tmp[i,j])/50 ),
                color = "black", fontsize = fontsize + 2 )
          else
            .gm.draw.text ( label = data.ana.tmp[i,j], x = x, y = y, xy_null = c( result[1]*abs(lwd*data.ana.tmp[i,j])/50, result[2]*abs(lwd*data.ana.tmp[i,j])/50 ),
                color = "black", alignment=c("right","bottom"), fontsize = fontsize + 2 )
        }
        else
          grid.lines ( x = x, y = y, arrow = arrow( length = unit( .005, "cm" ),type="open"),
                gp = gpar( fill="white", col="black", lty = l.type, lwd=lwd/5 ) )



    }

  popViewport(2)
  index = index + 1
  if(index > length(model))
    break
  }
  popViewport()

  TRUE
}
