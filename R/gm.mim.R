`gm.mim` <-
function(data,strategy=c("backwards","forwards","eh"),model=FALSE,chain=FALSE,options="")

# analyzes a data table for a graphical model
# parameter: the table or data frame
# returns: the model as a string
# output: the mim tests

{
  require(mimR,quietly=TRUE)

  strategy = match.arg(strategy)
  rejected = NULL
  accepted = NULL
  tests = NULL

  if(is.array(data))
  {
    require(epitools,quietly=TRUE)

    if(!length(dimnames(data)))
      for(i in length(dim(data)):1)
        dimnames(data)[[i]] = as.character(c(1:dim(data)[i]))

    data = data.frame(expand.table(data))
  }

  var.names = cbind(letters[1:dim(data)[2]],dimnames(data)[2][[1]])
  dimnames(data)[2][[1]] = letters[1:dim(data)[2]]

  for(i in 1:dim(data)[2])
    data[[names(data)[i]]] = as.factor(data[[names(data)[i]]])

  toMIM(as.gmData(data))

  if(strategy == "backwards")
  {
    if(chain != FALSE)
    {
      mim.cmd(paste(c("setblocks ",chain),collapse=""))
      mim.cmd("blockmode +")
      if(model != FALSE)
        mim.cmd(paste(c("brmodel ",model),collapse=""))
      else
      {
        mim.cmd("satmod")
        model = chain
      }
    }
    else if(model != FALSE)
      mim.cmd(paste(c("model ",model),collapse=""))
    else
    {
      mim.cmd("satmod")
      model = paste(letters[1:dim(data)[2]],collapse="")
    }
    mim.text = mim.cmd(paste(c("stepwise ",options),collapse=""))

    if(chain != FALSE)
    {
      accepted = mim.cmd("pr")
      accepted = accepted[7:length(accepted)]
      accepted[1:length(accepted)%%2 == 0] = "|"
      accepted = paste(accepted,collapse="")
    }
    else
      accepted = toString(mim.text[length(mim.text)])

    if(accepted != paste(letters[1:dim(data)[2]],collapse=""))
    {
      tests = mim.cmd("test")
      tests = tests[8:length(tests)]
    }
    if(chain != FALSE)
      mim.cmd("blockmode -")
  }
  else if(strategy == "forwards")
  {
    if(chain != FALSE)
    {
      mim.cmd(paste(c("setblocks ",chain),collapse=""))
      mim.cmd("blockmode +")
      if(model != FALSE)
        mim.cmd(paste(c("brmodel ",model),collapse=""))
      else
      {
        mim.cmd("maine")
        model = paste(letters[1:dim(data)[2]],collapse=",")
      }
    }
    else if(model != FALSE)
      mim.cmd(paste(c("Model ",model),collapse=""))
    else
    {
      mim.cmd("maine")
      model = paste(letters[1:dim(data)[2]],collapse=",")
    }
    mim.text = mim.cmd(paste(c("stepwise f",options),collapse=""))

    if(chain != FALSE)
    {
      accepted = mim.cmd("pr")
      accepted = accepted[7:length(accepted)]
      accepted[1:length(accepted)%%2 == 0] = "|"
      accepted = paste(accepted,collapse="")
    }
    else
      accepted = toString(mim.text[length(mim.text)])

    if(accepted != paste(letters[1:dim(data)[2]],collapse=""))
    {
      tests = mim.cmd("test")
      tests = tests[-(1:(grep("LR:",tests)-1))]
    }
    if(chain != FALSE)
      mim.cmd("blockmode -")
  }
  else if(strategy == "eh")
  {
    if(chain != FALSE)
    {
      mim.cmd(paste(c("setblocks ",chain),collapse=""))
      mim.cmd("blockmode +")
    }
    if(model == FALSE)
      model = paste(c(paste(letters[1:dim(data)[2]],collapse=",")," - ",paste(letters[1:dim(data)[2]],collapse="")),collapse="")
    mim.text = mim.cmd(paste(c("initsearch ",model,";startsearch ",options),collapse=""))
    accepted = mim.cmd("ehshow a")
    rejected = mim.cmd("ehshow r")
    tests = c(tests,accepted[5:length(accepted)][1:(length(accepted)-4)%%7 != 1])
    tests = c(tests,rejected[5:length(rejected)][1:(length(rejected)-4)%%7 != 1])
    accepted = accepted[5:length(accepted)][1:(length(accepted)-4)%%7 == 1]
    rejected = rejected[5:length(rejected)][1:(length(rejected)-4)%%7 == 1]

    if(chain != FALSE)
    {
      mim.cmd("blockmode -")
    }
  }
  else
    stop("No valid strategy!")

  if(length(tests) > 0)
  {
    tests = tests[1:length(tests)%%2 == 0]
    tests = as.numeric(tests)
    dim(tests) = c(length(tests)/3,3)
    if(length(rejected) > 0)
      dimnames(tests) = list(c(accepted,rejected),c("LR","df","p.value"))
    else
      dimnames(tests) = list(accepted,c("LR","df","p.value"))
  }

  p.list = list()
  if(chain == FALSE)
  for(ll in 1:length(accepted))
  {
    p.test = mim.cmd(paste(c("model ",accepted[ll],"; step ou"),collapse=""),look.nice=FALSE)
    p.test = p.test[-which(p.test == "+")]
    p.values = as.numeric(matrix(p.test[-(1:grep("P",p.test)[2])],ncol=4,byrow=TRUE)[,4])
    edges = matrix(p.test[-(1:grep("P",p.test)[2])],ncol=4,byrow=TRUE)[,1]
    p.matrix = matrix(NA,nrow=dim(data)[2],ncol=dim(data)[2])
    dimnames(p.matrix) = list(letters[1:dim(data)[2]],letters[1:dim(data)[2]])

    for(i in 1:length(edges))
    {
        e = strsplit(edges[i],"")[[1]]
        p.matrix[which(letters == e[2]),which(letters == e[3])] = p.values[i]
    }

      p.list[[accepted[ll]]] = p.matrix
  }

  list("accepted"=accepted,"rejected"=rejected,"base"=model,"strategy"=strategy,"tests against saturated"=tests,"p values"=p.list,"variable names"=var.names)
}

