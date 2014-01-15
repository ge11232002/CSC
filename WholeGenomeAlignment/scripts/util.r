
my.grep = function(patternList, x, combine="or", ...){
  result = logical(length(x))
  if (combine == "or"){
    for(pattern in patternList){
      result = grepl(pattern, x, ...) | result
    }
  }else{
    for(pattern in patternList){
      result = grepl(pattern, x, ...) & result
    }
  }
  return(result)
}

my.system = function(cmd, echo=TRUE, intern=FALSE, ...){
  if (echo){
    message(cmd)
  }
  res = system(cmd, intern=intern, ...)
  if (!intern){
    stopifnot(res == 0)
  }
  return(res)
}

