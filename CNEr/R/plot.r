CNEAnnotate = function(CNE, whichAssembly=c(1,2), chr, CNEstart, CNEend, windowSize, min_length){
  # This is the pipeline of doing the density plot
  # The windowSize is in kb.
  windowSize = windowSize * 1000
  CNElength = CNEend - CNEstart + 1
  pixel_width = 800
  if(CNElength <= pixel_width) {
    step_size = 1
  }else{
    step_size = as.integer(CNElength/pixel_width)
    if(step_size > windowSize/10)
      step_size = windowSize/10
    while(windowSize %% step_size){
      step_size = step_size - 1
    }
  }
  # make things easier
  if(windowSize %% 2 == 0)
    windowSize = windowSize - 1
  context_start = max(CNEstart - (windowSize-1)/2, 1)
  context_end = CNEend + (windowSize-1)/2 
  #win_nr_steps = windowSize / step_size
  #context_start = CNEstart - as.integer(((win_nr_steps-1)*step_size)/2+0.5)
  #if(context_start < 1)
  #  context_start = 1
  #context_end = CNEend + as.integer(((win_nr_steps-1)*step_size)/2+step_size+0.5)
  ranges = get_cne_ranges_in_region(CNE, whichAssembly, chr, context_start, context_end, min_length)
  # Implement get_cne_ranges_in_region_partitioned_by_other_chr later!!!
  ranges = reduce(ranges)
  covAll = coverage(ranges, width=context_end)
  runMeanAll = runmean(covAll, k=windowSize, "constant")
  resStart = max(CNEstart, (windowSize-1)/2+1)
  resEnd = min(CNEend, computeEnd-(windowSize-1)/2)
  resCoords = seq(resStart, resEnd, by=step_size)
  runMeanRes = runMeanAll[resCoords]*100
  return(list(resCoords=resCoords, runMeanRes=runMeanRes))
}


#calc_window_scores = function(CNEstart, CNEend, ranges, win_nr_steps, step_size){
#  ## Here the starts and ends are 1-based.
#  CNElength = CNEend - CNEstart + 1
#  win_size = win_nr_steps * step_size
#  offsetBlk = as.integer(((win_nr_steps-1)*step_size)/2+0.5)
#  context_start = CNEstart - offsetBlk
#  if(context_start < 1)
#    context_start = 1
#  context_end = CNEend + offsetBlk
#  context_size = context_end - context_start + 1
#  #nr_blocks = as.integer(context_size/step_size) + ifelse(context_size%%step_size, 1, 0)
#  #blk_scores = numeric(ifelse(nr_blocks>win_nr_steps, nr_blocks, win_nr_steps+1))
#
#  covAll = coverage(ranges, width=context_end)
#   
#  #runMeanAll = runmean(covAll, k=windowSize, "constant")
#  #resStart = max(CNEstart, (windowSize-1)/2+1)
#  #resEnd = min(CNEend, computeEnd-(windowSize-1)/2)
#  #height = runMeanAll[resStart:resEnd]*100
#}


plotCNE = function(listToPlot){
  for(ele in listToPlot){



}




