

gff_fn = "/export/data/CNEs/dm3/annotation/dmel-all-r5.51.gff"



## Implementation
transcript = list("gene_id"=c(), "type"=c(), "transcript_symbol"=c(), "chr"=c(),
                  "start"=c(), "end"=c(), "strand"=c(), "cds_start"=c(), 
                  "cds_end"=c(), "exons"=c(), NULL)

gene = list("gene_symbol"=c(), "chr"=c(), "start"=c(), "end"=c(), "strand"=c())



read_flybase_gff = function(gff_fn){
  gffFH = file(gff_fn, open="r")
  while(length(oneLine = readLines(gffFH, n=1)) > 0){
    if(grepl("^#", oneLine)){
      next
    }
    oneLineArray = strsplit(oneLine, split="\t")[[1]]
    names(oneLineArray) = c("chr", "source", "type", "start", "end", 
                            "score", "strand", "frame", "extras")
    if(is.na(oneLineArray["source"]) || oneLineArray["source"] != "FlyBase"){
      next
    }
    if(oneLineArray["type"] == "gene"){
      # Handle gene feature
      gene_id = grep("ID={\\w+|-};", oneLineArray["extras"], value=TRUE, perl=TRUE)

  }

}
