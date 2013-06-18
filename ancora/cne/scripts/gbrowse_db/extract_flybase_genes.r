selfDir = "~/Repositories/genomics"
selfScripts = list.files(path=selfDir, pattern='.*\\.r', full.names=TRUE, recursive=TRUE)
for(rs in selfScripts){message(rs);source(rs)}




## Implementation

# FlyBase
gff_fn = "/export/data/CNEs/dm3/annotation/dmel-all-r5.51.gff"
gff = read_flybase_gff(gff_fn)
gff$transcripts = validate_transcripts(gff$transcripts)
write_gbrowse_gff(gff$genes, gff$transcripts, source="FlyBase", output="flybase.gff")

# WormBase
gff_fn = "/export/data/CNEs/ce10/annotation/c_elegans.WS220.annotations.gff3"
genesID_fn = "/export/data/CNEs/ce10/annotation/geneIDs.WS220"

gff = read_wormbase_gff(gff_fn, genesID_fn)
gff$transcripts = validate_transcripts(gff$transcripts)

write_gbrowse_gff(gff$genes, gff$transcripts, source="WormBase", output="wormbase.gff")

