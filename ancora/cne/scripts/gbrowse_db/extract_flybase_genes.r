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
gff_fn = "/export/data/CNEs/ce4/annotation/c_elegans.WS170.annotations.gff2"
genesID_fn = "/export/data/CNEs/ce4/annotation/geneIDs.WS170"

gff = read_wormbase_gff(gff_fn)

