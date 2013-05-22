selfDir = "~/Repositories/genomics"
selfScripts = list.files(path=selfDir, pattern='.*\\.r', full.names=TRUE, recursive=TRUE)
for(rs in selfScripts){message(rs);source(rs)}

gff_fn = "/export/data/CNEs/dm3/annotation/dmel-all-r5.51.gff"



## Implementation

# FlyBase
gff = read_flybase_gff(gff_fn)
gff$transcripts = validate_transcripts(gff$transcripts)
write_gbrowse_gff(gff$genes, gff$transcripts, source="FlyBase", output="flybase.gff")



