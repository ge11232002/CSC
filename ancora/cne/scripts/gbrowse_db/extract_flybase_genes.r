selfDir = "~/Repositories/genomics"
selfScripts = list.files(path=selfDir, pattern='.*\\.r', full.names=TRUE, recursive=TRUE)
for(rs in selfScripts){message(rs);source(rs)}

gff_fn = "/export/data/CNEs/dm3/annotation/dmel-all-r5.51.gff"



## Implementation

foo = read_flybase_gff("debug.gff")
foo = read_flybase_gff_whole(gff_fn)

