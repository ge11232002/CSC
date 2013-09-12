# Bioconductor data package for JASPAR

## The sqlite file
The sqlite file can be generated from the mysql dump. 
The conversion script can be downloaded from [github](https://raw.github.com/gist/1287049/mysql2sqlite.sh).

First you need to load the mysql dump file into mysql.
```sh
mysql -u username -p -e 'create database jaspar'
mysql -u username -p jaspar < jaspar.sql
```
Then do the convertion,
Usage:
```sh
./mysql2sqlite.sh -u root -p jaspar | sqlite3 JASPAR_2010.sqlite
```


## The sites sequences
use DNAStringSetList to store all the sequences. 
There are thousands of sets of sequences there, so the speed might not be too bad to fetch from DNAStringSetList.


