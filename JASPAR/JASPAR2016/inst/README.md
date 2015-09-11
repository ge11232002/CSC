# Bioconductor data package for JASPAR2016

## The sqlite file
The sqlite file can be generated from the mysql dump.
The conversion script can be downloaded from [github](https://raw.github.com/gist/1287049/mysql2sqlite.sh).

First you need to load the mysql dump file into mysql.
```sh
mysql -u root -p -e 'create database jaspar2016'
mysql -u root -p jaspar2016 < JASPAR_2016_beta.20150910.sql
```
Then do the convertion,
Usage:
```sh
./mysql2sqlite.sh -u root -p jaspar2016 | sqlite3 JASPAR2016.sqlite
```


## The sites sequences
use DNAStringSetList to store all the sequences.
There are thousands of sets of sequences there, so the speed might not be too bad to fetch from DNAStringSetList.

