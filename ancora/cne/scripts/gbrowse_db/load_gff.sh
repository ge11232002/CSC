#!/bin/bash

# Prompt for MySQL username and password
echo -n "MySQL username: "
read USER
echo -n "MySQL password: "
stty -echo
read PASS
stty echo
echo ""         # force a carriage return to be output

#hg19
bp_bulk_load_gff.pl -u $USER -p $PASS -c -d gbrowse_gff_hg19 --maxfeature 1000000000 /export/data/CNEs/gff/hg19.gff

#mm10
bp_bulk_load_gff.pl -u $USER -p $PASS -c -d gbrowse_gff_mm10 --maxfeature 1000000000 /export/data/CNEs/gff/mm10.gff



# hg18
#bp_bulk_load_gff.pl -u $USER -p $PASS -c -d gbrowse_gff_hg18 --maxfeature 1000000000 hg18.gff
#
## mm8
#bp_bulk_load_gff.pl -u $USER -p $PASS -c -d gbrowse_gff_mm8 --maxfeature 1000000000 mm8.gff
#
## mm9
#bp_bulk_load_gff.pl -u $USER -p $PASS -c -d gbrowse_gff_mm9 --maxfeature 1000000000 mm9.gff
#
## danRer4
#bp_bulk_load_gff.pl -u $USER -p $PASS -c -d gbrowse_gff_danRer4 --maxfeature 1000000000 danRer4.gff
#
## danRer5
#bp_bulk_load_gff.pl -u $USER -p $PASS -c -d gbrowse_gff_danRer5 --maxfeature 1000000000 danRer5.gff
#
## tetNig1
#bp_bulk_load_gff.pl -u $USER -p $PASS -c -d gbrowse_gff_tetNig1 --maxfeature 1000000000 tetNig1.gff
#
## dm3
#bp_bulk_load_gff.pl -u $USER -p $PASS -c -d gbrowse_gff_dm3 --maxfeature 1000000000 dm3.gff
#
## ce4
#bp_bulk_load_gff.pl -u $USER -p $PASS -c -d gbrowse_gff_ce4 --maxfeature 1000000000 ce4.gff
