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
#bp_bulk_load_gff.pl -u $USER -p $PASS -c -d gbrowse_gff_hg19 --maxfeature 1000000000 /export/data/CNEs/gff/hg19.gff

#mm10
#bp_bulk_load_gff.pl -u $USER -p $PASS -c -d gbrowse_gff_mm10 --maxfeature 1000000000 /export/data/CNEs/gff/mm10.gff

#danRer7
#bp_bulk_load_gff.pl -u $USER -p $PASS -c -d gbrowse_gff_danRer7 --maxfeature 1000000000 /export/data/CNEs/gff/danRer7.gff

#tetNig2
#bp_bulk_load_gff.pl -u $USER -p $PASS -c -d gbrowse_gff_tetNig2 --maxfeature 1000000000 /export/data/CNEs/gff/tetNig2.gff

#equCab2
#bp_bulk_load_gff.pl -u $USER -p $PASS -c -d gbrowse_gff_equCab2 --maxfeature 1000000000 /export/data/CNEs/gff/equCab2.gff

#galGal4
#bp_bulk_load_gff.pl -u $USER -p $PASS -c -d gbrowse_gff_galGal4 --maxfeature 1000000000 /export/data/CNEs/gff/galGal4.gff

#dm3
#bp_bulk_load_gff.pl -u $USER -p $PASS -c -d gbrowse_gff_dm3 --maxfeature 1000000000 /export/data/CNEs/gff/dm3.gff

#canFam3
#bp_bulk_load_gff.pl -u $USER -p $PASS -c -d gbrowse_gff_canFam3 --maxfeature 1000000000 /export/data/CNEs/gff/canFam3.gff

#ce4
#bp_bulk_load_gff.pl -u $USER -p $PASS -c -d gbrowse_gff_ce4 --maxfeature 1000000000 /export/data/CNEs/gff/ce4.gff

#ce10
bp_bulk_load_gff.pl -u $USER -p $PASS -c -d gbrowse_gff_ce10 --maxfeature 1000000000 /export/data/CNEs/gff/ce10.gff








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
