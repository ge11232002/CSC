#!/bin/sh
#
# -- our name ---
#$ -N PWMrand
# -- request /bin/sh --
#$ -S /bin/sh
# -- run in the current working (submission) directory --
#$ -cwd
matlab -nojvm < test.m
