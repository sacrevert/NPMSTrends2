#!/bin/bash
filename='bHab.csv'
filelines=`cat $filename`
echo Start
OLDIFS=$IFS
IFS=$'\n'
for line in $filelines ; do
    qsub -v TEXTARG=$line npms2020_test2.job
done
IFS=$OLDIFS