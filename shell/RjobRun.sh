#!/bin/bash

# Take off the headerline from the line count
N=128 #Number of rows in all the possible combination

for ((i = 1; i<=$N ; i++))
  do
echo $i                          
R -q --save --args $i < "R/getDataFigure3andTable1.R" >> Rjob.out
done 


