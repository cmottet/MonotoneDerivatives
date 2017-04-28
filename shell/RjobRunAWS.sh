#!/bin/bash

# Take off the headerline from the line count
prog_dir="/home/ec2-user/MonotoneDerivatives"
while read line;
  do
  echo $line
  R -q --save --args $line < $prog_dir""/R/runMinMaxSyntheticTransformed_Hill_Horror_Pareto_GLP.R >> Rjob.out
done < $prog_dir""/data/log_Hill_data_rows_to_run.txt