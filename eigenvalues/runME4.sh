#!/bin/bash
export LANG=en_US.UTF-8
for N in 500 1000 1500 3000 5500 7000 10000 
do
  for E in $(seq -3 1 5)
  do
   echo "RUN:E=$E , scan_steps=$N"
   /usr/bin/time -f%E -o time_${N}_${E}.log ./ME4 -c "scan_steps=$N;E=math.exp($E*math.log(10))" 
  done
done 
