#!/bin/bash

for ((INDEX=0; INDEX<20000; INDEX++))
do
   qsub -v JOBID=$INDEX RunAll.sh
   sleep 20
done
