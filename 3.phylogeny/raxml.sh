#!/bin/bash/

d=$1 # directory to write tree	
f1=$2 # alignment file
cd $d
f2=$(echo $f | awk -F '/' '{print $NF}')
cp $f $f2
echo $f2
/home/bsgroup/soft/standard-RAxML-8.2.3/raxmlHPC-PTHREADS-SSE3 -f a -x 12345 -p 12345 -# 100 -m GTRGAMMA -T 16 -s $f2 -n $f2
rm $f2

