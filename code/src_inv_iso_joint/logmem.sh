#!/bin/bash

i=10000000
j=1
# file=top-output_SW.txt
# file=top-output_SW_interval_2_ckbd.txt
# file=top-output_SW_interval_2.txt
# file=top-output_SW_interval_2_voro100.txt
# file=top-output_SW_interval_2_voro200.txt
# file=top-output_SW_interval_2_voro50.txt
# file=top-output_SW_interval_2_voro300.txt
# file=top-output_SW_interval_2_voro150.txt
# file=top-output_SW_interval_2_voro200_iter6.txt
# file=top-output_SW_interval_2_voro100_iter6.txt
file=top-output_SW_interval_2_voro200_iter10.txt
while [ $j -lt $i ]; do
    echo "$(date '+%Y-%m-%d %H:%M:%S')" >> $file
    top -b -n1 | grep "DAzimSurfTomo*" >> $file
    ((j + 1))
    sleep 1
done