#!/usr/bin/bash

n=3
liste=("tota" "toa" "ta")
for ((i=0 ; $i - $n ; i++))
do
	echo "./lala" ${liste[i]}
done
