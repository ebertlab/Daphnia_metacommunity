#!/bin/bash

for i in {1..546}
do
	while [ `ps | grep julia | wc -l` -gt 10 ]
	do
		echo "Wait for it ... "
		sleep 30s
	done
	echo "Sites $i"
	julia simulated_ss.jl -s $i &
	julia simulated_ss.jl -s $i --withoutComp &
done