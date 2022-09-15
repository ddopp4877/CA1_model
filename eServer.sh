#!/bin/bash
echo "starting CA1 run" | sendmail dpd4k4@umsystem.edu

rm "out.txt"

if [ -d "output" ]; then
	rm -rf "output"
	mkdir "output"
fi

time mpirun -np 60 nrniv -mpi -python run_network.py &> out.txt && echo "CA1 done in $(($SECONDS/60)) minutes and $(($SECONDS - $SECONDS/60 * 60)) seconds" | sendmail dpd4k4@umsystem.edu&
