#!/bin/bash



if [ -d "output" ]; then
	rm -rf "output"
	mkdir "output"
fi

time mpirun -np 60 nrniv -mpi -python run_network.py
