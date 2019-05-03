#!/bin/bash



ibrun -np 1 root_init_arrays.out > task_1.txt &
ibrun -np 4 root_init_arrays.out > task_4.txt &
ibrun -np 16 root_init_arrays.out > task_16.txt & 

