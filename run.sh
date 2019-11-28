#PBS -N a5_strong 
#PBS -l walltime=01:58:00
#PBS -l nodes=2:ppn=32
#PBS -l mem=2GB
#PBS -S /bin/bash

cd /N/u/liushuj/BigRed2/a5

time aprun -n 64 ./heat 10000 10000 5 5 50.0 8 8 1e-3 
time aprun -n 64 ./heat 10000 10000 5 5 50.0 8 8 1e-3 
time aprun -n 49 ./heat 10000 10000 5 5 50.0 7 7 1e-3 
time aprun -n 49 ./heat 10000 10000 5 5 50.0 7 7 1e-3 
time aprun -n 36 ./heat 10000 10000 5 5 50.0 6 6 1e-3 
time aprun -n 36 ./heat 10000 10000 5 5 50.0 6 6 1e-3 
time aprun -n 25 ./heat 10000 10000 5 5 50.0 5 5 1e-3 
time aprun -n 25 ./heat 10000 10000 5 5 50.0 5 5 1e-3 
time aprun -n 16 ./heat 10000 10000 5 5 50.0 4 4 1e-3 
time aprun -n 16 ./heat 10000 10000 5 5 50.0 4 4 1e-3 
time aprun -n 9 ./heat 10000 10000 5 5 50.0 3 3 1e-3 
time aprun -n 9 ./heat 10000 10000 5 5 50.0 3 3 1e-3 
time aprun -n 4 ./heat 10000 10000 5 5 50.0 2 2 1e-3 
time aprun -n 4 ./heat 10000 10000 5 5 50.0 2 2 1e-3 
time aprun -n 1 ./heat 10000 10000 5 5 50.0 1 1 1e-3 
time aprun -n 1 ./heat 10000 10000 5 5 50.0 1 1 1e-3 
