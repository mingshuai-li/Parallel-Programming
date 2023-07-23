# MPI Heat Equation Solver, Efficient Programming of Multicore Processors and Supercomputers, SS23, TUM

## Group 7: Jingtian Zhao, Mingshuai Li, and Jinwen Pan

## 1. Build
With the default modules already loaded on the login node of SuperMUC-NG, just type `$ make` to compile the program. 

## 2. Run
We have merged the code of the four versions: blocking (default), non-blocking, hybrid blocking, and hybrid non-blocking, into a single version. After compilation, in addition to specifying the relevant experimental parameters in `test.dat`, you can also specify the desired execution mode within the file.

The program is run by executing the command `$ sbatch job.scp`. The necessary script is provided along with the code, and you can modify the corresponding parameters in the script to achieve the desired execution mode. All potential required settings are provided in the script.