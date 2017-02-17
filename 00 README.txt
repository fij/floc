floc (flocking): 

This code runs the 3d simulation (midpoint integration of coupled differential equations)
defined by a statistical physics based model of interacting self-driven particles. The 
equations of motion (in the model) are continuous in both space and time.

The start state can be both ordered and disordered. The model contains also
self-propelling (adjustment of speed) toward a preferred speed, repulsive radial interactions,
white noise and time lag.

I have done a significant amount of post-processing on the output data, mainly with Perl scripts and Gnuplot.

Compiling: gcc csim3d.c -Wall -lm -O3 -o csim3d

Detailed doc: please see inside csim3d.c and sim3d_lib.c

For the description of the arguments: run 'csim3d' without arguments (or with an incorrect number of arguments)

Test run: csim3d 123456 200 1 50 10 5 1 0.0 10 0.001 50 0 1.0 110 1 30 0 -1

Output line format: simulation time, efficiency (order parameter), direction of the vectorial sum of velocity vectors
