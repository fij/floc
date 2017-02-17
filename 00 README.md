floc (flocking): 3d simulation of interacting self-driven particles

Midpoint integration of the equations of motion that are continuous in both space and time

Start state ordered/disordered, self-propelling to preferred speed, repulsive radial interactions, noise, time lag

Compiling: gcc csim3d.c -Wall -lm -O3 -o csim3d

Detailed doc: please see inside csim3d.c and sim3d_lib.c

For the description of the arguments: run 'csim3d' without arguments (or with an incorrect number of arguments)

Test run: csim3d 123456 200 1 50 10 5 1 0.0 10 0.001 50 0 1.0 110 1 30 0 -1

Output line: simulation time, efficiency (order parameter), direction of the vectorial sum of velocity vectors