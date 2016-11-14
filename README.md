floc (flocking): 3d simulation of interacting self-driven particles

midpoint integration of the equations of motion that are continuous in space and time

start state ordered/disordered, self-propelling to preferred speed, repulsive radial interactions, noise, time lag

compiling: gcc csim3d.c -Wall -lm -O3 -o csim3d

detailed doc: please see inside csim3d.c and sim3d_lib.c

test run: csim3d 123456 200 1 50 10 5 1 0.0 10 0.001 50 0 1.0 110 1 30 0 -1

output line: simulation time, efficiency (order parameter), direction of the vectorial sum of velocity vectors (x, y, z components)
