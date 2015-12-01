# Fluid_Dynamics
Computational Physics Fall 2015 Project #2

This code is still definitely wrong.

I am not sure what's going on, but the code is having a temper tantrum about resetting values and I can't determine why that is.
I've tried resetting the all the values of xi and psi with each value update so the code starts again (and we're not doing math on matricies that we've already heavily editted with old values). It makes set to reset everything mathematically. My current method appears to be incorrect though because the residuals should not be nan for velocities over two since when we run this code on any individual set of these values with velocity greater than 2, we get a value. 

The alternate file is currently running to collect data, but I don't particularly trust it. I'll try something different in the morning.
