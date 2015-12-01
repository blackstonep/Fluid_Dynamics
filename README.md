# Fluid_Dynamics
Computational Physics Fall 2015 Project #2

Code solves 2D steady-state fluid dynamic problem using coupled elliptic partial differential equations with finite differencing.

Data collecting code (collecting_data.cc) works sufficiently well, though takes about 45 minutes to an hour with screen on on a computer using 3.5 GHz processor running OS 10.10.5 for iMacs for each w value. Tried to circumvent this by terminating sweeps if normalized residual value for Psi did not change after 100 sweeps, implying that it converged, which vaguely helps. Unclear how to make this shorter without evaluating less data points (therefore losing information).

Collected data in Google Doc spreadsheet: https://docs.google.com/spreadsheets/d/1uY0ZdyoeFviORP0PHMX97txD-4wUzvq_J03Ajajl34M/edit#gid=0&vpid=A1
