fluids: fluid_dynamics.cc
	-g++ -Wall -o $@ fluid_dynamics.cc

clean:
	rm -f fluids
