# Locate the header and the library for the `gsl` dependence.
# [] locate gsl_sf_coupling.h
# [] locate libgsl.a
#
# For example, I found those at my machine at the root folder
# PATH_TO_LIBRARY="${HOME}/Tools/root-git/localbuild/lib/"
# PATH_TO_GSL_FOLDER="${HOME}/Tools/root-git/localbuild/include/"
#
# export PATH_TO_LIBRARY="${HOME}/Tools/root-git/localbuild/lib/
# export PATH_TO_GSL_FOLDER="${HOME}/Tools/root-git/localbuild/include/"
# then, the `make` command will work
#
# make Example

CXXFLAGS = -std=c++11 -O2
GSL = -L$(PATH_TO_GSL_LIBRARY) -I$(PATH_TO_GSL_FOLDER)

%: %.cc
	g++ $(CXXFLAGS) $(GSL) -o $@ $< -lgsl
