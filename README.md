# Formalism-I
Support material for https://inspirehep.net/record/1642229

The repository contains `jupyter notebooks` with the implementation of the formalism and integration over the final state.
`Julia` kernel is used for the notebooks. It is supposed to be extended with the C++ library to use the formalism in the decay.

The folder `cpp-implementation` contains the c++ version of the formalism and the example of usage.
  - The file `FormalismI.cc` provides the functionality to compute the amplitude.
  - The file `Generator.cc` gives the set of functions to generate four-vectors
  - The file `Example.cc` shows the usage example: the amplitude and density is computed for a single event.
  - The file `DalitzPlot.cc` generate events according to the phase space and weight them with the density calculated using the formalism. The program outputs `s`,`t`, and the `weight`.
The program takes the number of events as the first argument.

### How to run notebooks
use [JuliaBox](https://juliabox.com), clone this repository and play around.

### How to compile and run the program
The Clebsch-Gordon coefficients are called from `GSL`-library.
Therefore, the program must be linked against `lgsl`.
It is already stored at your machine because it is used by `ROOT`.
  1. First, locate the header and the library for the `gsl` dependence.
```bash
locate gsl_sf_coupling.h
locate libgsl.a
```
You need to specify the path to folder `gsl/` where the headers are stored with the frag `-I`. The path to the library has to be given with the flag `-L`.

  2. Compile and link
```bash
g++ -std=c++11 -L${PATH_TO_GSL_LIBRARY} -I${PATH_TO_GSL_FOLDER} -o Example Example.cc -lgsl
```
Alternatively, export the variables `PATH_TO_GSL_LIBRARY` and `PATH_TO_GSL_FOLDER`
then use the `make`-command

  3. Run
```bash
./Example
```
