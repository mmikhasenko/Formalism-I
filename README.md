# Formalism-I
Support material for https://inspirehep.net/record/1642229

The repository contains `jupyter notebooks` with the implementation of the formalism and integration over the final state.
`Julia` kernel is used for the notebooks. It is supposed to be extended with the C++ library to use the formalism in the decay.

The folder `cpp-implementation` contains the c++ version of the formalism and the example of usage.
  - The file `FormalismI.cc` provides the functionality to compute the amplitude.
  - The file `Generator.cc` gives the set of functions to generate four-vectors
  - The file `Example.cc` shows the usage example: the amplitude and density is computed for a single event.

### How to run notebooks
use [JuliaBox](https://juliabox.com), clone this repository and play around.
