discrete-prolate-spheroidal-sequences
=====================================

Code to generate discrete prolate spheroidal sequences for use in
multitaper spectral estimation methods.

The generated sequences are stored in a HDF5 file.

The code appears to work but uses only a single core on my 16 core
machine. I should look into how to use LAPACK and BLAS properly.