sidm-nbody
==========
Monte Carlo N-body Simulation for Self-Interacting Dark Matter


## About this code

Self-Interacting Dark Matter (SIDM) is a hypothetical model for cold
dark matter in the Universe. A strong interaction between dark
matter particles introduce a different physics inside dark-matter
haloes, making the density profile cored, reduce the number of subhaloes,
and trigger gravothermal collapse.

This is an N-body simulation code with Direct Simulation Monte Carlo
scattering for self interaction, and some codes to analyse
gravothermal collapse of isolated haloes (see Koda & Shapiro 2011
[[1]]). The N-body simulation is based on GADGET 1.1 [[2]].

## Compile

See `nbody/Makefile` and run

```bash
$ make
```

in `nbody/`

MPI C compiler and FFTW2 is required; GSL is recomended.

## Run

Give the parameter file as the argument. A sample `parameter.txt` is
in `nbody` directory.

```bash
mpirun -n 2 ./sidm-gadget parameter.txt
```


## References
1. Koda J., Shapiro P. R., 2011, MNRAS, 415, 1125
2. Springel V., Yoshida N., White S. D. M., 2001, New Astronomy, 6, 79

[1]: http://adsabs.harvard.edu/abs/2011MNRAS.415.1125K 
[2]: http://adsabs.harvard.edu/abs/2001NewA....6...79S
