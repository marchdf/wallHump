2D NASA Wall-Mounted Hump Separated Flow
========================================

This presents simulation efforts for the SST and SST-DES
implementations in `Nalu <https://github.com/NaluCFD/Nalu>`_ using
NASA's `2D NASA Wall-Mounted Hump Separated Flow
<https://turbmodels.larc.nasa.gov/nasahump_val.html>`_. The setup,
grids, and NASA code results are all taken from that website.

The initial conditions are chosen to ensure the same flow
conditions. The inflow velocity is set to 34.6 m/s and the Mach number
is 0.1. Assuming atmospheric pressure, the density is therefore 1.185
kg/m^3. The Reynolds number is 936000. The bump chord length, c, is
0.42m. The viscosity is therefore 1.8398e-5 kg/ms. No pressure
gradient is imposed between the inlet and outlet (though the NASA
setup does indicate a small pressure drop across the domain). To
ensure as close a setup as the NASA test cases, no wall function is
used to model the SST wall boundary conditions and the BC for the SST
model are set according to the `NASA specifications
<https://turbmodels.larc.nasa.gov/flatplate_sst.html>`_: k_farfield =
9e-9 a^2_\infty = 9e-9 \frac{34.6^2}{0.1^2} = 0.00108 m^2/s^2,
omega_farfield = 1e-6 \frac{\rho_\infty a^2_\infty}{\mu_\infty} = 1e-6
\frac{1.185}{1.8398e-5} \frac{34.6^2}{0.1^2} = 7710.9 1/s. The
`turbulence intensity
<https://en.wikipedia.org/wiki/Turbulence_kinetic_energy>`_ used by
`NASA for the SST
<https://turbmodels.larc.nasa.gov/nasahump_val_sst.html>`_ is 0.077%,
this yields k_farfield = \frac{3}{2} (IU)^2 =
\frac{3}{2}*(\frac{0.077}{100} 34.6)^2 = 0.00106 m^2/s^2, which is
very similar to the value derived using the other method. Similarly,
according to the SST model, \frac{\mu_t}{\mu} = \frac{\rho a_1
k_farfield}{\mu max(a_1 \omega_farfield, \Omega F_2)} = \frac{\rho
k_farfield}{\mu \omega_farfield} = 0.009, which also matches `the NASA
SST number <https://turbmodels.larc.nasa.gov/nasahump_val_sst.html>`_.

## Using this repository
------------------------

A.  Generating the meshes

#. Get CGNS mesh from `the NASA website <https://turbmodels.larc.nasa.gov/nasahump_grids.html>`_
#. Use Pointwise to label the surfaces and set the BC
#. Use Shreyas' `near-distance-to-wall calculator <https://github.com/NaluCFD/NaluWindUtils>`_ to get the NDTW.

B. Running
.. code-block:: bash

   $ mpiexec -np 1 ./naluX -i wallHump.i

C. Post-processing::
.. code-block:: bash

   $ ./plot_validation.py

## RANS SST 
-----------

### Verification
~~~~~~~~~~~~~~~~

There is good agreement between Nalu's SST implementation
and
[NASA's SST implementation](https://turbmodels.larc.nasa.gov/flatplate_sst.html).

#### Convergence of skin friction coefficient at x = 0.75
<img src="./cf_0.75.png" alt="Cf" width="400">

#### Convergence of skin friction coefficient at x = 0.6321075
<img src="./cf_0.63.png" alt="Cf" width="400">

#### Convergence of skin friction coefficient at x = 0.8678025
<img src="./cf_0.86.png" alt="Cf" width="400">

#### Convergence of drag coefficient
<img src="./cd.png" alt="Cd" width="400">

#### Convergence of pressure drag coefficient
<img src="./cdp.png" alt="Cd" width="400">

#### Convergence of velocity drag coefficient
<img src="./cdv.png" alt="Cd" width="400">

#### Convergence of lift coefficient
<img src="./cl.png" alt="Cd" width="400">

#### Skin friction coefficient along the bump at t = 0.5 (1409x641 mesh)
<img src="./wall_cf.png" alt="wall_Cf" width="400">

#### Pressure coefficient along the bump at t = 0.5 (1409x641 mesh)
<img src="./wall_cp.png" alt="wall_Cf" width="400">

## Thanks
---------
Thanks to Shreyas Ananthan, Ganesh Vijayakumar, and Matt Barone for
their helpful insight and input.
