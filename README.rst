Wall-Mounted Hump Separated Flow
================================

This presents simulation efforts for the SST and SST-DES
implementations in `Nalu <https://github.com/NaluCFD/Nalu>`_ using
NASA's `2D NASA Wall-Mounted Hump Separated Flow
<https://turbmodels.larc.nasa.gov/nasahump_val.html>`_. The setup,
grids, and NASA code results are all taken from that website.

The initial conditions are chosen to ensure the same flow
conditions. The inflow velocity is set to :math:`34.6 \unitfrac{m}{s}`
and the Mach number is 0.1. Assuming atmospheric pressure, the density
is therefore :math:`1.185 \unitfrac{kg}{m^3}`. The Reynolds number is 936000. The bump
chord length, :math:`c`, is :math:`0.42\unit{m}`. The viscosity is
therefore :math:`1.8398 \times 10^{-5} \unitfrac{kg}{ms}`. No pressure gradient is
imposed between the inlet and outlet (though the NASA setup does
indicate a small pressure drop across the domain). To ensure as close
a setup as the NASA test cases, no wall function is used to model the
SST wall boundary conditions and the BC for the SST model are set
according to the `NASA specifications
<https://turbmodels.larc.nasa.gov/flatplate_sst.html>`_:

.. math::

   k_farfield = 9 \times 10^{-9} a_{\infty}^2 = 9 \times 10^{-9} \frac{34.6^2}{0.1^2} = 0.00108 \unitfrac{m^2}{s^2}

   \omega_farfield = 10^{-6} \frac{\rho_\infty a_{\infty}^2}{\mu_\infty} = 10^{-6} \frac{1.185}{1.8398 \times 10^{-5}} \frac{34.6^2}{0.1^2} = 7710.9 \unitfrac{1}{s}.


The `turbulence intensity
<https://en.wikipedia.org/wiki/Turbulence_kinetic_energy>`_ used by
`NASA for the SST
<https://turbmodels.larc.nasa.gov/nasahump_val_sst.html>`_ is 0.077%,
this yields

.. math::

   k_farfield = \frac{3}{2} (IU)^2 = \frac{3}{2}\left(\frac{0.077}{100} 34.6\right)^2 = 0.00106 \unitfrac{m^2}{s^2},

which is very similar to the value derived using the other method. Similarly,
according to the SST model,

.. math::

   \frac{\mu_t}{\mu} = \frac{\rho a_1 k_farfield}{\mu \max(a_1 \omega_farfield, \Omega F_2)} = \frac{\rho k_farfield}{\mu \omega_farfield} = 0.009,

which also matches `the NASA SST number <https://turbmodels.larc.nasa.gov/nasahump_val_sst.html>`_.

Using this repository
---------------------

A.  Generating the meshes

    1. Put the CGNS meshes from `the NASA website <https://turbmodels.larc.nasa.gov/nasahump_grids.html>`_ in each run directory
    #. Use Pointwise to execute the `pw_gen_exo.glf` scripts
    #. Run the `add_ndtw.sh` script to get the NDTW (uses Shreyas' `near-distance-to-wall calculator <https://github.com/NaluCFD/NaluWindUtils>`_).

B. Running

   .. code-block:: bash

		   $ mpiexec -np 1 ./naluX -i wallHump.i

C. Post-processing

   .. code-block:: bash

		   $ ./plot_validation.py

RANS SST
--------

Verification
~~~~~~~~~~~~

There is good agreement between Nalu's SST implementation and `NASA's
SST implementation
<https://turbmodels.larc.nasa.gov/flatplate_sst.html>`_.


Thanks
------
Thanks to Shreyas Ananthan, Ganesh Vijayakumar, and Matt Barone for
their helpful insight and input.
