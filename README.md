# ff_compute is a tool for the scattering form factors calculation.

This program calculates X-ray and neutron scattering form factors from lipid bilayer simulation trajectories.

The majority of trajectory formats are supported, and the trajectory processing is done in parallel.

# Installation

This code depends on the [Pteros](<https://yesint.github.io/pteros/>) molecular modeling library, so you'll have to install Pteros first. Pteros source code is [available](<https://github.com/yesint/pteros>) on GitHub and a detailed installation guide can be found [here](<https://yesint.github.io/pteros/install.html>). This code doesn't depend on Python bindings available in Pteros, so you can build Pteros with -DWITH_PYTHON="OFF" CMake option to disable Python support.
