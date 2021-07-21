# ff_compute is a tool for the scattering form factors calculation.

This program calculates X-ray and neutron scattering form factors from lipid bilayer simulation trajectories.

The majority of trajectory formats are supported, and the trajectory processing is done in parallel.

# Installation

This code depends on the [Pteros](<https://yesint.github.io/pteros/>) molecular modeling library, so you'll have to install Pteros first. Pteros source code is [available](<https://github.com/yesint/pteros>) on GitHub and a detailed installation guide can be found [here](<https://yesint.github.io/pteros/install.html>). This code doesn't depend on Python bindings available in Pteros, so you can build Pteros with -DWITH_PYTHON="OFF" CMake option to disable Python support.

Once Pteros is installed, you can proceed to install the ```ff_compute``` program. First, create ```/source``` and ```/build``` directories somewhere nice. Next, ```cd``` into the ```/source``` directory and clone the source code of ```ff_compute``` to it with:

```
git clone https://github.com/YevChern/ff_compute . 
```

Once you have the source code of ```ff_compute```, ```cd``` to the ```/build``` directory and compile ```ff_compute``` with:

```
cmake ../source
make -j
```

If you encounter an error that CMake can not find Pteros, check if you sourced ```pterosrc``` file after Pteros installation.

Now you're ready to use ```ff_compute``` for the scattering form factors calculation with lipid bilayers.

# Usage

To run ```ff_compute``` just type:

```
ff_compute -f some_structure.pdb some_trajectory.xtc <other_options>
```

You can see all the available options with ```-help```:

```
ff_compute -help
```

Options:

    -cutoff <nm>
	      Only atoms within this distance form the bilayer COM would be considered.
    -w_dens <1/nm^3>
	      Average number density of water.
    -w_dens_sqr <1/nm^6>
	      Average of water number density squared.
    -water_ind <file>
	      File with water atoms indices (starting with 0).
    -orig_ind <file>
      	File with bilayer atoms indices (starting with 0). Used for the bilayer COM calculation.
    -ions_ind <file>    (optional)
	      File with ions indices (starting with 0).
    -charges <file>     (optional)
      	File with partial charges. Should contain same number of inputs as atoms in the system. If not provided
	      all partial charges are considered to be zero.
    -exch_h <file>      (optional)
      	File with exchangeable hydrogens indices (starting with 0).
    -q_xray_file <file>
    -q_neutron_file <file>
      	File with X-ray and neutron q values.
    -q_start <1/nm>  default: 0.0
    -q_finish <1/nm> default: 1.0
    -q_step <1/nm>   default: 0.01
	      Alternatively, range of q values at which form factors will be calculated can be defined with q_start, 
      	q_finish and q_step. These options are overridden with the values provided in q_xray_file or q_neutron_file files.
    -out_pref <string>
      	Optional prefix for the output file names.


Trajectory processing options:

	General usage:
    -f filename1 filename2 ... <processing options>

Files:
	
    * Exactly one structure file (PDB or GRO)
      If not specified, topology TPR file or TNG trajectory file
      must be given instead.
    * Gromacs topology file (TPR)
      If structure file is also present only topology is read from this file.
      If structure file is not present the coordinates are also read.
    * One or more trajectory files (TRR, XTC, TNG or DCD).
      TNG files also contain the structure, so if no structure file
      is given the structure is read from the first TNG file.

    Files may appear in any order, but trajectory files will be processed
    in the order of their appearance.

Processing options:

    -path <string>
        optional path which will be prepended to all data files, default: empty string
    -b <value[suffix]>
        beginning of processing (starting frame or time), default: 0
    -e <value[suffix]>
        end of processing (end frame or time), default: -1 (up to the end)
    -skip <n>
        Process only each n'th frame, default: -1 (process each frame)
    -t0 <t[suffix]>
        Custom starting time, default: -1 (use value from first frame)
        Useful if trajectory does not contain time stamps
        or if the starting time is incorrect.
        If set and dt is not given sets dt to 1.0!
    -dt <t[suffix]>
        Cutom time step, default: -1 (use value from trajectory)
        Useful if trajectory does not contain time stamps.
        If set and start is not given sets start to 0.0!
    -log <n>
        Prints logging information on each n-th frame, default: -1 (no logging)
    -buffer <n>
        Number of frames, which are kept in memory, default: 10
        Only touch this if individual frames are very large.

Suffixes:

    All parameters marked as <value[suffix]> accept the following optional suffixes:
        (no suffix) - value is in frames
        fr - value is in frames
        t - value is time in picoseconds (value used as is)
        ps - value is time in picoseconds (value used as is)
        ns - value is time in nanoseconds (value multiplied by 10^3)
        us - value is time in microseconds (value multiplied by 10^6)
        ms - value is time in milliseconds (value multiplied by 10^9)
    Parameters marked as <t[suffix]> does not accept fr suffix.
    In this case no suffix means ps.
