# Computing X-ray scatting intensities from MD simulation data

Cameron F. Abrams -- cfa22@drexel.edu

This program uses VMD to compute X-ray scattering intensities from a
trajectory of MD data.  Periodic boundary conditions are assumed.  To
use, follow these steps:

1.  Compile the module in `src/` using

    `make libskdataspace.so`

    (this requires that swig is installed)

2.  Edit `demos/do_calc.tcl` as you like (default values will run the demo)

3.  Run using
  
    `vmd -dispdev text -e do_calc.tcl`

This will create a file "dump.sk" which is S vs k, where S is the
scattering intensity in arb. units and k is the magnitude of the
k-vector.

Some notes: This program assumes the atomic scattering lengths are
proportional to the atomic number.  Because of the periodic
boundaries, it is limited to k-vector magnitudes greater than k_min,
where

$$
k_{\rm min} = 2\pi\sqrt{1/L_x^2 + 1/L_y^2 + 1/L_z^2}
$$

where $L_x$, $L_y$, and $L_z$ are the box dimensions in the x, y, and z
directions, respectively.  E.g., a 40x40x40 angstrom box has a k_min
of about 0.27 angstroms$^{-1}$.
