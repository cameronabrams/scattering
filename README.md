June 4, 2014

Computing X-ray scatting intensities from MD simulation data:  Version 1

This program uses VMD to compute X-ray scattering intensities from a
trajectory of MD data.  Periodic boundary conditions are assumed.  To
use, must follow these steps:

1.  Compile the module using

    make libskdataspace.so

    (this requires swig)

2.  Edit cfa_sk.tcl as you like

3.  Run, using
  
    vmd -dispdev text -e cfa_sk.tcl

This will create a file "dump.sk" which is S vs k, where S is the
scattering intensity in arb. units and k is the magnitude of the
k-vector.

Some notes: This program assumes the atomic scattering lengths are
proportional to the atomic number.  Because of the periodic
boundaries, it is limited to k-vector magnitudes greater than k_min,
where

k_min = 2*pi*sqrt(1/L_x^2 + 1/L_y^2 + 1/L_z^2)

where L_x, L_y, and L_z are the box dimensions in the x, y, and z
directions, respectively.  E.g., a 40x40x40 angstrom box has a k_min
of about 0.27 \AA^{-1}. To get k_min down to 0.01 \AA^{-1}, we'd need
a box 1000x1000x1000 \AA.

cfa 2014    