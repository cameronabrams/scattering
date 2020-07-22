# example input file to compute S(k) from MD data
# using VMD
#
# (c) 2020 cameron f abrams
# cfa22@drexel.edu
#
if {[info exists env(CFACV_BASEDIR)]} {
    set CFACV_BASEDIR $env(CFACV_BASEDIR)
} else {
    set CFACV_BASEDIR ../src
}

source ${CFACV_BASEDIR}/skdataspace.tcl

# call below does it all.  Arguments are
# PSF file
# DCD file (contains trajectory)
# number of k-vector increments in the x, y, and z directions
# the discretization of S(k) along k (i.e., "delta-k") in \AA^{-1}

Tcl_doSk wat.psf wat_md.dcd 50 50 50 0.1

# the above command computes S(k) from a small MD simulation 
# of a water box with approx 6,000 atoms (40x40x40\AA).

exit

