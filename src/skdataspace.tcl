
load ${CFACV_BASEDIR}/libskdataspace.so libskdataspace

###########################################################
# Primitive data handling procedures
###########################################################

# ListToArray: allocates a new array and assigns it values
# from the Tcl list $l; returns pointer to new array.  List
# elements are treated as double-precision by default.
#
proc ListToArray {l} {
    set length [llength $l]
    set a [new_array $length]
    set i 0
    foreach item $l {
        array_setitem $a $i $item
        incr i 1
    }
    return $a
}

# ListToArray_Data:  Assigns elements of an existing 
# array pointed to by $a from the tcl list $l
#
proc ListToArray_Data { a l } {
    set i 0
    foreach item $l {
	array_setitem $a $i $item
	incr i 1
    }
    return $a
}

# intListToArray: like listToArray except elements of
# $l are treated as integers
#
proc intListToArray {l} {
    set length [llength $l]
    set a [new_arrayint $length]
    set i 0
    foreach item $l {
        arrayint_setitem $a $i $item
        incr i 1
    }
    return $a
}

# intListToArray_Data: list listToArray_Data except elements
# of $l are treated as integers
#
proc intListToArray_Data {a l} {
    set length [llength $l]
    set i 0
    foreach item $l {
        arrayint_setitem $a $i $item
        incr i 1
    }
    return $a
}

# ArrayToList:  Creates a new tcl list of length $n and
# assigns it elements of existing array $a
#
proc ArrayToList {a n} {
    set length $n
    set l {}
    for {set i 0} { $i < $length} {incr i} {
	lappend l [array_getitem $a $i]
    }
    return $l
}

# intArrayToList: like ArrayToList but treats
# elements of $a as integers
#
proc intArrayToList {a n} {
    set length $n
    set l {}
    for {set i 0} { $i < $length} {incr i} {
	lappend l [arrayint_getitem $a $i]
    }
    return $l
}

proc Tcl_SkDataSpace_WriteCoords { sds d clist } {
    set addr [SkDataSpace_getCoordAddr $sds $d]
    ListToArray_Data $addr $clist
    return $addr
}

proc Tcl_SkDataSpace_WriteAndScaleCoords { sds d clist } {
    set addr [Tcl_SkDataSpace_WriteCoords $sds $d $clist]
    SkDataSpace_scaleCoords $sds $d
    return $addr
}

proc Tcl_SkDataSpace_WriteCellSize { sds clist } {
    set addr [SkDataSpace_getCellSizeAddr $sds]
    ListToArray_Data $addr $clist
}

proc Tcl_SkDataSpace_WriteAtomicNumbers { sds clist } {
    set addr [SkDataSpace_getScatLengthAddr $sds]
    ListToArray_Data $addr $clist
}

proc Tcl_doSk { psf dcd nx ny nz dq } {
    mol new $psf
    mol addfile $dcd waitfor all autobonds off

    set base [molinfo top get id]
    set nframes [molinfo $base get numframes]
    puts "INFO/CFASK) $dcd has $nframes frames."
    set all [atomselect $base all]
    set natom [$all num]

    set sds [NewSkDataSpace $natom $nframes $nx $ny $nz]
    Tcl_SkDataSpace_WriteAtomicNumbers $sds [$all get atomicnumber]

    SkDataSpace_setdq $sds $dq

    for { set f 0 } { $f < $nframes } { incr f } {
	puts "INFO/CFASK) processing frame $f..."
	SkDataSpace_setFrame $sds $f
	set a [molinfo top get a frame $f]
	set b [molinfo top get b frame $f]
	set c [molinfo top get c frame $f]
	Tcl_SkDataSpace_WriteCellSize $sds [list $a $b $c]

	$all frame $f

	Tcl_SkDataSpace_WriteAndScaleCoords $sds 0 [$all get x]
	Tcl_SkDataSpace_WriteAndScaleCoords $sds 1 [$all get y]
	Tcl_SkDataSpace_WriteAndScaleCoords $sds 2 [$all get z]

	#SkDataSpace_report $sds

	SkDataSpace_updateSk $sds

    }

    SkDataSpace_outputSk $sds

    mol delete top
    FreeSkDataSpace $sds
}
