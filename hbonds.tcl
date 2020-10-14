# Find all h-bonds in the current frame using geometric criteria.
#
# H-bonds are detected between the atoms in *source* and *target*, which
# can act as both donors and acceptors. The two selections must be from
# the same molecule and should not share atoms. Water bridges connecting
# atoms in *source* and *target* are also supported with up to *-waters*
# water molecules.
#
# Note that this function uses the `[measure hbonds ...]` command from
# VMD internally to locate h-bonds.
#
# @param source atom selection
# @param target another atom selection
# @param -radius select atoms in *target* that are within *radius* (in
#        angstroms) of any atom in *source* for detecting h-bonds
# @param -dist Distance cutoff (in angstrom) between donor and acceptor
#        atoms. Defaults to 3 A
# @param -angle Minimum angle (in degrees) formed by the donor,
#        hydrogen, and acceptor atoms. Defaults to 120 degrees
# @param -waters Number of waters allowed in a water-bridge hydrogen
#        bond between *source* and *target*. Defaults to 2
# @return list of lists of the indices of residues (`resindex`) that
#         form h-bonds. Its length indicates the number of
#         h-bonds detected. The first and last elements of each sublist
#         correspond to the residues in *source* and *target*,
#         respectively, where up to *-waters* intermediate water
#         molecules may be included.
proc hbond_search {source target args} {
    array set opts [list -radius 5.0 -dist 3.0 -angle 120 -waters 2 {*}$args]
    set radius $opts(-radius)
    set max_distance $opts(-dist)
    set min_angle $opts(-angle)
    set max_waters $opts(-waters)

    set hbonds []
    if {[$source num] > 0 && [$target num] > 0} {
        set sel [atomselect top "within $radius of (index [$source get index])"]
        set hbond_table [_hbond_make_table $sel $max_distance $min_angle]
        set target [atomselect top "(index [$sel get index]) and (index [$target get index])"]
        set target [lsort -unique [$target get residue]]
        foreach resindex [lsort -unique [$source get residue]] {
            _hbond_walk $hbond_table $target $resindex {} [expr $max_waters + 2] hbonds
        }
    }
    return $hbonds
}

# Find all residues establishing h-bonds in the trajectory using
# geometric criteria.
#
# Refer to `hbond_search` for additional arguments.
#
# @param source atom selection
# @param target another atom selection
# @return two lists: (i) list of residues (`resindex`) in *target* that
#         establish h-bonds with any residue in *source* along the 
#         trajectory, and (ii) list of list of residues (`resindex`) in
#         *target* that establish h-bonds in the ith frame
proc hbond_traj {source target args} {
    set residues []
    set frame_hbonds []
    for {set f 0} {$f < [molinfo top get numframes]} {incr f} {
        animate goto $f
        set hbonds [hbond_search $source $target {*}$args]
        set hbond_list []
        foreach hbond $hbonds {
            set j [lindex $hbond end]
            lappend residues $j
            lappend hbond_list $j
        }
        lappend frame_hbonds $hbond_list
    }
    return [list $residues $frame_hbonds]
}

# Print out a table of h-bonds found in a trajectory.
#
# Table columns (j) correspond to the residues that formed h-bonds to be
# reported, and table rows (i) are the frames. Each table cell has 1 if
# the jth residue formed a h-bond in the ith frame, else 0. An
# additional row at the end will count the number of frames in which
# each residue formed a h-bond. The table is printed as comma-separated
# values. An example is the following:
#
#     frame,A:ARG14,A:TYR31,A:LEU104
#     0,0,1,1
#     1,1,1,0
#     2,1,1,0
#     3,1,1,1
#     4,0,1,1
#     ...
#     total,102,264,205
#
# @param output output channel (e.g., stdout)
# @param residues list of residues (`resindex`) to use as columns
# @param hbonds list of list of residues that formed h-bonds in the ith
#        frame. Its length indicates the number of frames/rows
proc hbond_traj_table {output residues hbonds} {
    set sorted_residues [lsort -real -unique $residues]
    set rescount [dict create]
    puts -nonewline $output "frame"
    foreach resindex $sorted_residues {
        set sel [atomselect top "residue $resindex"]
        set chain [lindex [$sel get chain] 0]
        set resname [lindex [$sel get resname] 0]
        set resid [lindex [$sel get resid] 0]
        $sel delete
        puts -nonewline $output ",$chain:$resname$resid"
    }
    puts $output ""
    set f 0
    foreach residues $hbonds {
        puts -nonewline $output "$f"
        foreach j $sorted_residues {
            set value [expr [lsearch $residues $j] > -1 ? 1 : 0]
            if {![dict exists $rescount $j]} {
                dict set rescount $j 0
            }
            dict set rescount $j [expr {[dict get $rescount $j] + $value}]
            puts -nonewline $output ",$value"
        }
        puts $output ""
        incr f
    }
    puts -nonewline $output "total"
    foreach resindex $sorted_residues {
        puts -nonewline $output ",[dict get $rescount $resindex]"
    }
    puts $output ""
}

proc _hbond_make_table {sel max_distance min_angle} {
    set hbonds [measure hbonds 3.0 [expr {180-$min_angle}] $sel]
    set hbond_table [dict create]
    for {set k 0} {$k < [llength [lindex $hbonds 0]]} {incr k} {
        set donor [atomselect top "index [lindex $hbonds 0 $k]"]
        set acceptor [atomselect top "index [lindex $hbonds 1 $k]"]
        set ri_donor [$donor get residue]
        set ri_acceptor [$acceptor get residue]

        foreach pair [list [list $ri_donor $ri_acceptor] [list $ri_acceptor $ri_donor]] {
            lassign $pair i j
            if {[dict exists $hbond_table $i]} {
                set hbond_list [dict get $hbond_table $i]
                set hbond_list [linsert $hbond_list end $j]
            } else {
                set hbond_list [list $j]
            }
            dict set hbond_table $i $hbond_list
        }

        $donor delete
        $acceptor delete
    }
    return $hbond_table
}

proc _hbond_walk {hbond_table target resindex history max_steps hbonds} {
    upvar $hbonds output
    lappend history $resindex
    if {[llength $history] > 1 && [lsearch $target $resindex] > -1} {
        lappend output $history
    } elseif {[llength $history] < $max_steps && [dict exists $hbond_table $resindex]} {
        foreach resindex [dict get $hbond_table $resindex] {
            if {[lsearch $history $resindex] == -1} {
                _hbond_walk $hbond_table $target $resindex $history $max_steps output
            }
        }
    }
}
