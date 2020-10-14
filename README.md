# VMDTools

A bunch of utility scripts for the VMD visualization software.

**Note**: below is given a summary of each script. However, refer to the
source file for documentation.

## Download

### Revision control

You can use [Git software][1] for revision control to download an
updatable version of this repository. If you don't have it, install it
by following the instructions at this [site][2]. Once you have it
working, go to a terminal and type:

    cd ../path-to-create-a-copy
    git clone https://github.com/franciscoadasme/vmdtools.git

This will generate a `vmdtools` folder at the current location, which is
under revision control. To update your local copy of the repository,
execute

    git pull

in your directory to fetch and merge remote changes.

### Manually

1. Click the **Code** green button in the repository webpage, and then
   **Download ZIP** inside the popup window.
2. Extract the files into a folder.

To update your local copy, you have to download it again and replace the
older files.

## Usage

Most VMD scripts are written in the Tcl language and they usually define
a bunch of functions that should be imported into the VMD session by
invoking the following:

```tcl
source <script>.tcl
```

Some scripts have command-line interface (indicated below), which can be
also run as follows:

```bash
vmd -dispdev none -e <script>.tcl -args ARG1 ARG2 ARG3 ...
```

The documentation of command-line arguments can be read in the header of
the script file.

## Scripts

### hbonds.tcl

Defines utility functions built upon VMD `[measure hbonds ...]` command
to detect h-bonds and water bridges in a structure and trajectory.
H-bonds are detected between two atom selections, which may be mediated
by water molecules. The residues (identified by their `resindex`)
involved in h-bonds are reported instead of individual atoms for better
trackability, as atoms that form h-bonds may change during the
trajectory, especially when waters are involved.

There are three main functions:

- `hbond_search` for finding h-bonds in a single structure
- `hbond_traj` for finding h-bonds in a trajectory
- `hbond_traj_table` for printing out a table of residues involved in
  h-bonds in a trajectory

The following is an example of finding h-bonds in a trajectory and
printing the h-bond table into a comma-separated value (.csv) file:

```tcl
mol new system.pdb
mol addfile traj.dcd
set sel [atomselect top "resname LIG and resid 1"]
set other [atomselect top "protein"]
lassign [hbond_traj $sel $other -waters 2] residues hbonds
set output [open hbonds.csv w]
hbond_traj_table $output $residues $hbonds
close $output
```

Note that `hbond_search` and `hbond_traj` accept optional arguments
(e.g., `-waters`) for controlling how h-bonds are detected. Refer to the
script file for more information.

## Contributors

Be sure to thank these contributors:

- [franciscoadasme](https://github.com/franciscoadasme) Francisco Adasme
  \- creator, maintainer
- [maurobedoya](http://github.com/maurobedoya) Mauricio Bedoya \-
  contributor

## License

Licensed under the MIT license, see the separate LICENSE file.

  [1]: http://en.wikipedia.org/wiki/Git_%28software%29
  [2]: https://help.github.com/articles/set-up-git