# Python Sander Stat Engine (pySS)

Script that retrieves information from AMBER.sander runs, runs user-customized calculations, and prints to a readable format.

Author: ethan.n.ho@gmail.com
Developed on Python 3.6.5, also tested on Python 2.7.10

## Changelog

*6/20/2018*
  * Started writing this script. Its purpose is to do the same thing as out_to_txt.py, but also import .dat files from VMD trajectories to get bond distances. To start, I will develop an object-oriented framework for file recognition and reading. As of now, this script should live in the calculation root directory (same location as batch_run.py). Sucessfully installed VMD as a module from https://github.com/Eigenstate/vmd-python. Docs for VMD module are found at http://www.ks.uiuc.edu/Research/vmd/current/ug/node160.html.
*7/5/2018*
  * Optimized for finding frames of trajectory files with lowest EPtot values. Big handling function convert_lmin_frames will peruse through calculation directories looking for MD trajectory files (NetCDF formatted out of Amber18.pmemd), then looks at the .out file with the same base file name to get EPtot values over the frames of the mdcrd file. Uses VMD to save these coordinates as an rst7 (ASCII) formatted restart file, which can then be converted into an Amber-friendly NetCDF-formatted input file using cpptraj. Updated to v0.1.1.
*7/5/2018*
  * Bug fixes. If .mdcrd file cannot be read as NetCDF, throws Exception that mdcrd file has 0 frames. Struct.loadmol now tries loading using netcdf plugin first, then tried crdbox. Updated to v0.1.2.
*7/12/2018*
  * Big update. In short, VMD API sucks and writes restart files in a non-sensical, non-standardized format. Instead of a .mdcrd > VMD > cpptraj > pmemd type workflow, I'll just remove the VMD portion and use the Amber18.cpptraj Python API (dubbed pytraj) to writes frames from the .mdcrd file. Removed VMD dependency for convert_lmin_frames, which resides in Struct.write_low_EPtot. VMD API dependency is not fully removed, for instance, from the Bond object. Updated to v0.2.0.
*7/18/2018*
  * I am now interested in trajectory analysis, tracking H-bonds between Cys69 and mEndoG's DNA substrate. Since I'm slowly ditching VMD support in favor of pytraj, I will need to more or less overhaul all the functions in this script. Distance tracking does look a lot easier on pytraj than it does on VMD, though.
