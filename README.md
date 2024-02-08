# LABFM
**The Local Anisotropic Basis Function Method**

**Dr Jack King**

*School of Engineering, University of Manchester*

This toy code was used for the simple convergence studies and stability analysis in
  - King, Lind & Nasar (2020) JCP 415:109549
  - Kind & Lind (2022) JCP 449:110760

Please cite the above if publishing based on this code.

Compiler options:
 - order = 2,3,4,5,6,7,...12 to set the order of LABFM
 - ob = 0,1 to use in-house SVD library or openblas respectively.

How to use:
 - Comment/uncomment routines in labfm.F90 to change type of node generation and type of convergence study, or other things.
 - Modify initial setup at bottom of labfm.F90 to change parameters (e.g. h/s).
 - Change compiler flag in analytic_functions.F90 to specify test functions.
 - Change compiler flag in moments.F90 to change type of basis function.
 - Use scripts in /plots with Octave (or Matlab) to generate plots.
 - If required, use code in vtk_conv to convert uv files (/data_out/uv/) to .vtu files in /paraview_files.



