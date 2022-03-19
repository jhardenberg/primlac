![stand with Ukraine](https://badgen.net/badge/stand%20with/UKRAINE/?color=0057B8&labelColor=FFD700)
# primlac
## Primitive equations solver for shallow lakes

Authors
- Jost von Hardenberg - ISAC-CNR - 3/2009 - 2/2010
- Maurice Rossi       - UPMC-CNRS    

Uses the Wu and Tsanis 1995 parametrization of vertical eddy viscosity

Arakawa C staggered grid                       
Vertical sigma coordinates
Second order interpolation/differencing 
Leapfrog + Adams-Moulton P-C timestepping                             

This code contains some very long lines, 
please use the 132 columns option of your compiler! 

     gfortran -ffixed-line-length-none lac16.F 
     ifort lac16.F -132

When the horizontal resolution NY is set to 1               
the preprocessor automatically removes not relevant lines 

BCs: no-stress+impermeability on the sides,
no-slip at the bottom and wind stress at the top 

