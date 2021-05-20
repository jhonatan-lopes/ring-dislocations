# Half-Space Ring Dislocations

MATLAB influence functions for ring dislocations in a half-space.

Each folder contains the functions required to calculate the influence functions for the stresses due to a ring dislocation in an axisymmetric half-space.

The stresses due to the dislocation can be calculated by:

```matlab
G = BX_KERNELS(rho,zeta,delta,a,mu,kap);
```

where:

* `G` - Stress structure (*e.g. contains G.rr, G.rz and G.zz for axial dislocation*)
* `BX` - Dislocation's identifier
* `rho` - Radial coordinates where stresses are evaluated
* `zeta` - Axial coordinates where stresses are evaluated
* `delta` - Depth at which the dislocation is located (*z-coordinate*)
* `a` - Dislocation's radius
* `mu` - Half-space's modulus of rigidity
* `kap` - Half-space's Kolosov's constant

The following dislocations (and functions) are available:

* `BZ` - Axial
* `BP` - Prismatic glide (*in cylindrical and cartesian*)
* `BQ` - Screw
* `BR` - Radial (*in four path cuts*)
* `BOU` - Boussinesq potential
* `KEL` - Kelvin potential

Each directory also contains a LH_INTEGRALS function to calculate the necessary [Lipschitz-Hankel integrals](https://github.com/jhonatan-lopes/LH_INTEGRALS). While it is redundant to have them all copied in each directory, it also makes it easier to fork only the directory that you need.

Please note that the radial dislocations are [path-cut dependent](https://www.sciencedirect.com/science/article/pii/S0020768308003648) and are offered in four "flavours" to match that. 

The complete elliptic integrals needed for the kernels are calculated via the built-in `ellipke` functions in MATLAB. While this is good for maintanability, these functions do not provide the best performance. They calculate the integrals symbolically and then convert them to numeric, which is quite slow. If performance is needed, another implementation can be used to calculate the integrals by modifying only the lines that contain the `ellipke` functions in the relevant influence function.

