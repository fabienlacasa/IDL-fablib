# IDL-fablib
A collection of personal IDL routines

Numerical integration
- tsumfast.pro : Trapezoidal summation of the area under a curve. Acceleration of the procedure TSUM from the ASTRON procedure library.
- tsum2D.pro : Trapezoidal summation of the volume under a surface. By making trapezoidal summation over each line then trapezoidal summation of the vector.
- trisum2D.pro : Tetrahedron summation of the volume under a surface.

Spherical harmonic analysis
- ell2nside.pro : Finds the smallest Healpix resolution such that ell<2*Nside
- spectrebin.pro : Computes the binned power spectrum of a map in healpix format
- couplingmatrix_cl : Computes the coupling matrix for power spectrum estimation in partial sky
- cl_masked.pro : Uses the pseudo-spectrum method to measure the power spectrum of a masked map. Calls couplingmatrix_cl.
- poor_inpainter.pro : Inpaints a map iteratively by filling a masked pixel with the average of neighbouring good pixels. Ok for small holes.
- threejzero.pro : Computes the Wigner 3J symbol with azimuthal paramethers m1=m2=m3=0

Miscellaneous
- derivlsq.pro : Derive a function in the presence of noise.
- randgenidl.pro : Generate pseudo random numbers uniformly distributed in [0,1]. Compared to the IDL builtin random number generator, it has a better period but is slower.
- spherbess.pro : Compute the spherical bessel function of order n. Or all orders up to n if keyword /hierarchy is present

