# darwin1-uvlf
Environment-dependent UV luminosity function from the DARWIN-1 simulation
(Kang et al., MNRAS, in press).

`uvlf(M_UV, z, delta, form='schechter'|'dpl')` returns Φ(M_UV) in cMpc⁻³ mag⁻¹ for
6 ≲ z ≲ 15 and −0.35 ≲ δ ≲ 0.35, where δ is the matter overdensity of a
NIRCam-pointing-sized cube (side ≈ 8–9 cMpc). Requires only numpy.

    python darwin1_uvlf.py 14 0.1      # table of Φ(M) at z = 14, δ = 0.1

Coefficients, 1σ uncertainties and the intrinsic scatter of the parameterisation
are given in `coefficients.csv` and in the docstring. Please cite the paper and
the Zenodo DOI if you use this code.
