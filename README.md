# darwin1-uvlf

Environment-dependent UV luminosity function from the DARWIN-1 simulation
(Kang et al., MNRAS, in revision).

`uvlf(M_UV, z, delta, form='schechter' | 'dpl')` returns Φ(M_UV) in cMpc⁻³ mag⁻¹ for
6 ≲ z ≲ 15 and −0.35 ≲ δ ≲ 0.35, where δ is the matter overdensity of a cube the size of one
NIRCam pointing (side ≈ 8–9 cMpc). Outside this range the function is an extrapolation.

```
python darwin1_uvlf.py 14 0.1      # table of Φ(M_UV) at z = 14, δ = 0.1
```

See `example.ipynb` for plots and for propagating the coefficient uncertainties.

Requirements: `numpy` (and `matplotlib` for the notebook).

## Files

| file | content |
|---|---|
| `darwin1_uvlf.py` | the parameterisation (Eq. 4 of the paper) with coefficients, 1σ uncertainties and intrinsic scatter |
| `coefficients.csv` | the same numbers in tabular form, with the submitted-version values for reference |
| `example.ipynb` | usage examples |

## Parameterisation

    alpha(z, δ)       = a0 + a1 z + a2 δ
    M*(z, δ)          = m0 + m1 z + m2 δ
    log10 φ*(z, δ)    = p0 + p1 z + p2 ln(1 + δ)        (φ* ∝ (1+δ)^p2)
    beta(z, δ)        = b0 + b1 z + b2 δ                 (double power law only)

Coefficients are posterior medians of a joint fit to the per-snapshot, per-δ parameters of the
DARWIN-1 sub-volume luminosity functions (50 walkers × 10,000 steps). The parameterisation
reproduces the stacked simulation luminosity functions with an rms of ≈ 0.2 dex; the intrinsic
scatter quoted in the files should be propagated when the function is used quantitatively.

## Citation

If you use this code, please cite the paper and the Zenodo record (DOI to be added at acceptance).
