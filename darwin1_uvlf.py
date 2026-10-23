#!/usr/bin/env python3
"""
darwin1_uvlf.py — environment-dependent UV luminosity function from the DARWIN-1 simulation
(Kang et al., MNRAS, in revision).  Minor comment 4: a small, dependency-free script for other users.

    phi = uvlf(M_UV, z, delta, form='schechter')     # cMpc^-3 mag^-1

Parameterisation (Eq. 4 of the paper; valid for 6 <~ z <~ 15 and -0.35 <~ delta <~ 0.35, where delta is
the matter overdensity of a cube of side ~ 8-9 cMpc, i.e. one NIRCam pointing):
    alpha(z, delta)    = a0 + a1 z + a2 delta
    M*(z, delta)       = m0 + m1 z + m2 delta
    log10 phi*(z,delta)= p0 + p1 z + p2 ln(1 + delta)          (natural log: phi* ~ (1+delta)^p2)
    beta(z, delta)     = b0 + b1 z + b2 delta                   (double power law only)
Coefficients: posterior medians of the joint fit with stage-1 uncertainties (COEF_*); the submitted-version
values are kept in COEF_*_V1 for reference.  Intrinsic scatter of the parameterisation about the stacked
simulation LFs: 0.10 (alpha), 0.47 mag (M*), 0.14 dex (log phi*) [Schechter].
Usage from the command line:   python darwin1_uvlf.py 14 0.1        -> prints phi(M) at z=14, delta=0.1
"""
import numpy as np, sys

# Schechter: rows alpha, Mstar, logphi; columns p0, p1 (z), p2 (delta | ln(1+delta))
COEF_SCH = np.array([[-2.2259, 0.0052, -0.1212],
                     [-22.8047, 0.2229, -1.9947],
                     [-3.6634, -0.0895, 1.2642]])
COEF_SCH_ERR = np.array([[0.0205, 0.0020, 0.0266],
                         [0.1002, 0.0099, 0.1299],
                         [0.0354, 0.0035, 0.0453]])
COEF_SCH_V1 = np.array([[-2.0766, -0.0078, 0.1004],
                        [-23.0587, 0.2432, -1.8685],
                        [-3.6724, -0.0849, 1.2385]])
# Double power law: rows alpha, Mstar, logphi, beta
COEF_DPL = np.array([[-2.3435, 0.0328, 0.1744],
                     [-21.7971, 0.2879, -0.6597],
                     [-3.0452, -0.0510, 2.3168],
                     [-3.2143, 0.0034, 0.1282]])
COEF_DPL_ERR = np.array([[0.0364, 0.0037, 0.0493],
                         [0.0603, 0.0060, 0.0831],
                         [0.0448, 0.0044, 0.0573],
                         [0.0742, 0.0074, 0.0982]])
COEF_DPL_V1 = np.array([[-2.1873, 0.0213, 0.2897],
                        [-21.8584, 0.2893, -0.5671],
                        [-3.0317, -0.0465, 2.0342],
                        [-2.9269, -0.0242, 0.4562]])


def params(z, delta, form='schechter', coef=None):
    """Return (alpha, Mstar, logphi[, beta]) at redshift z and overdensity delta."""
    c = coef if coef is not None else (COEF_SCH if form.startswith('s') else COEF_DPL)
    z = np.asarray(z, float); delta = np.asarray(delta, float)
    alpha = c[0, 0] + c[0, 1]*z + c[0, 2]*delta
    mstar = c[1, 0] + c[1, 1]*z + c[1, 2]*delta
    logphi = c[2, 0] + c[2, 1]*z + c[2, 2]*np.log1p(delta)
    if form.startswith('s'):
        return alpha, mstar, logphi
    beta = c[3, 0] + c[3, 1]*z + c[3, 2]*delta
    return alpha, mstar, logphi, beta


def uvlf(M, z, delta=0.0, form='schechter', coef=None):
    """UV luminosity function Phi(M_UV) [cMpc^-3 mag^-1] at redshift z for a region of overdensity delta."""
    M = np.asarray(M, float)
    if form.startswith('s'):
        a, ms, lp = params(z, delta, 'schechter', coef)
        t = 10**(-0.4*(M - ms))
        return 0.4*np.log(10)*10**lp*t**(a + 1)*np.exp(-t)
    a, ms, lp, b = params(z, delta, 'dpl', coef)
    dm = M - ms
    return (np.log(10)/2.5)*10**lp/(10**(0.4*(a + 1)*dm) + 10**(0.4*(b + 1)*dm))


if __name__ == '__main__':
    z = float(sys.argv[1]) if len(sys.argv) > 1 else 10.0
    d = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
    M = np.arange(-23, -15.9, 0.5)
    print(f'# DARWIN-1 UVLF at z = {z}, delta = {d}   (Schechter | DPL), cMpc^-3 mag^-1')
    for m, ps, pd in zip(M, uvlf(M, z, d, 'schechter'), uvlf(M, z, d, 'dpl')):
        print(f'{m:6.1f}  {ps:10.3e}  {pd:10.3e}')
