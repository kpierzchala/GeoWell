import math
import numpy as np
import CoolProp.CoolProp as cp


def c_heat(t, pressure, salt_mass_fraction):
    temperature = t - 273.15
    pbar = pressure / 1e5
    molar_mass_H2O = 18.015
    molar_mass_NaCl = 58.443
    n_NaCl = salt_mass_fraction * 1000 / molar_mass_NaCl
    n_H2O = 1000 / molar_mass_H2O
    xmol = n_NaCl / (n_NaCl + n_H2O)
    xmol1 = 1 - xmol
    xmol12 = xmol1 * xmol1

    q11 = -32.1724 + 0.0621255 * pbar
    q21 = -1.69513 - (4.52781e-4) * pbar - (6.04279e-8) * (pbar**2)
    q22 = 0.0612567 + (1.88082e-5) * pbar

    q1x1 = 47.9048 - (9.36994e-3) * pbar + (6.51059e-6) * (pbar**2)
    q2x1 = 0.241022 + (3.45087e-5) * pbar - (4.28356e-9) * (pbar**2)

    q10 = q1x1
    q20 = 1 - q21 * np.sqrt(q22)
    q12 = -q11 - q10
    q23 = q2x1 - q20 - q21 * np.sqrt(1 + q22)

    q1 = q10 + q11 * xmol1 + q12 * xmol12
    q2 = q20 + (q21 * np.sqrt(xmol + q22)) + q23 * xmol

    tstar_h = q1 + q2 * temperature

    c_heat = q2 * cp.PropsSI(
        "C", "T", tstar_h + 273.15, "P", pressure, "Water"
    )
    return c_heat

def density(t, pressure, salt_mass_fraction):
    T=t-273.15
    P=pressure/1e6
    S=salt_mass_fraction
    Rho_TP=1+0.000001*(-80*T-3.3*T**2+0.00175*T**3+489*P-2*P*T+0.016*T**2*P-0.000013*T**3*P-0.333*P**2-0.002*T*P**2)    
    Rho_TPS=(Rho_TP+S*(0.668+0.44*S+(0.000001*(300*S-2400*P*S+(T*(80+3*T-3300*S-13*P+47*P*S))))))*1000
    return Rho_TPS

def dyn_viscosity(t, p, s):
    rho_H2O = (cp.PropsSI("D", "T", t, "P", p, "Water")) / 1000
    p = p / 1e5
    m = (s * 1000) / 58.44
    d1 = 0.28853170e7
    d2 = -0.11072577e5
    d3 = -0.90834095e1
    d4 = 0.30925651e-1
    d5 = -1 * 0.27407100e-4
    d6 = -0.19283851e7
    d7 = 0.56216046e4
    d8 = 0.13827250e2
    d9 = -0.47609523e-1
    d10 = 0.35545041e-4
    a0 = -0.21319213
    a1 = 0.13651589e-2
    a2 = -0.12191756e-5
    b0 = 0.69161945e-1
    b1 = -0.27292263e-3
    b2 = 0.20852448e-6
    c0 = -0.25988855e-2
    c1 = 0.77989227e-5
    d = np.array([d1, d2, d3, d4, d5, d6, d7, d8, d9, d10])
    a = a0 + a1 * t + a2 * (t**2)
    b = b0 + b1 * t + b2 * (t**2)
    c = c0 + c1 * t
    ln_nr = a * m + b * (m**2) + c * (m**3)
    nr = np.exp(ln_nr)
    sum1 = 0
    for i in range(1, 6):
        term = d[i - 1] * t ** (i - 3)
        sum1 += term

    sum2 = 0
    for i in range(6, 11):
        term = d[i - 1] * rho_H2O * t ** (i - 8)
        sum2 += term
    ln_eta_H2O = sum1 + sum2
    n_H2O = np.exp(ln_eta_H2O)
    n_sol = n_H2O * nr
    return n_sol


def brine_compressibility(t, p, s):
    out = np.abs(
        15.5
        / (0.485 * p / 1.0e5 + 0.5415 * (s * 1000.0) - 966.6 * t - 156.28e3)
    )
    return out


def compressibility(tr, pr, s):
    p = pr / 6894.757
    t = (9.0 / 5.0) * (tr - 273.15) + 32.0
    dVwt = -1.0001e-2 + 1.33391e-4 * t + 5.50654e-7 * t**2.0
    dVwp = (
        -1.95301e-9 * p * t
        - 1.72834e-13 * p**2.0 * t
        - 3.58922e-7 * p
        - 2.25341e-10 * p**2.0
    )
    b = (1 - dVwp) * (1 + dVwt)
    return b


def heat_convection(t, p, s, d, w, lch, t_wall):
    g = 9.81
    Re = (density(t, p, s) * w * d) / dyn_viscosity(t, p, s)
    Pr = c_heat(t, p, s) * dyn_viscosity(t, p, s) / heat_conduction(t, p, s)

    if Re < 2000.0:
        Pr_wsf = (
            c_heat(t_wall, p, s)
            * dyn_viscosity(t_wall, p, s)
            / heat_conduction(t_wall, p, s)
        )
        Gr = (
            g
            * lch**3.0
            * brine_compressibility(t, p, s)
            * (np.abs(t_wall - t))
        ) / (dyn_viscosity(t, p, s) / density(t, p, s)) ** 2.0
        alfa = (
            (0.15 * Re**0.33 * Gr**0.1 * Pr**0.43 * (Pr / Pr_wsf) ** 0.25)
            * heat_conduction(t, p, s)
            / d
        )
    else:
        alfa = (0.0155 * Re**0.83 * Pr**0.5) * (heat_conduction(t, p, s) / d)

    return alfa


def heat_conduction(t, pr, s):
    p = pr / 1e6
    k_H20 = (
        ((7e-9) * (t**3))
        - ((1.5113e-5) * (t**2))
        + ((8.801e-3) * t)
        - (0.8624)
        + ((1.57e-6) * p * t)
    )
    salt_molecular_weight = 58.44
    water_molecular_weight = 18.015
    # Na_molecular_weight = 22.99
    # Cl_molecular_weight = 35.34
    salt_moles = (s * 1000) / salt_molecular_weight
    water_moles = 1000 / water_molecular_weight
    X_Na = salt_moles / (salt_moles + water_moles)
    X_Cl = salt_moles / (salt_moles + water_moles)
    a1Na = 0
    a2Na = 0
    a_Na = a1Na + a2Na * math.exp(-0.023 * (t - 273.15))
    a1Cl = -0.360439
    a2Cl = 0.006076
    a_Cl = a1Cl + a2Cl * math.exp(-0.023 * (t - 273.15))
    k_ions_solvent = a_Na * X_Na + a_Cl * X_Cl
    k = k_H20 + k_ions_solvent
    return k


def dp_flow(t, p, s, fi, w, l, d):

    Re = (density(t, p, s) * w * d) / dyn_viscosity(t, p, s)

    if Re == 0:
        beta = 0
    elif (Re > 0) and (Re < 3.0e3):
        beta = 64.0 / Re
    elif (Re >= 3.0e3) and (Re <= 1.0e4):
        beta = 0.3164 / Re**0.25
    elif Re > 1.0e4:
        beta = (0.221 / Re**0.237) + 0.0032

    dp = (fi * beta * (w**2.0) * density(t, p, s) * l) / (2.0 * d)

    return dp
