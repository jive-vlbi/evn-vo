import numpy as np
from numpy import sin, sqrt
from scipy.special import j1
from scipy.constants import pi, c
from scipy.optimize import root_scalar
from astropy import units as u

def primary_beam(theta, D, nu):
    x = pi * (D * nu / c) * sin(theta)
    return (2 * j1(x) / x) ** 2

def synthesized_beam(nu, B):
    return (c / nu) / B

def bandwidth_smearing(theta, nu, delta_nu, theta_HPBW):
    beta = (delta_nu / nu) * theta / theta_HPBW
    return 1 / sqrt(1 + beta ** 2)

def time_smearing(theta, tau, theta_HPBW):
    return 1 - 1.22e-9 * (theta / theta_HPBW) ** 2 * tau ** 2

def total(theta, nu, delta_nu, D, B, tau):
    theta_HPBW = synthesized_beam(nu, B)
    x = primary_beam(theta, D, nu)
    y = bandwidth_smearing(theta, nu, delta_nu, theta_HPBW)
    z = time_smearing(theta, tau, theta_HPBW)
    return x * y * z

def fov(nu, delta_nu, D, B, tau):
    f = lambda theta: total(theta, nu, delta_nu, D, B, tau) - 0.50

    theta_min = 1e-64
    theta_max = 1.22 * (c / nu) / D
    sol = root_scalar(f, bracket=(theta_min, theta_max))
    return sol.root

def resolution(nu, B):
    return synthesized_beam(nu, B)

if __name__ == "__main__":
    tau = 1
    D = 64
    B = 5000000
    nu = 1e9
    delta_nu = 1e6
    s_fov = fov(nu, delta_nu, D, B, tau) * u.rad
    s_resolution = resolution(nu, B) * u.rad
    print(s_fov.to_value(u.arcmin))
    print(s_resolution.to_value(u.arcsec), s_resolution.to_value(u.mas))
    print(total(s_fov.value, nu, delta_nu, D, B, tau))
