import astropy.units as u
import numpy as np

def flux_astro(espl):
    I_0 = 1e-12*u.W/u.m**2
    D = 300 * u.Mpc
    flux = np.power(10, espl/10) * I_0 * D**2
    return flux.to("Lsun cm2 cm-2").value