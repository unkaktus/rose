from astropy import units as u, constants as const

class GU():
    def __init__(self, v):
        self.v = v*u.Msun
        self.M = self._to(u.g)
        self.L = self._to(u.cm)
        self.T = self._to(u.s)
        self.Density = self.MLT(1, -3, 0)
        self.Energy = self.MLT(1, 2, -2)
        self.MagneticFluxDensity = self.MLT(1/2, -1/2, -1)
        self.Luminosity = self.Energy/self.T
    def _to(self, unit):
        if unit.is_equivalent(u.Msun):
            return self.v.to(unit)
        if unit.is_equivalent(u.cm):
            return (self.v*const.G/const.c**2).to(unit)
        if unit.is_equivalent(u.s):
            return (self.v*const.G/const.c**3).to(unit)
    def MLT(self, m, l, t):
        return self.M**m * self.L**l * self.T**t