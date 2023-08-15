
"""
CROSSOVER FUNCTIONS CALLED BY ZeroDmodel CLASS
"""


class Crossover:
    """
    Provides crossover mechanism for optional use in ZeroDModel

    Parameters
    ----------
    membrane_constant : float
        Cell geometric area divided by membrane thickness (cm)
    permeability_ox : float
        Permeability of oxidized species through membrane (cm^2/s)
    permeability_red : float
        Permeability of reduced species through membrane (cm^2/s)


    Notes
    -----

    """

    def __init__(self, membrane_constant: float, permeability_ox: float, permeability_red: float):
        """Initialize Crossover"""
        self.membrane_constant = membrane_constant
        self.p_ox = permeability_ox
        self.p_red = permeability_red

        if self.membrane_constant <= 0.0:
            raise ValueError("'membrane_constant' must be a non-zero, positive value")

        if self.p_ox < 0.0 or self.p_red < 0.0:
            raise ValueError("Permeabilities must be positive values")

    def crossover(self, c_ox_cls: float, c_red_cls: float, c_ox_ncls: float, c_red_ncls: float, timestep: float,
                  vol_cls: float, vol_ncls: float) -> tuple[float, float, float, float, float, float]:
        """
        Calculation of crossover species, considering permeabilities of oxidized/reduced species.

        Parameters
        ----------
        c_ox_cls : float
            CLS concentration of oxidized species (M).
        c_red_cls : float
            CLS concentration of reduced species (M).
        c_ox_ncls : float
            NCLS concentration of oxidized species (M).
        c_red_ncls : float
            NCLS concentration of reduced species (M).
        timestep : float
            Simulation time step (s).
        vol_cls : float
            Volume of CLS reservoir (L).
        vol_ncls : float
            Volume of NCLS reservoir (L).

        Returns
        -------
        concentration_ox_cls :
            Updated CLS concentration of oxidized species (M).
        concentration_red_cls :
            Updated CLS concentration of reduced species (M).
        concentration_ox_ncls :
            Updated NCLS concentration of oxidized species (M).
        concentration_red_ncls :
            Updated NCLS concentration of reduced species (M).
        delta_c_ox :
            Concentration difference (CLS-NCLS) of oxidized species (M).
        delta_c_red :
            Concentration difference (CLS-NCLS) of reduced species (M).

        """

        if self.p_ox == 0.0 and self.p_red == 0.0:
            return c_ox_cls, c_red_cls, c_ox_ncls, c_red_ncls, 0.0, 0.0

        # set volume term based on side that is source of crossing species
        # multiply by 1000 for liters to cm^3 conversion
        if c_red_cls < c_red_ncls:
            volume_red = vol_ncls * 1000
        else:
            volume_red = vol_cls * 1000

        if c_ox_cls < c_ox_ncls:
            volume_ox = vol_ncls * 1000
        else:
            volume_ox = vol_cls * 1000

        # driven force via concentration differences (mol/L)
        delta_c_ox = c_ox_cls - c_ox_ncls
        delta_c_red = c_red_cls - c_red_ncls

        # amount added/subtracted from concentrations, units of deltas are mol/L
        delta_red = (timestep * (self.p_red * self.membrane_constant / volume_red) * delta_c_red)
        delta_ox = (timestep * (self.p_ox * self.membrane_constant / volume_ox) * delta_c_ox)

        # update concentrations based on permeation
        concentration_red_cls = c_red_cls - delta_red
        concentration_red_ncls = c_red_ncls + delta_red

        concentration_ox_cls = c_ox_cls - delta_ox
        concentration_ox_ncls = c_ox_ncls + delta_ox

        return (concentration_ox_cls, concentration_red_cls,
                concentration_ox_ncls, concentration_red_ncls,
                delta_c_ox, delta_c_red)


if __name__ == '__main__':
    print('testing')
