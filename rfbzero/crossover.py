
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
            raise ValueError("'membrane_constant' must be > 0.0")

        if self.p_ox < 0.0 or self.p_red < 0.0:
            raise ValueError("'p_ox' and 'p_red' must be >= 0.0")

    def crossover(self, c_ox_cls: float, c_red_cls: float, c_ox_ncls: float, c_red_ncls: float, vol_cls: float,
                  vol_ncls: float, timestep: float) -> tuple[float, float, float, float, float, float]:
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
        vol_cls : float
            Volume of CLS reservoir (L).
        vol_ncls : float
            Volume of NCLS reservoir (L).
        timestep : float
            Time interval size (s).

        Returns
        -------
        c_ox_cls : float
            Updated CLS concentration of oxidized species (M).
        c_red_cls : float
            Updated CLS concentration of reduced species (M).
        c_ox_ncls : float
            Updated NCLS concentration of oxidized species (M).
        c_red_ncls : float
            Updated NCLS concentration of reduced species (M).
        delta_ox_mols : float
            Oxidized species crossing at given timestep (mol).
        delta_red_mols : float
            Reduced species crossing at given timestep (mol).

        """

        if self.p_ox == 0.0 and self.p_red == 0.0:
            return c_ox_cls, c_red_cls, c_ox_ncls, c_red_ncls, 0.0, 0.0

        """
        # old way
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
        
        # amount added/subtracted from concentrations, units of deltas are mol/L
        delta_red = (timestep * (self.p_red * self.membrane_constant / volume_red) * delta_c_red)
        delta_ox = (timestep * (self.p_ox * self.membrane_constant / volume_ox) * delta_c_ox)

        """

        # driving force from concentration differences (mol/L)
        c_ox_difference = c_ox_cls - c_ox_ncls
        c_red_difference = c_red_cls - c_red_ncls

        # amount added/subtracted from concentrations (mol), divided by 1000 for L to cm^3 conversion
        delta_ox_mols = (timestep * self.p_ox * self.membrane_constant * c_ox_difference) / 1000
        delta_red_mols = (timestep * self.p_red * self.membrane_constant * c_red_difference) / 1000

        # update concentrations (mol/L)
        c_ox_cls -= (delta_ox_mols / vol_cls)
        c_ox_ncls += (delta_ox_mols / vol_ncls)

        c_red_cls -= (delta_red_mols / vol_cls)
        c_red_ncls += (delta_red_mols / vol_ncls)

        return c_ox_cls, c_red_cls, c_ox_ncls, c_red_ncls, delta_ox_mols, delta_red_mols
