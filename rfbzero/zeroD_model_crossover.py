
"""
CROSSOVER FUNCTIONS CALLED BY ZeroDmodel CLASS
"""

# membrane_const is geo Area divide by membrane thickness
# could change that so just thickness is input, geo area already in model


class Crossover:
    """
    Provides crossover mechanism for RFBzero model

    Parameters
    ----------
    membrane_constant : float
        Cell geometric area divided by membrane thickness (cm)
    permeability_ox : float
        Permeability of oxidized species through membrane (cm^2/s)
    permeability_red : float
        Permeability of reduced species through membrane (cm^2/s)
    """

    def __init__(self, membrane_constant: float, permeability_ox: float, permeability_red: float):
        self.membrane_constant = membrane_constant
        self.p_ox = permeability_ox
        self.p_red = permeability_red

    def crossover(self, c_ox_cls: float, c_red_cls: float, c_ox_ncls: float, c_red_ncls: float,
                  timestep: float, vol_cls: float, vol_ncls: float) -> tuple[float, float, float, float, float, float]:
        # skip all this if no crossover is desired
        if self.p_ox == 0.0 and self.p_red == 0.0:
            return c_ox_cls, c_red_cls, c_ox_ncls, c_red_ncls, 0.0, 0.0

        # sets volume term based on side that is source of crossing species
        if c_red_cls < c_red_ncls:
            volume_red = vol_ncls * 1000  # converts liters to cm^3
        else:
            volume_red = vol_cls * 1000  # converts liters to cm^3

        if c_ox_cls < c_ox_ncls:
            volume_ox = vol_ncls * 1000  # converts liters to cm^3
        else:
            volume_ox = vol_cls * 1000  # converts liters to cm^3

        # amount added/subtracted from concentrations, units of deltas are mol/L
        delta_red = (timestep * (self.p_red * self.membrane_constant / volume_red)
                     * (c_red_cls - c_red_ncls))
        delta_ox = (timestep * (self.p_ox * self.membrane_constant / volume_ox)
                    * (c_ox_cls - c_ox_ncls))

        # update concentrations based on permeation
        concentration_red_cls = c_red_cls - delta_red
        concentration_red_ncls = c_red_ncls + delta_red

        concentration_ox_cls = c_ox_cls - delta_ox
        concentration_ox_ncls = c_ox_ncls + delta_ox

        return (concentration_ox_cls, concentration_red_cls,
                concentration_ox_ncls, concentration_red_ncls,
                delta_ox, delta_red)

