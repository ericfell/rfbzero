"""
Class to include redox-active crossover mechanisms.
"""


class Crossover:
    """
    Provides crossover mechanism for optional use in ZeroDModel.

    Parameters
    ----------
    membrane_thickness : float
        Membrane thickness (microns).
    permeability_ox : float
        Permeability of oxidized species through membrane (cm^2/s).
    permeability_red : float
        Permeability of reduced species through membrane (cm^2/s).


    Notes
    -----
    Crossover is only available for a symmetric cell, as it is indistinguishable from a first order degradation
    mechanism in a single reservoir of a full cell, when no interactions between crossing posolyte/negolyte species
    are taken into consideration. Current-driven crossover cannot currently be simulated in rfbzero.py.

    """

    def __init__(self, membrane_thickness: float, permeability_ox: float, permeability_red: float) -> None:
        self.membrane_thickness_cm = membrane_thickness / 10000  # convert microns to cm
        self.permeability_ox = permeability_ox
        self.permeability_red = permeability_red

        if self.membrane_thickness_cm <= 0.0:
            raise ValueError("'membrane_thickness' must be > 0.0")

        if self.permeability_ox < 0.0:
            raise ValueError("'permeability_ox' must be >= 0.0")

        if self.permeability_red < 0.0:
            raise ValueError("'permeability_red' must be >= 0.0")

        if self.permeability_ox == 0.0 and self.permeability_red == 0.0:
            raise ValueError("'permeability_ox' and 'permeability_red' cannot both be 0.0")

    def crossover(
            self,
            geometric_area: float,
            c_ox_cls: float,
            c_red_cls: float,
            c_ox_ncls: float,
            c_red_ncls: float,
            volume_cls: float,
            volume_ncls: float,
            time_step: float,
    ) -> tuple[float, float, float, float, float, float]:
        """
        Calculation of crossover species, considering permeabilities of oxidized/reduced species.

        Parameters
        ----------
        geometric_area : float
            Geometric area of cell (cm^2).
        c_ox_cls : float
            CLS concentration of oxidized species (M).
        c_red_cls : float
            CLS concentration of reduced species (M).
        c_ox_ncls : float
            NCLS concentration of oxidized species (M).
        c_red_ncls : float
            NCLS concentration of reduced species (M).
        volume_cls : float
            Volume of CLS reservoir (L).
        volume_ncls : float
            Volume of NCLS reservoir (L).
        time_step : float
            Time interval size (s).

        Returns
        -------
        delta_ox_cls : float
            Change in concentration of oxidized species (M) in the CLS.
        delta_red_cls : float
            Change in concentration of reduced species (M) in the CLS.
        delta_ox_ncls : float
            Change in concentration of oxidized species (M) in the NCLS.
        delta_red_ncls : float
            Change in concentration of reduced species (M) in the NCLS.
        crossed_ox_mols : float
            Oxidized species crossing at given time step (mol).
        crossed_red_mols : float
            Reduced species crossing at given time step (mol).

        """

        # Cell geometric area divided by membrane thickness (cm^2/cm)
        membrane_constant = geometric_area / self.membrane_thickness_cm

        # driving force from concentration differences (M)
        c_ox_difference = c_ox_cls - c_ox_ncls
        c_red_difference = c_red_cls - c_red_ncls

        # amount of species (mols) added/subtracted, divide by 1000 for L to cm^3 conversion
        crossed_ox_mols = time_step * self.permeability_ox * membrane_constant * (c_ox_difference / 1000)
        crossed_red_mols = time_step * self.permeability_red * membrane_constant * (c_red_difference / 1000)

        # changes in concentrations (M)
        delta_ox_cls = -(crossed_ox_mols / volume_cls)
        delta_ox_ncls = crossed_ox_mols / volume_ncls

        delta_red_cls = -(crossed_red_mols / volume_cls)
        delta_red_ncls = crossed_red_mols / volume_ncls

        return delta_ox_cls, delta_red_cls, delta_ox_ncls, delta_red_ncls, crossed_ox_mols, crossed_red_mols
