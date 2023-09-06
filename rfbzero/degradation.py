
"""
DEGRADATION FUNCTIONS CALLED BY ZeroDmodel CLASS
"""

from abc import ABC, abstractmethod


class DegradationMechanism(ABC):
    """
    Base class to be overridden by specific degradation mechanisms.
    """
    @abstractmethod
    def degrade(self, c_ox: float, c_red: float, timestep: float) -> tuple[float, float]:
        """Applies (possible) degradation mechanisms to oxidized/reduced species at each timestep"""
        raise NotImplementedError


class ChemicalDegradation(DegradationMechanism):
    """
    Provides an N-th order chemical degradation mechanism.

    Parameters
    ----------
    rate_order : int
        Rate order for chemical degradation reaction.
    rate_constant : float
        N-th order rate constant of chemical oxidation (unit is rate order-dependent).
    species : str
        Species ('red' or 'ox') undergoing chemical degradation.

    """

    def __init__(self, rate_order: int, rate_constant: float, species: str = 'red'):
        self.rate_order = rate_order
        self.rate_constant = rate_constant
        self.species = species

        if not isinstance(self.rate_order, int) or self.rate_order < 0:
            raise ValueError("'rate_order' must be an integer")

        if self.rate_constant <= 0.0:
            raise ValueError("'rate_constant' must be > 0.0")

        if self.species not in ['red', 'ox']:
            raise ValueError(f"Invalid value: {self.species}. 'species' options: 'red', 'ox' ")

    def degrade(self, c_ox: float, c_red: float, timestep: float) -> tuple[float, float]:
        """
        N-th order chemical degradation of
        N*[redox-active species] --> [redox-inactive species]

        Parameters
        ----------
        c_ox : float
            Concentration of oxidized species (M).
        c_red : float
            Concentration of reduced species (M).
        timestep : float
            Time interval size (s).

        Returns
        -------
        c_ox : float
            (Possibly) updated concentration of oxidized species (M).
        c_red : float
            (Possibly) updated concentration of reduced species (M).

        """

        if self.species == 'red':
            c_red -= (timestep * self.rate_constant * (c_red**self.rate_order))
        else:
            c_ox -= (timestep * self.rate_constant * (c_ox**self.rate_order))

        return c_ox, c_red


class AutoOxidation(DegradationMechanism):
    """
    Provides a 1st order auto-oxidation mechanism, (red --> ox) with no loss of active material. Can be thought of as a
    chemical oxidation of the redox-active, balanced by the reduction of some species not of interest to the model i.e.,
    water reduction (HER). Typically, this could occur in a low-potential negolyte and be considered a self-discharge.

    Parameters
    ----------
    rate_constant : float
        First order rate of auto-oxidation (1/s).

    """
    def __init__(self, rate_constant: float):
        self.rate_constant = rate_constant
        if self.rate_constant <= 0.0:
            raise ValueError("'rate_constant' must be > 0.0")

    def degrade(self, c_ox: float, c_red: float, timestep: float) -> tuple[float, float]:
        """
        Assumes first order process: red --> ox

        Parameters
        ----------
        c_ox : float
            Concentration of oxidized species (M).
        c_red : float
            Concentration of reduced species (M).
        timestep : float
            Time interval size (s).

        Returns
        -------
        c_ox : float
            Updated concentration of oxidized species (M).
        c_red : float
            Updated concentration of reduced species (M).

        """

        delta_concentration = timestep * self.rate_constant * c_red

        c_ox += delta_concentration
        c_red -= delta_concentration
        return c_ox, c_red


class AutoReduction(DegradationMechanism):
    """
    Provides a 1st order auto-reduction mechanism, (ox --> red) with no loss of active material. Can be thought of as a
    chemical reduction of the redox-active, balanced by the oxidation of some species not of interest to the model i.e.,
    water oxidation (OER). Typically, this could occur in a high-potential posolyte and be considered a self-discharge.

    Parameters
    ----------
    rate_constant : float
        First order rate of auto-oxidation (1/s).

    """
    def __init__(self, rate_constant: float):
        self.rate_constant = rate_constant
        if self.rate_constant <= 0.0:
            raise ValueError("'rate_constant' must be > 0.0")

    def degrade(self, c_ox: float, c_red: float, timestep: float) -> tuple[float, float]:
        """
        Assumes first order process: ox --> red

        Parameters
        ----------
        c_ox : float
            Concentration of oxidized species (M).
        c_red : float
            Concentration of reduced species (M).
        timestep : float
            Time interval size (s).

        Returns
        -------
        c_ox : float
            Updated concentration of oxidized species (M).
        c_red : float
            Updated concentration of reduced species (M).

        """

        delta_concentration = timestep * self.rate_constant * c_ox

        c_ox -= delta_concentration
        c_red += delta_concentration
        return c_ox, c_red


class MultiDegradationMechanism(DegradationMechanism):
    """
    Provides option to input multiple degradation mechanisms available within the DegradationMechanism abstract class.
    Allows for different and/or multiple mechanisms to be applied to reduced and/or oxidized species.

    Parameters
    ----------
    mechanisms_list : list[DegradationMechanism]
        List of degradation mechanism subclass instances.

    """
    def __init__(self, mechanisms_list: list[DegradationMechanism]):
        self.mechanisms_list = mechanisms_list
        for mechanism in self.mechanisms_list:
            if not isinstance(mechanism, DegradationMechanism):
                raise ValueError(f"Mechanism {mechanism} is not of type DegradationMechanism")

    def degrade(self, c_ox: float, c_red: float, timestep: float) -> tuple[float, float]:
        """
        Input of different degradation mechanisms.

        Parameters
        ----------
        c_ox : float
            Concentration of oxidized species (M).
        c_red : float
            Concentration of reduced species (M).
        timestep : float
            Time interval size (s).

        Returns
        -------
        c_ox : float
            (Possibly) updated concentration of oxidized species (M).
        c_red : float
            (Possibly) updated concentration of reduced species (M).

        """

        for mechanism in self.mechanisms_list:
            c_ox, c_red = mechanism.degrade(c_ox, c_red, timestep)
        return c_ox, c_red


"""
# to-dos
def auto_red_test(c_ox, c_red, timestep, extra_ratio=1, rate=0):
    # assumes first order process: ox --> red
    delta_concentration = timestep*rate*c_ox*extra_ratio
    concentration_red = c_red + delta_concentration
    concentration_ox = c_ox - delta_concentration
    return concentration_ox, concentration_red


# def oxygen_oxidation(concentration_o2, c_ox, c_red, timestep, rate=0):
    

# def potential_dependent(potential, c_ox, c_red, timestep, rate=0):
"""
