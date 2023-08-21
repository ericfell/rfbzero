
"""
DEGRADATION FUNCTIONS CALLED BY ZeroDmodel CLASS
"""

from abc import ABC, abstractmethod


class DegradationMechanism(ABC):
    @abstractmethod
    def degrade(self, c_ox: float, c_red: float, timestep: float) -> tuple[float, float]:
        raise NotImplementedError


class ChemicalDegradation(DegradationMechanism):
    """
    Subclass for a chemical degradation mechanism.

    Parameters
    ----------
    rate_order : int
        Rate order for chemical degradation.
    rate : float
        Nth order rate of chemical oxidation (unit is rate-dependent).
    species : str
        Species ('red' or 'ox') undergoing chemical degradation.

    """

    def __init__(self, rate_order: int, rate: float, species: str = 'red'):
        self.rate_order = rate_order
        self.rate = rate
        self.species = species

        if not isinstance(self.rate_order, int):
            raise ValueError("Rate order must be an integer")

        if self.rate <= 0.0:
            raise ValueError("Rate must be a non-zero, positive value")

        if self.species not in ['red', 'ox']:
            raise ValueError("'species' options: 'red', 'ox' ")

    def degrade(self, c_ox: float, c_red: float, timestep: float) -> tuple[float, float]:
        """
        nth order chemical degradation of n[species] --> [redox-inactive species]

        Parameters
        ----------
        c_ox : float
            Concentration of oxidized species (M).
        c_red : float
            Concentration of reduced species (M).
        timestep : float
            Simulation time step (s).

        Returns
        -------
        concentration_ox : float
            Updated concentration of oxidized species (M).
        concentration_red : float
            Updated concentration of reduced species (M).

        """

        if self.species == 'red':
            concentration_red = c_red - (timestep * self.rate * (c_red**self.rate_order))
            return c_ox, concentration_red
        else:
            concentration_ox = c_ox - (timestep * self.rate * (c_ox**self.rate_order))
            return concentration_ox, c_red


class AutoOxidation(DegradationMechanism):
    """
    Subclass for an auto-oxidation mechanism,
    (red --> ox) with no loss of active material.

    Parameters
    ----------
    rate : float
        First order rate of auto-oxidation (1/s).

    """
    def __init__(self, rate: float):
        self.rate = rate
        if self.rate <= 0.0:
            raise ValueError("Rate must be a non-zero, positive value")

    def degrade(self, c_ox: float, c_red: float, timestep: float) -> tuple[float, float]:
        """assumes first order process: red --> ox"""

        delta_concentration = timestep * self.rate * c_red

        concentration_red = c_red - delta_concentration
        concentration_ox = c_ox + delta_concentration
        return concentration_ox, concentration_red


class AutoReduction(DegradationMechanism):
    """
    Subclass for an auto-reduction mechanism,
    (ox --> red) with no loss of active material.

    Parameters
    ----------
    rate : float
        First order rate of auto-oxidation (1/s).

    """
    def __init__(self, rate: float):
        self.rate = rate
        if self.rate <= 0.0:
            raise ValueError("Rate must be a non-zero, positive value")

    def degrade(self, c_ox: float, c_red: float, timestep: float) -> tuple[float, float]:
        """assumes first order process: ox --> red"""

        delta_concentration = timestep * self.rate * c_ox

        concentration_red = c_red + delta_concentration
        concentration_ox = c_ox - delta_concentration
        return concentration_ox, concentration_red


class MultiDegradationMechanism(DegradationMechanism):
    """
    Subclass for input of multiple degradation mechanisms.

    Parameters
    ----------
    mechanisms : list
        List of degradation mechanism subclass instances.

    """
    def __init__(self, mechanisms: list[DegradationMechanism]):
        self.mechanisms = mechanisms
        for i in self.mechanisms:
            if not isinstance(i, DegradationMechanism):
                raise ValueError("Mechanism is not of type DegradationMechanism")

    def degrade(self, c_ox: float, c_red: float, timestep: float) -> tuple[float, float]:

        for mechanism in self.mechanisms:
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

if __name__ == '__main__':
    print('testing')
