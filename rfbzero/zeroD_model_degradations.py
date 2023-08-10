
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
    to-do
    """

    def __init__(self, rate_order: int, rate: float, species: str = 'red'):#, reservoir: str = 'both'):
        self.rate_order = rate_order
        self.rate = rate
        self.species = species
        #self.reservoir = reservoir
        #if self.reservoir not in ['both', 'cls', 'ncls']:
        #    raise ValueError("Options: 'both', 'cls', or 'ncls' ")

    def degrade(self, c_ox: float, c_red: float, timestep: float) -> tuple[float, float]:
        """
        nth order chemical degradation of n[species] --> [redox-inactive species]

        Parameters
        ----------
        c_ox
        c_red
        timestep
        is_cls

        Returns
        -------

        """
        # checks if degradation is applied to ClS, NCLS, or both reservoirs
        #if (self.reservoir == 'cls' and not is_cls) or (self.reservoir == 'ncls' and is_cls):
        #    return c_ox, c_red

        if self.species == 'red':
            concentration_red = c_red - (timestep * self.rate * (c_red**self.rate_order))
            return c_ox, concentration_red
        else:
            concentration_ox = c_ox - (timestep * self.rate * (c_ox**self.rate_order))
            return concentration_ox, c_red


class AutoOxidation(DegradationMechanism):
    """
    to-dos
    """
    def __init__(self, rate: float):#, reservoir: str = 'both'):
        self.rate = rate
        #self.reservoir = reservoir
        #if self.reservoir not in ['both', 'cls', 'ncls']:
        #    raise ValueError("Options: 'both', 'cls', or 'ncls' ")

    def degrade(self, c_ox: float, c_red: float, timestep: float) -> tuple[float, float]:
        """assumes first order process: red --> ox"""
        # checks if degradation is applied to ClS, NCLS, or both reservoirs

        #if (self.reservoir == 'cls' and not is_cls) or (self.reservoir == 'ncls' and is_cls):
        #    return c_ox, c_red

        delta_concentration = timestep * self.rate * c_red

        concentration_red = c_red - delta_concentration
        concentration_ox = c_ox + delta_concentration
        return concentration_ox, concentration_red


class AutoReduction(DegradationMechanism):
    """
    to-dos
    """
    def __init__(self, rate: float):#, reservoir: str = 'both'):
        self.rate = rate
        #self.reservoir = reservoir
        #if self.reservoir not in ['both', 'cls', 'ncls']:
        #    raise ValueError("Options: 'both', 'cls', or 'ncls' ")

    def degrade(self, c_ox: float, c_red: float, timestep: float) -> tuple[float, float]:
        """assumes first order process: ox --> red"""
        # checks if degradation is applied to ClS, NCLS, or both reservoirs

        #if (self.reservoir == 'cls' and not is_cls) or (self.reservoir == 'ncls' and is_cls):
        #    return c_ox, c_red

        delta_concentration = timestep * self.rate * c_ox

        concentration_red = c_red + delta_concentration
        concentration_ox = c_ox - delta_concentration
        return concentration_ox, concentration_red


class MultiDegradationMechanism(DegradationMechanism):
    """
    to-dos
    """
    def __init__(self, mechanisms: list[DegradationMechanism]):
        self.mechanisms = mechanisms

    def degrade(self, c_ox: float, c_red: float, timestep: float) -> tuple[float, float]:

        for mechanism in self.mechanisms:
            c_ox, c_red = mechanism.degrade(c_ox, c_red, timestep)#, is_cls)

        return c_ox, c_red


"""
# to-dos
def auto_red_test(conc_ox_t, conc_red_t, timestep, extra_ratio=1, rate=0):
    # assumes first order process: ox --> red
    delta_conc = timestep*rate*conc_ox_t*extra_ratio
    concentration_red = conc_red_t + delta_conc
    concentration_ox = conc_ox_t - delta_conc
    return concentration_ox, concentration_red


# def oxygen_oxidation(conc_o2, conc_ox_t, conc_red_t, timestep, rate=0):
    

# def potential_dependent(potential, conc_ox_t, conc_red_t, timestep, rate=0):
"""

if __name__ == '__main__':
    print('testing')
