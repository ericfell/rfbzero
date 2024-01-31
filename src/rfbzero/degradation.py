"""
Classes for defining capacity fade mechanisms.
"""

from abc import ABC, abstractmethod
from functools import reduce
from operator import ior


class DegradationMechanism(ABC):
    """Abstract base class to be implemented by specific degradation mechanisms."""

    def __init__(self, **c_products: float):
        self.c_products = c_products

    @abstractmethod
    def degrade(self, c_ox: float, c_red: float, time_step: float) -> tuple[float, float]:
        """Applies desired degradation mechanisms to oxidized/reduced species at each time step."""
        raise NotImplementedError


class ChemicalDegradationOxidized(DegradationMechanism):
    """
    Provides an N-th order chemical degradation mechanism for an oxidized species.

    Parameters
    ----------
    rate_order : int
        Rate order for chemical degradation reaction.
    rate_constant : float
        N-th order rate constant of chemical degradation of oxidized species (unit is rate order-dependent).

    """

    def __init__(self, rate_order: int, rate_constant: float) -> None:
        super().__init__()
        self.rate_order = rate_order
        self.rate_constant = rate_constant

        if self.rate_order < 0:
            raise ValueError("'rate_order' must be >= 0")

        if self.rate_constant <= 0.0:
            raise ValueError("'rate_constant' must be > 0.0")

    def degrade(self, c_ox: float, c_red: float, time_step: float) -> tuple[float, float]:
        """
        Applies N-th order chemical degradation of N*[redox-active species] --> [redox-inactive species] at each
        time step. Returns updated concentrations; concentration may be unchanged if species does not degrade.

        Parameters
        ----------
        c_ox : float
            Concentration of oxidized species (M).
        c_red : float
            Concentration of reduced species (M).
        time_step : float
            Time interval size (s).

        Returns
        -------
        delta_ox : float
            Change in concentration of oxidized species (M).
        delta_red : float
            Change in concentration of reduced species (M). This will always be zero.

        """

        delta_ox = -(time_step * self.rate_constant * (c_ox ** self.rate_order))
        delta_red = 0.0
        return delta_ox, delta_red


class ChemicalDegradationReduced(DegradationMechanism):
    """
    Provides an N-th order chemical degradation mechanism for a reduced species.

    Parameters
    ----------
    rate_order : int
        Rate order for chemical degradation reaction.
    rate_constant : float
        N-th order rate constant of chemical degradation of reduced species (unit is rate order-dependent).

    """

    def __init__(self, rate_order: int, rate_constant: float) -> None:
        super().__init__()
        self.rate_order = rate_order
        self.rate_constant = rate_constant

        if self.rate_order < 0:
            raise ValueError("'rate_order' must be >= 0")

        if self.rate_constant <= 0.0:
            raise ValueError("'rate_constant' must be > 0.0")

    def degrade(self, c_ox: float, c_red: float, time_step: float) -> tuple[float, float]:
        """
        Applies N-th order chemical degradation of N*[redox-active species] --> [redox-inactive species] at each
        time step. Returns updated concentrations; concentration may be unchanged if species does not degrade.

        Parameters
        ----------
        c_ox : float
            Concentration of oxidized species (M).
        c_red : float
            Concentration of reduced species (M).
        time_step : float
            Time interval size (s).

        Returns
        -------
        delta_ox : float
            Change in concentration of oxidized species (M). This will always be zero.
        delta_red : float
            Change in concentration of reduced species (M).

        """

        delta_ox = 0.0
        delta_red = -(time_step * self.rate_constant * (c_red ** self.rate_order))
        return delta_ox, delta_red


class AutoOxidation(DegradationMechanism):
    """
    Provides a 1st order auto-oxidation mechanism, (red --> ox) with no loss of active material. This can be thought
    of as a chemical oxidation of the redox-active, balanced by an oxidant e.g., water splitting (HER). This could
    occur in a low-potential negolyte and be considered a self-discharge. If it is desired for the concentration of
    the oxidant to affect the chemical oxidation rate, the initial oxidant concentration and the stoichiometric factor
    i.e., red + n*oxidant --> ox + ..., can be input. This could simulate the effect of H2 leaving the system.

    Parameters
    ----------
    rate_constant : float
        First order rate of auto-oxidation (1/s).
    c_oxidant : float
        Initial concentration of oxidant, defaults to 0.0 (M).
    oxidant_stoich : int
        Number of oxidants involved in the chemical oxidation of the redox-active species, defaults to 0.

    """

    def __init__(self, rate_constant: float, c_oxidant: float = 0.0, oxidant_stoich: int = 0) -> None:
        super().__init__()
        self.rate_constant = rate_constant
        self.c_oxidant = c_oxidant
        self.oxidant_stoich = oxidant_stoich

        if self.rate_constant <= 0.0:
            raise ValueError("'rate_constant' must be > 0.0")

        if self.c_oxidant < 0.0:
            raise ValueError("'c_oxidant' must be >= 0.0")

        if self.oxidant_stoich < 0:
            raise ValueError("'oxidant_stoich' must be >= 0")

        if (self.oxidant_stoich > 0 and c_oxidant == 0.0) or (self.oxidant_stoich == 0 and c_oxidant > 0.0):
            raise ValueError("'c_oxidant' and 'oxidant_stoich' must both be zero, or both be positive")

    def degrade(self, c_ox: float, c_red: float, time_step: float) -> tuple[float, float]:
        """
        Applies an auto-oxidation mechanism to oxidized/reduced species at each time step.
        Defaults to first order process: red --> ox.

        Parameters
        ----------
        c_ox : float
            Concentration of oxidized species (M).
        c_red : float
            Concentration of reduced species (M).
        time_step : float
            Time interval size (s).

        Returns
        -------
        delta_ox : float
            Change in concentration of oxidized species (M).
        delta_red : float
            Change in concentration of reduced species (M).

        """

        delta_concentration = time_step * self.rate_constant * c_red * (self.c_oxidant ** self.oxidant_stoich)

        delta_ox = delta_concentration
        delta_red = -delta_concentration
        self.c_oxidant -= delta_concentration * self.oxidant_stoich
        self.c_oxidant = max(self.c_oxidant, 0.0)

        return delta_ox, delta_red


class AutoReduction(DegradationMechanism):
    """
    Provides a 1st order auto-reduction mechanism, (ox --> red) with no loss of active material. This can be thought
    of as a chemical reduction of the redox-active, balanced by a reductant e.g., water splitting (OER). This could
    occur in a high-potential posolyte and be considered a self-discharge. If it is desired for the concentration of
    the reductant to affect the chemical reduction rate, the initial reductant concentration and the stoichiometric
    factor i.e., ox + n*reductant --> red + ..., can be input. This could simulate the effect of O2 leaving the system.

    Parameters
    ----------
    rate_constant : float
        First order rate of auto-reduction (1/s).
    c_reductant : float
        Initial concentration of reductant, defaults to 0.0 (M).
    reductant_stoich : int
        Number of reductants involved in the chemical reduction of the redox-active species, defaults to 0.

    """
    def __init__(self, rate_constant: float, c_reductant: float = 0.0, reductant_stoich: int = 0) -> None:
        super().__init__()
        self.rate_constant = rate_constant
        self.c_reductant = c_reductant
        self.reductant_stoich = reductant_stoich

        if self.rate_constant <= 0.0:
            raise ValueError("'rate_constant' must be > 0.0")

        if self.c_reductant < 0.0:
            raise ValueError("'c_reductant' must be >= 0.0")

        if self.reductant_stoich < 0:
            raise ValueError("'reductant_stoich' must be >= 0")

        if (self.reductant_stoich > 0 and c_reductant == 0.0) or (self.reductant_stoich == 0 and c_reductant > 0.0):
            raise ValueError("'c_reductant' and 'reductant_stoich' must both be zero, or both be positive")

    def degrade(self, c_ox: float, c_red: float, time_step: float) -> tuple[float, float]:
        """
        Applies an auto-reduction mechanism to oxidized/reduced species at each time step.
        Defaults to first order process: ox --> red.

        Parameters
        ----------
        c_ox : float
            Concentration of oxidized species (M).
        c_red : float
            Concentration of reduced species (M).
        time_step : float
            Time interval size (s).

        Returns
        -------
        delta_ox : float
            Change in concentration of oxidized species (M).
        delta_red : float
            Change in concentration of reduced species (M).

        """

        delta_concentration = time_step * self.rate_constant * c_ox * (self.c_reductant ** self.reductant_stoich)

        delta_ox = -delta_concentration
        delta_red = delta_concentration
        self.c_reductant -= delta_concentration * self.reductant_stoich
        self.c_reductant = max(self.c_reductant, 0.0)

        return delta_ox, delta_red


class Dimerization(DegradationMechanism):
    """
    Provides a reversible dimerization mechanism: ox + red <--> dimer.

    Parameters
    ----------
    forward_rate_constant : float
        Second order rate constant for forward reaction (1/(M s)).
    backward_rate_constant : float
        First order rate constant for backward reaction (1/s).
    c_dimer : float
        Initial concentration of dimer, defaults to 0.0 (M).

    """

    def __init__(self, forward_rate_constant: float, backward_rate_constant: float, c_dimer: float = 0.0) -> None:
        super().__init__(c_dimer=c_dimer)
        self.forward_rate_constant = forward_rate_constant
        self.backward_rate_constant = backward_rate_constant

        if self.forward_rate_constant <= 0.0:
            raise ValueError("'forward_rate_constant' must be > 0.0")

        if self.backward_rate_constant <= 0.0:
            raise ValueError("'backward_rate_constant' must be > 0.0")

        if c_dimer < 0.0:
            raise ValueError("'c_dimer' must be >= 0.0")

    def degrade(self, c_ox: float, c_red: float, time_step: float) -> tuple[float, float]:
        """
        Applies a reversible dimerization mechanism to oxidized/reduced species at each time step.
        Returns updated concentrations.

        Parameters
        ----------
        c_ox : float
            Concentration of oxidized species (M).
        c_red : float
            Concentration of reduced species (M).
        time_step : float
            Time interval size (s).

        Returns
        -------
        delta_ox : float
            Change in concentration of oxidized species (M).
        delta_red : float
            Change in concentration of reduced species (M).

        """

        delta_concentration = time_step * (
                (self.forward_rate_constant * c_ox * c_red) - (self.backward_rate_constant * self.c_products['c_dimer'])
        )

        delta_ox = -delta_concentration
        delta_red = -delta_concentration
        self.c_products['c_dimer'] += delta_concentration

        return delta_ox, delta_red


class MultiDegradationMechanism(DegradationMechanism):
    """
    Provides usage of multiple degradation mechanisms that implement the DegradationMechanism abstract base class.
    Allows for different and/or multiple mechanisms to be applied to reduced and/or oxidized species. Degradation
    mechanisms are applied in the same order as the input list.

    Parameters
    ----------
    mechanisms : list[DegradationMechanism]
        List of degradation mechanism subclass instances. Note that multiple instances of the same mechanism type cannot
        be used in the same MultiDegradationMechanism.

    """
    def __init__(self, mechanisms: list[DegradationMechanism]) -> None:
        self.mechanisms = mechanisms

        c_products: dict[str, float] = {}
        for mechanism in self.mechanisms:
            if not isinstance(mechanism, DegradationMechanism):
                raise ValueError(f"Mechanism {mechanism} is not of type DegradationMechanism")
            # Compute the union of c_products across all mechanisms
            c_products |= mechanism.c_products

        super().__init__(**c_products)

    def degrade(self, c_ox: float, c_red: float, time_step: float) -> tuple[float, float]:
        """
        Applies multiple degradation mechanisms to oxidized/reduced species at each time step.
        Concentration may be unchanged if species does not degrade.

        Parameters
        ----------
        c_ox : float
            Concentration of oxidized species (M).
        c_red : float
            Concentration of reduced species (M).
        time_step : float
            Time interval size (s).

        Returns
        -------
        delta_ox : float
            Change in concentration of oxidized species (M).
        delta_red : float
            Change in concentration of reduced species (M).

        """

        deltas = [mechanism.degrade(c_ox, c_red, time_step) for mechanism in self.mechanisms]
        delta_ox, delta_red = (sum(x) for x in zip(*deltas))

        # Update the c_products dictionary to be the union of all mechanisms' dictionaries
        self.c_products = reduce(ior, [mechanism.c_products for mechanism in self.mechanisms], {})

        return delta_ox, delta_red
