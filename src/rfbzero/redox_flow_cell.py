"""
Class for cell setup and declaring electrolyte parameters.
"""

from math import log

import scipy.constants as spc

from .degradation import DegradationMechanism
from .crossover import Crossover

# Faraday constant (C/mol)
F = spc.value('Faraday constant')

# Molar gas constant (J/K/mol)
R = spc.R


class ZeroDModel:
    """
    Zero dimensional model for redox flow battery (RFB) cycling [1].

    Parameters
    ----------
    volume_cls : float
        Volume of capacity-limiting side (CLS) reservoir (L).
    volume_ncls : float
        Volume of non-capacity-limiting side (NCLS) reservoir (L).
    c_ox_cls : float
        CLS initial concentration of oxidized species (M).
    c_red_cls : float
        CLS initial concentration of reduced species (M).
    c_ox_ncls : float
        NCLS initial concentration of oxidized species (M).
    c_red_ncls : float
        NCLS initial concentration of reduced species (M).
    ocv_50_soc : float
        Cell voltage (formal potentials E_+ - E_-) (V).
        If voltage > 0 then it's a Full cell.
        If voltage = 0 then it's a Symmetric cell.
    resistance : float
        Cell ohmic resistance (ohms).
    k_0_cls : float
        Electrochemical rate constant, CLS redox couple (cm/s).
    k_0_ncls : float
        Electrochemical rate constant, NCLS redox couple (cm/s).
    alpha_cls : float
        Charge transfer coefficient of CLS redox couple, dimensionless.
        Default is 0.5, which is standard in electrochemistry.
    alpha_ncls : float
        Charge transfer coefficient of NCLS redox couple, dimensionless.
        Default is 0.5, which is standard in electrochemistry.
    geometric_area : float
        Geometric area of cell (cm^2).
        Default is 5.0, a typical lab-scale cell size.
    cls_negolyte : bool
        True if negolyte is the CLS, False if posolyte is the CLS.
        Default is True.
    time_step : float
        Simulation time step (s).
        Default is 0.01, providing adequate balance of accuracy vs compute time.
    k_mt : float
        Mass transport coefficient (cm/s).
        Default is 0.8, as used in [1].
    roughness_factor : float
        Roughness factor, dimensionless.
        Total surface area divided by geometric surface area.
        Default is 26.0, as used in [1].
    num_electrons_cls : int
        Number of electrons transferred per active species molecule in the CLS.
        Default is 1.
    num_electrons_ncls : int
        Number of electrons transferred per active species molecule in the NCLS.
        Default is 1.
    temperature : float
        Temperature of battery and electrolytes (K).
        Default is 298 K (25C).

    Notes
    -----
    Most equations are adapted from [1]. If ZeroDModel has been significant to your research please cite the paper.

    [1] Modak, S.; Kwabi, D. G. A Zero-Dimensional Model for Electrochemical Behavior and Capacity Retention
    in Organic Flow Cells, Journal of The Electrochemical Society, 168, 2021, 080528. DOI 10.1149/1945-7111/ac1c1f

    """

    def __init__(
            self,
            volume_cls: float,
            volume_ncls: float,
            c_ox_cls: float,
            c_red_cls: float,
            c_ox_ncls: float,
            c_red_ncls: float,
            ocv_50_soc: float,
            resistance: float,
            k_0_cls: float,
            k_0_ncls: float,
            alpha_cls: float = 0.5,
            alpha_ncls: float = 0.5,
            geometric_area: float = 5.0,
            cls_negolyte: bool = True,
            time_step: float = 0.01,
            k_mt: float = 0.8,
            roughness_factor: float = 26.0,
            num_electrons_cls: int = 1,
            num_electrons_ncls: int = 1,
            temperature: float = 298.0,
    ) -> None:
        self.volume_cls = volume_cls
        self.volume_ncls = volume_ncls
        self.c_ox_cls = c_ox_cls
        self.c_red_cls = c_red_cls
        self.c_ox_ncls = c_ox_ncls
        self.c_red_ncls = c_red_ncls
        self.ocv_50_soc = ocv_50_soc
        self.resistance = resistance
        self.k_0_cls = k_0_cls
        self.k_0_ncls = k_0_ncls
        self.alpha_cls = alpha_cls
        self.alpha_ncls = alpha_ncls
        self.geometric_area = geometric_area
        self.cls_negolyte = cls_negolyte
        self.time_step = time_step
        self.k_mt = k_mt
        self.const_i_ex = F * roughness_factor * self.geometric_area
        self.num_electrons_cls = num_electrons_cls
        self.num_electrons_ncls = num_electrons_ncls

        self.prev_c_ox_cls = self.c_ox_cls
        self.prev_c_red_cls = self.c_red_cls
        self.prev_c_ox_ncls = self.c_ox_ncls
        self.prev_c_red_ncls = self.c_red_ncls
        self.crossed_ox_mols = 0.0
        self.crossed_red_mols = 0.0

        self.nernst_const = (R * temperature) / F

        if not isinstance(self.num_electrons_cls, int) or not isinstance(self.num_electrons_ncls, int):
            raise ValueError("'num_electrons_cls' and 'num_electrons_ncls' must be non-zero integers")

        for key, value in {'volume_cls': self.volume_cls, 'volume_ncls': self.volume_ncls, 'c_ox_cls': self.c_ox_cls,
                           'c_red_cls': self.c_red_cls, 'c_ox_ncls': self.c_ox_ncls, 'c_red_ncls': self.c_red_ncls,
                           'ocv_50_soc': self.ocv_50_soc, 'resistance': self.resistance, 'k_0_cls': self.k_0_cls,
                           'k_0_ncls': self.k_0_ncls, 'geometric_area': self.geometric_area,
                           'time_step': self.time_step, 'k_mt': self.k_mt, 'const_i_ex': self.const_i_ex,
                           'temperature': temperature}.items():

            if key not in ['ocv_50_soc', 'resistance',
                           'c_ox_cls', 'c_red_cls',
                           'c_ox_ncls', 'c_red_ncls'] and value <= 0.0:
                raise ValueError(f"'{key}' must be > 0.0")

            if value < 0.0:
                raise ValueError(f"'{key}' must be >= 0.0")

        if not 0.0 < self.alpha_cls < 1.0 or not 0.0 < self.alpha_ncls < 1.0:
            raise ValueError("Alpha parameters must be between 0.0 and 1.0")

        if self.num_electrons_cls < 1:
            raise ValueError("'num_electrons_cls' must be >= 1")

        if self.num_electrons_ncls < 1:
            raise ValueError("'num_electrons_ncls' must be >= 1")

        if self.ocv_50_soc == 0.0 and self.volume_cls >= self.volume_ncls:
            raise ValueError("'volume_cls' must be < 'volume_ncls' in a symmetric cell ('ocv_50_soc' = 0.0)")

        if self.ocv_50_soc == 0.0 and self.num_electrons_cls != self.num_electrons_ncls:
            raise ValueError("'num_electrons_cls' and 'num_electrons_ncls' must be equal (same species) "
                             "in a symmetric cell ('ocv_50_soc' = 0.0)")

        self.init_cls_capacity = self.volume_cls * self.num_electrons_cls * (self.c_ox_cls + self.c_red_cls)
        init_ncls_capacity = self.volume_ncls * self.num_electrons_ncls * (self.c_ox_ncls + self.c_red_ncls)
        if self.init_cls_capacity >= init_ncls_capacity:
            raise ValueError("Initial capacity of CLS must be less than initial capacity of NCLS")

        if self.time_step >= 1.0:
            print("WARNING: 'time_step' >= 1 second will result in very coarse data.\
                  \nzero-D model approaches theory as time step decreases.")

    def _exchange_current(self) -> tuple[float, float]:
        """
        Calculates exchange current (i_0) of redox couples in the CLS and NCLS.

        Returns
        -------
        i_0_cls : float
            Exchange current of CLS redox couple at a given time step (A).
        i_0_ncls : float
            Exchange current of NCLS redox couple at a given time step (A).

        """
        # division by 1000 for conversion from L to cm^3
        i_0_cls = (self.num_electrons_cls * self.const_i_ex * self.k_0_cls * (self.c_red_cls ** self.alpha_cls)
                   * (self.c_ox_cls ** (1 - self.alpha_cls)) * 0.001)
        i_0_ncls = (self.num_electrons_ncls * self.const_i_ex * self.k_0_ncls * (self.c_red_ncls ** self.alpha_ncls)
                    * (self.c_ox_ncls ** (1 - self.alpha_ncls)) * 0.001)
        return i_0_cls, i_0_ncls

    def _limiting_current(self, c_lim: float) -> float:
        """
        Calculates limiting current (i_lim) for a single reservoir.
        This is equation 6 of [1].

        """
        # div by 1000 for conversion from L to cm^3
        return F * self.k_mt * c_lim * self.geometric_area * 0.001

    def _limiting_concentration(self, charge: bool) -> tuple[float, float]:
        """
        Selects limiting concentration and calculates limiting current for CLS and NCLS.
        Multiplies by number of electrons transferred per molecule, for the given species.

        Parameters
        ----------
        charge : bool
            True if charging, False if discharging.

        Returns
        -------
        i_lim_cls : float
            Limiting current of CLS redox couple at a given time step (A).
        i_lim_ncls : float
            Limiting current of NCLS redox couple at a given time step (A).

        """
        if self.cls_negolyte == charge:
            i_lim_cls = self._limiting_current(self.c_ox_cls) * self.num_electrons_cls
            i_lim_ncls = self._limiting_current(self.c_red_ncls) * self.num_electrons_ncls
        else:
            i_lim_cls = self._limiting_current(self.c_red_cls) * self.num_electrons_cls
            i_lim_ncls = self._limiting_current(self.c_ox_ncls) * self.num_electrons_ncls

        return i_lim_cls, i_lim_ncls

    def _activation_overpotential(self, current: float, i_0_cls: float, i_0_ncls: float) -> float:
        """
        Calculates overall cell activation overpotential.
        This is equation 4 of [1].

        Parameters
        ----------
        current : float
            Instantaneous current flowing (A). Positive if charging, negative if discharging.
        i_0_cls : float
            Exchange current of CLS redox couple at a given time step (A).
        i_0_ncls : float
            Exchange current of NCLS redox couple at a given time step (A).

        Returns
        -------
        n_act : float
            Combined (CLS+NCLS) activation overpotential (V).

        """

        z_cls = abs(current) / (2 * i_0_cls)
        z_ncls = abs(current) / (2 * i_0_ncls)
        n_act = self.nernst_const * ((log(z_cls + ((z_cls ** 2) + 1) ** 0.5) / self.num_electrons_cls)
                                     + (log(z_ncls + ((z_ncls ** 2) + 1) ** 0.5) / self.num_electrons_ncls))
        return n_act

    def _negative_concentrations(self) -> bool:
        """Return True if any concentration is negative."""
        return any(x < 0.0 for x in [self.c_ox_cls, self.c_red_cls, self.c_ox_ncls, self.c_red_ncls])

    def __mass_transport_overpotential(self, current: float, i_lim_cls: float, i_lim_ncls: float) -> float:
        """
        Calculates overall cell mass transport overpotential.
        This is equation 8 of [1].

        Parameters
        ----------
        current : float
            Instantaneous current flowing (A). Positive if charging, negative if discharging.
         i_lim_cls : float
            Limiting current of CLS redox couple at a given time step (A).
        i_lim_ncls : float
            Limiting current of NCLS redox couple at a given time step (A).

        Returns
        -------
        n_mt : float
            Combined (CLS+NCLS) mass transport overpotential (V).

        """

        if self._negative_concentrations():
            raise ValueError('Negative concentration detected')

        c_tot_cls = self.c_red_cls + self.c_ox_cls
        c_tot_ncls = self.c_red_ncls + self.c_ox_ncls

        charge = current >= 0.0
        if self.cls_negolyte == charge:
            c1_cls = self.c_ox_cls
            c2_cls = self.c_red_cls
            c1_ncls = self.c_red_ncls
            c2_ncls = self.c_ox_ncls
        else:
            c1_cls = self.c_red_cls
            c2_cls = self.c_ox_cls
            c1_ncls = self.c_ox_ncls
            c2_ncls = self.c_red_ncls

        i = abs(current)

        n_mt = self.nernst_const * (
                (log(1 - ((c_tot_cls * i) / ((c1_cls * i_lim_cls) + (c2_cls * i)))) / self.num_electrons_cls)
                + (log(1 - ((c_tot_ncls * i) / ((c1_ncls * i_lim_ncls) + (c2_ncls * i)))) / self.num_electrons_ncls)
        )
        return n_mt * -1

    def _total_overpotential(self, current: float, i_lim_cls: float, i_lim_ncls: float) -> tuple[float, float, float]:
        """
        Calculates total cell overpotential.
        This is the sum of overpotentials of equation 2 in [1].

        Parameters
        ----------
        current : float
            Instantaneous current flowing (A). Positive if charging, negative if discharging.
        i_lim_cls : float
            Limiting current of CLS redox couple at a given time step (A).
        i_lim_ncls : float
            Limiting current of NCLS redox couple at a given time step (A).

        Returns
        -------
        total_overpotential : float
            Total cell overpotential (V).
        n_act : float
            Total activation overpotential (V).
        n_mt : float
            Total mass transport overpotential (V).

        """

        i_0_cls, i_0_ncls = self._exchange_current()
        # calculate overpotentials
        n_ohmic = abs(current) * self.resistance
        n_act = self._activation_overpotential(current, i_0_cls, i_0_ncls)
        n_mt = self.__mass_transport_overpotential(current, i_lim_cls, i_lim_ncls)

        total_overpotential = n_ohmic + n_act + n_mt

        return total_overpotential, n_act, n_mt

    def _open_circuit_voltage(self) -> float:
        """
        Nernstian calculation of the cell open circuit voltage.
        This is equivalent to equation 3 of [1].

        Returns
        -------
        ocv : float
            Cell open circuit voltage (V).

        """

        if self._negative_concentrations():
            raise ValueError('Negative concentration detected')

        direction = 1 if self.cls_negolyte else -1

        ocv = (self.ocv_50_soc
               + direction
               * (((self.nernst_const / self.num_electrons_cls) * log(self.c_red_cls / self.c_ox_cls))
                  + ((self.nernst_const / self.num_electrons_ncls) * log(self.c_ox_ncls / self.c_red_ncls))))
        return ocv

    @staticmethod
    def _cell_voltage(ocv: float, total_overpotential: float, charge: bool) -> float:
        """If charging, add overpotentials to OCV, else subtract them."""
        return ocv + total_overpotential if charge else ocv - total_overpotential

    def _coulomb_counter(
            self,
            current: float,
            cls_degradation: DegradationMechanism = None,
            ncls_degradation: DegradationMechanism = None,
            cross_over: Crossover = None,
    ) -> tuple[dict[str, float], dict[str, float]]:
        """
        Updates all species' concentrations at each time step. Contributions from faradaic current, (optional)
        degradation mechanisms, and (optional) crossover mechanism.

        Parameters
        ----------
        current : float
            Instantaneous current flowing (A). Positive if charging, negative if discharging.
        cls_degradation : DegradationMechanism, optional
            Degradation mechanism for CLS.
        ncls_degradation: DegradationMechanism, optional
            Degradation mechanism for NCLS.
        cross_over : Crossover, optional
            Crossover class instance.

        Returns
        -------
        c_products_cls : dict[str, float]
            Updated concentrations of all CLS product species (M).
        c_products_ncls : dict[str, float]
            Updated concentrations of all NCLS product species (M).

        """

        # Change in concentration from coulomb counting based solely on current
        direction = 1 if self.cls_negolyte else -1
        delta_cls = ((self.time_step * current) / (F * self.num_electrons_cls * self.volume_cls)) * direction
        delta_ncls = ((self.time_step * current) / (F * self.num_electrons_ncls * self.volume_ncls)) * direction

        self.prev_c_ox_cls = self.c_ox_cls
        self.prev_c_red_cls = self.c_red_cls
        self.prev_c_ox_ncls = self.c_ox_ncls
        self.prev_c_red_ncls = self.c_red_ncls

        # update CLS and NCLS concentrations
        new_c_ox_cls = self.c_ox_cls - delta_cls
        new_c_red_cls = self.c_red_cls + delta_cls
        new_c_ox_ncls = self.c_ox_ncls + delta_ncls
        new_c_red_ncls = self.c_red_ncls - delta_ncls

        # for no crossover situation
        crossed_ox_mols = 0.0
        crossed_red_mols = 0.0

        # Coulomb counting from optional degradation and/or crossover mechanisms
        c_products_cls = {}
        if cls_degradation is not None:
            delta_ox_cls, delta_red_cls = cls_degradation.degrade(self.c_ox_cls, self.c_red_cls,
                                                                  self.time_step)
            new_c_ox_cls += delta_ox_cls
            new_c_red_cls += delta_red_cls
            c_products_cls = cls_degradation.c_products

        c_products_ncls = {}
        if ncls_degradation is not None:
            delta_ox_ncls, delta_red_ncls = ncls_degradation.degrade(self.c_ox_ncls, self.c_red_ncls,
                                                                     self.time_step)
            new_c_ox_ncls += delta_ox_ncls
            new_c_red_ncls += delta_red_ncls
            c_products_ncls = ncls_degradation.c_products

        if cross_over is not None:
            delta_ox_cls, delta_red_cls, delta_ox_ncls, delta_red_ncls, crossed_ox_mols, crossed_red_mols = \
                cross_over.crossover(self.geometric_area, self.c_ox_cls, self.c_red_cls, self.c_ox_ncls,
                                     self.c_red_ncls, self.volume_cls, self.volume_ncls, self.time_step)
            new_c_ox_cls += delta_ox_cls
            new_c_red_cls += delta_red_cls
            new_c_ox_ncls += delta_ox_ncls
            new_c_red_ncls += delta_red_ncls

        # Update new concentrations to self
        self.c_ox_cls = new_c_ox_cls
        self.c_red_cls = new_c_red_cls
        self.c_ox_ncls = new_c_ox_ncls
        self.c_red_ncls = new_c_red_ncls

        self.crossed_ox_mols = crossed_ox_mols
        self.crossed_red_mols = crossed_red_mols

        return c_products_cls, c_products_ncls

    def _revert_concentrations(self) -> None:
        """Resets concentrations to previous value if a (invalid) negative concentration is calculated."""
        self.c_ox_cls = self.prev_c_ox_cls
        self.c_red_cls = self.prev_c_red_cls
        self.c_ox_ncls = self.prev_c_ox_ncls
        self.c_red_ncls = self.prev_c_red_ncls
