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

# TODO make temperature a variable
TEMPERATURE = 298  # Kelvins, for S.T.P.
NERNST_CONST = (R * TEMPERATURE) / F


class ZeroDModel:
    """
    Zero dimensional model for redox flow battery (RFB) cycling [1].

    Parameters
    ----------
    cls_volume : float
        Volume of capacity-limiting side (CLS) reservoir (L).
    ncls_volume : float
        Volume of non-capacity-limiting side (NCLS) reservoir (L).
    cls_start_c_ox : float
        CLS initial concentration of oxidized species (M).
    cls_start_c_red : float
        CLS initial concentration of reduced species (M).
    ncls_start_c_ox : float
        NCLS initial concentration of oxidized species (M).
    ncls_start_c_red : float
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
    time_increment : float
        Simulation time step (s).
        Default is 0.01, providing adequate balance of accuracy vs compute time.
    k_mt : float
        Mass transport coefficient (cm/s).
        Default is 0.8, as used in [1].
    roughness_factor : float
        Roughness factor, dimensionless.
        Total surface area divided by geometric surface area.
        Default is 26.0, as used in [1].
    n_cls : int
        Number of electrons transferred per active species molecule in the CLS.
    n_ncls : int
        Number of electrons transferred per active species molecule in the NCLS.


    Notes
    -----
    Most equations are adapted from [1]. If ZeroDModel has been significant to your research please cite the paper.

    [1] Modak, S.; Kwabi, D. G. A Zero-Dimensional Model for Electrochemical Behavior and Capacity Retention
    in Organic Flow Cells, Journal of The Electrochemical Society, 168, 2021, 080528. DOI 10.1149/1945-7111/ac1c1f

    """

    def __init__(
            self,
            cls_volume: float,
            ncls_volume: float,
            cls_start_c_ox: float,
            cls_start_c_red: float,
            ncls_start_c_ox: float,
            ncls_start_c_red: float,
            ocv_50_soc: float,
            resistance: float,
            k_0_cls: float,
            k_0_ncls: float,
            alpha_cls: float = 0.5,
            alpha_ncls: float = 0.5,
            geometric_area: float = 5.0,
            cls_negolyte: bool = True,
            time_increment: float = 0.01,
            k_mt: float = 0.8,
            roughness_factor: float = 26.0,
            n_cls: int = 1,
            n_ncls: int = 1
    ) -> None:
        self.cls_volume = cls_volume
        self.ncls_volume = ncls_volume
        self.c_ox_cls = cls_start_c_ox
        self.c_red_cls = cls_start_c_red
        self.c_ox_ncls = ncls_start_c_ox
        self.c_red_ncls = ncls_start_c_red
        self.ocv_50_soc = ocv_50_soc
        self.resistance = resistance
        self.k_0_cls = k_0_cls
        self.k_0_ncls = k_0_ncls
        self.alpha_cls = alpha_cls
        self.alpha_ncls = alpha_ncls
        self.geometric_area = geometric_area
        self.cls_negolyte = cls_negolyte
        self.time_increment = time_increment
        self.k_mt = k_mt
        self.const_i_ex = F * roughness_factor * self.geometric_area
        self.n_cls = n_cls
        self.n_ncls = n_ncls

        self.prev_c_ox_cls = self.c_ox_cls
        self.prev_c_red_cls = self.c_red_cls
        self.prev_c_ox_ncls = self.c_ox_ncls
        self.prev_c_red_ncls = self.c_red_ncls
        self.delta_ox = 0.0
        self.delta_red = 0.0

        for key, value in {'cls_volume': self.cls_volume, 'ncls_volume': self.ncls_volume, 'k_0_cls': self.k_0_cls,
                           'cls_start_c_ox': self.c_ox_cls, 'cls_start_c_red': self.c_red_cls,
                           'ncls_start_c_ox': self.c_ox_ncls, 'ncls_start_c_red': self.c_red_ncls,
                           'k_0_ncls': self.k_0_ncls, 'geometric_area': self.geometric_area,
                           'time_increment': self.time_increment, 'k_mt': self.k_mt, 'const_i_ex': self.const_i_ex,
                           'ocv_50_soc': self.ocv_50_soc, 'resistance': self.resistance, 'n_cls': self.n_cls,
                           'n_ncls': self.n_ncls}.items():

            if key not in ['ocv_50_soc', 'resistance',
                           'cls_start_c_ox', 'cls_start_c_red',
                           'ncls_start_c_ox', 'ncls_start_c_red'] and value <= 0.0:
                raise ValueError(f"'{key}' must be > 0.0")

            if value < 0.0:
                raise ValueError(f"'{key}' must be >= 0.0")

        if not 0.0 < self.alpha_cls < 1.0 or not 0.0 < self.alpha_ncls < 1.0:
            raise ValueError("Alpha parameters must be between 0.0 and 1.0")

        if not isinstance(self.n_cls, int) or not isinstance(self.n_ncls, int):
            raise ValueError("'n_cls' and 'n_ncls' must be integers")

        if self.ocv_50_soc == 0.0 and self.cls_volume >= self.ncls_volume:
            raise ValueError("'cls_volume' must be < 'ncls_volume' in a symmetric cell")

        if self.ocv_50_soc == 0.0 and self.n_cls != self.n_ncls:
            raise ValueError("Symmetric cell (0 volt OCV) requires 'n_cls' and 'n_ncls' to be equal (same species)")

        self.init_cls_capacity = self.cls_volume * self.n_cls * (self.c_ox_cls + self.c_red_cls)
        init_ncls_capacity = self.ncls_volume * self.n_ncls * (self.c_ox_ncls + self.c_red_ncls)
        if self.init_cls_capacity >= init_ncls_capacity:
            raise ValueError("Initial capacity of CLS must be less than initial capacity of NCLS")

        if self.time_increment >= 1.0:
            print("WARNING: 'time_increment' >= 1 second will result in very coarse data.\
                  \nzero-D model approaches theory as timestep decreases.")

    def _exchange_current(self) -> tuple[float, float]:
        """
        Calculates exchange current (i_0) of redox couples in the CLS and NCLS.
        Value returned is in Amps.

        Returns
        -------
        i_0_cls : float
            Exchange current of CLS redox couple at a given timestep (A).
        i_0_ncls : float
            Exchange current of NCLS redox couple at a given timestep (A).

        """
        # division by 1000 for conversion from L to cm^3
        i_0_cls = (self.n_cls * self.const_i_ex * self.k_0_cls * (self.c_red_cls ** self.alpha_cls)
                   * (self.c_ox_cls ** (1 - self.alpha_cls)) * 0.001)
        i_0_ncls = (self.n_ncls * self.const_i_ex * self.k_0_ncls * (self.c_red_ncls ** self.alpha_ncls)
                    * (self.c_ox_ncls ** (1 - self.alpha_ncls)) * 0.001)
        return i_0_cls, i_0_ncls

    def _limiting_current(self, c_lim: float) -> float:
        """
        Calculates limiting current (i_lim) for a single reservoir.
        Value returned is in Amps. This is equation 6 of [1].

        """
        # div by 1000 for conversion from L to cm^3
        return F * self.k_mt * c_lim * self.geometric_area * 0.001

    def limiting_concentration(self, charge: bool) -> tuple[float, float]:
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
            Limiting current of CLS redox couple at a given timestep (A).
        i_lim_ncls : float
            Limiting current of NCLS redox couple at a given timestep (A).

        """
        if self.cls_negolyte == charge:
            i_lim_cls = self._limiting_current(self.c_ox_cls) * self.n_cls
            i_lim_ncls = self._limiting_current(self.c_red_ncls) * self.n_ncls
        else:
            i_lim_cls = self._limiting_current(self.c_red_cls) * self.n_cls
            i_lim_ncls = self._limiting_current(self.c_ox_ncls) * self.n_ncls

        return i_lim_cls, i_lim_ncls

    def _activation_overpotential(self, current: float, i_0_cls: float, i_0_ncls: float) -> float:
        """
        Calculates overall cell activation overpotential.
        This is equation 4 of [1].

        Parameters
        ----------
        current : float
            Instantaneous current flowing (A).
        i_0_cls : float
            Exchange current of CLS redox couple at a given timestep (A).
        i_0_ncls : float
            Exchange current of NCLS redox couple at a given timestep (A).

        Returns
        -------
        n_act : float
            Combined (CLS+NCLS) activation overpotential (V).

        """

        z_cls = abs(current) / (2 * i_0_cls)
        z_ncls = abs(current) / (2 * i_0_ncls)
        n_act = NERNST_CONST * ((log(z_cls + ((z_cls ** 2) + 1) ** 0.5) / self.n_cls)
                                + (log(z_ncls + ((z_ncls ** 2) + 1) ** 0.5) / self.n_ncls))
        return n_act

    def negative_concentrations(self) -> bool:
        """Return True if any concentration is negative."""
        return any(x < 0.0 for x in [self.c_ox_cls, self.c_red_cls, self.c_ox_ncls, self.c_red_ncls])

    def _mass_transport_overpotential(self, current: float, i_lim_cls: float, i_lim_ncls: float) -> float:
        """
        Calculates overall cell mass transport overpotential.
        This is equation 8 of [1].

        Parameters
        ----------
        current : float
            Instantaneous current flowing (A). Positive if charging, negative if discharging.
         i_lim_cls : float
            Limiting current of CLS redox couple at a given timestep (A).
        i_lim_ncls : float
            Limiting current of NCLS redox couple at a given timestep (A).

        Returns
        -------
        n_mt : float
            Combined (CLS+NCLS) mass transport overpotential (V).

        """

        if self.negative_concentrations():
            raise ValueError('Negative concentration detected')

        c_tot_cls = self.c_red_cls + self.c_ox_cls
        c_tot_ncls = self.c_red_ncls + self.c_ox_ncls

        charge = current >= 0
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

        n_mt = NERNST_CONST * ((log(1 - ((c_tot_cls * i) / ((c1_cls * i_lim_cls) + (c2_cls * i)))) / self.n_cls)
                               + (log(1 - ((c_tot_ncls * i) / ((c1_ncls * i_lim_ncls) + (c2_ncls * i)))) / self.n_ncls))
        return n_mt * -1

    def total_overpotential(self, current: float, i_lim_cls: float, i_lim_ncls: float) -> tuple[float, float, float]:
        """
        Calculates total cell overpotential.
        This is the sum of overpotentials of equation 2 in [1].

        Parameters
        ----------
        current : float
            Instantaneous current flowing (A). Positive if charging, negative if discharging.
        i_lim_cls : float
            Limiting current of CLS redox couple at a given timestep (A).
        i_lim_ncls : float
            Limiting current of NCLS redox couple at a given timestep (A).

        Returns
        -------
        n_loss : float
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
        n_mt = self._mass_transport_overpotential(current, i_lim_cls, i_lim_ncls)

        n_loss = n_ohmic + n_act + n_mt

        return n_loss, n_act, n_mt

    def open_circuit_voltage(self) -> float:
        """
        Nernstian calculation of the cell open circuit voltage.
        This is equivalent to equation 3 of [1].

        Returns
        -------
        ocv : float
            Cell open circuit voltage (V).

        """

        if self.negative_concentrations():
            raise ValueError('Negative concentration detected')

        direction = 1 if self.cls_negolyte else -1

        ocv = (self.ocv_50_soc
               + direction * (((NERNST_CONST / self.n_cls) * log(self.c_red_cls / self.c_ox_cls))
                              + ((NERNST_CONST / self.n_ncls) * log(self.c_ox_ncls / self.c_red_ncls))))
        return ocv

    @staticmethod
    def cell_voltage(ocv: float, losses: float, charge: bool) -> float:
        """If charging, add overpotentials to OCV, else subtract them."""
        return ocv + losses if charge else ocv - losses

    def coulomb_counter(
            self,
            current: float,
            cls_degradation: DegradationMechanism = None,
            ncls_degradation: DegradationMechanism = None,
            cross_over: Crossover = None
    ) -> None:
        """
        Updates all species' concentrations at each timestep. Contributions from faradaic current, (optional)
        degradation mechanisms, and (optional) crossover mechanism.

        Parameters
        ----------
        current : float
            Instantaneous current flowing (A).
        cls_degradation : DegradationMechanism, optional
            Degradation mechanism for CLS.
        ncls_degradation: DegradationMechanism, optional
            Degradation mechanism for NCLS.
        cross_over : Crossover, optional
            Crossover class instance.

        Returns
        -------
        delta_ox : float
            Concentration difference (CLS-NCLS) of oxidized species (M).
        delta_red : float
            Concentration difference (CLS-NCLS) of reduced species (M).

        """

        # Change in concentration from coulomb counting based solely on current
        direction = 1 if self.cls_negolyte else -1
        delta_cls = ((self.time_increment * current) / (F * self.n_cls * self.cls_volume)) * direction
        delta_ncls = ((self.time_increment * current) / (F * self.n_ncls * self.ncls_volume)) * direction

        # update CLS and NCLS concentrations
        c_ox_cls = self.c_ox_cls - delta_cls
        c_red_cls = self.c_red_cls + delta_cls
        c_ox_ncls = self.c_ox_ncls + delta_ncls
        c_red_ncls = self.c_red_ncls - delta_ncls

        # for no crossover situation
        delta_ox = 0.0
        delta_red = 0.0

        # Coulomb counting from optional degradation and/or crossover mechanisms
        if cls_degradation is not None:
            c_ox_cls, c_red_cls = cls_degradation.degrade(c_ox_cls, c_red_cls, self.time_increment)

        if ncls_degradation is not None:
            c_ox_ncls, c_red_ncls = ncls_degradation.degrade(c_ox_ncls, c_red_ncls, self.time_increment)

        if cross_over is not None:
            (c_ox_cls, c_red_cls, c_ox_ncls, c_red_ncls, delta_ox,
             delta_red) = cross_over.crossover(self.geometric_area, c_ox_cls, c_red_cls, c_ox_ncls, c_red_ncls,
                                               self.cls_volume, self.ncls_volume, self.time_increment)
        # update concentrations to self
        self.prev_c_ox_cls = self.c_ox_cls
        self.prev_c_red_cls = self.c_red_cls
        self.prev_c_ox_ncls = self.c_ox_ncls
        self.prev_c_red_ncls = self.c_red_ncls

        self.c_ox_cls = c_ox_cls
        self.c_red_cls = c_red_cls
        self.c_ox_ncls = c_ox_ncls
        self.c_red_ncls = c_red_ncls

        self.delta_ox = delta_ox
        self.delta_red = delta_red

    def revert_concentrations(self) -> None:
        """Resets concentrations to previous value if a (invalid) negative concentration is calculated."""
        self.c_ox_cls = self.prev_c_ox_cls
        self.c_red_cls = self.prev_c_red_cls
        self.c_ox_ncls = self.prev_c_ox_ncls
        self.c_red_ncls = self.prev_c_red_ncls
