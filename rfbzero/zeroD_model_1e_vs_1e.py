

from math import log
import scipy.constants as spc
from zeroD_model_degradations import DegradationMechanism
from zeroD_model_crossover import Crossover


# Faraday constant (C/mol)
F = spc.value('Faraday constant')

# Molar gas constant (J/K/mol)
R = spc.R

# make these parameters at some point?
TEMPERATURE = 298  # Kelvins, for S.T.P.
NERNST_CONST = (R * TEMPERATURE) / F  # should have n_electrons input option


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
    init_ocv : float
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
    alpha_ncls : float
        Charge transfer coefficient of NCLS redox couple, dimensionless.
    geometric_area : float
        Geometric area of cell (cm^2).
    cls_negolyte : bool
        If True, negolyte is the CLS.
    time_increment : float
        Simulation time step (s).
    k_mt : float
        Mass transport coefficient (cm/s).
    roughness_factor : float
        Roughness factor, dimensionless.
        Surface area divided by geometric surface area.


    Notes
    -----
    All equations are adapted from [1]. If ZeroDModel has been
    significant to your research please cite the paper.

    [1] Modak, S.; Kwabi, D. G. A Zero-Dimensional Model for Electrochemical
    Behavior and Capacity Retention in Organic Flow Cells, Journal of The
    Electrochemical Society, 168, 2021, 080528.
    """

    def __init__(self, cls_volume: float, ncls_volume: float, cls_start_c_ox: float, cls_start_c_red: float,
                 ncls_start_c_ox: float, ncls_start_c_red: float, init_ocv: float, resistance: float, k_0_cls: float,
                 k_0_ncls: float, alpha_cls: float = 0.5, alpha_ncls: float = 0.5, geometric_area: float = 5.0,
                 cls_negolyte: bool = True, time_increment: float = 0.01, k_mt: float = 0.8,
                 roughness_factor: float = 26.0) -> None:
        """Inits ZeroDModel"""
        self.cls_volume = cls_volume
        self.ncls_volume = ncls_volume
        self.c_ox_cls = cls_start_c_ox
        self.c_red_cls = cls_start_c_red
        self.c_ox_ncls = ncls_start_c_ox
        self.c_red_ncls = ncls_start_c_red
        self.init_ocv = init_ocv
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

    # option for just measuring capacity over time, doesn't need to make all arrays?

    @staticmethod
    def current_direction(charge: bool) -> int:
        """Make current positive for charge, negative for discharge"""
        return 1 if charge else -1

    def _exchange_current(self) -> tuple[float, float]:
        """
        Calculates exchange current (i_0) of redox couples in the CLS and NCLS.
        Value returned is in Amps.


        Returns
        -------
        i_0_cls : float
            Exchange current of CLS redox couple
            at a given timestep (A).
        i_0_ncls : float
            Exchange current of NCLS redox couple
            at a given timestep (A)

        """
        # division by 1000 for conversion from mol/L to mol/cm^3
        i_0_cls = (self.const_i_ex * self.k_0_cls * (self.c_red_cls ** self.alpha_cls)
                   * (self.c_ox_cls ** (1 - self.alpha_cls)) * 0.001)
        i_0_ncls = (self.const_i_ex * self.k_0_ncls * (self.c_red_ncls ** self.alpha_ncls)
                    * (self.c_ox_ncls ** (1 - self.alpha_ncls)) * 0.001)
        return i_0_cls, i_0_ncls

    def _limiting_current(self, c_lim: float) -> float:
        """Calculates limiting current (i_lim) for a single reservoir.
        Value returned is in Amps.
        This is equation 6 of [1].
        """
        # div by 1000 for conversion from mol/L to mol/cm^3
        # will require n electrons param
        return F * self.k_mt * c_lim * self.geometric_area * 0.001

    def limiting_concentration(self, charge: bool) -> tuple[float, float]:
        """Selects limiting concentration and calculates limiting current for CLS and NCLS."""
        if (self.cls_negolyte and charge) or (not self.cls_negolyte and not charge):
            i_lim_cls = self._limiting_current(self.c_ox_cls)
            i_lim_ncls = self._limiting_current(self.c_red_ncls)
        else:
            i_lim_cls = self._limiting_current(self.c_red_cls)
            i_lim_ncls = self._limiting_current(self.c_ox_ncls)

        return i_lim_cls, i_lim_ncls

    @staticmethod
    def _activation_overpotential(current: float, i_0_cls: float, i_0_ncls: float) -> float:
        """
        This is equation 4 of [1].
        Parameters
        ----------
        current
        i_0_cls
        i_0_ncls

        Returns
        -------

        """

        z_cls = abs(current) / (2 * i_0_cls)
        z_ncls = abs(current) / (2 * i_0_ncls)
        n_act = NERNST_CONST * (log(z_ncls + ((z_ncls**2) + 1)**0.5) + log(z_cls + ((z_cls**2) + 1)**0.5))
        return n_act

    def negative_concentrations(self) -> bool:
        """Return True if any concentration is negative"""
        return any(x < 0.0 for x in [self.c_ox_cls, self.c_red_cls, self.c_ox_ncls, self.c_red_ncls])

    def _mass_transport_overpotential(self, charge: bool, current: float, i_lim_cls: float, i_lim_ncls: float) -> float:
        """
        This is equation 8 of [1].

        Parameters
        ----------
        charge
        current
        i_lim_cls
        i_lim_ncls

        Returns
        -------

        """
        # Raise ValueError if a negative concentration is detected
        if self.negative_concentrations():
            raise ValueError('Negative concentration detected')

        c_tot_cls = self.c_red_cls + self.c_ox_cls
        c_tot_ncls = self.c_red_ncls + self.c_ox_ncls

        i = abs(current)

        if (self.cls_negolyte and charge) or (not self.cls_negolyte and not charge):
            n_mt = NERNST_CONST * log((1 - ((c_tot_cls * i) / ((self.c_red_cls * i_lim_cls) + (self.c_ox_cls * i))))
                                      * (1 - ((c_tot_ncls * i) / ((self.c_ox_ncls * i_lim_ncls)
                                                                  + (self.c_red_ncls * i)))))
        else:
            n_mt = NERNST_CONST * log(((1 - ((c_tot_cls * i) / ((self.c_ox_cls * i_lim_cls) + (self.c_red_cls * i))))
                                       * (1 - ((c_tot_ncls * i) / ((self.c_red_ncls * i_lim_ncls)
                                                                   + (self.c_ox_ncls * i))))))
        return n_mt

    def total_overpotential(self, current: float, charge: bool,
                            i_lim_cls: float, i_lim_ncls: float) -> tuple[float, float, float]:
        """
        This is the overpotentials of equation 2 in [1].

        Parameters
        ----------
        current
        charge
        i_lim_cls
        i_lim_ncls

        Returns
        -------

        """

        i_0_cls, i_0_ncls = self._exchange_current()
        # calculate ohmic, activation, mass transport overpotentials
        n_ohmic = abs(current)*self.resistance
        n_act = self._activation_overpotential(current, i_0_cls, i_0_ncls)
        n_mt = self._mass_transport_overpotential(charge, current, i_lim_cls, i_lim_ncls)

        n_loss = n_ohmic + n_act + n_mt

        return n_loss, n_act, n_mt

    def open_circuit_voltage(self) -> float:
        """
        Nernstian calculation of cell's open circuit voltage.
        This is equivalent to equation 3 of [1].

        Returns
        -------

        """

        # Raise ValueError if a negative concentration is detected
        if self.negative_concentrations():
            raise ValueError('Negative concentration detected')

        # will need n_electrons input
        # CLS is negolyte
        if self.cls_negolyte:
            ocv = (self.init_ocv
                   + (NERNST_CONST * log(self.c_red_cls / self.c_ox_cls))
                   + (NERNST_CONST * log(self.c_ox_ncls / self.c_red_ncls)))

        # CLS is posolyte
        else:
            ocv = (self.init_ocv
                   - (NERNST_CONST * log(self.c_red_cls / self.c_ox_cls))
                   - (NERNST_CONST * log(self.c_ox_ncls / self.c_red_ncls)))
        return ocv

    @staticmethod
    def cell_voltage(ocv: float, losses: float, charge: bool) -> float:
        """If charging, add overpotentials to OCV, else subtract them."""
        return ocv + losses if charge else ocv - losses

    def coulomb_counter(self, current: float,
                        cls_degradation: DegradationMechanism = None,
                        ncls_degradation: DegradationMechanism = None,
                        crossover_params: Crossover = None) -> tuple[float, float]:

        # Coulomb counting based solely on current
        direction = 1 if self.cls_negolyte else -1
        delta_cls = ((self.time_increment * current) / (F * self.cls_volume)) * direction
        delta_ncls = ((self.time_increment * current) / (F * self.ncls_volume)) * direction

        # update CLS and NCLS concentrations
        c_ox_cls = self.c_ox_cls - delta_cls
        c_red_cls = self.c_red_cls + delta_cls
        c_ox_ncls = self.c_ox_ncls + delta_ncls
        c_red_ncls = self.c_red_ncls - delta_ncls

        # for no crossover situation
        delta_ox = 0.0
        delta_red = 0.0

        # Coulomb counting from optional degradation/crossover mechanisms

        if cls_degradation is not None:
            c_ox_cls, c_red_cls = cls_degradation.degrade(c_ox_cls, c_red_cls, self.time_increment)

        if ncls_degradation is not None:
            c_ox_ncls, c_red_ncls = ncls_degradation.degrade(c_ox_ncls, c_red_ncls, self.time_increment)

        if crossover_params is not None:
            (c_ox_cls, c_red_cls, c_ox_ncls, c_red_ncls, delta_ox,
             delta_red) = crossover_params.crossover(c_ox_cls, c_red_cls, c_ox_ncls, c_red_ncls, self.time_increment,
                                                     self.cls_volume, self.ncls_volume)
        # update concentrationss to self
        self.c_ox_cls = c_ox_cls
        self.c_red_cls = c_red_cls
        self.c_ox_ncls = c_ox_ncls
        self.c_red_ncls = c_red_ncls

        return delta_ox, delta_red

    @staticmethod
    def state_of_charge(c_ox_cls, c_red_cls, c_ox_ncls, c_red_ncls) -> tuple[float, float]:
        """Calculate state-of-charge in each reservoir"""
        soc_cls = (c_red_cls / (c_ox_cls + c_red_cls)) * 100
        soc_ncls = (c_red_ncls / (c_ox_ncls + c_red_ncls)) * 100
        return soc_cls, soc_ncls

    # is below proper *args unpacking naming style?
    def cv_current_solver(self, current: float, *data: float) -> float:
        (cell_v, ocv, charge, i_lim_cls, i_lim_ncls) = data
        # curr has sign but total_overpotential makes it always positive
        loss_solve, _, _ = self.total_overpotential(current, charge, i_lim_cls,
                                                    i_lim_ncls)
        # returns what solver will try to minimize
        return cell_v - ocv - loss_solve if charge else cell_v - ocv + loss_solve


if __name__ == '__main__':
    print('testing')           
        