from abc import ABC, abstractmethod

from scipy.optimize import fsolve

from redox_flow_cell import ZeroDModel
from degradation import DegradationMechanism
from crossover import Crossover


class CyclingProtocolResults:
    """
    Returns a data object of the simulation results.

    Parameters
    ----------
    size : int
        Total number of time steps in desired simulation

    """

    def __init__(self, size: int, charge_first: bool):
        self.size = size
        self.charge_first = charge_first

        self.step = 0
        self.cycles = 0

        self.current_profile = [0.0] * size
        self.c_ox_cls_profile = [0.0] * size
        self.c_red_cls_profile = [0.0] * size
        self.c_ox_ncls_profile = [0.0] * size
        self.c_red_ncls_profile = [0.0] * size
        self.cell_v_profile = [0.0] * size
        self.soc_profile_cls = [0.0] * size
        self.soc_profile_ncls = [0.0] * size
        self.ocv_profile = [0.0] * size
        self.times = [0.0] * size
        self.act_profile = [0.0] * size
        self.mt_profile = [0.0] * size
        self.loss_profile = [0.0] * size
        self.del_ox = [0.0] * size
        self.del_red = [0.0] * size

        # total cycles is unknown at start, thus size is undetermined
        self.cycle_capacity = []
        self.cycle_time = []

        self.time_charge = []
        self.time_discharge = []
        self.charge_capacity = []
        self.discharge_capacity = []

    def finalize(self):
        self._state_of_charge()
        self._structure_data()

    def _state_of_charge(self):
        """Calculate state-of-charge in each reservoir"""
        for i, (cls_ox, cls_red, ncls_ox, ncls_red) in enumerate(zip(
                self.c_ox_cls_profile, self.c_red_cls_profile, self.c_ox_ncls_profile, self.c_red_ncls_profile
        )):

            if cls_ox + cls_red == 0.0 or ncls_ox + ncls_red == 0.0:
                break

            soc_cls = (cls_red / (cls_ox + cls_red)) * 100
            soc_ncls = (ncls_red / (ncls_ox + ncls_red)) * 100
            self.soc_profile_cls[i] = soc_cls
            self.soc_profile_ncls[i] = soc_ncls

    def _structure_data(self):
        """Create separate charge/discharge cycle capacities and times from model outputs"""
        charge_index = int(not self.charge_first)
        self.time_charge = self.cycle_time[charge_index::2]
        self.charge_capacity = self.cycle_capacity[charge_index::2]

        discharge_index = int(self.charge_first)
        self.time_discharge = self.cycle_time[discharge_index::2]
        self.discharge_capacity = self.cycle_capacity[discharge_index::2]


class CyclingProtocol(ABC):
    """
    Base class to be overridden by specific cycling protocol choice.

    Parameters
    ----------
    current : float
        Instantaneous current flowing (A).
    charge_first : bool
        True if CLS charges first, False if CLS discharges first.

    """

    def __init__(self, current: float, charge_first: bool):
        self.current = current
        self.charge = charge_first
        self.charge_first = charge_first

        if self.current <= 0.0:
            raise ValueError("'current' must be > 0.0")

        if not isinstance(self.charge, bool):
            raise ValueError("'charge_first' must be a boolean")

        self.cell_model = None
        self.cls_degradation = None
        self.ncls_degradation = None
        self.cross_over = None
        self.results = None

    @abstractmethod
    def run(
            self,
            duration: int,
            cell_model: ZeroDModel,
            degradation: DegradationMechanism = None,
            cls_degradation: DegradationMechanism = None,
            ncls_degradation: DegradationMechanism = None,
            cross_over: Crossover = None,
    ) -> CyclingProtocolResults:
        """Applies a cycling protocol and (optional) degradation mechanisms to a cell model"""
        raise NotImplementedError

    def _init_protocol(
            self,
            duration: int,
            cell_model: ZeroDModel,
            degradation: DegradationMechanism,
            cls_degradation: DegradationMechanism,
            ncls_degradation: DegradationMechanism,
            cross_over: Crossover
    ) -> None:
        if degradation is not None and (cls_degradation is not None or ncls_degradation is not None):
            raise ValueError("Cannot specify both 'degradation' and '(n)cls_degradation'")

        if cell_model.negative_concentrations():
            raise ValueError('Negative concentration detected')

        if degradation is not None:
            cls_degradation = degradation
            ncls_degradation = degradation

        self.cell_model = cell_model
        self.cls_degradation = cls_degradation
        self.ncls_degradation = ncls_degradation
        self.cross_over = cross_over

        # initialize data results object to be sent to user
        self.results = CyclingProtocolResults(int(duration / cell_model.time_increment), self.charge_first)

    def current_direction(self) -> int:
        """Make current positive for charge, negative for discharge"""
        return 1 if self.charge else -1


class ConstantCurrent(CyclingProtocol):
    """
    Provides a constant current (CC) cycling method.

    Parameters
    ----------
    voltage_limit_charge : float
        Voltage above which cell will switch to discharge (V).
    voltage_limit_discharge : float
        Voltage below which cell will switch to charge (V).
    current : float
        Instantaneous current flowing (A).
    charge_first : bool
        True if CLS charges first, False if CLS discharges first.

    """

    def __init__(
            self,
            voltage_limit_charge: float,
            voltage_limit_discharge: float,
            current: float,
            charge_first: bool = True,
    ):
        super().__init__(current, charge_first)
        self.voltage_limit_charge = voltage_limit_charge
        self.voltage_limit_discharge = voltage_limit_discharge

    def run(
            self,
            duration: int,
            cell_model: ZeroDModel,
            degradation: DegradationMechanism = None,
            cls_degradation: DegradationMechanism = None,
            ncls_degradation: DegradationMechanism = None,
            cross_over: Crossover = None
    ) -> CyclingProtocolResults:
        """

        Parameters
        ----------
        duration : int
            Simulation time (s).
        cell_model : ZeroDModel
            Defined cell parameters for simulating.
        degradation : DegradationMechanism, optional
            Degradation mechanism(s) applied to CLS and NCLS.
        cls_degradation : DegradationMechanism, optional
            Degradation mechanism(s) applied to CLS.
        ncls_degradation : DegradationMechanism, optional
            Degradation mechanism(s) applied to NCLS.
        cross_over : Crossover, optional
            Crossover mechanism applied to cell.

        Returns
        -------
        results : CyclingProtocolResults
            Object with results from simulation

        """

        self._init_protocol(duration, cell_model, degradation, cls_degradation, ncls_degradation, cross_over)

        if not self.voltage_limit_discharge < cell_model.init_ocv < self.voltage_limit_charge:
            raise ValueError("Ensure that 'voltage_limit_discharge' < 'init_ocv' < 'voltage_limit_charge'")

        i_lim_cls_t, i_lim_ncls_t = cell_model.limiting_concentration(self.charge)

        if not self._is_current_in_limits(i_lim_cls_t, i_lim_ncls_t):
            raise ValueError("Desired current > limiting current, cell can't run")

        if not self._is_voltage_in_limits(i_lim_cls_t, i_lim_ncls_t):
            raise ValueError("Desired current too high, overpotentials place cell voltage outside voltage limits")

        print(f"{duration} sec of cycling, time steps: {cell_model.time_increment} sec")

        capacity = 0.0
        end_simulation = False

        while not end_simulation and self.results.step < self.results.size:
            capacity, i_lim_cls_t, i_lim_ncls_t, end_simulation = self._cc_mode_cycle(capacity, i_lim_cls_t, i_lim_ncls_t)

        self.results.finalize()
        return self.results

    def _is_current_in_limits(self, i_lim_cls_t: float, i_lim_ncls_t: float) -> bool:
        return self.current < min(i_lim_cls_t, i_lim_ncls_t)

    def _is_voltage_in_limits(self, i_lim_cls_t, i_lim_ncls_t):
        i = self.current_direction() * self.current
        losses, *_ = self.cell_model.total_overpotential(i, self.charge, i_lim_cls_t, i_lim_ncls_t)
        ocv = self.cell_model.open_circuit_voltage()
        cell_v = self.cell_model.cell_voltage(ocv, losses, self.charge)
        return self.voltage_limit_discharge < cell_v < self.voltage_limit_charge

    def _cc_mode_cycle(
            self,
            capacity: float,
            i_lim_cls_t: float,
            i_lim_ncls_t: float
    ) -> tuple[float, float, float, bool]:
        end_simulation = False

        # set current
        i = self.current_direction() * self.current

        # calculate species' concentrations
        self.cell_model.coulomb_counter(i, self.cls_degradation, self.ncls_degradation, self.cross_over)

        # edge case where the voltage limits are never reached, i.e. straight CC cycling
        if self.cell_model.negative_concentrations():
            end_simulation = self._end_cycle(capacity)
            if not end_simulation:
                # set self back to previous, valid, concentration value
                self.cell_model.revert_concentrations()
                capacity = 0.0
                i_lim_cls_t, i_lim_ncls_t = self.cell_model.limiting_concentration(self.charge)
            return capacity, i_lim_cls_t, i_lim_ncls_t, end_simulation

        # calculate overpotentials and resulting cell voltage
        losses, n_act, n_mt = self.cell_model.total_overpotential(i, self.charge, i_lim_cls_t, i_lim_ncls_t)
        ocv = self.cell_model.open_circuit_voltage()
        cell_v = self.cell_model.cell_voltage(ocv, losses, self.charge)

        # check if V limit is reached?
        if not self.voltage_limit_discharge < cell_v < self.voltage_limit_charge:
            capacity, i_lim_cls_t, i_lim_ncls_t, end_simulation = self._handle_cc_mode_voltage_limit_reached(capacity, i_lim_cls_t, i_lim_ncls_t)
            if end_simulation:
                return capacity, i_lim_cls_t, i_lim_ncls_t, end_simulation

        # update capacity
        capacity += abs(i * self.cell_model.time_increment)

        # update concentrations
        self._record_common_results(i, cell_v, ocv)
        self.results.act_profile[self.results.step] = n_act
        self.results.mt_profile[self.results.step] = n_mt
        self.results.loss_profile[self.results.step] = losses

        self.results.step += 1

        return capacity, i_lim_cls_t, i_lim_ncls_t, end_simulation

    # Returns true if the simulation should end
    def _end_cycle(self, capacity: float) -> bool:
        self.results.cycle_capacity.append(capacity)
        self.results.cycle_time.append(self.results.step * self.cell_model.time_increment)
        self.results.cycles += 1

        # end the simulation if capacity nears zero
        if capacity < 1.0 and self.results.cycles > 2:
            print(f"Simulation stopped after {self.results.step} time steps, due to capacity < 1 coulomb")
            return True

        # switch charge to discharge or vice-versa
        self.charge = not self.charge

        return False

    def _handle_cc_mode_voltage_limit_reached(self, capacity, i_lim_cls_t, i_lim_ncls_t):
        # at some point maybe record capacity here too, so you know capacity due to CC and due to CV?
        end_simulation = self._end_cycle(capacity)
        if not end_simulation:
            capacity = 0.0
            i_lim_cls_t, i_lim_ncls_t = self.cell_model.limiting_concentration(self.charge)

        return capacity, i_lim_cls_t, i_lim_ncls_t, end_simulation

    def _record_common_results(self, i: float, cell_v: float, ocv: float) -> None:
        step = self.results.step

        self.results.current_profile[step] = i
        self.results.cell_v_profile[step] = cell_v
        self.results.ocv_profile[step] = ocv

        self.results.c_ox_cls_profile[step] = self.cell_model.c_ox_cls
        self.results.c_red_cls_profile[step] = self.cell_model.c_red_cls
        self.results.c_ox_ncls_profile[step] = self.cell_model.c_ox_ncls
        self.results.c_red_ncls_profile[step] = self.cell_model.c_red_ncls

        self.results.del_ox[step] = self.cell_model.delta_ox
        self.results.del_red[step] = self.cell_model.delta_red

        self.results.times[step] = self.cell_model.time_increment * (step + 1)


class ConstantCurrentConstantVoltage(ConstantCurrent):
    """
    Provides a constant current constant voltage (CCCV) cycling method which, in the limit of a high current demanded of
    a cell that it cannot maintain, becomes a constant voltage (CV) cycling method.

    Parameters
    ----------
    voltage_limit_charge : float
        Voltage above which cell will switch to CV mode, charging (V).
    voltage_limit_discharge : float
        Voltage below which cell will switch to CV mode, discharging (V).
    current_cutoff_charge : float
        Current below which CV charging will switch to CC portion of CCCV discharge (A).
    current_cutoff_discharge : float
        Current above which CV discharging will switch to CC portion of CCCV charge (A).
    current : float
        Instantaneous current flowing (A).
    charge_first : bool
        True if CLS charges first, False if CLS discharges first.

    """

    def __init__(
            self,
            voltage_limit_charge: float,
            voltage_limit_discharge: float,
            current_cutoff_charge: float,
            current_cutoff_discharge: float,
            current: float,
            charge_first: bool = True,
    ):
        super().__init__(voltage_limit_charge, voltage_limit_discharge, current, charge_first)
        self.current_cutoff_charge = current_cutoff_charge
        self.current_cutoff_discharge = current_cutoff_discharge
        self.cc_mode = True

        if self.current_cutoff_discharge >= 0.0 or self.current_cutoff_charge <= 0.0:
            raise ValueError("Ensure 'current_cutoff_discharge' < 0.0, 'current_cutoff_charge' > 0.0")

    def run(
            self,
            duration: int,
            cell_model: ZeroDModel,
            degradation: DegradationMechanism = None,
            cls_degradation: DegradationMechanism = None,
            ncls_degradation: DegradationMechanism = None,
            cross_over: Crossover = None
    ) -> CyclingProtocolResults:
        """

        Parameters
        ----------
        duration : int
            Simulation time (s).
        cell_model : ZeroDModel
            Defined cell parameters for simulating.
        degradation : DegradationMechanism, optional
            Degradation mechanism(s) applied to CLS and NCLS.
        cls_degradation : DegradationMechanism, optional
            Degradation mechanism(s) applied to CLS.
        ncls_degradation : DegradationMechanism, optional
            Degradation mechanism(s) applied to NCLS.
        cross_over : Crossover, optional
            Crossover mechanism applied to cell.

        Returns
        -------
        results : CyclingProtocolResults
            Object with results from simulation

        """
        self._init_protocol(duration, cell_model, degradation, cls_degradation, ncls_degradation, cross_over)

        print(f"{duration} sec of cycling, time steps: {cell_model.time_increment} sec")

        i_lim_cls_t, i_lim_ncls_t = cell_model.limiting_concentration(self.charge)

        # check if cell needs to go straight to CV
        self.cc_mode = True  # tries to start in cc mode first
        if not self._is_current_in_limits(i_lim_cls_t, i_lim_ncls_t):
            print("Skip to CV cycling due to high current")
            self.cc_mode = False
        elif not self._is_voltage_in_limits(i_lim_cls_t, i_lim_ncls_t):
            print("Skip to CV cycling due to high overpotential")
            self.cc_mode = False

        capacity = 0.0
        i_cv = 0.0 # This will be set to an initial guess of the current before being used
        end_simulation = False

        while not end_simulation and self.results.step < self.results.size:
            if self.current >= min(i_lim_cls_t, i_lim_ncls_t):
                print('trouble')  # todo
                self.cc_mode = False

            if self.cc_mode:
                capacity, i_lim_cls_t, i_lim_ncls_t, end_simulation = self._cc_mode_cycle(capacity, i_lim_cls_t, i_lim_ncls_t)
            else:
                capacity, i_cv, i_lim_cls_t, i_lim_ncls_t, end_simulation = self._cv_mode_cycle(capacity, i_cv, i_lim_cls_t, i_lim_ncls_t)

        self.results.finalize()
        return self.results

    def _handle_cc_mode_voltage_limit_reached(self, capacity, i_lim_cls_t, i_lim_ncls_t):
        self.cc_mode = False
        return capacity, i_lim_cls_t, i_lim_ncls_t, False

    # all cycling is now constant voltage
    def _cv_mode_cycle(
            self,
            capacity: float,
            i_cv: float,
            i_lim_cls_t: float,
            i_lim_ncls_t: float,
    ) -> tuple[float, float, float, float, bool]:
        end_simulation = False

        if i_cv and abs(i_cv) >= min(i_lim_cls_t, i_lim_ncls_t):
            print('i_cv > i_lim, cell stopped')
            end_simulation = True
            return capacity, i_cv, i_lim_cls_t, i_lim_ncls_t, end_simulation

        # calculate the first current value for guess
        if not i_cv:
            i_cv = self.current_direction() * self.current  # this can be too big but solver can handle

        cell_v = self.voltage_limit_charge if self.charge else self.voltage_limit_discharge
        ocv = self.cell_model.open_circuit_voltage()

        # adapting the solver's guess to the updated current
        i_cv = self._get_min_current(i_cv, cell_v, ocv, i_lim_cls_t, i_lim_ncls_t)

        # if current is below cutoffs, record cycle data and switch to charge/discharge
        if (self.charge and i_cv <= self.current_cutoff_charge) or \
                (not self.charge and i_cv >= self.current_cutoff_discharge):
            end_simulation = self._end_cycle(capacity)
            if not end_simulation:
                # set self back to previous, valid, concentration value
                self.cell_model.revert_concentrations()
                capacity = 0.0
                i_lim_cls_t, i_lim_ncls_t = self.cell_model.limiting_concentration(self.charge)
                self.cc_mode = True

            return capacity, i_cv, i_lim_cls_t, i_lim_ncls_t, end_simulation

        self.cell_model.coulomb_counter(i_cv, self.cls_degradation, self.ncls_degradation, self.cross_over)

        # check if any reactant remains
        if self.cell_model.negative_concentrations():
            end_simulation = self._end_cycle(capacity)
            if not end_simulation:
                # TODO: duplicate code of current_cutoff if statement, try to refactor this
                # set self back to previous, valid, concentration value
                self.cell_model.revert_concentrations()
                capacity = 0.0
                i_lim_cls_t, i_lim_ncls_t = self.cell_model.limiting_concentration(self.charge)
                self.cc_mode = True

            return capacity, i_cv, i_lim_cls_t, i_lim_ncls_t, end_simulation

        # update capacity
        capacity += abs(i_cv * self.cell_model.time_increment)

        # update concentrations
        self._record_common_results(i_cv, cell_v, ocv)
        self.results.step += 1

        return capacity, i_cv, i_lim_cls_t, i_lim_ncls_t, end_simulation

    def _get_min_current(self, i_guess: float, cell_v: float, ocv: float, i_lim_cls: float, i_lim_ncls: float) -> float:
        """
        Method wrapper to solve for current during constant voltage cycling.
        Attempts to minimize the difference of voltage, OCV, and losses (function of current).

        Parameters
        ----------
        i_guess : float
            Initial guess for root of solver. Current at constant voltage (A).
        cell_v : float
            Cell voltage (V).
        ocv : float
            Cell open circuit voltage (V).
        i_lim_cls : float
            Limiting current of CLS redox couple
            at a given timestep (A).
        i_lim_ncls : float
            Limiting current of NCLS redox couple
            at a given timestep (A).

        Returns
        -------
        min_current : float
            Solved current at given timestep of constant voltage cycling (A).
        """
        def solver(current: float) -> float:
            """Numerical solver for current during constant voltage cycling"""
            loss_solve, *_ = self.cell_model.total_overpotential(current, self.charge, i_lim_cls, i_lim_ncls)
            return cell_v - ocv - self.current_direction() * loss_solve

        min_current, *_ = fsolve(solver, i_guess, xtol=1e-5)
        return min_current
