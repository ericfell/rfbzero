from abc import ABC, abstractmethod

import numpy as np
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
    capacity_only : bool
        True if only cycle times and capacities will be returned,
        False if all time series profiles will also be returned.

    """

    def __init__(self, size: int, capacity_only: bool):
        if not capacity_only:
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

    def state_of_charge(self):
        """Calculate state-of-charge in each reservoir"""
        for i, (cls_ox, cls_red, ncls_ox, ncls_red) in enumerate(zip(self.c_ox_cls_profile, self.c_red_cls_profile,
                                                                     self.c_ox_ncls_profile, self.c_red_ncls_profile)):

            if cls_ox + cls_red == 0.0 or ncls_ox + ncls_red == 0.0:
                break

            soc_cls = (cls_red / (cls_ox + cls_red)) * 100
            soc_ncls = (ncls_red / (ncls_ox + ncls_red)) * 100
            self.soc_profile_cls[i] = soc_cls
            self.soc_profile_ncls[i] = soc_ncls

    def structure_data(self, charge_first: bool):
        """Create separate charge/discharge cycle capacities and times from model outputs"""
        if charge_first:
            self.time_charge = self.cycle_time[::2]
            self.time_discharge = self.cycle_time[1::2]
            self.charge_capacity = self.cycle_capacity[::2]
            self.discharge_capacity = self.cycle_capacity[1::2]
        else:
            self.time_charge = self.cycle_time[1::2]
            self.time_discharge = self.cycle_time[::2]
            self.charge_capacity = self.cycle_capacity[1::2]
            self.discharge_capacity = self.cycle_capacity[::2]


class CyclingProtocol(ABC):
    """
    Base class to be overridden by specific cycling protocol choice.

    Parameters
    ----------
    current : float
        Instantaneous current flowing (A).
    charge_first : bool
        True if CLS charges first, False if CLS discharges first.
    capacity_only : bool
        True if only cycle times and capacities will be returned,
        False if all time series profiles will also be returned.

    """

    def __init__(self, current: float, charge_first: bool, capacity_only: bool = True):
        self.current = current
        self.charge = charge_first
        self.charge_first = charge_first
        self.capacity_only = capacity_only

        if self.current <= 0.0:
            raise ValueError("'current must be > 0.0")

        if not isinstance(self.charge, bool):
            raise ValueError("'charge_first' must be a boolean")

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

    @staticmethod
    def _validate_model(cell_model: ZeroDModel,
                        degradation: DegradationMechanism = None,
                        cls_degradation: DegradationMechanism = None,
                        ncls_degradation: DegradationMechanism = None):

        if degradation is not None and (cls_degradation is not None or ncls_degradation is not None):
            raise ValueError("Cannot specify both 'degradation' and '(n)cls_degradation'")

        # Raise ValueError if user inputs a negative concentration
        if cell_model.negative_concentrations():
            raise ValueError('Negative concentration detected')

    def current_direction(self) -> int:
        """Make current positive for charge, negative for discharge"""
        return 1 if self.charge else -1


class ConstantCurrent(CyclingProtocol):
    """
    Provides a constant current (CC) cycling method.

    Parameters
    ----------
    voltage_cutoff_charge : float
        Voltage above which cell will switch to discharge (V).
    voltage_cutoff_discharge : float
        Voltage below which cell will switch to charge (V).
    current : float
        Instantaneous current flowing (A).
    charge_first : bool
        True if CLS charges first, False if CLS discharges first.
    capacity_only : bool
        True if only cycle times and capacities will be returned,
        False if all time series profiles will also be returned.

    """

    def __init__(self, voltage_cutoff_charge: float, voltage_cutoff_discharge: float, current: float,
                 charge_first: bool = True, capacity_only: bool = True):
        self.voltage_cutoff_charge = voltage_cutoff_charge
        self.voltage_cutoff_discharge = voltage_cutoff_discharge
        super().__init__(current, charge_first, capacity_only)

    def run(self, duration: int, cell_model: ZeroDModel,
            degradation: DegradationMechanism = None,
            cls_degradation: DegradationMechanism = None,
            ncls_degradation: DegradationMechanism = None,
            cross_over: Crossover = None) -> CyclingProtocolResults:
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

        self._validate_model(cell_model, degradation, cls_degradation, ncls_degradation)

        if degradation is not None:
            cls_degradation = degradation
            ncls_degradation = degradation

        if not self.voltage_cutoff_discharge < cell_model.init_ocv < self.voltage_cutoff_charge:
            raise ValueError("Ensure that 'voltage_cutoff_discharge' < 'init_ocv' < 'voltage_cutoff_charge'")

        print(f"{duration} sec of cycling, time steps: {cell_model.time_increment} sec")

        # initialize data results object to be sent to user
        length_data = int(duration / cell_model.time_increment)
        """
        if int(duration // cell_model.time_increment) != length_data:
            print("WARNING: 'duration' not evenly dividable by 'time_increment',\
             \nexperiment duration will not be exactly as requested")
        """

        results = CyclingProtocolResults(length_data, self.capacity_only)

        count = 0
        capacity = 0.0

        i_lim_cls_t, i_lim_ncls_t = cell_model.limiting_concentration(self.charge)

        if self.current >= i_lim_cls_t or self.current >= i_lim_ncls_t:
            raise ValueError("Desired current > limiting current, cell can't run")

        # assign + current to charge, - current to discharge
        i = self.current_direction() * self.current

        losses, _, _ = cell_model.total_overpotential(i, self.charge, i_lim_cls_t, i_lim_ncls_t)
        ocv = cell_model.open_circuit_voltage()
        cell_v = cell_model.cell_voltage(ocv, losses, self.charge)

        if cell_v <= self.voltage_cutoff_discharge or cell_v >= self.voltage_cutoff_charge:
            raise ValueError("Desired current too high, overpotentials place cell voltage outside voltage cutoffs")

        while count != length_data:
            # set current
            i = self.current_direction() * self.current

            # record temporary values of concentrations for all species
            c_ox_cls_temp = cell_model.c_ox_cls
            c_red_cls_temp = cell_model.c_red_cls
            c_ox_ncls_temp = cell_model.c_ox_ncls
            c_red_ncls_temp = cell_model.c_red_ncls

            delta_ox, delta_red = cell_model.coulomb_counter(i, cls_degradation, ncls_degradation, cross_over)

            # EDGE CASE where voltage limits never reached i.e. straight CC cycling until concentration runs out
            if cell_model.negative_concentrations():
                # record capacity here
                results.cycle_capacity.append(capacity)
                results.cycle_time.append(count * cell_model.time_increment)

                #  Break out of loop if capacity approaches zero
                if capacity < 1.0 and len(results.cycle_capacity) > 2:
                    print(f"Simulation stopped after {count} time steps, due to capacity < 1 coulomb")
                    break

                capacity = 0.0

                # switch charge to discharge or vice-versa
                self.charge = not self.charge

                # set self back to previous, valid, concentration value
                cell_model.c_ox_cls = c_ox_cls_temp
                cell_model.c_red_cls = c_red_cls_temp
                cell_model.c_ox_ncls = c_ox_ncls_temp
                cell_model.c_red_ncls = c_red_ncls_temp

                # set limiting current for next cycle, with previous allowable concentrations
                i_lim_cls_t, i_lim_ncls_t = cell_model.limiting_concentration(self.charge)
                continue

            # calculate overpotentials and resulting cell voltage
            losses, n_act, n_mt = cell_model.total_overpotential(i, self.charge, i_lim_cls_t, i_lim_ncls_t)
            ocv = cell_model.open_circuit_voltage()
            cell_v = cell_model.cell_voltage(ocv, losses, self.charge)

            # did it hit voltage limit ?
            if cell_v >= self.voltage_cutoff_charge or cell_v <= self.voltage_cutoff_discharge:
                # record cycle capacity, cycle time
                results.cycle_capacity.append(capacity)
                results.cycle_time.append(count * cell_model.time_increment)

                if capacity < 1.0 and count > 2:
                    print(f"Simulation stopped after {count} time steps, due to capacity < 1 coulomb")
                    break

                capacity = 0.0

                # switch charge to discharge or vice-versa
                self.charge = not self.charge

                # set limiting current for next cycle
                i_lim_cls_t, i_lim_ncls_t = cell_model.limiting_concentration(self.charge)
                continue

            # update capacity
            capacity += abs(i * cell_model.time_increment)

            # update concentrations
            if not self.capacity_only:
                results.current_profile[count] = i
                results.c_ox_cls_profile[count] = cell_model.c_ox_cls
                results.c_red_cls_profile[count] = cell_model.c_red_cls
                results.c_ox_ncls_profile[count] = cell_model.c_ox_ncls
                results.c_red_ncls_profile[count] = cell_model.c_red_ncls

                results.del_ox[count] = delta_ox
                results.del_red[count] = delta_red

                results.cell_v_profile[count] = cell_v
                results.ocv_profile[count] = ocv
                results.act_profile[count] = n_act
                results.mt_profile[count] = n_mt
                results.loss_profile[count] = losses

                # record time of model data point
                results.times[count] = (count * cell_model.time_increment) + cell_model.time_increment

            count += 1

        # now calculating SOC of cls and ncls
        if not self.capacity_only:
            results.state_of_charge()
        # structures data into individual charge and discharge cycle times and capacities
        results.structure_data(self.charge_first)
        return results


class ConstantCurrentConstantVoltage(CyclingProtocol):
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
    capacity_only : bool
        True if only cycle times and capacities will be returned,
        False if all time series profiles will also be returned.

    """

    def __init__(self, voltage_limit_charge: float, voltage_limit_discharge: float, current_cutoff_charge: float,
                 current_cutoff_discharge: float, current: float, charge_first: bool = True, capacity_only: bool = True):
        self.voltage_limit_charge = voltage_limit_charge
        self.voltage_limit_discharge = voltage_limit_discharge
        self.current_cutoff_charge = current_cutoff_charge
        self.current_cutoff_discharge = current_cutoff_discharge
        super().__init__(current, charge_first, capacity_only)

        if self.current_cutoff_discharge >= 0.0 or self.current_cutoff_charge <= 0.0:
            raise ValueError("Ensure 'current_cutoff_discharge' < 0.0, 'current_cutoff_charge' > 0.0")

    def run(self, duration: int, cell_model: ZeroDModel,
            degradation: DegradationMechanism = None,
            cls_degradation: DegradationMechanism = None,
            ncls_degradation: DegradationMechanism = None,
            cross_over: Crossover = None) -> CyclingProtocolResults:
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

        self._validate_model(cell_model, degradation, cls_degradation, ncls_degradation)

        if degradation is not None:
            cls_degradation = degradation
            ncls_degradation = degradation

        if not self.voltage_limit_discharge < cell_model.init_ocv < self.voltage_limit_charge:
            raise ValueError("Ensure that 'voltage_limit_discharge' < 'init_ocv' < 'voltage_limit_charge'")

        print(f"{duration} sec of cycling, time steps: {cell_model.time_increment} sec")

        cc_mode = True  # tries to start in cc mode first
        i_first = True  # True if first current value of the cycle

        # initialize data object to be sent to user
        length_data = int(duration / cell_model.time_increment)

        results = CyclingProtocolResults(length_data, self.capacity_only)

        count = 0
        capacity = 0.0

        i_lim_cls_t, i_lim_ncls_t = cell_model.limiting_concentration(self.charge)
        i = self.current_direction() * self.current

        # current surpasses transport limitations, must switch to constant voltage
        if self.current >= i_lim_cls_t or self.current >= i_lim_ncls_t:
            cc_mode = False
            print("Current >= limiting current, forced to do CV cycling")
        else:
            losses, n_act, n_mt = cell_model.total_overpotential(i, self.charge, i_lim_cls_t, i_lim_ncls_t)
            ocv = cell_model.open_circuit_voltage()
            cell_v = cell_model.cell_voltage(ocv, losses, self.charge)

            if cell_v >= self.voltage_limit_charge or cell_v <= self.voltage_limit_discharge:
                cc_mode = False
                print("Overpotential would put cell outside voltage window, forced to do CV cycling")

        while count != length_data:
            # record temporary values of concentrations for all species
            c_ox_cls_temp = cell_model.c_ox_cls
            c_red_cls_temp = cell_model.c_red_cls
            c_ox_ncls_temp = cell_model.c_ox_ncls
            c_red_ncls_temp = cell_model.c_red_ncls

            if self.current >= i_lim_cls_t or self.current >= i_lim_ncls_t:
                print('trouble')
                cc_mode = False

            # check if in CC or CV mode
            if cc_mode:
                # set current
                i = self.current_direction() * self.current

                # calculate species' concentrations
                delta_ox, delta_red = cell_model.coulomb_counter(i, cls_degradation, ncls_degradation, cross_over)

                # EDGE CASE where voltage limits never reached i.e. straight CC cycling
                if cell_model.negative_concentrations():
                    # record cycle capacity, cycle time
                    results.cycle_capacity.append(capacity)
                    results.cycle_time.append(count * cell_model.time_increment)

                    #  Break out of loop if capacity near zero
                    if capacity < 1.0 and len(results.cycle_capacity) > 2:
                        print(f"Simulation stopped after {count} time steps, due to capacity < 1 coulomb")
                        break

                    capacity = 0.0

                    # switch charge to discharge or vice-versa
                    self.charge = not self.charge

                    # set self back to previous, valid, concentration value
                    cell_model.c_ox_cls = c_ox_cls_temp
                    cell_model.c_red_cls = c_red_cls_temp
                    cell_model.c_ox_ncls = c_ox_ncls_temp
                    cell_model.c_red_ncls = c_red_ncls_temp

                    # set limiting current for next cycle
                    i_lim_cls_t, i_lim_ncls_t = cell_model.limiting_concentration(self.charge)
                    continue

                # calculate overpotentials and resulting cell voltage
                losses, n_act, n_mt = cell_model.total_overpotential(i, self.charge, i_lim_cls_t, i_lim_ncls_t)
                ocv = cell_model.open_circuit_voltage()
                cell_v = cell_model.cell_voltage(ocv, losses, self.charge)

                # update capacity
                capacity += abs(i * cell_model.time_increment)

                # check if V limit is reached?
                if cell_v >= self.voltage_limit_charge or cell_v <= self.voltage_limit_discharge:
                    cc_mode = not cc_mode
                    # at some point maybe record capacity here too, so you know capacity due to CC and due to CV?
                    continue

                # update concentrations
                if not self.capacity_only:
                    results.current_profile[count] = i  # CC section
                    results.c_ox_cls_profile[count] = cell_model.c_ox_cls
                    results.c_red_cls_profile[count] = cell_model.c_red_cls
                    results.c_ox_ncls_profile[count] = cell_model.c_ox_ncls
                    results.c_red_ncls_profile[count] = cell_model.c_red_ncls

                    results.del_ox[count] = delta_ox
                    results.del_red[count] = delta_red

                    results.cell_v_profile[count] = cell_v
                    results.ocv_profile[count] = ocv
                    results.act_profile[count] = n_act
                    results.mt_profile[count] = n_mt
                    results.loss_profile[count] = losses

                    # times
                    results.times[count] = (count * cell_model.time_increment) + cell_model.time_increment

                count += 1

            else:  # now we're in CV mode

                # record temporary values of concentrations for all species
                c_ox_cls_temp = cell_model.c_ox_cls
                c_red_cls_temp = cell_model.c_red_cls
                c_ox_ncls_temp = cell_model.c_ox_ncls
                c_red_ncls_temp = cell_model.c_red_ncls

                # all cycling is now constant voltage
                cell_v = self.voltage_limit_charge if self.charge else self.voltage_limit_discharge
                ocv = cell_model.open_circuit_voltage()

                # adapting the solver's guess to the updated current
                # calculate the first current value for guess
                if i_first:
                    i_guess = self.current_direction() * self.current  # this can be too big but solver can handle
                    i_first = False
                else:
                    i_guess = i_cv

                    if abs(i_cv) >= i_lim_cls_t or abs(i_cv) >= i_lim_ncls_t:
                        print('i_cv > i_lim, cell stopped')
                        break

                i_cv = self._get_min_current(cell_model, i_guess, cell_v, ocv, self.charge, i_lim_cls_t, i_lim_ncls_t)

                # if current is below cutoffs, record cycle data and switch to charge/discharge
                if (self.charge and i_cv <= self.current_cutoff_charge) or \
                        (not self.charge and i_cv >= self.current_cutoff_discharge):

                    results.cycle_capacity.append(capacity)
                    results.cycle_time.append(count * cell_model.time_increment)

                    # Break out of full simulation if capacity nears zero
                    if capacity < 1.0 and len(results.cycle_capacity) > 2:
                        print(f"Simulation stopped after {count} time steps, due to capacity < 1 coulomb")
                        break

                    capacity = 0.0

                    # switch charge to discharge or vice-versa
                    self.charge = not self.charge
                    cc_mode = True
                    i_first = True

                    # set self back to previous, valid, concentration value
                    cell_model.c_ox_cls = c_ox_cls_temp
                    cell_model.c_red_cls = c_red_cls_temp
                    cell_model.c_ox_ncls = c_ox_ncls_temp
                    cell_model.c_red_ncls = c_red_ncls_temp

                    # set limiting current for next cycle
                    i_lim_cls_t, i_lim_ncls_t = cell_model.limiting_concentration(self.charge)
                    continue

                delta_ox, delta_red = cell_model.coulomb_counter(i_cv, cls_degradation, ncls_degradation, cross_over)

                # check if any reactant remains
                if cell_model.negative_concentrations():
                    results.cycle_capacity.append(capacity)
                    results.cycle_time.append(count * cell_model.time_increment)

                    # Break out of loop if capacity nears zero
                    if capacity < 1.0 and len(results.cycle_capacity) > 2:
                        print(f"Simulation stopped after {count} time steps, due to capacity < 1 coulomb")
                        break

                    capacity = 0.0

                    # switch charge to discharge or vice-versa
                    self.charge = not self.charge
                    cc_mode = True
                    i_first = True

                    # set self back to previous, valid, concentration value
                    cell_model.c_ox_cls = c_ox_cls_temp
                    cell_model.c_red_cls = c_red_cls_temp
                    cell_model.c_ox_ncls = c_ox_ncls_temp
                    cell_model.c_red_ncls = c_red_ncls_temp

                    # set limiting current for next cycle
                    i_lim_cls_t, i_lim_ncls_t = cell_model.limiting_concentration(self.charge)
                    continue

                # update capacity
                capacity += abs(i_cv * cell_model.time_increment)

                # update concentrations
                if not self.capacity_only:
                    results.current_profile[count] = i_cv
                    results.c_ox_cls_profile[count] = cell_model.c_ox_cls
                    results.c_red_cls_profile[count] = cell_model.c_red_cls
                    results.c_ox_ncls_profile[count] = cell_model.c_ox_ncls
                    results.c_red_ncls_profile[count] = cell_model.c_red_ncls

                    results.del_ox[count] = delta_ox
                    results.del_red[count] = delta_red

                    results.cell_v_profile[count] = cell_v
                    results.ocv_profile[count] = ocv
                    # these are undefined in CV mode
                    results.act_profile[count] = 0.0
                    results.mt_profile[count] = 0.0
                    results.loss_profile[count] = 0.0

                    # times
                    results.times[count] = (count * cell_model.time_increment) + cell_model.time_increment

                count += 1

        # now calculating SOC of cls and ncls
        if not self.capacity_only:
            results.state_of_charge()
        # structures data into individual charge and discharge cycle times and capacities
        results.structure_data(self.charge_first)
        return results

    @staticmethod
    def _get_min_current(cell_model: ZeroDModel, i_guess: float, cell_v: float, ocv: float, charge: bool,
                         i_lim_cls: float, i_lim_ncls: float) -> float:
        """
        Method wrapper to solve for current during constant voltage cycling.
        Attempts to minimize the difference of voltage, OCV, and losses (function of current).

        Parameters
        ----------
        cell_model : ZeroDModel
            Defined cell parameters for simulating.
        i_guess : float
            Initial guess for root of solver. Current at constant voltage (A).
        cell_v : float
            Cell voltage (V).
        ocv : float
            Cell open circuit voltage (V).
        charge : bool
            Positive if charging, negative if discharging.
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
            loss_solve, *_ = cell_model.total_overpotential(current, charge, i_lim_cls, i_lim_ncls)
            return cell_v - ocv - loss_solve if charge else cell_v - ocv + loss_solve

        min_current, *_ = fsolve(solver, np.array([i_guess]), xtol=1e-5)
        return min_current
