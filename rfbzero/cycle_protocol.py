

from abc import ABC, abstractmethod
from scipy.optimize import fsolve
import numpy as np
from zeroD_model_1e_vs_1e import ZeroDModel
from zeroD_model_degradations import DegradationMechanism
from zeroD_model_crossover import Crossover


class CyclingProtocol(ABC):
    @abstractmethod
    def run(self, duration: int, cell_model: ZeroDModel,
            degradation: DegradationMechanism = None,
            cls_degradation: DegradationMechanism = None,
            ncls_degradation: DegradationMechanism = None,
            crossover_params: Crossover = None,
            ) -> object:
        raise NotImplementedError


class ConstantCurrent(CyclingProtocol):
    """
    Subclass for constant current (CC) cycling method.

    Parameters
    ----------
    voltage_cutoff_charge : float
        Voltage above which cell will switch to discharge (V).
    voltage_cutoff_discharge : float
        Voltage below which cell will switch to charge (V).
    current : float
        Instantaneous current flowing (A).
    charge_first : bool
        If True, charges the cell first.

    """

    def __init__(self, voltage_cutoff_charge: float, voltage_cutoff_discharge: float,
                 current: float, charge_first: bool = True):
        self.voltage_cutoff_charge = voltage_cutoff_charge
        self.voltage_cutoff_discharge = voltage_cutoff_discharge
        self.current = current
        self.charge = charge_first

        if self.current < 0.0:
            raise ValueError("Current must be a positive value")

    def run(self, duration: int, cell_model: ZeroDModel,
            degradation: DegradationMechanism = None,
            cls_degradation: DegradationMechanism = None,
            ncls_degradation: DegradationMechanism = None,
            crossover_params: Crossover = None) -> object:
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
        crossover_params : Crossover, optional
            Crossover mechanism applied to cell.

        Returns
        -------

        """

        if degradation is not None and (cls_degradation is not None or ncls_degradation is not None):
            raise ValueError("Cannot specify both 'degradation' and '(n)cls_degradation'")

        if degradation is not None:
            cls_degradation = degradation
            ncls_degradation = degradation

        # Raise ValueError if user inputs a negative concentration
        if cell_model.negative_concentrations():
            raise ValueError('Negative concentration detected')

        (current_profile, c_ox_cls_profile, c_red_cls_profile, cell_v_profile, soc_profile_cls, ocv_profile,
         c_ox_ncls_profile, c_red_ncls_profile, soc_profile_ncls, cycle_capacity, cycle_time, act_profile,
         mt_profile, loss_profile) = [], [], [], [], [], [], [], [], [], [], [], [], [], []
        # testing
        del_ox, del_red = [], []

        # record temporary values of concentrations for all species
        c_ox_cls_temp = cell_model.c_ox_cls
        c_red_cls_temp = cell_model.c_red_cls
        c_ox_ncls_temp = cell_model.c_ox_ncls
        c_red_ncls_temp = cell_model.c_red_ncls

        length_data = int(duration / cell_model.time_increment)
        times = [x*cell_model.time_increment for x in range(1, length_data + 1)]
        print(f"{duration} sec of cycling, timesteps: {cell_model.time_increment} sec")

        count = 0
        cap = 0.0
        cap_low = False
        # initialized in case simulation has to stop due to no more available capacity
        final_count = length_data

        i_lim_cls_t, i_lim_ncls_t = cell_model.limiting_concentration(self.charge)

        if self.current >= i_lim_cls_t or self.current >= i_lim_ncls_t:
            raise ValueError("Desired current > limiting current, cell can't run")

        # assign + current to charge, - current to discharge
        i = cell_model.current_direction(self.charge) * self.current

        losses, _, _ = cell_model.total_overpotential(i, self.charge, i_lim_cls_t, i_lim_ncls_t)
        ocv = cell_model.open_circuit_voltage()
        cell_v = cell_model.cell_voltage(ocv, losses, self.charge)

        if cell_v <= self.voltage_cutoff_discharge or cell_v >= self.voltage_cutoff_charge:
            raise ValueError("Desired current too high, overpotentials place cell voltage outside voltage cutoffs")

        while count != length_data:
            # set current
            i = cell_model.current_direction(self.charge) * self.current

            delta_ox, delta_red = cell_model.coulomb_counter(i, cls_degradation, ncls_degradation, crossover_params)

            # EDGE CASE where voltage limits never reached i.e. straight CC cycling until concentration runs out
            if cell_model.negative_concentrations():
                # record capacity here
                cycle_capacity.append(cap)

                #  Break out of loop if capacity approaches zero
                if cap < 1.0 and len(cycle_capacity) > 2:
                    print(str(count) + 'count')
                    print('Simulation stopped, capacity < 1 coulomb')
                    final_count = count
                    cap_low = True
                    break
                ##############
                cap = 0.0
                # record cycle time
                cycle_time.append(count * cell_model.time_increment)
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
            else:
                pass
            # calculate overpotentials and resulting cell voltage
            losses, n_act, n_mt = cell_model.total_overpotential(i, self.charge, i_lim_cls_t, i_lim_ncls_t)
            ocv = cell_model.open_circuit_voltage()
            cell_v = cell_model.cell_voltage(ocv, losses, self.charge)

            # did it hit voltage limit ?
            if cell_v >= self.voltage_cutoff_charge or cell_v <= self.voltage_cutoff_discharge:
                # record capacity here
                cycle_capacity.append(cap)
                #####################
                if cap < 1.0 and len(cycle_capacity) > 2:
                    print(str(count) + 'count')
                    print('Simulation stopped, capacity < 1 coulomb')
                    final_count = count
                    cap_low = True
                    break
                ##################
                cap = 0.0
                # record cycle time
                cycle_time.append(count * cell_model.time_increment)
                # switch charge to discharge or vice-versa
                self.charge = not self.charge

                # set limiting current for next cycle
                i_lim_cls_t, i_lim_ncls_t = cell_model.limiting_concentration(self.charge)
                continue
            else:
                pass

            # update capacity
            cap += abs(i * cell_model.time_increment)
            # update concentrations

            current_profile.append(i)
            c_ox_cls_profile.append(cell_model.c_ox_cls)
            c_red_cls_profile.append(cell_model.c_red_cls)
            c_ox_ncls_profile.append(cell_model.c_ox_ncls)
            c_red_ncls_profile.append(cell_model.c_red_ncls)

            del_ox.append(delta_ox)
            del_red.append(delta_red)

            cell_v_profile.append(cell_v)
            ocv_profile.append(ocv)

            act_profile.append(n_act)
            mt_profile.append(n_mt)
            loss_profile.append(losses)
            count += 1
        # adjusts time points if capacity decreased past set point
        if cap_low:
            times = times[:final_count + 1]

        # now calculating SOC of cls and ncls
        # can definitely be written better
        for a,b,c,d in zip(c_ox_cls_profile, c_red_cls_profile, c_ox_ncls_profile, c_red_ncls_profile):
            c, n = cell_model.state_of_charge(a, b, c, d)
            soc_profile_cls.append(c)
            soc_profile_ncls.append(n)

        return (current_profile, c_ox_cls_profile, c_red_cls_profile, c_ox_ncls_profile, c_red_ncls_profile,
                cell_v_profile, soc_profile_cls, soc_profile_ncls, ocv_profile, cycle_capacity, cycle_time, times,
                act_profile, mt_profile, loss_profile, del_ox, del_red)


class ConstantCurrentConstantVoltage(CyclingProtocol):
    """
    Subclass for constant current constant voltage (CCCV) cycling method.

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
        If True, charges the cell first.

    """

    def __init__(self, voltage_limit_charge: float, voltage_limit_discharge: float,
                 current_cutoff_charge: float, current_cutoff_discharge: float, current: float,
                 charge_first: bool = True):
        self.voltage_limit_charge = voltage_limit_charge
        self.voltage_limit_discharge = voltage_limit_discharge
        self.current_cutoff_charge = current_cutoff_charge
        self.current_cutoff_discharge = current_cutoff_discharge
        self.current = current
        self.charge = charge_first

        if self.current_cutoff_discharge > 0 or self.current_cutoff_charge < 0:
            raise ValueError("Invalid (dis)charge current cutoff")

    def run(self, duration: int, cell_model: ZeroDModel,
            degradation: DegradationMechanism = None,
            cls_degradation: DegradationMechanism = None,
            ncls_degradation: DegradationMechanism = None,
            crossover_params: Crossover = None) -> object:
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
        crossover_params : Crossover, optional
            Crossover mechanism applied to cell.

        Returns
        -------

        """

        if degradation is not None and (cls_degradation is not None or ncls_degradation is not None):
            raise ValueError("Cannot specify both 'degradation' and '(n)cls_degradation'")

        if degradation is not None:
            cls_degradation = degradation
            ncls_degradation = degradation

        # Raise ValueError if user inputs a negative concentration
        if cell_model.negative_concentrations():
            raise ValueError('Negative concentration detected')

        (c_ox_cls_profile, c_red_cls_profile, c_ox_ncls_profile, c_red_ncls_profile, cycle_capacity,
         current_profile, cell_v_profile, soc_profile_cls, ocv_profile, soc_profile_ncls, cycle_time, act_profile,
         mt_profile, loss_profile) = [], [], [], [], [], [], [], [], [], [], [], [], [], []
        del_ox, del_red = [], []

        # record temporary values of concentrations for all species
        c_ox_cls_temp = cell_model.c_ox_cls
        c_red_cls_temp = cell_model.c_red_cls
        c_ox_ncls_temp = cell_model.c_ox_ncls
        c_red_ncls_temp = cell_model.c_red_ncls

        cc_mode = True
        i_first = True
        cv_only = False

        length_data = int(duration / cell_model.time_increment)
        times = [x * cell_model.time_increment for x in range(1, length_data + 1)]
        print(f"{duration} sec of cycling, timesteps: {cell_model.time_increment} sec")

        count = 0
        cap = 0.0
        cap_low = False
        final_count = length_data

        i_lim_cls_t, i_lim_ncls_t = cell_model.limiting_concentration(self.charge)

        i = cell_model.current_direction(self.charge) * self.current
        #  check if cell needs to go straight to CV
        if self.current >= i_lim_cls_t or self.current >= i_lim_ncls_t:
            cc_mode = False
            cv_only = True
            print("Goes straight to CV cycling")
        else:
            losses, n_act, n_mt = cell_model.total_overpotential(i, self.charge, i_lim_cls_t, i_lim_ncls_t)
            ocv = cell_model.open_circuit_voltage()
            cell_v = cell_model.cell_voltage(ocv, losses, self.charge)

            if cell_v >= self.voltage_limit_charge or cell_v <= self.voltage_limit_discharge:
                cc_mode = False
                cv_only = True
                print("Has now switched to CV cycling")
            else:
                pass
        ##############
        while count != length_data:
            # check if in CC or CV mode
            if cc_mode:
                # set current
                i = cell_model.current_direction(self.charge) * self.current

                # calculate species' concentrations
                delta_ox, delta_red = cell_model.coulomb_counter(i, cls_degradation, ncls_degradation, crossover_params)

                # EDGE CASE where voltage limits never reached i.e straight CC cycling
                if cell_model.negative_concentrations():
                    # record capacity here
                    cycle_capacity.append(cap)

                    #  Break out of loop if capacity near zero
                    if cap < 1.0 and len(cycle_capacity) > 2:
                        print(str(count) + 'count')
                        print('Simulation stopped, capacity < 1 coulomb')
                        final_count = count
                        cap_low = True
                        break

                    #####################
                    cap = 0.0
                    # record cycle time
                    cycle_time.append(count * cell_model.time_increment)
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
                else:
                    pass

                # calculate overpotentials and resulting cell voltage
                losses, n_act, n_mt = cell_model.total_overpotential(i, self.charge, i_lim_cls_t, i_lim_ncls_t)
                ocv = cell_model.open_circuit_voltage()
                cell_v = cell_model.cell_voltage(ocv, losses, self.charge)

                # update capacity
                cap += abs(i * cell_model.time_increment)

                # check if V limit is reached?
                if cell_v >= self.voltage_limit_charge or cell_v <= self.voltage_limit_discharge:
                    cc_mode = not cc_mode
                    # at some point maybe record cap here too, so you know cap due to CC and due to CV?
                    continue
                else:
                    pass

                # update concentrations
                current_profile.append(i)  # CC section
                c_ox_cls_profile.append(cell_model.c_ox_cls)
                c_red_cls_profile.append(cell_model.c_red_cls)
                c_ox_ncls_profile.append(cell_model.c_ox_ncls)
                c_red_ncls_profile.append(cell_model.c_red_ncls)

                del_ox.append(delta_ox)
                del_red.append(delta_red)

                cell_v_profile.append(cell_v)
                ocv_profile.append(ocv)
                act_profile.append(n_act)
                mt_profile.append(n_mt)
                loss_profile.append(losses)
                count += 1
            ##########################################################
            else:  # now we're in CV mode

                # record temporary values of concentrations for all species
                c_ox_cls_temp = cell_model.c_ox_cls
                c_red_cls_temp = cell_model.c_red_cls
                c_ox_ncls_temp = cell_model.c_ox_ncls
                c_red_ncls_temp = cell_model.c_red_ncls

                # all cycling is now constant voltage
                cell_v = self.voltage_limit_charge if self.charge else self.voltage_limit_discharge
                ocv = cell_model.open_circuit_voltage()

                # variables to be passed into numerical solver to calculate current
                data = (cell_v, ocv, self.charge, i_lim_cls_t, i_lim_ncls_t)
                # adapting the solver's guess to the updated current
                # calculate the first current value for guess
                if i_first:
                    if cv_only:  # the case where you're straight CV cycling, and an initial current guess is needed
                        i_guess = cell_model.current_direction(self.charge) * self.current  # ?? think about this more
                    else:
                        i_guess = i
                        i_first = not i_first
                else:
                    pass
                ########################
                # testing
                '''
                if charge:
                    min_bracket = current_cutoff_charge*-1
                    max_bracket = ZeroDModel.current_direction(charge)*current
                else:
                    max_bracket = current_cutoff_discharge*-1
                    min_bracket = ZeroDModel.current_direction(charge)*current

                if i_first:
                    if charge:
                        max_bracket = 10
                    else:
                        min_bracket = -10

                i_cv = brentq(self.cv_current_solver, min_bracket, max_bracket, args=data)
                '''

                # consider different solver that ensure max allowable current upon switch to CV ?
                i_cv = float(fsolve(cell_model.cv_current_solver, np.array([i_guess]), args=data)[0])
                i_guess = i_cv

                # check if current is below cutoffs
                if (self.charge and i_cv <= self.current_cutoff_charge) or \
                        (not self.charge and i_cv >= self.current_cutoff_discharge):
                    # CV part of cycle has now ended, record capacity data

                    cycle_capacity.append(cap)

                    # Break out of full simulation if capacity nears zero
                    if cap < 1.0 and len(cycle_capacity) > 2:
                        print(str(count) + 'count')
                        print('Simulation stopped, capacity < 1 coulomb')
                        final_count = count
                        cap_low = True
                        break

                    ############################
                    cap = 0.0
                    # record cycle time
                    cycle_time.append(count * cell_model.time_increment)
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
                else:
                    pass

                delta_ox, delta_red = cell_model.coulomb_counter(i_cv, cls_degradation, ncls_degradation,
                                                                 crossover_params)

                # check if any reactant remains
                if cell_model.negative_concentrations():
                    cycle_capacity.append(cap)

                    # Break out of loop if capacity nears zero
                    if cap < 1.0 and len(cycle_capacity) > 2:
                        print(str(count) + 'count')
                        print('Simulation stopped, capacity < 1 coulomb')
                        final_count = count
                        cap_low = True
                        break

                    cap = 0.0
                    # record cycle time
                    cycle_time.append(count * cell_model.time_increment)
                    # switch charge to discharge or vice-versa
                    self.charge = not self.charge
                    cc_mode = True
                    i_first = True  # not sure if this would be needed

                    # set self back to previous, valid, concentration value
                    cell_model.c_ox_cls = c_ox_cls_temp
                    cell_model.c_red_cls = c_red_cls_temp
                    cell_model.c_ox_ncls = c_ox_ncls_temp
                    cell_model.c_red_ncls = c_red_ncls_temp

                    # set limiting current for next cycle
                    i_lim_cls_t, i_lim_ncls_t = cell_model.limiting_concentration(self.charge)
                    continue
                else:
                    pass

                # update capacity
                cap += abs(i_cv * cell_model.time_increment)
                # update concentrations

                current_profile.append(i_cv)
                c_ox_cls_profile.append(cell_model.c_ox_cls)
                c_red_cls_profile.append(cell_model.c_red_cls)
                c_ox_ncls_profile.append(cell_model.c_ox_ncls)
                c_red_ncls_profile.append(cell_model.c_red_ncls)

                # test
                del_ox.append(delta_ox)
                del_red.append(delta_red)

                cell_v_profile.append(cell_v)
                ocv_profile.append(ocv)
                act_profile.append(n_act)
                mt_profile.append(n_mt)
                loss_profile.append(losses)
                count += 1

        if cap_low:
            times = times[:final_count + 1]

        # now calculating SOC of cls and ncls
        # can defs make this better
        for a,b,c,d in zip(c_ox_cls_profile, c_red_cls_profile, c_ox_ncls_profile, c_red_ncls_profile):
            c, n = cell_model.state_of_charge(a, b, c, d)
            soc_profile_cls.append(c)
            soc_profile_ncls.append(n)

        return (current_profile, c_ox_cls_profile, c_red_cls_profile, c_ox_ncls_profile, c_red_ncls_profile,
                cell_v_profile, soc_profile_cls, soc_profile_ncls, ocv_profile, cycle_capacity, cycle_time, times,
                act_profile, mt_profile, loss_profile, del_ox, del_red)


if __name__ == '__main__':
    print('testing')

