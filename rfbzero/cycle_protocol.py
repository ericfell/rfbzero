

from abc import ABC, abstractmethod
from scipy.optimize import fsolve
from zeroD_model_1e_vs_1e import ZeroDModel
from zeroD_model_degradations import DegradationMechanism
from zeroD_model_crossover import Crossover


class CyclingProtocol(ABC):
    @abstractmethod
    def run(self, duration: int, cell_model: ZeroDModel,
            degradation: DegradationMechanism = None,
            crossover_params: Crossover = None,
            ) -> object:
        raise NotImplementedError


class ConstantCurrent(CyclingProtocol):

    def __init__(self, voltage_cutoff_charge: float, voltage_cutoff_discharge: float,
                 current: float, charge_first: bool = True):
        self.voltage_cutoff_charge = voltage_cutoff_charge
        self.voltage_cutoff_discharge = voltage_cutoff_discharge
        self.current = current
        self.charge_first = charge_first

    def run(self, duration: int, cell_model: ZeroDModel,
            degradation: DegradationMechanism = None,
            crossover_params: Crossover = None) -> object:

        (cls_start_c_ox, cls_start_c_red,
         ncls_start_c_ox, ncls_start_c_red) = cell_model.starting_concentrations()

        # Raise ValueError if user inputs a negative concentration
        if cell_model._negative_concentrations(cls_start_c_ox, cls_start_c_red,
                                               ncls_start_c_ox, ncls_start_c_red):
            raise ValueError('Negative concentration detected')

        (current_profile, c_ox_cls_profile, c_red_cls_profile, cell_v_profile, soc_profile_cls, ocv_profile,
         c_ox_ncls_profile, c_red_ncls_profile, soc_profile_ncls, cycle_capacity, cycle_time, act_profile,
         mt_profile, loss_profile) = [] ,[] ,[] ,[] ,[] ,[] ,[] ,[] ,[] ,[] ,[] ,[] ,[] ,[]
        # testing
        del_ox ,del_red = [] ,[]
        # set starting concentrations for all species
        conc_ox_now_CLS = cls_start_c_ox
        conc_red_now_CLS = cls_start_c_red
        conc_ox_now_NCLS = ncls_start_c_ox
        conc_red_now_NCLS = ncls_start_c_red
        #
        charge = self.charge_first  # ??
        #times = cell_model.times #self.times
        length_data = int(duration / cell_model.time_increment)
        times = [x*cell_model.time_increment for x in range(1, length_data + 1)]
        print(f"{duration} sec of cycling, timesteps: {cell_model.time_increment} sec")


        count = 0
        cap = 0.0
        cap_low = False
        # initialized in case simulation has to stop due to no more capacity
        final_count = length_data

        i_lim_cls_t, i_lim_ncls_t = cell_model.limiting_reactant_selector(charge, conc_ox_now_CLS, conc_red_now_CLS,
                                                                    conc_ox_now_NCLS, conc_red_now_NCLS)

        if (self.current >= (i_lim_cls_t * 2)) or (self.current >= (i_lim_ncls_t * 2)):
            raise ValueError("Desired current > limiting current, cell can't run")

        # assign + current to charge, - current to discharge
        i = cell_model._current_direction(charge) * self.current

        losses ,_ ,_ = cell_model.v_losses(i, charge, conc_ox_now_CLS, conc_red_now_CLS, conc_ox_now_NCLS, conc_red_now_NCLS,
                                   i_lim_cls_t, i_lim_ncls_t)
        #
        OCV = cell_model.nernst_OCV_full(conc_ox_now_CLS, conc_red_now_CLS, conc_ox_now_NCLS, conc_red_now_NCLS)

        cell_V = cell_model._cell_voltage(OCV, losses, charge)

        if (cell_V < self.voltage_cutoff_discharge) or (cell_V > self.voltage_cutoff_charge):
            raise ValueError("Desired current too high, overpotentials place cell voltage outside voltage cutoffs")

        while count != length_data:
            # set current
            i = cell_model._current_direction(charge) * self.current

            # need to do this for CCCV method below too
            (concentration_ox_CLS,
             concentration_red_CLS,
             concentration_ox_NCLS,
             concentration_red_NCLS,
             delta_ox, delta_red) = cell_model.coulomb_counter(i, conc_ox_now_CLS, conc_red_now_CLS, conc_ox_now_NCLS,
                                                         conc_red_now_NCLS, degradation, crossover_params)

            # EDGE CASE where voltage limits never reached i.e straight CC cycling until concentration runs out
            if cell_model._negative_concentrations(concentration_ox_CLS, concentration_red_CLS,
                                                   concentration_ox_NCLS, concentration_red_NCLS):
                # record capacity here
                cycle_capacity.append(cap)
                ############### Break out of loop if cpacity approaches zero
                if (cap < 1.0) and (len(cycle_capacity) > 2):
                    print(str(count) + 'count')
                    print('Simulation stopped, capacity < 1 coulomb')
                    final_count = count
                    # break out of while loop
                    count = length_data # not needed?
                    cap_low = True
                    break
                ##############
                cap = 0.0
                # record cycle time
                cycle_time.append(count * cell_model.time_increment)
                # switch charge to discharge or vice-versa
                charge = not charge

                # set limiting current for next cycle, with previous allowable concentrations
                i_lim_cls_t, i_lim_ncls_t = cell_model.limiting_reactant_selector(charge, conc_ox_now_CLS, conc_red_now_CLS,
                                                                            conc_ox_now_NCLS, conc_red_now_NCLS)
                continue
            else:
                pass
            # calculate overpotentials and resulting cell voltage
            losses, n_act, n_mt = cell_model.v_losses(i, charge, concentration_ox_CLS, concentration_red_CLS,
                                                concentration_ox_NCLS, concentration_red_NCLS, i_lim_cls_t, i_lim_ncls_t)

            OCV = cell_model.nernst_OCV_full(concentration_ox_CLS, concentration_red_CLS, concentration_ox_NCLS,
                                       concentration_red_NCLS)

            cell_V = cell_model._cell_voltage(OCV, losses, charge)

            # did it hit voltage limit ?
            if (cell_V >= self.voltage_cutoff_charge) or (cell_V <= self.voltage_cutoff_discharge):
                # record capacity here
                cycle_capacity.append(cap)
                #####################
                if (cap < 1.0) and (len(cycle_capacity) > 2):
                    print(str(count) + 'count')
                    print('Simulation stopped, capacity < 1 coulomb')
                    final_count = count
                    # to break out of while loop
                    count = length_data # not needed?
                    cap_low = True
                    break
                ##################
                cap = 0.0
                # record cycle time
                cycle_time.append(count * cell_model.time_increment)
                # switch charge to discharge or viceversa
                charge = not charge

                # set limiting current for next cycle
                i_lim_cls_t, i_lim_ncls_t = cell_model.limiting_reactant_selector(charge, conc_ox_now_CLS, conc_red_now_CLS,
                                                                            conc_ox_now_NCLS, conc_red_now_NCLS)
                continue
            else:
                pass

            # calculate SOC (local, in this case) ** for CLS make this a function, can be other file
            soc_CLS, soc_NCLS = cell_model._soc(concentration_ox_CLS, concentration_red_CLS, concentration_ox_NCLS,
                                                concentration_red_NCLS)
            # update capacity
            cap += abs( i * cell_model.time_increment)
            # update concentrations
            conc_ox_now_CLS = concentration_ox_CLS
            conc_red_now_CLS = concentration_red_CLS
            conc_ox_now_NCLS = concentration_ox_NCLS
            conc_red_now_NCLS = concentration_red_NCLS
            #
            current_profile.append(i)
            c_ox_cls_profile.append(concentration_ox_CLS)
            c_red_cls_profile.append(concentration_red_CLS)
            c_ox_ncls_profile.append(concentration_ox_NCLS)
            c_red_ncls_profile.append(concentration_red_NCLS)

            # test
            del_ox.append(delta_ox)
            del_red.append(delta_red)

            cell_v_profile.append(cell_V)
            ocv_profile.append(OCV)
            soc_profile_cls.append(soc_CLS)
            soc_profile_ncls.append(soc_NCLS)
            act_profile.append(n_act)
            mt_profile.append(n_mt)
            loss_profile.append(losses)
            count += 1
        # adjusts time points if capacity decreased past set point
        if cap_low:
            times = times[:final_count + 1]
        return (current_profile, c_ox_cls_profile, c_red_cls_profile, c_ox_ncls_profile, c_red_ncls_profile,
                cell_v_profile, soc_profile_cls, soc_profile_ncls, ocv_profile, cycle_capacity, cycle_time, times,
                act_profile, mt_profile, loss_profile, del_ox, del_red)


class ConstantCurrentConstantVoltage(CyclingProtocol):

    def __init__(self, voltage_limit_charge: float, voltage_limit_discharge: float,
                 current_cutoff_charge: float, current_cutoff_discharge: float, current: float,
                 charge_first: bool = True):
        self.voltage_limit_charge = voltage_limit_charge
        self.voltage_limit_discharge = voltage_limit_discharge
        self.current_cutoff_charge = current_cutoff_charge
        self.current_cutoff_discharge = current_cutoff_discharge
        self.current = current
        self.charge_first = charge_first

    def run(self,  duration: int, cell_model: ZeroDModel,
            degradation: DegradationMechanism = None,
            crossover_params: Crossover = None) -> object:

        (cls_start_c_ox, cls_start_c_red,
         ncls_start_c_ox, ncls_start_c_red) = cell_model.starting_concentrations()

        # to dos: if full CV, then should be appending voltage to overpotential

        # Raise ValueError if user inputs a negative concentration
        if cell_model._negative_concentrations(cls_start_c_ox, cls_start_c_red,
                                               ncls_start_c_ox, ncls_start_c_red):
            raise ValueError('Negative concentration detected')

        (conc_ox_CLS_profile, conc_red_CLS_profile, conc_ox_NCLS_profile, conc_red_NCLS_profile, cycle_capacity,
         current_profile, cell_V_profile, soc_profile_CLS, ocv_profile, soc_profile_NCLS, cycle_time, act_profile,
         mt_profile, loss_profile) = [], [], [], [], [], [], [], [], [], [], [], [], [], []
        del_ox, del_red = [], []
        # set starting concentrations for all species
        # simplify this ?????
        conc_ox_now_CLS = cls_start_c_ox
        conc_red_now_CLS = cls_start_c_red
        conc_ox_now_NCLS = ncls_start_c_ox
        conc_red_now_NCLS = ncls_start_c_red

        assert self.current_cutoff_discharge < 0, "invalid discharge current cutoff"
        #
        charge = self.charge_first  # True
        CC_mode = True
        i_first = True
        CV_only = False
        #times = cell_model.times #self.times
        length_data = int(duration / cell_model.time_increment)
        times = [x * cell_model.time_increment for x in range(1, length_data + 1)]
        print(f"{duration} sec of cycling, timesteps: {cell_model.time_increment} sec")

        count = 0
        cap = 0.0
        cap_low = False
        final_count = length_data

        i_lim_cls_t, i_lim_ncls_t = cell_model.limiting_reactant_selector(charge, conc_ox_now_CLS, conc_red_now_CLS,
                                                                    conc_ox_now_NCLS, conc_red_now_NCLS)

        i = cell_model._current_direction(charge) * self.current
        ##### check if need to go straight to CV
        # allow option to say 0 current for straight CV?
        if (self.current >= (i_lim_cls_t * 2)) or (self.current >= (i_lim_ncls_t * 2)):
            CC_mode = False
            CV_only = True
            print("Goes straight to CV cycling")
        else:
            losses, n_act, n_mt = cell_model.v_losses(i, charge, conc_ox_now_CLS, conc_red_now_CLS, conc_ox_now_NCLS,
                                                conc_red_now_NCLS, i_lim_cls_t, i_lim_ncls_t)
            ## OCV due to CLS and NCLS
            OCV = cell_model.nernst_OCV_full(conc_ox_now_CLS, conc_red_now_CLS, conc_ox_now_NCLS, conc_red_now_NCLS)
            cell_V = cell_model._cell_voltage(OCV, losses, charge)
            if (cell_V >= self.voltage_limit_charge) or (cell_V <= self.voltage_limit_discharge):
                CC_mode = False
                CV_only = True
                print("Has now switched to CV cycling")
            else:
                pass
        ##############
        while count != length_data:
            # check if in CC or CV mode
            if CC_mode:
                # set current
                i = cell_model._current_direction(charge) * self.current

                # calculate species' concentrations
                (concentration_ox_CLS,
                 concentration_red_CLS,
                 concentration_ox_NCLS,
                 concentration_red_NCLS,
                 delta_ox, delta_red) = cell_model.coulomb_counter(i, conc_ox_now_CLS, conc_red_now_CLS, conc_ox_now_NCLS,
                                                             conc_red_now_NCLS, degradation, crossover_params)

                # EDGE CASE where voltage limits never reached i.e straight CC cycling
                if cell_model._negative_concentrations(concentration_ox_CLS, concentration_red_CLS,
                                                       concentration_ox_NCLS, concentration_red_NCLS):
                    # record capacity here
                    cycle_capacity.append(cap)

                    ############ Break out of loop if capacity near zero
                    if (cap < 1.0) and (len(cycle_capacity) > 2):
                        print(str(count) + 'count')
                        print('Simulation stopped, capacity < 1 coulomb')
                        final_count = count
                        # count = self.length_data # not needed?
                        cap_low = True
                        break

                    #####################
                    cap = 0.0
                    # record cycle time
                    cycle_time.append(count * cell_model.time_increment)
                    # switch charge to discharge or viceversa
                    charge = not charge

                    # set limiting current for next cycle
                    i_lim_cls_t, i_lim_ncls_t = cell_model.limiting_reactant_selector(charge, conc_ox_now_CLS,
                                                                                conc_red_now_CLS, conc_ox_now_NCLS,
                                                                                conc_red_now_NCLS)
                    continue
                else:
                    pass

                # calculate overpotentials and resulting cell voltage
                losses, n_act, n_mt = cell_model.v_losses(i, charge, concentration_ox_CLS, concentration_red_CLS,
                                                    concentration_ox_NCLS, concentration_red_NCLS, i_lim_cls_t,
                                                    i_lim_ncls_t)

                # OCV due to CLS and NCLS
                OCV = cell_model.nernst_OCV_full(concentration_ox_CLS, concentration_red_CLS, concentration_ox_NCLS,
                                           concentration_red_NCLS)

                cell_V = cell_model._cell_voltage(OCV, losses, charge)

                # calculate SOC (local, not global if there's loss mechanisms)
                soc_CLS, soc_NCLS = cell_model._soc(concentration_ox_CLS, concentration_red_CLS, concentration_ox_NCLS,
                                                    concentration_red_NCLS)
                # update capacity
                cap += abs(i * cell_model.time_increment)

                # check if V limit is reached?
                if (cell_V >= self.voltage_limit_charge) or (cell_V <= self.voltage_limit_discharge):
                    CC_mode = not CC_mode
                    # at some point maybe record cap here too so you know cap due to CC and due to CV?
                    continue
                else:
                    pass

                # update concentrations
                conc_ox_now_CLS = concentration_ox_CLS
                conc_red_now_CLS = concentration_red_CLS
                conc_ox_now_NCLS = concentration_ox_NCLS
                conc_red_now_NCLS = concentration_red_NCLS
                ##
                current_profile.append(i)  # CC section
                conc_ox_CLS_profile.append(concentration_ox_CLS)
                conc_red_CLS_profile.append(concentration_red_CLS)
                conc_ox_NCLS_profile.append(concentration_ox_NCLS)
                conc_red_NCLS_profile.append(concentration_red_NCLS)

                # test
                del_ox.append(delta_ox)
                del_red.append(delta_red)

                cell_V_profile.append(cell_V)
                ocv_profile.append(OCV)
                soc_profile_CLS.append(soc_CLS)
                soc_profile_NCLS.append(soc_NCLS)
                act_profile.append(n_act)
                mt_profile.append(n_mt)
                loss_profile.append(losses)
                count += 1
            ##########################################################
            else:  # now we're in CV mode

                # all constant voltage here
                cell_V = self.voltage_limit_charge if charge else self.voltage_limit_discharge

                OCV = cell_model.nernst_OCV_full(conc_ox_now_CLS, conc_red_now_CLS, conc_ox_now_NCLS, conc_red_now_NCLS)

                data = (cell_V, OCV, charge, conc_ox_now_CLS, conc_red_now_CLS, conc_ox_now_NCLS, conc_red_now_NCLS,
                        i_lim_cls_t, i_lim_ncls_t)
                # adapting the solver's guess to the updated current
                # calculate the first current value for guess
                if i_first:
                    if CV_only:  # the case where you're straight CV cycling, and an initial current guess is needed
                        i_guess = cell_model._current_direction(charge) * self.current  # ??
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
                    max_bracket = ZeroDModel._current_direction(charge)*current
                else:
                    max_bracket = current_cutoff_discharge*-1
                    min_bracket = ZeroDModel._current_direction(charge)*current

                if i_first:
                    if charge:
                        max_bracket = 10
                    else:
                        min_bracket = -10

                i_CV = brentq(self.cv_current_solver, min_bracket, max_bracket, args=data)
                '''
                ######
                ######### end testing ############
                # consider different solver that ensure max allowable current upon switch to CV ?
                i_CV = float(fsolve(cell_model.cv_current_solver, i_guess, args=data)[0])
                i_guess = i_CV

                # check if current is below cutoffs
                if (charge and (i_CV <= self.current_cutoff_charge)) or (not charge and (i_CV >= self.current_cutoff_discharge)):
                    # CV part of cycle has now ended, record capacity data

                    cycle_capacity.append(cap)
                    ############ Break out of full simulation if capacity nears zero

                    if (cap < 1.0) and (len(cycle_capacity) > 2):
                        print(str(count) + 'count')
                        print('Simulation stopped, capacity < 1 coulomb')
                        final_count = count
                        # count = self.length_data # not needed?
                        cap_low = True
                        break

                    ############################
                    cap = 0.0
                    # record cycle time
                    cycle_time.append(count * cell_model.time_increment)
                    # switch charge to discharge or viceversa
                    charge = not charge
                    CC_mode = True
                    i_first = True

                    # set limiting current for next cycle
                    i_lim_cls_t, i_lim_ncls_t = cell_model.limiting_reactant_selector(charge, conc_ox_now_CLS,
                                                                                conc_red_now_CLS, conc_ox_now_NCLS,
                                                                                conc_red_now_NCLS)
                    continue
                else:
                    pass

                # testing new. should i_CV be i_guesss?
                (concentration_ox_CLS,
                 concentration_red_CLS,
                 concentration_ox_NCLS,
                 concentration_red_NCLS,
                 delta_ox, delta_red) = cell_model.coulomb_counter(i_CV, conc_ox_now_CLS, conc_red_now_CLS, conc_ox_now_NCLS,
                                                             conc_red_now_NCLS, degradation, crossover_params)

                # check if any reactant remains
                if cell_model._negative_concentrations(concentration_ox_CLS, concentration_red_CLS,
                                                       concentration_ox_NCLS, concentration_red_NCLS):

                    cycle_capacity.append(cap)
                    ############### Break out of loop if capacity nears zero

                    if (cap < 1.0) and (len(cycle_capacity) > 2):
                        print(str(count) + 'count')
                        print('Simulation stopped, capacity < 1 coulomb')
                        final_count = count
                        # count = self.length_data # not needed?
                        cap_low = True
                        break

                    ###########
                    cap = 0.0
                    # record cycle time
                    cycle_time.append(count * cell_model.time_increment)
                    # switch charge to discharge or viceversa
                    charge = not charge
                    CC_mode = True
                    i_first = True  # not sure if this would be needed

                    # set limiting current for next cycle
                    i_lim_cls_t, i_lim_ncls_t = cell_model.limiting_reactant_selector(charge, conc_ox_now_CLS,
                                                                                conc_red_now_CLS, conc_ox_now_NCLS,
                                                                                conc_red_now_NCLS)
                    continue
                else:
                    pass

                # calculate SOC (local, in this case)
                soc_CLS, soc_NCLS = cell_model._soc(concentration_ox_CLS, concentration_red_CLS,
                                                    concentration_ox_NCLS, concentration_red_NCLS)
                # update capacity
                cap += abs(i_CV * cell_model.time_increment)
                # update concentrations
                conc_ox_now_CLS = concentration_ox_CLS
                conc_red_now_CLS = concentration_red_CLS
                conc_ox_now_NCLS = concentration_ox_NCLS
                conc_red_now_NCLS = concentration_red_NCLS
                ##
                current_profile.append(i_CV)
                conc_ox_CLS_profile.append(concentration_ox_CLS)
                conc_red_CLS_profile.append(concentration_red_CLS)
                conc_ox_NCLS_profile.append(concentration_ox_NCLS)
                conc_red_NCLS_profile.append(concentration_red_NCLS)

                # test
                del_ox.append(delta_ox)
                del_red.append(delta_red)

                cell_V_profile.append(cell_V)
                ocv_profile.append(OCV)
                soc_profile_CLS.append(soc_CLS)
                soc_profile_NCLS.append(soc_NCLS)
                act_profile.append(n_act)
                mt_profile.append(n_mt)
                loss_profile.append(losses)
                count += 1

        if cap_low:
            times = times[:final_count + 1]
        return (current_profile, conc_ox_CLS_profile, conc_red_CLS_profile, conc_ox_NCLS_profile, conc_red_NCLS_profile,
                cell_V_profile, soc_profile_CLS, soc_profile_NCLS, ocv_profile, cycle_capacity, cycle_time, times,
                act_profile, mt_profile, loss_profile, del_ox, del_red)








