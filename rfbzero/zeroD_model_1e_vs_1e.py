

from math import log
import scipy.constants as spc
from scipy.optimize import fsolve
from zeroD_model_degradations import DegradationMechanism #degradation_mechanism,
from zeroD_model_crossover import Crossover


# Faraday constant (C/mol)
F = spc.value('Faraday constant')

# Molar gas constant, J/K/mol
R = spc.R

# make these parameters at some point?
TEMPERATURE = 298  # Kelvins, for S.T.P.
NERNST_CONST = (R*TEMPERATURE) / F  # should have n_electrons input option


class ZeroDModel:
    """
    Zero dimensional model for redox flow battery (RFB) cycling [1].

    Parameters
    ----------
    geometric_area : float
        Geometric area of cell (cm^2).
    resistance : float
        Cell ohmic resistance (ohms).
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
    duration : int
        Amount of experiment time to simulate (s).
    time_increment : float
        Simulation time step (s).
    init_ocv : float
        Cell voltage (formal potentials E_+ - E_-) (V).
        If voltage > 0 then it's a Full cell.
        If voltage = 0 then it's a Symmetric cell.
    k_mt : float
        Mass transport coefficient (cm/s).
    roughness_factor : float
        Roughness factor, dimensionless.
        Surface area divided by geometric surface area.
    k_0_cls : float
        Electrochemical rate constant, CLS redox couple (cm/s).
    k_0_ncls : float
        Electrochemical rate constant, NCLS redox couple (cm/s).
    alpha_cls : float
        Charge transfer coefficient of CLS redox couple, dimensionless.
    alpha_ncls : float
        Charge transfer coefficient of NCLS redox couple, dimensionless.
    cls_negolyte : bool
        If True, negolyte is the CLS.
    mechanism_list : list, optional
        List of degradation mechanisms to include in simulation.
    ###### mechanism_params : dict, optional
    #####    Parameters for mechanisms specified in `mechanism_list`.
    crossover_params : list, optional
        List of parameters for crossover mechanism.


    Notes
    -----
    All equations are adapted from [1]. If ZeroDModel has been
    significant to your research please cite the paper.

    [1] Modak, S.; Kwabi, D. G. A Zero-Dimensional Model for Electrochemical
    Behavior and Capacity Retention in Organic Flow Cells, Journal of The
    Electrochemical Society, 168, 2021, 080528.
    """

    def __init__(self, geometric_area: float, resistance: float, cls_volume: float, ncls_volume: float,
                 cls_start_c_ox: float, cls_start_c_red: float, ncls_start_c_ox: float, ncls_start_c_red: float,
                 duration: int, time_increment: float, init_ocv: float, k_mt: float, roughness_factor: float,
                 k_0_cls: float, k_0_ncls: float, alpha_cls: float, alpha_ncls: float, cls_negolyte: bool = True,
                 #mechanism_list: DegradationMechanism = None, mechanism_params: dict = None, crossover_params: list = None) -> None:
                 mechanism_list: DegradationMechanism = None, crossover_params: Crossover = None) -> None:
        """Inits ZeroDModel"""

        self.geometric_area = geometric_area
        self.resistance = resistance
        self.cls_volume = cls_volume
        self.ncls_volume = ncls_volume
        self.cls_start_c_ox = cls_start_c_ox
        self.cls_start_c_red = cls_start_c_red
        self.ncls_start_c_ox = ncls_start_c_ox
        self.ncls_start_c_red = ncls_start_c_red
        self.duration = duration
        self.time_increment = time_increment
        self.init_ocv = init_ocv
        self.k_mt = k_mt
        self.k_0_cls = k_0_cls
        self.k_0_ncls = k_0_ncls
        self.alpha_cls = alpha_cls
        self.alpha_ncls = alpha_ncls
        self.cls_negolyte = cls_negolyte
        self.mechanism_list = mechanism_list  # need to rename
        # self.mechanism_params = mechanism_params
        self.crossover_params = crossover_params
        self.const_i_ex = F * roughness_factor * self.geometric_area
        self.length_data = int(self.duration / self.time_increment)
        self.times = [x*self.time_increment for x in range(1, self.length_data + 1)]
        print(f"{self.duration} s of cycling, {self.time_increment} s timesteps,"
              f"{self.length_data} total datapoints")

    # option for just measuring capacity over time, doesn't need to make all arrays?

    @staticmethod
    def _current_direction(charge: bool) -> int:
        """Make current positive for charge, negative for discharge"""
        return 1 if charge else -1

    def i_exchange_current(self, c_ox_cls: float, c_red_cls: float, c_ox_ncls: float,
                           c_red_ncls: float) -> tuple[float, float]:
        """
        Calculates exchange current of redox couples in the CLS and NCLS.

        Parameters
        ----------
        c_ox_cls: float
            Concentration of oxidized species in CLS
             at a given timestep (M).
        c_red_cls: float
            Concentration of reduced species in CLS
             at a given timestep (M).
        c_ox_ncls: float
            Concentration of oxidized species in NCLS
             at a given timestep (M).
        c_red_ncls: float
            Concentration of reduced species in NCLS
             at a given timestep (M).

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
        i_0_cls = (self.const_i_ex * self.k_0_cls * (c_red_cls ** self.alpha_cls)
                   * (c_ox_cls ** (1 - self.alpha_cls)) * 0.001)
        i_0_ncls = (self.const_i_ex * self.k_0_ncls * (c_red_ncls ** self.alpha_ncls)
                    * (c_ox_ncls ** (1 - self.alpha_ncls)) * 0.001)
        return i_0_cls, i_0_ncls

    def i_limiting(self, c_lim: float) -> float:
        """Calculates limiting current for a single reservoir.
        This is equation 6 of [1].
        """
        # div by 1000 for conversion from mol/L to mol/cm^3
        # will require n electrons param
        return F * self.k_mt * c_lim * self.geometric_area * 0.001

    def limiting_reactant_selector(self, charge: bool, c_ox_cls: float, c_red_cls: float, c_ox_ncls: float,
                                   c_red_ncls: float) -> tuple[float, float]:
        """Selects limiting concentration and calculates limiting current for CLS and NCLS."""
        if (self.cls_negolyte and charge) or (not self.cls_negolyte and not charge):
            i_lim_cls = self.i_limiting(c_ox_cls)
            i_lim_ncls = self.i_limiting(c_red_ncls)
        else:
            i_lim_cls = self.i_limiting(c_red_cls)
            i_lim_ncls = self.i_limiting(c_ox_ncls)

        return i_lim_cls, i_lim_ncls

    @staticmethod
    def n_activation(current: float, i_0_cls: float, i_0_ncls: float) -> float:
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

    @staticmethod
    def _negative_concentrations(c_ox_cls: float, c_red_cls: float, c_ox_ncls: float, c_red_ncls: float) -> bool:
        """Return True if any concentration is negative"""
        return any(x < 0.0 for x in [c_ox_cls, c_red_cls, c_ox_ncls, c_red_ncls])

    def n_mass_transport(self, charge: bool, current: float, c_ox_cls: float, c_red_cls: float,
                         c_ox_ncls: float, c_red_ncls: float, i_lim_cls: float, i_lim_ncls: float) -> float:
        """
        This is equation 8 of [1].

        Parameters
        ----------
        charge
        current
        c_ox_cls
        c_red_cls
        c_ox_ncls
        c_red_ncls
        i_lim_cls
        i_lim_ncls

        Returns
        -------

        """
        # Raise ValueError if a negative concentration is detected
        if ZeroDModel._negative_concentrations(c_ox_cls, c_red_cls, c_ox_ncls, c_red_ncls):
            raise ValueError('Negative concentration detected')

        c_tot_cls = c_red_cls + c_ox_cls
        c_tot_ncls = c_red_ncls + c_ox_ncls

        i = abs(current)

        if (self.cls_negolyte and charge) or (not self.cls_negolyte and not charge):
            n_mt = NERNST_CONST * log((1 - ((c_tot_cls * i) / ((c_red_cls * i_lim_cls) + (c_ox_cls * i))))
                                      * (1 - ((c_tot_ncls * i) / ((c_ox_ncls * i_lim_ncls) + (c_red_ncls * i)))))
        else:
            n_mt = NERNST_CONST * log(((1 - ((c_tot_cls * i) / ((c_ox_cls * i_lim_cls) + (c_red_cls * i))))
                                       * (1 - ((c_tot_ncls * i) / ((c_red_ncls * i_lim_ncls) + (c_ox_ncls * i))))))
        return n_mt

    def v_losses(self, current: float, charge: bool, c_ox_cls: float, c_red_cls: float, c_ox_ncls: float,
                 c_red_ncls: float, i_lim_cls: float, i_lim_ncls: float) -> tuple[float, float, float]:
        """
        This is the overpotentials of equation 2 in [1].

        Parameters
        ----------
        current
        charge
        c_ox_cls
        c_red_cls
        c_ox_ncls
        c_red_ncls
        i_lim_cls
        i_lim_ncls

        Returns
        -------

        """

        i_0_cls, i_0_ncls = self.i_exchange_current(c_ox_cls, c_red_cls, c_ox_ncls, c_red_ncls)
        # calculate ohmic, activation, mass transport overpotentials
        n_ohmic = abs(current)*self.resistance
        n_act = self.n_activation(current, i_0_cls, i_0_ncls)
        n_mt = self.n_mass_transport(charge, current, c_ox_cls, c_red_cls, c_ox_ncls, c_red_ncls, i_lim_cls, i_lim_ncls)

        n_loss = n_ohmic + n_act + n_mt

        return n_loss, n_act, n_mt

    def nernst_OCV_full(self, c_ox_cls: float, c_red_cls: float, c_ox_ncls: float, c_red_ncls: float) -> float:
        """
        This is equivalent to equation 3 of [1].

        Parameters
        ----------
        c_ox_cls
        c_red_cls
        c_ox_ncls
        c_red_ncls

        Returns
        -------

        """

        # Raise ValueError if a negative concentration is detected
        if ZeroDModel._negative_concentrations(c_ox_cls, c_red_cls, c_ox_ncls, c_red_ncls):
            raise ValueError('Negative concentration detected')

        # will need n_electrons input
        # CLS is negolyte
        if self.cls_negolyte:
            ocv = (self.init_ocv
                   + (NERNST_CONST * log(c_red_cls / c_ox_cls)) + (NERNST_CONST * log(c_ox_ncls / c_red_ncls)))

        # CLS is posolyte
        else:
            ocv = (self.init_ocv
                   - (NERNST_CONST * log(c_red_cls / c_ox_cls)) - (NERNST_CONST * log(c_ox_ncls / c_red_ncls)))
        return ocv

    @staticmethod
    def _cell_voltage(ocv: float, losses: float, charge: bool) -> float:
        """If charging, add overpotentials to OCV, else subtract them."""
        return ocv + losses if charge else ocv - losses

    def coulomb_counter(self, current: float, c_ox_cls: float, c_red_cls: float, c_ox_ncls: float,
                        c_red_ncls: float) -> tuple[float, float, float, float, float, float]:

        direction = 1 if self.cls_negolyte else -1
        delta_cls = ((self.time_increment * current) / (F * self.cls_volume)) * direction
        delta_ncls = ((self.time_increment * current) / (F * self.ncls_volume)) * direction

        # update CLS concentrations
        c_ox_cls = c_ox_cls - delta_cls
        c_red_cls = c_red_cls + delta_cls
        # update NCLS concentrations
        c_ox_ncls = c_ox_ncls + delta_ncls
        c_red_ncls = c_red_ncls - delta_ncls

        # for no crossover situation
        delta_ox = 0.0
        delta_red = 0.0

        # Now consider effects of user-input degradations and/or crossover

        # no degradation / no crossover
        if (self.mechanism_list is None) and (self.crossover_params is None):
            pass

        # degradation / no crossover
        elif (self.mechanism_list is not None) and (self.crossover_params is None):
            # possible CLS degradation
            """
            c_ox_cls, c_red_cls = degradation_mechanism(c_ox_cls, c_red_cls, self.time_increment,
                                                        *self.mechanism_list, **self.mechanism_params)
            # possible NCLS degradation
            c_ox_ncls, c_red_ncls = degradation_mechanism(c_ox_ncls, c_red_ncls, self.time_increment,
                                                          *self.mechanism_list, **self.mechanism_params)
            """
            # testing abstract method class
            c_ox_cls, c_red_cls = self.mechanism_list.degrade(c_ox_cls, c_red_cls, self.time_increment)
            c_ox_ncls, c_red_ncls = self.mechanism_list.degrade(c_ox_ncls, c_red_ncls, self.time_increment)

        # no degradation / crossover
        elif (self.mechanism_list is None) and (self.crossover_params is not None):
            """
            (c_ox_cls, c_red_cls, c_ox_ncls, c_red_ncls, delta_ox,
             delta_red) = crossover(c_ox_cls, c_red_cls, c_ox_ncls, c_red_ncls, self.time_increment,
                                    self.cls_volume, self.ncls_volume, *self.crossover_params)
            """

            (c_ox_cls, c_red_cls, c_ox_ncls, c_red_ncls, delta_ox,
             delta_red) = self.crossover_params.crossover(c_ox_cls, c_red_cls, c_ox_ncls, c_red_ncls,
                                                          self.time_increment, self.cls_volume, self.ncls_volume)

        # degradation AND crossover
        else:
            """
            # CLS degradation
            c_ox_cls, c_red_cls = degradation_mechanism(c_ox_cls, c_red_cls, self.time_increment,
                                                        *self.mechanism_list, **self.mechanism_params)
            # NCLS degradation
            c_ox_ncls, c_red_ncls = degradation_mechanism(c_ox_ncls, c_red_ncls, self.time_increment,
                                                          *self.mechanism_list, **self.mechanism_params)
            """
            # testing abstract method class
            c_ox_cls, c_red_cls = self.mechanism_list.degrade(c_ox_cls, c_red_cls, self.time_increment)
            c_ox_ncls, c_red_ncls = self.mechanism_list.degrade(c_ox_ncls, c_red_ncls, self.time_increment)

            # crossover
            """
            (c_ox_cls, c_red_cls, c_ox_ncls, c_red_ncls, delta_ox,
             delta_red) = crossover(c_ox_cls, c_red_cls, c_ox_ncls, c_red_ncls, self.time_increment,
                                    self.cls_volume, self.ncls_volume, *self.crossover_params)
            """
            (c_ox_cls, c_red_cls, c_ox_ncls, c_red_ncls, delta_ox,
             delta_red) = self.crossover_params.crossover(c_ox_cls, c_red_cls, c_ox_ncls, c_red_ncls,
                                                          self.time_increment, self.cls_volume, self.ncls_volume)

        return c_ox_cls, c_red_cls, c_ox_ncls, c_red_ncls, delta_ox, delta_red

    @staticmethod
    def _soc(c_ox_cls: float, c_red_cls: float, c_ox_ncls: float, c_red_ncls: float) -> tuple[float, float]:
        """Calculate state-of-charge in each reservoir"""
        soc_cls = (c_red_cls / (c_ox_cls + c_red_cls)) * 100
        # this could be defined differently i.e., cell vs reservoir definition of SOC
        soc_ncls = (c_red_ncls / (c_ox_ncls + c_red_ncls)) * 100
        return soc_cls, soc_ncls

    # is below proper *args unpacking naming style?
    def cv_current_solver(self, current: float, *data: float) -> float:
        (cell_v, ocv, charge, c_ox_cls, c_red_cls, c_ox_ncls, c_red_ncls, i_lim_cls, i_lim_ncls) = data
        # curr has sign but v_losses makes it always positive
        loss_solve, _, _ = self.v_losses(current, charge, c_ox_cls, c_red_cls, c_ox_ncls, c_red_ncls, i_lim_cls,
                                         i_lim_ncls)
        # returns what solver will try to minimize
        return cell_v - ocv - loss_solve if charge else cell_v - ocv + loss_solve


    ##### Section: Cycling models

    def cc_experiment(self, voltage_cutoff_charge: float, voltage_cutoff_discharge: float,
                      current: float, charge_first: bool) -> object: # edit the return type soon

        # Raise ValueError if user inputs a negative concentration
        if ZeroDModel._negative_concentrations(self.cls_start_c_ox, self.cls_start_c_red,
                                               self.ncls_start_c_ox, self.ncls_start_c_red):
            raise ValueError('Negative concentration detected')

        (current_profile, c_ox_cls_profile, c_red_cls_profile, cell_v_profile, soc_profile_cls, ocv_profile,
         c_ox_ncls_profile, c_red_ncls_profile, soc_profile_ncls, cycle_capacity, cycle_time, act_profile,
         mt_profile, loss_profile) = [],[],[],[],[],[],[],[],[],[],[],[],[],[]
        # testing
        del_ox,del_red = [],[]
        # set starting concentrations for all species
        conc_ox_now_CLS = self.cls_start_c_ox
        conc_red_now_CLS = self.cls_start_c_red
        conc_ox_now_NCLS = self.ncls_start_c_ox
        conc_red_now_NCLS = self.ncls_start_c_red
        # 
        charge = charge_first #??
        times = self.times

        count = 0
        cap = 0.0
        cap_low = False
        # initialized in case simulation has to stop due to no more capacity
        final_count = self.length_data
        
        i_lim_cls_t, i_lim_ncls_t = self.limiting_reactant_selector(charge, conc_ox_now_CLS, conc_red_now_CLS,
                                                                    conc_ox_now_NCLS, conc_red_now_NCLS)

        if (current >= (i_lim_cls_t * 2)) or (current >= (i_lim_ncls_t * 2)):
            raise ValueError("Desired current > limiting current, cell can't run")

        # assign + current to charge, - current to discharge
        i = ZeroDModel._current_direction(charge) * current
        
        losses,_,_ = self.v_losses(i, charge, conc_ox_now_CLS, conc_red_now_CLS, conc_ox_now_NCLS, conc_red_now_NCLS,
                                   i_lim_cls_t, i_lim_ncls_t)
        #
        OCV = self.nernst_OCV_full(conc_ox_now_CLS, conc_red_now_CLS, conc_ox_now_NCLS, conc_red_now_NCLS)
    
        cell_V = ZeroDModel._cell_voltage(OCV, losses, charge)

        if (cell_V < voltage_cutoff_discharge) or (cell_V > voltage_cutoff_charge):
            raise ValueError("Desired current too high, overpotentials place cell voltage outside voltage cutoffs")
        
        while count != self.length_data:
            # set current
            i = ZeroDModel._current_direction(charge) * current

            # need to do this for CCCV method below too
            (concentration_ox_CLS,
             concentration_red_CLS,
             concentration_ox_NCLS,
             concentration_red_NCLS,
             delta_ox, delta_red) = self.coulomb_counter(i, conc_ox_now_CLS, conc_red_now_CLS, conc_ox_now_NCLS,
                                                         conc_red_now_NCLS)

            # EDGE CASE where voltage limits never reached i.e straight CC cycling until concentration runs out
            if ZeroDModel._negative_concentrations(concentration_ox_CLS, concentration_red_CLS,
                                                   concentration_ox_NCLS, concentration_red_NCLS):
                # record capacity here
                cycle_capacity.append(cap)
                ############### Break out of loop if cpacity approaches zero
                if (cap < 1.0) and (len(cycle_capacity) > 2):
                    print(str(count) + 'count')
                    print('Simulation stopped, capacity < 1 coulomb')
                    final_count = count
                    # break out of while loop
                    count = self.length_data # not needed?
                    cap_low = True
                    break
                ##############
                cap = 0.0
                # record cycle time
                cycle_time.append(count*self.time_increment)
                # switch charge to discharge or vice-versa
                charge = not charge
    
                # set limiting current for next cycle, with previous allowable concentrations
                i_lim_cls_t, i_lim_ncls_t = self.limiting_reactant_selector(charge, conc_ox_now_CLS, conc_red_now_CLS,
                                                                            conc_ox_now_NCLS, conc_red_now_NCLS)
                continue
            else:
                pass
            # calculate overpotentials and resulting cell voltage
            losses, n_act, n_mt = self.v_losses(i, charge, concentration_ox_CLS, concentration_red_CLS,
                                                concentration_ox_NCLS, concentration_red_NCLS, i_lim_cls_t, i_lim_ncls_t)
        
            OCV = self.nernst_OCV_full(concentration_ox_CLS, concentration_red_CLS, concentration_ox_NCLS,
                                       concentration_red_NCLS)
    
            cell_V = ZeroDModel._cell_voltage(OCV, losses, charge)
    
            # did it hit voltage limit ?
            if (cell_V >= voltage_cutoff_charge) or (cell_V <= voltage_cutoff_discharge):
                # record capacity here
                cycle_capacity.append(cap)
                #####################
                if (cap < 1.0) and (len(cycle_capacity) > 2):
                    print(str(count) + 'count')
                    print('Simulation stopped, capacity < 1 coulomb')
                    final_count = count
                    # to break out of while loop
                    count = self.length_data # not needed?
                    cap_low = True
                    break
                ##################    
                cap = 0.0
                # record cycle time
                cycle_time.append(count*self.time_increment)
                # switch charge to discharge or viceversa
                charge = not charge
    
                # set limiting current for next cycle
                i_lim_cls_t, i_lim_ncls_t = self.limiting_reactant_selector(charge, conc_ox_now_CLS, conc_red_now_CLS,
                                                                            conc_ox_now_NCLS, conc_red_now_NCLS)
                continue
            else:
                pass

            # calculate SOC (local, in this case) ** for CLS make this a function, can be other file
            soc_CLS, soc_NCLS = ZeroDModel._soc(concentration_ox_CLS, concentration_red_CLS, concentration_ox_NCLS,
                                                concentration_red_NCLS)
            # update capacity
            cap += abs(i*self.time_increment)
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
    ##########################################################################
    ##########################################################################

    def CCCV_experiment(self, voltage_limit_charge: float, voltage_limit_discharge: float, current_cutoff_charge: float,
                        current_cutoff_discharge: float, current: float, charge_first: bool) -> object:

        # to dos: if full CV, then should be appending voltage to overpotential

        # Raise ValueError if user inputs a negative concentration
        if ZeroDModel._negative_concentrations(self.cls_start_c_ox, self.cls_start_c_red,
                                               self.ncls_start_c_ox, self.ncls_start_c_red):
            raise ValueError('Negative concentration detected')
        
        (conc_ox_CLS_profile, conc_red_CLS_profile, conc_ox_NCLS_profile, conc_red_NCLS_profile, cycle_capacity,
         current_profile, cell_V_profile, soc_profile_CLS, ocv_profile, soc_profile_NCLS, cycle_time, act_profile,
         mt_profile, loss_profile) = [],[],[],[],[],[],[],[],[],[],[],[],[],[]
        del_ox,del_red = [],[]
        # set starting concentrations for all species
        # simplify this ?????
        conc_ox_now_CLS = self.cls_start_c_ox
        conc_red_now_CLS = self.cls_start_c_red
        conc_ox_now_NCLS = self.ncls_start_c_ox
        conc_red_now_NCLS = self.ncls_start_c_red

        assert current_cutoff_discharge < 0, "invalid discharge current cutoff"
        #
        charge = charge_first #True
        CC_mode = True
        i_first = True
        CV_only = False
        times = self.times

        count = 0
        cap = 0.0
        cap_low = False
        final_count = self.length_data
        
        i_lim_cls_t, i_lim_ncls_t = self.limiting_reactant_selector(charge, conc_ox_now_CLS, conc_red_now_CLS,
                                                                    conc_ox_now_NCLS, conc_red_now_NCLS)
     
        i = ZeroDModel._current_direction(charge) * current
        ##### check if need to go straight to CV
        # allow option to say 0 current for straight CV?
        if (current >= (i_lim_cls_t*2)) or (current >= (i_lim_ncls_t*2)):
            CC_mode = False
            CV_only = True
            print("Goes straight to CV cycling")
        else:
            losses,n_act,n_mt = self.v_losses(i, charge, conc_ox_now_CLS, conc_red_now_CLS, conc_ox_now_NCLS,
                                              conc_red_now_NCLS, i_lim_cls_t, i_lim_ncls_t)
            ## OCV due to CLS and NCLS
            OCV = self.nernst_OCV_full(conc_ox_now_CLS, conc_red_now_CLS, conc_ox_now_NCLS, conc_red_now_NCLS)
            cell_V = ZeroDModel._cell_voltage(OCV, losses, charge)
            if (cell_V >= voltage_limit_charge) or (cell_V <= voltage_limit_discharge):
                CC_mode = False
                CV_only = True
                print("Has now switched to CV cycling")
            else:
                pass
        ##############
        while count != self.length_data:
            # check if in CC or CV mode
            if CC_mode: 
                # set current
                i = ZeroDModel._current_direction(charge) * current

                # calculate species' concentrations
                (concentration_ox_CLS,
                 concentration_red_CLS,
                 concentration_ox_NCLS,
                 concentration_red_NCLS,
                 delta_ox, delta_red) = self.coulomb_counter(i, conc_ox_now_CLS, conc_red_now_CLS, conc_ox_now_NCLS,
                                                             conc_red_now_NCLS)

                # EDGE CASE where voltage limits never reached i.e straight CC cycling
                if ZeroDModel._negative_concentrations(concentration_ox_CLS, concentration_red_CLS,
                                                       concentration_ox_NCLS, concentration_red_NCLS):
                    # record capacity here
                    cycle_capacity.append(cap)

                    ############ Break out of loop if capacity near zero
                    if (cap < 1.0) and (len(cycle_capacity) > 2):
                        print(str(count) + 'count')
                        print('Simulation stopped, capacity < 1 coulomb')
                        final_count = count
                        #count = self.length_data # not needed?
                        cap_low = True
                        break
                    
                    #####################
                    cap = 0.0
                    # record cycle time
                    cycle_time.append(count*self.time_increment)
                    # switch charge to discharge or viceversa
                    charge = not charge
                    
                    # set limiting current for next cycle
                    i_lim_cls_t, i_lim_ncls_t = self.limiting_reactant_selector(charge, conc_ox_now_CLS,
                                                                                conc_red_now_CLS, conc_ox_now_NCLS,
                                                                                conc_red_now_NCLS)
                    continue
                else:
                    pass
    
                # calculate overpotentials and resulting cell voltage
                losses, n_act, n_mt = self.v_losses(i, charge, concentration_ox_CLS, concentration_red_CLS,
                                                    concentration_ox_NCLS, concentration_red_NCLS, i_lim_cls_t,
                                                    i_lim_ncls_t)
    
                # OCV due to CLS and NCLS
                OCV = self.nernst_OCV_full(concentration_ox_CLS, concentration_red_CLS, concentration_ox_NCLS,
                                           concentration_red_NCLS)
                
                cell_V = ZeroDModel._cell_voltage(OCV, losses, charge)
    
                # calculate SOC (local, not global if there's loss mechanisms)
                soc_CLS, soc_NCLS = ZeroDModel._soc(concentration_ox_CLS, concentration_red_CLS, concentration_ox_NCLS,
                                                    concentration_red_NCLS)
                # update capacity
                cap += abs(i*self.time_increment)
    
                #check if V limit is reached? 
                if (cell_V >= voltage_limit_charge) or (cell_V <= voltage_limit_discharge):
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
                current_profile.append(i) # CC section
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
            else: #now we're in CV mode
    
                # all constant voltage here
                cell_V = voltage_limit_charge if charge else voltage_limit_discharge

                OCV = self.nernst_OCV_full(conc_ox_now_CLS, conc_red_now_CLS, conc_ox_now_NCLS, conc_red_now_NCLS)
    
                data = (cell_V, OCV, charge, conc_ox_now_CLS, conc_red_now_CLS, conc_ox_now_NCLS, conc_red_now_NCLS,
                        i_lim_cls_t, i_lim_ncls_t)
                # adapting the solver's guess to the updated current
                # calculate the first current value for guess
                if i_first:
                    if CV_only: # the case where you're straight CV cycling, and an initial current guess is needed
                        i_guess = ZeroDModel._current_direction(charge) * current #??
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
                i_CV = float(fsolve(self.cv_current_solver, i_guess, args=data)[0])
                i_guess = i_CV

                # check if current is below cutoffs
                if (charge and (i_CV <= current_cutoff_charge)) or (not charge and (i_CV >= current_cutoff_discharge)):
                    # CV part of cycle has now ended, record capacity data

                    cycle_capacity.append(cap)
                    ############ Break out of full simulation if capacity nears zero
                    
                    if (cap < 1.0) and (len(cycle_capacity) > 2):
                        print(str(count) + 'count')
                        print('Simulation stopped, capacity < 1 coulomb')
                        final_count = count
                        #count = self.length_data # not needed?
                        cap_low = True
                        break
                    
                    ############################
                    cap = 0.0
                    # record cycle time
                    cycle_time.append(count*self.time_increment)
                    # switch charge to discharge or viceversa
                    charge = not charge
                    CC_mode = True 
                    i_first = True
    
                    # set limiting current for next cycle
                    i_lim_cls_t, i_lim_ncls_t = self.limiting_reactant_selector(charge, conc_ox_now_CLS,
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
                 delta_ox, delta_red) = self.coulomb_counter(i_CV, conc_ox_now_CLS, conc_red_now_CLS, conc_ox_now_NCLS,
                                                             conc_red_now_NCLS)
    
                # check if any reactant remains
                if ZeroDModel._negative_concentrations(concentration_ox_CLS, concentration_red_CLS,
                                                       concentration_ox_NCLS, concentration_red_NCLS):
                    
                    cycle_capacity.append(cap)
                    ############### Break out of loop if capacity nears zero
                    
                    if (cap < 1.0) and (len(cycle_capacity) > 2):
                        print(str(count) + 'count')
                        print('Simulation stopped, capacity < 1 coulomb')
                        final_count = count
                        #count = self.length_data # not needed?
                        cap_low = True
                        break
                    
                    ###########
                    cap = 0.0
                    # record cycle time
                    cycle_time.append(count*self.time_increment)
                    # switch charge to discharge or viceversa
                    charge = not charge
                    CC_mode = True
                    i_first = True # not sure if this would be needed
    
                    # set limiting current for next cycle
                    i_lim_cls_t, i_lim_ncls_t = self.limiting_reactant_selector(charge, conc_ox_now_CLS,
                                                                                conc_red_now_CLS, conc_ox_now_NCLS,
                                                                                conc_red_now_NCLS)
                    continue
                else:
                    pass
    
                # calculate SOC (local, in this case)
                soc_CLS, soc_NCLS = ZeroDModel._soc(concentration_ox_CLS, concentration_red_CLS,
                                                    concentration_ox_NCLS, concentration_red_NCLS)
                # update capacity
                cap += abs(i_CV*self.time_increment)
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


if __name__ == '__main__':
    print('testing')           
        