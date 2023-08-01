

from math import log
import scipy.constants as spc
from scipy.optimize import fsolve
from zeroD_model_degradations import degradation_mechanism
from zeroD_model_crossover import crossover_mechanism

F = spc.physical_constants['Faraday constant'][0]  # or just input 96485.3321 ? same for R
R = spc.R
# make these parameters at some point?
TEMPERATURE = 298  # kelvin
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
    CLS_volume : float
        Volume of capacity-limiting side (CLS) reservoir (L).
    NCLS_volume : float
        Volume of non-capacity-limiting side (NCLS) reservoir (L).
    CLS_start_conc_ox : float
        CLS initial concentration of oxidized species (M).
    CLS_start_conc_red : float
        CLS initial concentration of reduced species (M).
    NCLS_start_conc_ox : float
        NCLS initial concentration of oxidized species (M).
    NCLS_start_conc_red : float
        NCLS initial concentration of reduced species (M).
    duration : int
        Amount of real time to simulate (s).
    time_increment : float
        Simulation time step (s).
    standard_E : float
        Cell voltage (formal potentials E_+ - E_-) (V).
        If voltage > 0 then it's a Full cell.
        If voltage = 0 then it's a Symmetric cell.
    k_mt : float
        Mass transport coefficient (cm/s).
    roughness_factor : float
        Roughness factor, dimensionless.
        Surface area divided by geometric surface area.
    k_0_CLS : float
        Electrochemical rate constant, CLS redox couple (cm/s).
    k_0_NCLS : float
        Electrochemical rate constant, NCLS redox couple (cm/s).
    alpha_CLS : float
        Charge transfer coefficient of CLS redox couple, dimensionless.
    alpha_NCLS : float
        Charge transfer coefficient of NCLS redox couple, dimensionless.
    CLS_negolyte : bool
        If True, negolyte is the CLS.
    mechanism_list : list, optional
        List of degradation mechanisms to include in simulation.
    mechanism_params : dict, optional
        Parameters for mechanisms specified in `mechanism_list`.
    crossover_list : list, optional
        List of crossover mechanisms to include in simulation.
    crossover_params : dict, optional
        Parameters for mechanisms specified in `crossover_list`.


    Notes
    -----
    All equations are adapted from [1]. If ZeroDModel has been
    significant to your research please cite the paper.

    [1] Modak, S.; Kwabi, D. G. A Zero-Dimensional Model for Electrochemical
    Behavior and Capacity Retention in Organic Flow Cells, Journal of The
    Electrochemical Society, 168, 2021, 080528.
    """

    @staticmethod
    def current_direction(charge: bool) -> int:
        """Makes current positive for charge, negative for discharge"""
        return 1 if charge else -1

    @staticmethod
    def SOC(conc_ox_CLS: float, conc_red_CLS: float, conc_ox_NCLS: float, conc_red_NCLS: float) -> tuple[float, float]:
        """Calculate state-of-charge in each reservoir"""
        soc_CLS = (conc_red_CLS / (conc_ox_CLS + conc_red_CLS)) * 100
        # this could be defined differently i.e., cell vs reservoir definition of SOC
        soc_NCLS = (conc_red_NCLS / (conc_ox_NCLS + conc_red_NCLS)) * 100
        return soc_CLS, soc_NCLS

    @staticmethod
    def negative_concentrations(conc_ox_CLS: float, conc_red_CLS: float, conc_ox_NCLS: float, conc_red_NCLS: float) -> bool:
        """Return True if any concentration is negative"""
        return any(x < 0.0 for x in [conc_ox_CLS, conc_red_CLS, conc_ox_NCLS, conc_red_NCLS])

    def __init__(self, geometric_area: float, resistance: float, CLS_volume: float, NCLS_volume: float,
                 CLS_start_conc_ox: float, CLS_start_conc_red: float, NCLS_start_conc_ox: float,
                 NCLS_start_conc_red: float, duration: int, time_increment: float,
                 standard_E: float, k_mt: float, roughness_factor: float, k_0_CLS: float, k_0_NCLS: float,
                 alpha_CLS: float, alpha_NCLS: float, CLS_negolyte: bool = True, mechanism_list: list = None, mechanism_params: dict = None,
                 crossover_list: list = None, crossover_params: dict = None) -> None:
        """Inits ZeroDModel"""

        self.geometric_area = geometric_area
        self.resistance = resistance
        self.CLS_volume = CLS_volume
        self.NCLS_volume = NCLS_volume
        self.CLS_start_conc_ox = CLS_start_conc_ox
        self.CLS_start_conc_red = CLS_start_conc_red
        self.NCLS_start_conc_ox = NCLS_start_conc_ox
        self.NCLS_start_conc_red = NCLS_start_conc_red
        self.duration = duration
        self.time_increment = time_increment
        self.standard_E = standard_E
        self.k_mt = k_mt
        self.k_0_CLS = k_0_CLS
        self.k_0_NCLS = k_0_NCLS
        self.alpha_CLS = alpha_CLS
        self.alpha_NCLS = alpha_NCLS
        self.CLS_negolyte = CLS_negolyte
        self.mechanism_list = mechanism_list
        self.mechanism_params = mechanism_params
        self.crossover_list = crossover_list
        self.crossover_params = crossover_params
        self.const_i_ex = F*roughness_factor*self.geometric_area
        self.length_data = int(self.duration / self.time_increment)
        self.times = [x*self.time_increment for x in range(1, self.length_data + 1)]
        print(f"{self.duration} s of cycling, {self.time_increment} s timesteps,"
              f"{self.length_data} total datapoints")

    # option for just measuring capacity over time, doesnt need to make all arrays?
    # should all four conc be a list passed into everything?

    def i_exchange_current(self, c_ox_CLS: float, c_red_CLS: float,
                           c_ox_NCLS: float, c_red_NCLS: float) -> tuple[float, float]:
        """
        Calculates exchange current of redox couples in the CLS and NCLS.

        Parameters
        ----------
        c_ox_CLS: float
            Concentration of oxidized species in CLS
             at a given timestep (M).
        c_red_CLS: float
            Concentration of reduced species in CLS
             at a given timestep (M).
        c_ox_NCLS: float
            Concentration of oxidized species in NCLS
             at a given timestep (M).
        c_red_NCLS: float
            Concentration of reduced species in NCLS
             at a given timestep (M).

        Returns
        -------
        i_ex_CLS : float
            Exchange current of CLS redox couple
            at a given timestep (A).
        i_ex_NCLS : float
            Exchange current of NCLS redox couple
            at a given timestep (A)

        """
        # division by 1000 for conversion from mol/L to mol/cm^3
        i_ex_CLS = (self.const_i_ex * self.k_0_CLS * (c_red_CLS**self.alpha_CLS) * (c_ox_CLS**(1 - self.alpha_CLS)) * 0.001)
        i_ex_NCLS = (self.const_i_ex * self.k_0_NCLS * (c_red_NCLS**self.alpha_NCLS) * (c_ox_NCLS**(1 - self.alpha_NCLS)) * 0.001)
        return i_ex_CLS, i_ex_NCLS

    def i_limiting(self, c_lim: float) -> float:
        """Calculates limiting current for a single reservoir.
        This is equation 6 of [1].
        """
        # div by 1000 for conversion from mol/L to mol/cm^3
        # requires n electrons param
        return F * self.k_mt * c_lim * self.geometric_area * 0.001

    def limiting_reactant_selector(self, charge: bool, conc_ox_now_CLS: float, conc_red_now_CLS: float,
                                   conc_ox_now_NCLS: float, conc_red_now_NCLS: float) -> tuple[float, float]:
        """Selects limiting concentration and calculates limiting current for CLS and NCLS."""
        if (self.CLS_negolyte and charge) or (not self.CLS_negolyte and not charge):
            i_lim_cls_t = self.i_limiting(conc_ox_now_CLS)
            i_lim_ncls_t = self.i_limiting(conc_red_now_NCLS)
        else:
            i_lim_cls_t = self.i_limiting(conc_red_now_CLS)
            i_lim_ncls_t = self.i_limiting(conc_ox_now_NCLS)

        return i_lim_cls_t, i_lim_ncls_t

    def n_activation(self, current: float, i_0_cls: float, i_0_ncls: float) -> float:
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

        z_cls = abs(current) / (2*i_0_cls)
        z_ncls = abs(current) / (2*i_0_ncls)
        n_act = NERNST_CONST*(log(z_ncls + ((z_ncls**2) + 1)**0.5)
                              + log(z_cls + ((z_cls**2) + 1)**0.5))
        return n_act

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

        assert c_red_cls > 0, "c_red_cls is less than 0"
        assert c_ox_cls > 0, "c_ox_cls is less than 0"
        assert c_red_ncls > 0, "c_red_ncls is less than 0"
        assert c_ox_ncls > 0, "c_ox_ncls is less than 0"

        c_tot_cls = c_red_cls + c_ox_cls
        c_tot_ncls = c_red_ncls + c_ox_ncls

        current = abs(current)

        if (self.CLS_negolyte and charge) or (not self.CLS_negolyte and not charge):
            n_mt = NERNST_CONST * log((1 - ((c_tot_cls * current)
                                            / ((c_red_cls * i_lim_cls)
                                               + (c_ox_cls * current))))
                                      * (1 - ((c_tot_ncls * current)
                                              / ((c_ox_ncls * i_lim_ncls)
                                                 + (c_red_ncls * current)))))
        else:
            n_mt = NERNST_CONST * log(((1 - ((c_tot_cls * current)
                                             / ((c_ox_cls * i_lim_cls)
                                                + (c_red_cls * current))))
                                       * (1 - ((c_tot_ncls * current)
                                               / ((c_red_ncls * i_lim_ncls)
                                                  + (c_ox_ncls * current))))))
        return n_mt

    def v_losses(self, current: float, charge: bool, c_ox_cls: float, c_red_cls: float,
                 c_ox_ncls: float, c_red_ncls: float, i_lim_cls: float, i_lim_ncls: float) -> tuple[float, float, float]:
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

    def nernst_OCV_full(self, conc_ox_CLS: float, conc_red_CLS: float, conc_ox_NCLS: float, conc_red_NCLS: float) -> float:
        """
        This is equivalent to equation 3 of [1].

        Parameters
        ----------
        conc_ox_CLS
        conc_red_CLS
        conc_ox_NCLS
        conc_red_NCLS

        Returns
        -------

        """
        # will need n_electrons input
        assert conc_red_CLS > 0, "CLS [red] is less than 0, nernst_OCV_full call"
        assert conc_ox_CLS > 0, "CLS [ox] is less than 0, nernst_OCV_full call"
        assert conc_red_NCLS > 0, "NCLS [red] is less than 0, nernst_OCV_full call"
        assert conc_ox_NCLS > 0, "NCLS [ox] is less than 0, nernst_OCV_full call"
     
        if self.CLS_negolyte:
            OCV = (self.standard_E 
                   + (NERNST_CONST*log(conc_red_CLS / conc_ox_CLS)) 
                   + (NERNST_CONST*log(conc_ox_NCLS / conc_red_NCLS)))
        else: # CLS is posolyte
            OCV = (self.standard_E 
                   - (NERNST_CONST*log(conc_red_CLS / conc_ox_CLS)) 
                   - (NERNST_CONST*log(conc_ox_NCLS / conc_red_NCLS)))
        return OCV

    def cell_voltage(self, OCV: float, losses: float, charge: bool) -> float:
        """If charging, add overpotentials to OCV, else subtract them."""
        return OCV + losses if charge else OCV - losses

    def coulomb_counter(self, current: float, conc_ox_CLS: float, conc_red_CLS: float,
                        conc_ox_NCLS: float, conc_red_NCLS: float) -> tuple[float, float, float, float, float, float]:

        direction = 1 if self.CLS_negolyte else -1
        delta_CLS = ((self.time_increment * current) / (F * self.CLS_volume)) * direction
        delta_NCLS = ((self.time_increment * current) / (F * self.NCLS_volume)) * direction

        # update CLS concentrations
        conc_ox_CLS = conc_ox_CLS - delta_CLS
        conc_red_CLS = conc_red_CLS + delta_CLS
        # update NCLS concentrations
        conc_ox_NCLS = conc_ox_NCLS + delta_NCLS
        conc_red_NCLS = conc_red_NCLS - delta_NCLS

        # for no crossover situation
        delta_ox = 0.0
        delta_red = 0.0

        # allow for degradations and/or crossover
        if (self.mechanism_list is None) and (self.crossover_list is None):  # no degrade/ no crossover
            pass

        elif (self.mechanism_list is not None) and (self.crossover_list is None): # degrade/ no crossover
            # possible CLS degradation
            conc_ox_CLS, conc_red_CLS = degradation_mechanism(conc_ox_CLS, conc_red_CLS, self.time_increment,
                *self.mechanism_list, **self.mechanism_params)
            # possible NCLS degradation
            conc_ox_NCLS, conc_red_NCLS = degradation_mechanism(conc_ox_NCLS, conc_red_NCLS, self.time_increment,
                                                              *self.mechanism_list, **self.mechanism_params)

        elif (self.mechanism_list is None) and (self.crossover_list is not None): # no degrade/ crossover
            (conc_ox_CLS, conc_red_CLS, conc_ox_NCLS, conc_red_NCLS, delta_ox,
             delta_red) = crossover_mechanism(conc_ox_CLS, conc_red_CLS, conc_ox_NCLS, conc_red_NCLS,
                                              self.time_increment, *self.crossover_list, **self.crossover_params)

        else: # degrade AND crossover
            # CLS degradation
            conc_ox_CLS, conc_red_CLS = degradation_mechanism(conc_ox_CLS, conc_red_CLS, self.time_increment,
                                                              *self.mechanism_list, **self.mechanism_params)
            # NCLS degradation
            conc_ox_NCLS, conc_red_NCLS = degradation_mechanism(conc_ox_NCLS, conc_red_NCLS, self.time_increment,
                                                                *self.mechanism_list, **self.mechanism_params)

            # crossover
            (conc_ox_CLS, conc_red_CLS, conc_ox_NCLS, conc_red_NCLS, delta_ox,
             delta_red) = crossover_mechanism(conc_ox_CLS, conc_red_CLS, conc_ox_NCLS, conc_red_NCLS,
                                              self.time_increment, *self.crossover_list, **self.crossover_params)

        return conc_ox_CLS, conc_red_CLS, conc_ox_NCLS, conc_red_NCLS, delta_ox, delta_red


    # is below proper *args unpacking naming style?
    def CV_i_solver(self, current: float, *data: float) -> float:
        (cell_V, OCV, charge, c_ox_cls, c_red_cls, c_ox_ncls, c_red_ncls,
         i_lim_cls, i_lim_ncls) = data
        # curr has sign but v_losses makes it always positive
        loss_solve,_,_ = self.v_losses(current, charge, c_ox_cls, c_red_cls,
                                       c_ox_ncls, c_red_ncls, i_lim_cls, i_lim_ncls)
        # returns what solver will try to minimize
        return cell_V - OCV - loss_solve if charge else cell_V - OCV + loss_solve


    ##### Section: Cycling models
    ###########################################################################

    def CC_experiment(self, voltage_cutoff_charge: float, voltage_cutoff_discharge: float,
                      current: float, charge_first: bool) -> object: # edit the return type soon

        (current_profile, conc_ox_CLS_profile, conc_red_CLS_profile, cell_V_profile, 
         soc_profile_CLS, ocv_profile, conc_ox_NCLS_profile, conc_red_NCLS_profile, 
         soc_profile_NCLS, cycle_capacity, cycle_time, act_profile, mt_profile, 
         loss_profile) = [],[],[],[],[],[],[],[],[],[],[],[],[],[]
        # testing
        del_ox,del_red = [],[]
        # set starting concentrations for all species
        conc_ox_now_CLS = self.CLS_start_conc_ox
        conc_red_now_CLS = self.CLS_start_conc_red
        conc_ox_now_NCLS = self.NCLS_start_conc_ox
        conc_red_now_NCLS = self.NCLS_start_conc_red
        # 
        charge = charge_first #??
        times = self.times

        count = 0
        cap = 0.0
        cap_low = False
        # initialized in case simulation has to stop due to no more capacity
        final_count = self.length_data
        
        i_lim_cls_t, i_lim_ncls_t = self.limiting_reactant_selector(charge,
                                                                    conc_ox_now_CLS, conc_red_now_CLS,
                                                                    conc_ox_now_NCLS, conc_red_now_NCLS)
        # why the *2 ??
        assert current < (i_lim_cls_t*2), "Cannot maintain desired current (i> CLS limiting current)"
        assert current < (i_lim_ncls_t*2), "Cannot maintain desired current (i> NCLS limiting current)"

        # assign + current to charge, - current to discharge
        i = ZeroDModel.current_direction(charge) * current
        
        losses,_,_ = self.v_losses(i, charge, conc_ox_now_CLS, conc_red_now_CLS,
                                    conc_ox_now_NCLS, conc_red_now_NCLS,
                                    i_lim_cls_t, i_lim_ncls_t)
        #
        OCV = self.nernst_OCV_full(conc_ox_now_CLS, conc_red_now_CLS, 
                                   conc_ox_now_NCLS, conc_red_now_NCLS)
    
        cell_V = self.cell_voltage(OCV, losses, charge)
        assert voltage_cutoff_discharge <= cell_V <= voltage_cutoff_charge, "won't cycle, overpotential is outside voltage cutoffs"
        
        while count != self.length_data:
            # set current
            i = ZeroDModel.current_direction(charge) * current

            # need to do this for CCCV method below too
            (concentration_ox_CLS,
             concentration_red_CLS,
             concentration_ox_NCLS,
             concentration_red_NCLS,
             delta_ox, delta_red) = self.coulomb_counter(i, conc_ox_now_CLS,
                                                         conc_red_now_CLS, conc_ox_now_NCLS, conc_red_now_NCLS)

            # EDGE CASE where voltage limits never reached i.e straight CC cycling until concentration runs out
            if ZeroDModel.negative_concentrations(concentration_ox_CLS, concentration_red_CLS,
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
                i_lim_cls_t, i_lim_ncls_t = self.limiting_reactant_selector(
                    charge, conc_ox_now_CLS, conc_red_now_CLS, conc_ox_now_NCLS, 
                    conc_red_now_NCLS)
                continue
            else:
                pass
            # calculate overpotentials and resulting cell voltage
            losses, n_act, n_mt = self.v_losses(i, charge, concentration_ox_CLS,
                                concentration_red_CLS, concentration_ox_NCLS,
                                concentration_red_NCLS, i_lim_cls_t, i_lim_ncls_t)
        
            OCV = self.nernst_OCV_full(concentration_ox_CLS, concentration_red_CLS, 
                                       concentration_ox_NCLS, concentration_red_NCLS)
    
            cell_V = self.cell_voltage(OCV, losses, charge)
    
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
                i_lim_cls_t, i_lim_ncls_t = self.limiting_reactant_selector(
                    charge, conc_ox_now_CLS, conc_red_now_CLS, 
                    conc_ox_now_NCLS, conc_red_now_NCLS)
                continue
            else:
                pass

            # calculate SOC (local, in this case) ** for CLS make this a function, can be other file
            soc_CLS, soc_NCLS = ZeroDModel.SOC(concentration_ox_CLS, concentration_red_CLS,
                                               concentration_ox_NCLS, concentration_red_NCLS)
            # update capacity
            cap += abs(i*self.time_increment)
            # update concentrations
            conc_ox_now_CLS = concentration_ox_CLS
            conc_red_now_CLS = concentration_red_CLS
            conc_ox_now_NCLS = concentration_ox_NCLS
            conc_red_now_NCLS = concentration_red_NCLS
            #
            current_profile.append(i)
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
        # adjusts time points if capacity decreased past set point
        if cap_low:
            times = times[:final_count + 1]
        return (current_profile, conc_ox_CLS_profile, conc_red_CLS_profile, 
                conc_ox_NCLS_profile, conc_red_NCLS_profile, cell_V_profile, 
                soc_profile_CLS, soc_profile_NCLS, ocv_profile, cycle_capacity, 
                cycle_time, times, act_profile, mt_profile, loss_profile, del_ox, del_red)
    ##########################################################################
    ##########################################################################

    def CCCV_experiment(self, voltage_limit_charge: float, voltage_limit_discharge: float,
                        current_cutoff_charge: float, current_cutoff_discharge: float,
                        current: float, charge_first: bool) -> object:

        # to dos: if full CV, then should be appending voltage to overpotential
        
        (conc_ox_CLS_profile, conc_red_CLS_profile, conc_ox_NCLS_profile, 
         conc_red_NCLS_profile, cycle_capacity, current_profile, cell_V_profile, 
         soc_profile_CLS, ocv_profile, soc_profile_NCLS, cycle_time, act_profile, 
         mt_profile, loss_profile) = [],[],[],[],[],[],[],[],[],[],[],[],[],[]
        del_ox,del_red = [],[]
        # set starting concentrations for all species
        # simplify this ?????
        conc_ox_now_CLS = self.CLS_start_conc_ox
        conc_red_now_CLS = self.CLS_start_conc_red
        conc_ox_now_NCLS = self.NCLS_start_conc_ox
        conc_red_now_NCLS = self.NCLS_start_conc_red
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
        
        i_lim_cls_t, i_lim_ncls_t = self.limiting_reactant_selector(charge, conc_ox_now_CLS,
                                                                    conc_red_now_CLS, conc_ox_now_NCLS, conc_red_now_NCLS)
     
        i = ZeroDModel.current_direction(charge) * current
        ##### check if need to go straight to CV
        # allow option to say 0 current for straight CV?
        if (current >= (i_lim_cls_t*2)) or (current >= (i_lim_ncls_t*2)):
            CC_mode = False
            CV_only = True
            print("Goes straight to CV cycling")
        else:
            losses,n_act,n_mt = self.v_losses(i, charge, conc_ox_now_CLS,
                                   conc_red_now_CLS, conc_ox_now_NCLS,
                                   conc_red_now_NCLS, i_lim_cls_t, i_lim_ncls_t)
            ## OCV due to CLS and NCLS
            OCV = self.nernst_OCV_full(conc_ox_now_CLS, conc_red_now_CLS, 
                                       conc_ox_now_NCLS, conc_red_now_NCLS)
            cell_V = self.cell_voltage(OCV, losses, charge)
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
                i = ZeroDModel.current_direction(charge) * current

                # calculate species' concentrations
                (concentration_ox_CLS,
                 concentration_red_CLS,
                 concentration_ox_NCLS,
                 concentration_red_NCLS,
                 delta_ox, delta_red) = self.coulomb_counter(i, conc_ox_now_CLS,
                                                             conc_red_now_CLS, conc_ox_now_NCLS, conc_red_now_NCLS)

                # EDGE CASE where voltage limits never reached i.e straight CC cycling
                if ZeroDModel.negative_concentrations(concentration_ox_CLS, concentration_red_CLS,
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
                    i_lim_cls_t, i_lim_ncls_t = self.limiting_reactant_selector(
                        charge, conc_ox_now_CLS, conc_red_now_CLS, conc_ox_now_NCLS, 
                        conc_red_now_NCLS)
                    continue
                else:
                    pass
    
                # calculate overpotentials and resulting cell voltage
                losses, n_act, n_mt = self.v_losses(i, charge, concentration_ox_CLS,
                                   concentration_red_CLS, concentration_ox_NCLS,
                                   concentration_red_NCLS, i_lim_cls_t, i_lim_ncls_t)
    
                # OCV due to CLS and NCLS
                OCV = self.nernst_OCV_full(concentration_ox_CLS, concentration_red_CLS, 
                                      concentration_ox_NCLS, concentration_red_NCLS)
                
                cell_V = self.cell_voltage(OCV, losses, charge)
    
                # calculate SOC (local, not global if there's loss mechanisms)
                soc_CLS, soc_NCLS = ZeroDModel.SOC(concentration_ox_CLS, concentration_red_CLS,
                                                   concentration_ox_NCLS, concentration_red_NCLS)
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

                OCV = self.nernst_OCV_full(conc_ox_now_CLS, conc_red_now_CLS, 
                                           conc_ox_now_NCLS, conc_red_now_NCLS)
    
                data = (cell_V, OCV, charge, conc_ox_now_CLS, conc_red_now_CLS,
                        conc_ox_now_NCLS, conc_red_now_NCLS, i_lim_cls_t, i_lim_ncls_t)
                # adapting the solver's guess to the updated current
                # calculate the first current value for guess
                if i_first:
                    if CV_only: # the case where you're straight CV cycling, and an initial current guess is needed
                        i_guess = ZeroDModel.current_direction(charge) * current #??
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
                    max_bracket = self.current_direction(charge)*current
                else:
                    max_bracket = current_cutoff_discharge*-1
                    min_bracket = self.current_direction(charge)*current
                
                if i_first:
                    if charge:
                        max_bracket = 10
                    else:
                        min_bracket = -10
                
                i_CV = brentq(self.CV_i_solver, min_bracket, max_bracket, args=data)
                '''
                ######
                ######### end testing ############
                # consider different solver that ensure max allowable current upon switch to CV ?
                i_CV = float(fsolve(self.CV_i_solver, i_guess, args=data)[0])
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
                    i_lim_cls_t, i_lim_ncls_t = self.limiting_reactant_selector(
                        charge, conc_ox_now_CLS, conc_red_now_CLS, 
                        conc_ox_now_NCLS, conc_red_now_NCLS) 
                    continue
                else:
                    pass

                # testing new. should i_CV be i_guesss?
                (concentration_ox_CLS,
                 concentration_red_CLS,
                 concentration_ox_NCLS,
                 concentration_red_NCLS,
                 delta_ox, delta_red) = self.coulomb_counter(i_CV, conc_ox_now_CLS,
                                                             conc_red_now_CLS, conc_ox_now_NCLS, conc_red_now_NCLS)
    
                # check if any reactant remains
                if ZeroDModel.negative_concentrations(concentration_ox_CLS, concentration_red_CLS,
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
                    i_lim_cls_t, i_lim_ncls_t = self.limiting_reactant_selector(
                        charge, conc_ox_now_CLS, conc_red_now_CLS, conc_ox_now_NCLS, 
                        conc_red_now_NCLS)         
                    continue
                else:
                    pass
    
                # calculate SOC (local, in this case)
                soc_CLS, soc_NCLS = ZeroDModel.SOC(concentration_ox_CLS, concentration_red_CLS,
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
        return (current_profile, conc_ox_CLS_profile, conc_red_CLS_profile, 
                conc_ox_NCLS_profile, conc_red_NCLS_profile, cell_V_profile, 
                soc_profile_CLS, soc_profile_NCLS, ocv_profile, cycle_capacity, 
                cycle_time, times, act_profile, mt_profile, loss_profile, del_ox, del_red)
##############################################################################
if __name__ == '__main__':
    print('testing')           
        