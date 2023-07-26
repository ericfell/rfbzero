

from math import log, isclose
import scipy.constants as spc
from scipy.optimize import fsolve #, brentq
#from model_code.zeroD_model_degradations import degradation_mechanism
#from model_code.zeroD_model_crossover import crossover_mechanism
from zeroD_model_degradations import degradation_mechanism
from zeroD_model_crossover import crossover_mechanism

F = spc.physical_constants['Faraday constant'][0]
R = spc.R
# make these parameters at some point?
TEMPERATURE = 298  # kelvin
NERNST_CONST = (R*TEMPERATURE) / F

# concentration cutoff
CONC_CUTOFF = 0.0#1e-17#1e-11  # should be function of timestep coulombs perhaps???


class ZeroDModelSingleVsSingle:
    """
    Model adapted from:

    Modak, S.; Kwabi, D. G. A Zero-Dimensional Model for Electrochemical
    Behavior and Capacity Retention in Organic Flow Cells, Journal of The
    Electrochemical Society, 168, 2021, 080528.
    """

    def __init__(self, geometric_area, resistance, CLS_volume, NCLS_volume,
                 CLS_start_conc_ox, CLS_start_conc_red, NCLS_start_conc_ox,
                 NCLS_start_conc_red, CLS_negolyte, duration, time_increment,
                 standard_E, k_mt, roughness_factor, k_0_CLS, k_0_NCLS,
                 alpha_CLS, alpha_NCLS, mechanism_list=[], mechanism_params={},
                 crossover_list=[], crossover_params={}): # new crssver items
        """Parameters to define the  setup"""

        self.geometric_area = geometric_area   # geometric_area of cell (cm**2)
        self.resistance = resistance           # cell resistance (ohms)
        self.CLS_volume = CLS_volume           # volume of CLS (L)
        self.NCLS_volume = NCLS_volume         # volume of NCLS (L)
        self.CLS_start_conc_ox = CLS_start_conc_ox    # start [ox] in CLS (M)
        self.CLS_start_conc_red = CLS_start_conc_red  # start [red] in CLS (M)
        self.NCLS_start_conc_ox = NCLS_start_conc_ox  # start [ox] in NCLS (M)
        self.NCLS_start_conc_red = NCLS_start_conc_red  # start [ox] in CLS (M)
        self.CLS = True
        self.NCLS = False
        self.CLS_negolyte = True  # True: negolyte = CLS, False: posolyte = CLS
        self.duration = duration  # time of experiment to simulate (sec)
        self.time_increment = time_increment  # time steps for simulation (sec)
        self.standard_E = standard_E  # OCV of cell with equal concentrations?
        self.k_mt = k_mt                    # mass transport coefficient (cm/s)
        self.k_0_CLS = k_0_CLS  # echem rate constant, CLS chemistry (cm/s)
        self.k_0_NCLS = k_0_NCLS  # echem rate constant, NCLS chemistry (cm/s)
        self.alpha_CLS = alpha_CLS      # charge transfer coefficient of CLS
        self.alpha_NCLS = alpha_NCLS    # charge transfer coefficient of NCLS
        self.mechanism_list = mechanism_list  # list of degradation functions
        self.mechanism_params = mechanism_params  # parameters for functions
        self.crossover_list = crossover_list
        self.crossover_params = crossover_params
        self.const_i_ex = F*roughness_factor*self.geometric_area

    def experiment_time(self):
        steps = int((self.duration / self.time_increment)) + 1
        times = [0.0 + x*self.time_increment for x in range(steps)]
        return times

    def current_direction(self, charge: bool) -> int:
        return 1 if charge else -1
    ###############################################################
    # Section: Overpotential functions

    def i_exchange_current(self, k_0, c_red, c_ox, alpha):
        # div by 1000 for conversion from mol/L to mol/cm^3
        i_ex = (self.const_i_ex * k_0 * (c_red**alpha)
                * (c_ox**(1-alpha)) * 0.001)
        return i_ex

    def i_limiting(self, c_lim):
        # div by 1000 for conversion from mol/L to mol/cm^3
        return F*self.k_mt*c_lim*self.geometric_area*0.001

    def limiting_reactant_selecter(self, charge,
                                   conc_ox_now_CLS, conc_red_now_CLS,
                                   conc_ox_now_NCLS, conc_red_now_NCLS):
        if self.CLS_negolyte:  # CLS is negolyte
            if charge:
                i_lim_cls_t = self.i_limiting(conc_ox_now_CLS)
                i_lim_ncls_t = self.i_limiting(conc_red_now_NCLS)
            else:  # discharge first
                i_lim_cls_t = self.i_limiting(conc_red_now_CLS)
                i_lim_ncls_t = self.i_limiting(conc_ox_now_NCLS)

        else:  # CLS is posolyte
            if charge:
                i_lim_cls_t = self.i_limiting(conc_red_now_CLS)
                i_lim_ncls_t = self.i_limiting(conc_ox_now_NCLS)
            else:  # discharge first
                i_lim_cls_t = self.i_limiting(conc_ox_now_CLS)
                i_lim_ncls_t = self.i_limiting(conc_red_now_NCLS)

        return i_lim_cls_t, i_lim_ncls_t

    def n_activation(self, current, i_0_cls, i_0_ncls):
        z_cls = abs(current) / (2*i_0_cls)
        z_ncls = abs(current) / (2*i_0_ncls)
        n_act = NERNST_CONST*(log(z_ncls + ((z_ncls**2) + 1)**0.5)
                              + log(z_cls + ((z_cls**2) + 1)**0.5))
        return n_act

    def n_mass_transport(self, charge, current, c_red_cls, c_ox_cls,
                         c_red_ncls, c_ox_ncls, i_lim_cls, i_lim_ncls):

        assert c_red_cls > 0, "c_red_cls is less than 0"
        assert c_ox_cls > 0, "c_ox_cls is less than 0"
        assert c_red_ncls > 0, "c_red_ncls is less than 0"
        assert c_ox_ncls > 0, "c_ox_ncls is less than 0"

        c_tot_cls = c_red_cls + c_ox_cls
        c_tot_ncls = c_red_ncls + c_ox_ncls

        if self.CLS_negolyte:
            if charge:
                n_mt = NERNST_CONST*log((1 - ((c_tot_cls*current)
                                              / ((c_red_cls*i_lim_cls)
                                                 + (c_ox_cls*current))))
                                        * (1 - ((c_tot_ncls*current)
                                                / ((c_ox_ncls*i_lim_ncls)
                                                   + (c_red_ncls*current)))))
            else:  # discharging
                # keep the abs() for now, but it's only for the strange
                # situation where capacity is balanced and sides are ambiguous
                n_mt = NERNST_CONST*log(abs((1 - ((c_tot_cls*current)
                                                  / ((c_ox_cls*i_lim_cls*-1)
                                                     + (c_red_cls*current))))
                                            * (1 - ((c_tot_ncls*current)
                                                    / ((c_red_ncls*i_lim_ncls*-1)
                                                       + (c_ox_ncls*current))))))
        else:  # CLS is posolyte
            if charge:
                n_mt = NERNST_CONST*log((1 - ((c_tot_cls*current)
                                              / ((c_ox_cls*i_lim_cls)
                                                 + (c_red_cls*current))))
                                        * (1 - ((c_tot_ncls*current)
                                                / ((c_red_ncls*i_lim_ncls)
                                                   + (c_ox_ncls*current)))))
            else:  # discharging
                n_mt = NERNST_CONST*log((1 - ((c_tot_cls*current)
                                              / ((c_red_cls*i_lim_cls*-1)
                                                 + (c_ox_cls*current))))
                                        * (1 - ((c_tot_ncls*current)
                                                / ((c_ox_ncls*i_lim_ncls*-1)
                                                   + (c_red_ncls*current)))))
        return n_mt

    def v_losses(self, current, charge, c_red_cls, c_ox_cls, 
                 c_red_ncls, c_ox_ncls, i_lim_cls, i_lim_ncls):
        
        i_0_cls = self.i_exchange_current(self.k_0_CLS, c_red_cls, c_ox_cls, 
                                          self.alpha_CLS) 
        i_0_ncls = self.i_exchange_current(self.k_0_NCLS, c_red_ncls, c_ox_ncls, 
                                           self.alpha_NCLS)
        #
        n_ohmic = current*self.resistance
        n_act = self.n_activation(current, i_0_cls, i_0_ncls)
        n_mt = self.n_mass_transport(charge, current, c_red_cls, c_ox_cls, 
                                c_red_ncls, c_ox_ncls, i_lim_cls, i_lim_ncls)
        
        if charge:
            n_loss = n_act + n_mt + n_ohmic
        else: # discharge
            n_loss = n_act + n_mt - n_ohmic 
        return n_loss, n_act, n_mt
    ##################################################################
    def nernst_OCV_full(self, conc_ox_CLS, conc_red_CLS, conc_ox_NCLS, conc_red_NCLS):
        
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
    ############################################################# 
    def cell_voltage(self, OCV, losses, charge):
        return OCV + losses if charge else OCV - losses
    #############################################################
    def coulomb_counter(self, current, volume, concentration_0_ox, 
                        concentration_0_red, CLS):

        direction = 1 if self.CLS_negolyte else -1
        delta = ((self.time_increment*current) / (F*volume))*direction

        if CLS: # this is the CLS
            concentration_ox = concentration_0_ox - delta 
            concentration_red = concentration_0_red + delta  
        else: # this is the NCLS
            concentration_ox = concentration_0_ox + delta  
            concentration_red = concentration_0_red - delta 
        
        # now apply possible degradations
        # this assumes symmetric cell, must adapt into above for full cell
        concentration_ox, concentration_red = degradation_mechanism(
                concentration_ox, concentration_red, self.time_increment,
                *self.mechanism_list, **self.mechanism_params)
        
        return concentration_ox, concentration_red    
    ############################################################# 
    def CV_i_solver(self, curr, *data):
        (cell_V, OCV, charge, c_red_cls, c_ox_cls, c_red_ncls, c_ox_ncls, 
         i_lim_cls, i_lim_ncls) = data
        loss_solve,_,_ = self.v_losses(curr, charge, c_red_cls, c_ox_cls, 
                                       c_red_ncls, c_ox_ncls, i_lim_cls, i_lim_ncls)
        
        return cell_V - OCV - loss_solve if charge else cell_V - OCV + loss_solve
    ###########################################################################
    ###########################################################################
    ##### Section: Cycling models
    ###########################################################################
    def CC_experiment(self, voltage_cutoff_charge, voltage_cutoff_discharge, 
                      current, charge_first):

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
        times = self.experiment_time()
        datapoints = len(times) - 1
        print(f"{self.duration} sec of cycling, {self.time_increment} s timesteps, \
              {datapoints} total datapoints")
        count = 0
        cap = 0.0
        cap_low = False
        
        i_lim_cls_t, i_lim_ncls_t = self.limiting_reactant_selecter(charge, 
                                        conc_ox_now_CLS, conc_red_now_CLS, 
                                        conc_ox_now_NCLS, conc_red_now_NCLS)
    
        #assert current <= (i_lim_cls_t*2), "Cannot maintain desired current (i> CLS limiting current)"
        #assert current <= (i_lim_ncls_t*2), "Cannot maintain desired current (i> NCLS limiting current)"
        assert isclose(current, (i_lim_cls_t*2)) or (current < (i_lim_cls_t*2)), \
            "Cannot maintain desired current (i> CLS limiting current)"
        assert isclose(current, (i_lim_ncls_t*2)) or (current < (i_lim_ncls_t*2)), \
            "Cannot maintain desired current (i> NCLS limiting current)"
        
        i = self.current_direction(charge)*current
        
        losses,_,_ = self.v_losses(i, charge, conc_red_now_CLS, conc_ox_now_CLS, 
                                    conc_red_now_NCLS, conc_ox_now_NCLS, 
                                    i_lim_cls_t, i_lim_ncls_t)
        #
        OCV = self.nernst_OCV_full(conc_ox_now_CLS, conc_red_now_CLS, 
                                   conc_ox_now_NCLS, conc_red_now_NCLS)
    
        cell_V = self.cell_voltage(OCV, losses, charge)
        #assert voltage_cutoff_discharge <= cell_V <= voltage_cutoff_charge, "won't cycle, overpotential is outside voltage cutoffs"
        assert isclose(voltage_cutoff_discharge, cell_V) or (voltage_cutoff_discharge < cell_V), \
                "won't cycle, overpotential is outside voltage cutoffs"
        assert isclose(voltage_cutoff_charge, cell_V) or (voltage_cutoff_charge > cell_V), \
                "won't cycle, overpotential is outside voltage cutoffs"
        
        while count != datapoints:
            # set current
            i = self.current_direction(charge)*current
            
            # determine concentrations by coulomb counting, CLS
            concentration_ox_CLS, concentration_red_CLS = self.coulomb_counter(
                i, self.CLS_volume, conc_ox_now_CLS, conc_red_now_CLS, self.CLS)
            # determine concentrations by coulomb counting, NCLS
            concentration_ox_NCLS, concentration_red_NCLS = self.coulomb_counter(
                i, self.NCLS_volume, conc_ox_now_NCLS, conc_red_now_NCLS, self.NCLS)
            
            ################# testing ######################
            # put crossover stuff here
            (concentration_ox_CLS,
             concentration_red_CLS, 
             concentration_ox_NCLS,
             concentration_red_NCLS,
             delta_ox, delta_red) = crossover_mechanism(
                 concentration_ox_CLS, concentration_red_CLS, 
                 concentration_ox_NCLS, concentration_red_NCLS, 
                 self.time_increment, *self.crossover_list, **self.crossover_params)
            ################# end of testing ############################
            
            # EDGE CASE where voltage limits never reached i.e straight CC cycling until concentration runs out
            if ((concentration_ox_CLS < CONC_CUTOFF) or (concentration_red_CLS < CONC_CUTOFF) 
                    or (concentration_ox_NCLS < CONC_CUTOFF) or (concentration_red_NCLS < CONC_CUTOFF)):
                # record capacity here
                cycle_capacity.append(cap)
                ############### Break out of loop if cpacity approaches zero
                if (cap < 1.0) and (len(cycle_capacity) > 2):
                    print(str(count) + 'count')
                    final_count = count
                    count = datapoints
                    cap_low = True
                    break
                ##############
                cap = 0.0
                # record cycle time
                cycle_time.append(count*self.time_increment)
                # switch charge to discharge or viceversa
                charge = not charge
    
                # set limiting current for next cycle
                i_lim_cls_t, i_lim_ncls_t = self.limiting_reactant_selecter(
                    charge, conc_ox_now_CLS, conc_red_now_CLS, conc_ox_now_NCLS, 
                    conc_red_now_NCLS)
                continue
            else:
                pass
            # calculate overpotentials and resulting cell voltage
            losses, n_act, n_mt = self.v_losses(i, charge, concentration_red_CLS, 
                                concentration_ox_CLS, concentration_red_NCLS, 
                                concentration_ox_NCLS, i_lim_cls_t, i_lim_ncls_t)
        
            OCV = self.nernst_OCV_full(concentration_ox_CLS, concentration_red_CLS, 
                                       concentration_ox_NCLS, concentration_red_NCLS)
    
            cell_V = self.cell_voltage(OCV, losses, charge)
    
            # did it hit voltage limit ?
            #if (cell_V >= voltage_cutoff_charge) or (cell_V <= voltage_cutoff_discharge):
            if (isclose(cell_V, voltage_cutoff_charge) or (cell_V > voltage_cutoff_charge)) or (isclose(cell_V, voltage_cutoff_discharge) or (cell_V < voltage_cutoff_discharge)):
                # record capacity here
                cycle_capacity.append(cap)
                #####################
                if (cap < 1.0) and (len(cycle_capacity) > 2):
                    print(str(count) + 'count')
                    final_count = count
                    count = datapoints
                    cap_low = True
                    break
                ##################    
                cap = 0.0
                # record cycle time
                cycle_time.append(count*self.time_increment)
                # switch charge to discharge or viceversa
                charge = not charge
    
                # set limiting current for next cycle
                i_lim_cls_t, i_lim_ncls_t = self.limiting_reactant_selecter(
                    charge, conc_ox_now_CLS, conc_red_now_CLS, 
                    conc_ox_now_NCLS, conc_red_now_NCLS)
                continue
            else:
                pass
    
            # assert these so that an initial overpotential, when switched to charge/discharge, will still let you cycle
            #assert voltage_cutoff_discharge <= (OCV - losses), "won't cycle, overpotential is outside lower voltage cutoff"
            #assert voltage_cutoff_charge >= (OCV + losses), "won't cycle, overpotential is outside upper voltage cutoff"
            assert isclose(voltage_cutoff_discharge, (OCV - losses)) or (voltage_cutoff_discharge < (OCV - losses)), \
                    "won't cycle, overpotential is outside lower voltage cutoff"
            assert isclose(voltage_cutoff_charge, (OCV + losses)) or (voltage_cutoff_charge > (OCV + losses)), \
                    "won't cycle, overpotential is outside upper voltage cutoff"
            
            # calculate SOC (local, in this case) ** for CLS
            soc_CLS = (concentration_red_CLS / (concentration_ox_CLS 
                                                + concentration_red_CLS))*100
            soc_NCLS = (concentration_red_NCLS / (concentration_ox_NCLS 
                                                  + concentration_red_NCLS))*100
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
    def CCCV_experiment(self, voltage_limit_charge, voltage_limit_discharge, 
                        current_cutoff_charge, current_cutoff_discharge, 
                        current, charge_first):

        # to dos: if full CV, then should be appending voltage to overpotential
        # data if wanting to plot that too
        
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
        charge = True 
        CC_mode = True
        i_first = True
        CV_only = False
        times = self.experiment_time()
        datapoints = len(times) - 1
        print(f"{self.duration}s of cycling, {self.time_increment}s timesteps, \
              {datapoints} total datapoints")
        count = 0
        cap = 0.0
        cap_low = False
        
        i_lim_cls_t, i_lim_ncls_t = self.limiting_reactant_selecter(charge, conc_ox_now_CLS,
                                   conc_red_now_CLS, conc_ox_now_NCLS, conc_red_now_NCLS)
     
        i = self.current_direction(charge)*current
        ##### check if need to go straight to CV
        # *** if going back to below line, typo, second thing should be for ncls, not cls
        #if (current >= (i_lim_cls_t*2)) or (current >= (i_lim_cls_t*2)): # allow option to say 0 current for straight CV?
        if (isclose(current, (i_lim_cls_t*2)) or (current > (i_lim_cls_t*2))) or (isclose(current, (i_lim_ncls_t*2)) or (current > (i_lim_ncls_t*2))):
            CC_mode = False
            CV_only = True
            print("Goes straight to CV cycling")
        else:
            losses,n_act,n_mt = self.v_losses(i, charge, conc_red_now_CLS, 
                                   conc_ox_now_CLS, conc_red_now_NCLS, 
                                   conc_ox_now_NCLS, i_lim_cls_t, i_lim_ncls_t)
            ## OCV due to CLS and NCLS
            OCV = self.nernst_OCV_full(conc_ox_now_CLS, conc_red_now_CLS, 
                                       conc_ox_now_NCLS, conc_red_now_NCLS)
            cell_V = self.cell_voltage(OCV, losses, charge)
            #if (cell_V >= voltage_limit_charge) or (cell_V <= voltage_limit_discharge):
            if (isclose(cell_V, voltage_limit_charge) or (cell_V > voltage_limit_charge)) or (isclose(cell_V, voltage_limit_discharge) or (cell_V < voltage_limit_discharge)):
                CC_mode = False
                CV_only = True
                print("Goes straight to CV cycling")
            else:
                pass
        ##############
        while count != datapoints:
            # check if in CC or CV mode
            if CC_mode: 
                # set current
                i = self.current_direction(charge)*current
    
                # determine concentrations by coulomb counting, CLS
                concentration_ox_CLS, concentration_red_CLS = self.coulomb_counter(i, 
                self.CLS_volume, conc_ox_now_CLS, conc_red_now_CLS, self.CLS)
                # determine concentrations by coulomb counting, NCLS
                concentration_ox_NCLS, concentration_red_NCLS = self.coulomb_counter(i, 
                self.NCLS_volume, conc_ox_now_NCLS, conc_red_now_NCLS, self.NCLS)
    
    
                ################# testing ######################
                # put crossover stuff here
                (concentration_ox_CLS,
                 concentration_red_CLS, 
                 concentration_ox_NCLS,
                 concentration_red_NCLS,
                 delta_ox, delta_red) = crossover_mechanism(
                     concentration_ox_CLS, concentration_red_CLS, 
                     concentration_ox_NCLS, concentration_red_NCLS, 
                     self.time_increment, *self.crossover_list, **self.crossover_params)
                ################# end of testing ############################
    
    
    
                # EDGE CASE where voltage limits never reached i.e straight CC cycling
                if ((concentration_ox_CLS < CONC_CUTOFF) or (concentration_red_CLS < CONC_CUTOFF) or
                    (concentration_ox_NCLS < CONC_CUTOFF) or (concentration_red_NCLS < CONC_CUTOFF)):
                    # record capacity here
                    cycle_capacity.append(cap)
                    ############ Break out of loop if capacity near zero
                    
                    if (cap < 1.0) and (len(cycle_capacity) > 2):
                        print(str(count) + 'count')
                        final_count = count
                        count = datapoints
                        cap_low = True
                        break
                    
                    #####################
                    cap = 0.0
                    # record cycle time
                    cycle_time.append(count*self.time_increment)
                    # switch charge to discharge or viceversa
                    charge = not charge
                    
                    # set limiting current for next cycle
                    i_lim_cls_t, i_lim_ncls_t = self.limiting_reactant_selecter(
                        charge, conc_ox_now_CLS, conc_red_now_CLS, conc_ox_now_NCLS, 
                        conc_red_now_NCLS)
                    continue
                else:
                    pass
    
                # calculate overpotentials and resulting cell voltage
                losses, n_act, n_mt = self.v_losses(i, charge, concentration_red_CLS, 
                                   concentration_ox_CLS, concentration_red_NCLS, 
                                   concentration_ox_NCLS, i_lim_cls_t, i_lim_ncls_t)
    
                # OCV due to CLS and NCLS
                OCV = self.nernst_OCV_full(concentration_ox_CLS, concentration_red_CLS, 
                                      concentration_ox_NCLS, concentration_red_NCLS)
                
                cell_V = self.cell_voltage(OCV, losses, charge)
    
                # calculate SOC (local, not global if there's loss mechanisms)
                soc_CLS = (concentration_red_CLS / (concentration_ox_CLS 
                                                    + concentration_red_CLS))*100
                soc_NCLS = (concentration_red_NCLS / (concentration_ox_NCLS 
                                                      + concentration_red_NCLS))*100
                ##############################
                # update capacity
                cap += abs(i*self.time_increment)
    
                #check if V limit is reached? 
                #if (cell_V >= voltage_limit_charge) or (cell_V <= voltage_limit_discharge):
                if (isclose(cell_V, voltage_limit_charge) or (cell_V > voltage_limit_charge)) or (isclose(cell_V, voltage_limit_discharge) or (cell_V < voltage_limit_discharge)):
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
                if charge:
                    cell_V = voltage_limit_charge
                else:
                    cell_V = voltage_limit_discharge
    
                OCV = self.nernst_OCV_full(conc_ox_now_CLS, conc_red_now_CLS, 
                                           conc_ox_now_NCLS, conc_red_now_NCLS)
    
                data = (cell_V, OCV, charge, conc_red_now_CLS, conc_ox_now_CLS, 
                        conc_red_now_NCLS, conc_ox_now_NCLS, i_lim_cls_t, i_lim_ncls_t)
                # adapting the solver's guess to the updated current
                # calculate the first current value for guess
                if i_first:
                    if CV_only: # the case where you're straight CV cycling, and an initial current guess is needed at start
                        i_guess = self.current_direction(charge)*current #??
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
                i_CV = fsolve(self.CV_i_solver, i_guess, args=data)[0]
                i_guess = i_CV
    
                #if ((charge and (i_CV <= current_cutoff_charge)) or 
                #    (not charge and (i_CV >= current_cutoff_discharge))):
                if (charge and (isclose(i_CV, current_cutoff_charge) or (i_CV < current_cutoff_charge))) or (not charge and (isclose(i_CV, current_cutoff_discharge) or (i_CV > current_cutoff_discharge))):
                
                    # record capacity here
                    cycle_capacity.append(cap)
                    ############ Break out of loop if capacity nears zero
                    
                    if (cap < 1.0) and (len(cycle_capacity) > 2):
                        print(str(count) + 'count')
                        final_count = count
                        count = datapoints
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
                    i_lim_cls_t, i_lim_ncls_t = self.limiting_reactant_selecter(
                        charge, conc_ox_now_CLS, conc_red_now_CLS, 
                        conc_ox_now_NCLS, conc_red_now_NCLS) 
                    continue
                else:
                    pass
    
                # determine concentrations by coulomb counting, CLS
                concentration_ox_CLS, concentration_red_CLS = self.coulomb_counter(
                    i_CV, self.CLS_volume, conc_ox_now_CLS, conc_red_now_CLS, self.CLS)
                # determine concentrations by coulomb counting, NCLS
                concentration_ox_NCLS, concentration_red_NCLS = self.coulomb_counter(
                    i_CV, self.NCLS_volume, conc_ox_now_NCLS, conc_red_now_NCLS, self.NCLS)
    
    
                ################# testing ######################
                # put crossover stuff here
                (concentration_ox_CLS,
                 concentration_red_CLS, 
                 concentration_ox_NCLS,
                 concentration_red_NCLS,
                 delta_ox, delta_red) = crossover_mechanism(
                     concentration_ox_CLS, concentration_red_CLS, 
                     concentration_ox_NCLS, concentration_red_NCLS, 
                     self.time_increment, *self.crossover_list, **self.crossover_params)
                ################# end of testing ############################
    
    
                # check if any reactant remains
                if ((concentration_ox_CLS < CONC_CUTOFF) or (concentration_red_CLS < CONC_CUTOFF) or
                    (concentration_ox_NCLS < CONC_CUTOFF) or (concentration_red_NCLS < CONC_CUTOFF)):
                    
                    cycle_capacity.append(cap)
                    ############### Break out of loop if capacity nears zero
                    
                    if (cap < 1.0) and (len(cycle_capacity) > 2):
                        print(str(count) + 'count')
                        final_count = count
                        count = datapoints
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
                    i_lim_cls_t, i_lim_ncls_t = self.limiting_reactant_selecter(
                        charge, conc_ox_now_CLS, conc_red_now_CLS, conc_ox_now_NCLS, 
                        conc_red_now_NCLS)         
                    continue
                else:
                    pass
    
                # calculate SOC (local, in this case)
                soc_CLS = (concentration_red_CLS / (concentration_ox_CLS 
                                                    + concentration_red_CLS))*100
                soc_NCLS = (concentration_red_NCLS / (concentration_ox_NCLS 
                                                      + concentration_red_NCLS))*100
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
        