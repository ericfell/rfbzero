from zeroD_model_1e_vs_1e import ZeroDModelSingleVsSingle as single_e

CLS_start_conc_ox = 0.01
CLS_start_conc_red = 0.01
NCLS_start_conc_ox = 0.01
NCLS_start_conc_red = 0.01
area = 5
CLS_vol = 0.020
NCLS_vol = 0.040
CLS_nego = False
#duration = 1e4
t_step = 0.1
rough = 26
#
voltage_limit_charge = 0.4
voltage_limit_discharge = -0.4
current = 0.2
kmt = 0.8

resistance = 0.5
k_species =  2.2e-3
duration = 1000

setup = single_e(area, resistance, CLS_vol, NCLS_vol, CLS_start_conc_ox, CLS_start_conc_red, NCLS_start_conc_ox,
                 NCLS_start_conc_red, CLS_nego, duration, t_step, 0.0, kmt, rough, k_species, k_species, 0.5, 0.5)

(current_profile, conc_ox_CLS_profile, conc_red_CLS_profile, conc_ox_NCLS_profile, conc_red_NCLS_profile,
 cell_V_profile, soc_profile_CLS, soc_profile_NCLS, ocv_profile, cycle_capacity, cycle_time, times,  act_profile,
 mt_profile, loss_profile, del_ox, del_red) = setup.CCCV_experiment(voltage_limit_charge, voltage_limit_discharge, 0.005, -0.005,
                                                   current, True)

print(len(current_profile))
print(cycle_capacity)