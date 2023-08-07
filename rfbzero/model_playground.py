from zeroD_model_1e_vs_1e import ZeroDModel as single_e
from zeroD_model_degradations import chemical_degradation
import matplotlib.pyplot as plt
from zeroD_model_degradations import ChemicalDegradation, AutoOxidation

CLS_start_conc_ox = 0.01
CLS_start_conc_red = 0.01
NCLS_start_conc_ox = 0.01
NCLS_start_conc_red = 0.01
area = 5.0
CLS_vol = 0.01
NCLS_vol = 0.050
#CLS_nego = False
t_step = 0.01
E_redox = 1.0
rough = 26
#
voltage_limit_charge = 1.4#0.4
voltage_limit_discharge = 0.6#-0.4
current = 0.3
kmt = 0.8

resistance = 1.0
k_species =  2.2e-3
duration = 5000

#mechanism_list = [chemical_degradation] #<< this is imported, probably should change
mechanism_params = {'red_degrade': ['red', 1, 1e-5]}



# for crossover
membrane_thickness = 183/10000 # cm, nafion 117
membrane_constant = area / membrane_thickness
p_ox = 1.0e-6 # cm^2/s
p_red = 1.0e-6 # cm^2/s

crossover_params = [membrane_constant, p_ox, p_red]

## testing of abstract method classes

#test_fade = ChemicalDegradation(rate_order=1, rate=9e-5, species='red')
test_fade = AutoOxidation(rate=9e-5)
mechanism_list = test_fade
#print(test_fade)
###############################

# setup cycling procedure
setup = single_e(area, resistance, CLS_vol, NCLS_vol, CLS_start_conc_ox, CLS_start_conc_red, NCLS_start_conc_ox,
                 NCLS_start_conc_red, duration, t_step, E_redox, kmt, rough, k_species, k_species, 0.5, 0.5,True,
                 mechanism_list=mechanism_list, mechanism_params=mechanism_params, crossover_params=None)#crossover_params)

(current_profile, conc_ox_CLS_profile, conc_red_CLS_profile, conc_ox_NCLS_profile, conc_red_NCLS_profile,
 cell_V_profile, soc_profile_CLS, soc_profile_NCLS, ocv_profile, cycle_capacity, cycle_time, times,  act_profile,
 mt_profile, loss_profile, del_ox, del_red,
) = setup.cc_experiment(voltage_limit_charge, voltage_limit_discharge, current, True)
#) = setup.CCCV_experiment(voltage_limit_charge, voltage_limit_discharge, 0.005, -0.005,current, True)

#print(len(current_profile))
print(cycle_capacity[:5])


def structure_data(x, y):
    time_charge = x[::2]
    time_discharge = x[1::2]
    cap_charge = y[::2]
    cap_discharge = y[1::2]

    if len(time_charge) < len(cap_charge):
        cap_charge = cap_charge[:-1]
    if len(time_discharge) < len(cap_discharge):
        cap_discharge = cap_discharge[:-1]
    return time_charge, time_discharge, cap_charge, cap_discharge

time_charge, time_discharge, cap_charge, cap_discharge = structure_data(
        cycle_time, cycle_capacity)

g = 0
fig,ax = plt.subplots(nrows=5,ncols=1,sharex=True)
ax[0].plot(time_discharge, cap_discharge, 'bo--')
ax[0].plot(time_charge, cap_charge, 'ro--')
ax[1].plot(times[g:], current_profile, 'r')
ax[2].plot(times[g:], cell_V_profile, 'b')
ax[3].plot(times[g:], ocv_profile, 'k')
ax[4].plot(times[g:], loss_profile, 'k')

ax[4].plot(times[g:], act_profile, 'r')
ax[4].plot(times[g:], mt_profile, 'b')

ax[0].set_ylabel('Capacity')
ax[1].set_ylabel('Current')
ax[2].set_ylabel('Voltage')
ax[3].set_ylabel('OCV')
ax[4].set_ylabel('Overpotential')
plt.show()

# plot out current and voltage cycyles