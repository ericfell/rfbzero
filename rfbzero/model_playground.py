from zeroD_model_1e_vs_1e import ZeroDModel as battery
import matplotlib.pyplot as plt
from zeroD_model_degradations import ChemicalDegradation, AutoOxidation, AutoReduction, MultiDegradationMechanism
from zeroD_model_crossover import Crossover
from cycle_protocol import ConstantCurrent, ConstantCurrentConstantVoltage

CLS_start_conc_ox = 0.01
CLS_start_conc_red = 0.01
NCLS_start_conc_ox = 0.01
NCLS_start_conc_red = 0.01
area = 5.0
CLS_vol = 0.01
NCLS_vol = 0.05
E_redox = 1.0
#
voltage_limit_charge = 1.4 #0.4
voltage_limit_discharge = 0.6 #-0.4
current = 0.3

resistance = 1.0
k_species = 2.2e-3


# for crossover
membrane_thickness = 183/10000 # cm, nafion 117
membrane_c = area / membrane_thickness
p_ox = 1.0e-5 # cm^2/s
p_red = 1.0e-5 # cm^2/s
crossover_f = Crossover(membrane_constant=membrane_c, permeability_ox=p_ox, permeability_red=p_red)

## testing of abstract method classes


# define the battery design parameters
setup = battery(CLS_vol, NCLS_vol,
                CLS_start_conc_ox, CLS_start_conc_red,
                NCLS_start_conc_ox, NCLS_start_conc_red,
                E_redox, resistance,
                k_species, k_species)

# define degradation mechanisms
test_f1 = ChemicalDegradation(rate_order=1, rate=10e-5, species='red')
test_f2 = AutoOxidation(rate=30e-5)
test_f3 = ChemicalDegradation(rate_order=1, rate=10e-5, species='red')
test_f4 = MultiDegradationMechanism([test_f1, test_f2]) # maybe have multi do *args

# define cycling protocol and run based on defined cell and optional degradations

bbb = ConstantCurrent(voltage_cutoff_charge=voltage_limit_charge,
                      voltage_cutoff_discharge=voltage_limit_discharge,
                      current=current)
"""
bbb = ConstantCurrentConstantVoltage(voltage_limit_charge=voltage_limit_charge,
                                     voltage_limit_discharge=voltage_limit_discharge,
                                     current_cutoff_charge=0.005, current_cutoff_discharge=-0.005,
                                     current=current)
"""
# run based on defined cell and optional degradations
(current_profile, conc_ox_CLS_profile, conc_red_CLS_profile, conc_ox_NCLS_profile, conc_red_NCLS_profile,
 cell_V_profile, soc_profile_CLS, soc_profile_NCLS, ocv_profile, cycle_capacity, cycle_time, times, act_profile,
 mt_profile, loss_profile, del_ox, del_red,
) = bbb.run(cell_model=setup,
            cls_degradation=test_f1,
            #degradation=test_f2,
            #ncls_degradation=test_f2,
            #crossover_params=crossover_f,
            duration=5000)


##### PLOTTING BELOW ################
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
#ax[4].plot(times[g:], loss_profile, 'k')
#ax[4].plot(times[g:], act_profile, 'r')
#ax[4].plot(times[g:], mt_profile, 'b')
ax[4].plot(times[g:], soc_profile_CLS, 'r')
ax[4].plot(times[g:], soc_profile_NCLS, 'k')

ax[0].set_ylabel('Capacity')
ax[1].set_ylabel('Current')
ax[2].set_ylabel('Voltage')
ax[3].set_ylabel('OCV')
#ax[4].set_ylabel('Overpotential')
ax[4].set_ylabel('SOC (%)')
plt.show()

# plot out current and voltage cycyles