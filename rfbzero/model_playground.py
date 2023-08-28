
import matplotlib.pyplot as plt
import time
from scipy import stats
import numpy as np

from redox_flow_cell import ZeroDModel
from degradation import ChemicalDegradation, AutoOxidation, AutoReduction, MultiDegradationMechanism
from crossover import Crossover
from experiment import ConstantCurrent, ConstantCurrentConstantVoltage


def fade_rate(capacity, time_capacity):
    y = np.log(capacity)
    slope, _, r, *_ = stats.linregress(time_capacity, y)
    return slope*-100, r**2


c_ox_cls_start = 0.01
c_red_cls_start = 0.01
c_ox_ncls_start = 0.01
c_red_ncls_start = 0.01
area = 5.0
CLS_vol = 0.005
NCLS_vol = 0.010
E_redox = 0.0
#
voltage_limit_charge = 0.2
voltage_limit_discharge = -0.2
current = 0.05

resistance = 0.1
k_species = 5e-3  # AQDS nature paper says 7.3e-3

# for crossover
membrane_thickness = 183 / 10000  # cm, nafion 117
membrane_c = area / membrane_thickness
p_ox = 1.0e-6  # cm^2/s
p_red = 1.0e-6  # cm^2/s
##############################
# working on timing here
"""
start_time = time.time()

# define the battery design parameters
cell = ZeroDModel(CLS_vol, NCLS_vol, c_ox_cls_start, c_red_cls_start, c_ox_ncls_start, c_red_ncls_start,
                  E_redox, resistance, k_species, k_species, n_cls=2, n_ncls=2)

# define degradation mechanisms
test_f1 = ChemicalDegradation(rate_order=1, rate_constant=10e-5, species='red')
# test_f2 = AutoOxidation(rate_constant=30e-5)
# test_f3 = ChemicalDegradation(rate_order=1, rate_constant=10e-5, species='red')
# test_f4 = MultiDegradationMechanism([test_f1, test_f2])
# crossover mechanism
crossover_f = Crossover(membrane_constant=membrane_c, permeability_ox=p_ox, permeability_red=p_red)

# define cycling protocol and run based on defined cell and optional degradations
run_CC = False

if run_CC:
    protocol = ConstantCurrent(voltage_cutoff_charge=voltage_limit_charge,
                               voltage_cutoff_discharge=voltage_limit_discharge,
                               current=current)
else:
    protocol = ConstantCurrentConstantVoltage(voltage_limit_charge=voltage_limit_charge,
                                              voltage_limit_discharge=voltage_limit_discharge,
                                              current_cutoff_charge=0.005,
                                              current_cutoff_discharge=-0.005,
                                              current=current)
# putting it all together
all_results = protocol.run(cell_model=cell,
                           cls_degradation=test_f1,
                           #degradation=test_f1,
                           #ncls_degradation=test_f2,
                           cross_over=crossover_f,
                           duration=10000)

# print([attr for attr in dir(all_results) if not attr.startswith('__')])
print(f"Execution time: {(time.time() - start_time):.2f}")
print(all_results.cycle_capacity[:5])


g = 0
fig, ax = plt.subplots(nrows=4, ncols=1, sharex=True)
ax[0].plot(time_discharge, cap_discharge, 'bo--')
ax[0].plot(time_charge, cap_charge, 'ro--')
ax[1].plot(all_results.times[g:], all_results.current_profile, 'r')
ax[2].plot(all_results.times[g:], all_results.cell_v_profile, 'b')
ax[3].plot(all_results.times[g:], all_results.ocv_profile, 'k')
# ax[4].plot(times[g:], loss_profile, 'k')
# ax[4].plot(times[g:], act_profile, 'r')
# ax[4].plot(times[g:], mt_profile, 'b')
# ax[4].plot(all_results.times[g:], all_results.soc_profile_cls, 'r')
# ax[4].plot(all_results.times[g:], all_results.soc_profile_ncls, 'k')
plt.show()

ax[0].set_ylabel('Capacity')
ax[1].set_ylabel('Current')
ax[2].set_ylabel('Voltage')
ax[3].set_ylabel('OCV')
# ax[4].set_ylabel('Overpotential')
# ax[4].set_ylabel('SOC (%)')


# plot out current and voltage cycyles
fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True)
ax[0].plot(all_results.times, all_results.del_ox, 'b--')
ax[0].plot(all_results.times, all_results.del_red, 'r--')

ax[1].plot(all_results.times[g:], all_results.soc_profile_cls, 'r')
ax[1].plot(all_results.times[g:], all_results.soc_profile_ncls, 'k')
ax[1].set_xlabel('Time (s)')
ax[0].set_ylabel("Mols of species crossed")  # (r'$_{CLS} - C_{NCLS}$ (M)')
ax[1].set_ylabel('SOC (%)')
plt.show()
"""


#param_try = [1.0e-6, 1.0e-7, 1.0e-8]


c1 = ChemicalDegradation(rate_order=1, rate_constant=1e-7, species='red')
c2 = ChemicalDegradation(rate_order=1, rate_constant=1e-7, species='ox')
d1 = MultiDegradationMechanism([c1, c2])

param_try = [c1, c2, d1]
labs = ["red: 1e-7 1/s", "ox: 1e-7 1/s", "red=ox: 1e-7 1/s"]
#labs = ["P_red: 1e-7", "P_ox: 1e-7", "P_red=P_ox: 1e-7"]

p1 = Crossover(membrane_constant=membrane_c, permeability_ox=1e-7, permeability_red=0.0)
p2 = Crossover(membrane_constant=membrane_c, permeability_ox=0.0, permeability_red=1.0e-7)
p3 = Crossover(membrane_constant=membrane_c, permeability_ox=1e-7, permeability_red=1.0e-7)

#param_try = [p1, p2, p3]

fig, ax = plt.subplots()
start_time = time.time()
for idx, z in enumerate(param_try):

    # define cell
    cell = ZeroDModel(CLS_vol, NCLS_vol, c_ox_cls_start, c_red_cls_start, c_ox_ncls_start, c_red_ncls_start,
                      E_redox, resistance, k_species, k_species, n_cls=2, n_ncls=2)

    # define degradations
    #crossover_f = Crossover(membrane_constant=membrane_c, permeability_ox=0.0, permeability_red=1.0e-7)
    #chem_deg = ChemicalDegradation(rate_order=1, rate_constant=1e-7, species='red')
    #chem_deg2 = ChemicalDegradation(rate_order=1, rate_constant=1e-7, species='ox')
    #deg = MultiDegradationMechanism([chem_deg, chem_deg2])

    # define cycling protocol
    protocol = ConstantCurrentConstantVoltage(voltage_limit_charge=voltage_limit_charge,
                                              voltage_limit_discharge=voltage_limit_discharge,
                                              current_cutoff_charge=0.005,
                                              current_cutoff_discharge=-0.005,
                                              current=current)#, charge_first=False)
    # simulate
    all_results = protocol.run(cell_model=cell,
                               #cls_degradation=test_f1,
                               degradation=z,
                               # ncls_degradation=test_f2,
                               #cross_over=z,
                               duration=30000)
    #print(all_results.cycle_capacity[:5])

    t_dis = [i/86400 for i in all_results.time_discharge]

    fade, r2 = fade_rate(all_results.discharge_capacity[1:], t_dis[1:])
    print(f"Fade rate: {fade:.3f}%/day, r^2: {r2:.3f}")
    # plot
    ax.plot(t_dis[1:], all_results.discharge_capacity[1:], 'o--', label=labs[idx]) # label="p_ox= " + str(pox))
    #ax[idx].plot(time_charge, cap_charge, 'ro--')
print(f"Execution time: {(time.time() - start_time):.2f}")
ax.set_xlabel('Time (day)')
ax.set_ylabel('Discharge Capacity (C)')
ax.legend()
plt.show()

#ax[0].set_ylabel('Capacity')
#ax[].set_ylabel('Current')
#ax[2].set_ylabel('Voltage')
#ax[3].set_ylabel('OCV')
