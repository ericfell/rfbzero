
import matplotlib.pyplot as plt
import time
from scipy import stats
import numpy as np

from redox_flow_cell import ZeroDModel
from degradation import ChemicalDegradation, AutoOxidation, AutoReduction, MultiDegradationMechanism, AutoReductionO2Release
from crossover import Crossover
from experiment import ConstantCurrent, ConstantCurrentConstantVoltage


def fade_rate(capacity, time_capacity):
    y = np.log(capacity)
    slope, _, r, *_ = stats.linregress(time_capacity, y)
    return slope*-100, r**2


c_ox_cls_start = 0.1
c_red_cls_start = 0.1
c_ox_ncls_start = 0.1
c_red_ncls_start = 0.1
area = 5.0
CLS_vol = 0.006
NCLS_vol = 0.05

E_redox = 0.0
voltage_limit_charge = 0.2
voltage_limit_discharge = -0.2
current = 0.5

resistance = 1.5
k_species = 5e-3  # AQDS nature paper says 7.3e-3

# for crossover
membrane_thickness = 183 / 10000  # cm, nafion 117 # 183, 120, 50, 25
membrane_c = area / membrane_thickness
p_ox = 1.0e-6  # cm^2/s
p_red = 1.0e-6  # cm^2/s
##############################
# working on timing here
"""
start_time = time.time()

# define the battery design parameters
cell = ZeroDModel(CLS_vol, NCLS_vol, c_ox_cls_start, c_red_cls_start, c_ox_ncls_start, c_red_ncls_start,
                  E_redox, resistance, k_species, k_species, n_cls=1, n_ncls=1)

# define degradation mechanisms
test_f1 = ChemicalDegradation(rate_order=1, rate_constant=1e-7, species='red')
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
                                              current=current, capacity_only=False)#, charge_first=False)
# putting it all together
all_results = protocol.run(cell_model=cell,
                           #cls_degradation=test_f1,
                           #degradation=test_f1,
                           #ncls_degradation=test_f2,
                           #cross_over=crossover_f,
                           duration=1000)

# print([attr for attr in dir(all_results) if not attr.startswith('__')])
print(f"Execution time: {(time.time() - start_time):.2f}")
print(all_results.cycle_capacity[:5])

t_dis = [i/86400 for i in all_results.time_discharge]
fade, r2 = fade_rate(all_results.discharge_capacity, t_dis)
print(f"Fade rate: {fade:.3f}%/day, r^2: {r2:.3f}")

fig, ax = plt.subplots(nrows=5, ncols=1, sharex=True)
#ax[0].plot(all_results.time_discharge, all_results.discharge_capacity, 'bo')
ax[0].plot(all_results.times, all_results.current_profile, 'r')
ax[1].plot(all_results.times, all_results.cell_v_profile, 'b')
ax[2].plot(all_results.times, all_results.ocv_profile, 'k')
ax[3].plot(all_results.times, all_results.act_profile, 'r')
ax[4].plot(all_results.times, all_results.mt_profile, 'k')
ax[0].set_ylabel('Current (A)')
ax[1].set_ylabel('Voltage (V)')
ax[2].set_ylabel('OCV (V)')
ax[3].set_ylabel(r'$\eta_{act}$ (V)')
ax[4].set_ylabel(r'$\eta_{mt}$ (V)')
ax[4].set_xlabel('Time (s)')
#ax[0].set_ylim(-0.4, 0.4)
ax[1].set_ylim(-0.22, 0.22)
ax[2].set_ylim(-0.22, 0.22)
#ax[2].set_xlabel('Time (s)')
ax[3].set_ylim(0, 0.04)
ax[4].set_ylim(0, 0.04)
plt.show()
"""


"""
ax[0].plot(all_results.time_discharge, all_results.discharge_capacity, 'bo--')
ax[0].plot(all_results.time_charge, all_results.charge_capacity, 'ro--')
ax[1].plot(all_results.times, all_results.current_profile, 'r')
ax[2].plot(all_results.times, all_results.cell_v_profile, 'b')
ax[3].plot(all_results.times, all_results.ocv_profile, 'k')
# ax[4].plot(times[g:], loss_profile, 'k')
# ax[4].plot(times[g:], act_profile, 'r')
# ax[4].plot(times[g:], mt_profile, 'b')
# ax[4].plot(all_results.times[g:], all_results.soc_profile_cls, 'r')
# ax[4].plot(all_results.times[g:], all_results.soc_profile_ncls, 'k')
#plt.show()

ax[0].set_ylabel('Capacity')
ax[1].set_ylabel('Current')
ax[2].set_ylabel('Voltage')
ax[3].set_ylabel('OCV')
# ax[4].set_ylabel('Overpotential')
# ax[4].set_ylabel('SOC (%)')
"""

"""
# plot out current and voltage cycyles
fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True)
ax[0].plot(all_results.times, all_results.del_ox, 'b--')
ax[0].plot(all_results.times, all_results.del_red, 'r--')

ax[1].plot(all_results.times, all_results.soc_profile_cls, 'r')
ax[1].plot(all_results.times, all_results.soc_profile_ncls, 'k')
ax[1].set_xlabel('Time (s)')
ax[0].set_ylabel("Mols of species crossed")  # (r'$_{CLS} - C_{NCLS}$ (M)')
ax[1].set_ylabel('SOC (%)')
plt.show()
"""

"""
#param_try = [1.0e-6, 1.0e-7, 1.0e-8]

c1 = ChemicalDegradation(rate_order=1, rate_constant=1e-7, species='red')
#c2 = ChemicalDegradation(rate_order=1, rate_constant=1e-7, species='ox')
#d1 = MultiDegradationMechanism([c1, c2])

#param_try = [c1, c2, d1]
labs = ["red: 1e-7 1/s", "ox: 1e-7 1/s", "red=ox: 1e-7 1/s"]
#labs = ["P_red: 1e-7", "P_ox: 1e-7", "P_red=P_ox: 1e-7"]

p1 = Crossover(membrane_constant=membrane_c, permeability_ox=1e-7, permeability_red=0.0)
p2 = Crossover(membrane_constant=membrane_c, permeability_ox=0.0, permeability_red=1.0e-7)
p3 = Crossover(membrane_constant=membrane_c, permeability_ox=1e-7, permeability_red=1.0e-7)

#param_try = [p1, p2, p3]
"""



#c1 = ChemicalDegradation(rate_order=1, rate_constant=1e-8, species='red')
#c2 = ChemicalDegradation(rate_order=1, rate_constant=2e-8, species='red')
#c3 = ChemicalDegradation(rate_order=1, rate_constant=3e-8, species='red')
#c4 = ChemicalDegradation(rate_order=1, rate_constant=4e-8, species='red')

#c1 = None #AutoReduction(rate_constant=1e-8)
#c2 = AutoReduction(rate_constant=1e-4)
c2 = AutoReduction(rate_constant=4e-5)
#c2 = AutoReductionO2Release(rate_constant=5e-4, release_factor=5e-5)
c3 = AutoReductionO2Release(rate_constant=4e-5, release_factor=1e-5)
c4 = AutoReductionO2Release(rate_constant=4e-5, release_factor=5e-5)
#c3 = AutoReduction(rate_constant=1e-6)
#c4 = AutoReduction(rate_constant=3e-5)
p_try = [
         AutoReduction(rate_constant=5e-5),
         AutoReduction(rate_constant=3e-4),
         AutoReduction(rate_constant=8e-4),
         ]

#p_try = [0.032, 0.016, 0.008, 0.006]#[c1, c2, c4]# c2, c3, c4]
#labs = ["auto 6e-5", "auto 8e-5", "auto 1e-4", "auto 1.2e-4", "auto 1.4e-4"]
#labs = ["32 mL", "16 mL", "8 mL", "6 mL"]#["0 1/s", "5e-4 1/s", "1e-3 1/s"]
labs = ["auto 5e-5", "auto 3e-4", "auto 8e-4"]
cols = ['r', 'b', 'orange']
#membrane_thickness = 183 / 10000  # cm, nafion 117 # 183, 120, 50, 25
#mem_thick = [25 / 10000, 50 / 10000, 125 / 10000, 183 / 10000]
#mem_cons = [5.0 / i for i in mem_thick]
#pox = [8.3e-9, 6.7e-9, 1.6e-9, 2.3e-9]
#p_try = mem_cons
#labs = ["NR211", "NR212", "N115", "N117"]
#res = [0.034, 0.036, 0.06, 0.075]

fig, ax = plt.subplots(nrows=2,ncols=1)#,sharex=True)
start_time = time.time()
#for idx, (z,p,r) in enumerate(zip(p_try, pox, res)):
theory = CLS_vol*c_ox_cls_start*2*96485.33
for idx, z in enumerate(p_try):
    # change resistance for diff membranes?
    cell = ZeroDModel(CLS_vol, 0.011, c_ox_cls_start, c_red_cls_start, c_ox_ncls_start, c_red_ncls_start,
                      E_redox, resistance, k_species, k_species, n_cls=1, n_ncls=1, time_increment=0.1)

    # define cycling protocol
    protocol = ConstantCurrentConstantVoltage(voltage_limit_charge=voltage_limit_charge,
                                              voltage_limit_discharge=voltage_limit_discharge,
                                              current_cutoff_charge=0.005,
                                              current_cutoff_discharge=-0.005,
                                              current=current, capacity_only=False)#, charge_first=False)
    # simulate
    all_results = protocol.run(cell_model=cell,
                               #cls_degradation=z,
                               degradation=z,
                               #ncls_degradation=AutoReduction(rate_constant=3e-5),
                               #cross_over=Crossover(membrane_constant=z, permeability_ox=p*10, permeability_red=p),
                               duration=40000)

    t_dis = [i/86400 for i in all_results.time_discharge]
    t_ch = [i / 86400 for i in all_results.time_charge]
    t_time = [i / 86400 for i in all_results.times]

    fade, r2 = fade_rate(all_results.discharge_capacity, t_dis)
    print(f"Fade rate: {fade:.3f}%/day, r^2: {r2:.3f}")
    # plot
    #ax[0].plot(t_dis, [(i/theory)*100 for i in all_results.discharge_capacity], 'o--', color=cols[idx], label=labs[idx])  # label="p_ox= " + str(pox))
    ax[0].plot(t_dis, all_results.discharge_capacity, 'o--', color=cols[idx], label=labs[idx])
    #ax[1].plot(t_time, all_results.current_profile, color=cols[idx])
    ax[1].plot([100 - i for i in all_results.soc_profile_ncls], [100 - i for i in all_results.soc_profile_cls], color=cols[idx])
    #ax.plot(t_ch, all_results.charge_capacity, 'rv--', label='Charge,' + labs[idx])
print(f"Execution time: {(time.time() - start_time):.2f}")

ax[0].set_xlabel('Time (days)')
ax[0].set_ylabel('Discharge Capacity (C)')
ax[0].axhline(theory, color='k', linestyle='--')
#ax.set_yticklabels([])
#ax[1].set_ylabel('Current (A)')
ax[0].set_ylim(bottom=0)
#ax.set_xticklabels([])
ax[0].legend()

ax[1].set_ylabel('CLS (% oxidized)')
ax[1].set_xlabel('NCLS (% oxidized)')
ax[1].set_ylim(bottom=0,top=100)
ax[1].set_xlim(left=0,right=100)
plt.show()





"""
# define the battery design parameters
cell = ZeroDModel(cls_volume=0.005,     # L
                  ncls_volume=0.05,     # L
                  cls_start_c_ox=0.1,   # M
                  cls_start_c_red=0.1,  # M
                  ncls_start_c_ox=0.1,  # M
                  ncls_start_c_red=0.1, # M
                  init_ocv=0.0,         # V
                  resistance=1.0,       # ohms
                  k_0_cls=1e-3,         # cm/s
                  k_0_ncls=1e-3,        # cm/s
                  n_cls=1,              # electrons
                  n_ncls=1,             # electrons
                  )

# define degradation mechanisms
deg = ChemicalDegradation(rate_order=1,
                          rate_constant=1e-7,  # 1/s
                          species='red',
                          )

# define cycling protocol
protocol = ConstantCurrent(voltage_cutoff_charge=0.2,       # V
                           voltage_cutoff_discharge=-0.2,   # V
                           current=0.2,                     # A
                           )

# putting it all together
all_results = protocol.run(cell_model=cell,
                           degradation=deg,
                           duration=1000,   # cycle time to simulate (s)
                           )
"""