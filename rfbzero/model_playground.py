from zeroD_model_1e_vs_1e import ZeroDModel
import matplotlib.pyplot as plt
from zeroD_model_degradations import ChemicalDegradation, AutoOxidation, AutoReduction, MultiDegradationMechanism
from zeroD_model_crossover import Crossover
from cycle_protocol import ConstantCurrent, ConstantCurrentConstantVoltage

CLS_start_conc_ox = 0.01
CLS_start_conc_red = 0.01
NCLS_start_conc_ox = 0.01
NCLS_start_conc_red = 0.01
area = 5.0
CLS_vol = 0.005
NCLS_vol = 0.020
E_redox = 1.0
#
voltage_limit_charge = 1.6
voltage_limit_discharge = 0.4
current = 0.2

resistance = 1.0
k_species = 2.2e-3

# for crossover
membrane_thickness = 183 / 10000  # cm, nafion 117
membrane_c = area / membrane_thickness
p_ox = 1.0e-8  # cm^2/s
p_red = 1.0e-8  # cm^2/s
##############################

# define the battery design parameters
cell = ZeroDModel(CLS_vol, NCLS_vol, CLS_start_conc_ox, CLS_start_conc_red, NCLS_start_conc_ox, NCLS_start_conc_red,
                  E_redox, resistance, k_species, k_species, n_cls=2, n_ncls=1)

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
                           # cls_degradation=test_f1,
                           #degradation=test_f1,
                           # #ncls_degradation=test_f2,
                           #crossover_params=crossover_f,
                           duration=5000)

# print([attr for attr in dir(all_results) if not attr.startswith('__')])
print(all_results.cycle_capacity[:5])


def structure_data(x, y):
    t_charge = x[::2]
    t_discharge = x[1::2]
    capacity_charge = y[::2]
    capacity_discharge = y[1::2]

    if len(t_charge) < len(capacity_charge):
        capacity_charge = capacity_charge[:-1]
    if len(t_discharge) < len(capacity_discharge):
        capacity_discharge = capacity_discharge[:-1]
    return t_charge, t_discharge, capacity_charge, capacity_discharge


time_charge, time_discharge, cap_charge, cap_discharge = structure_data(all_results.cycle_time,
                                                                        all_results.cycle_capacity)

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

"""
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
# if __name__ == '__main__':
#    print('testing')
