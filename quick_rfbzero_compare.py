
from src.rfbzero.redox_flow_cell import ZeroDModel
from src.rfbzero.degradation import ChemicalDegradation, AutoOxidation, AutoReduction, MultiDegradationMechanism, AutoReductionO2Release
# from src.rfbzero.crossover import Crossover
from src.rfbzero.experiment import ConstantCurrent, ConstantCurrentConstantVoltage


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
k_species = 5e-3

# for crossover
# membrane_thickness = 183 / 10000  # cm, nafion 117
# membrane_c = area / membrane_thickness
# p_ox = 1.0e-6  # cm^2/s
# p_red = 1.0e-6  # cm^2/s
##############################

# define the battery design parameters
cell = ZeroDModel(cls_volume=0.005,     # L
                  ncls_volume=0.05,     # L
                  cls_start_c_ox=0.01,   # M
                  cls_start_c_red=0.01,  # M
                  ncls_start_c_ox=0.01,  # M
                  ncls_start_c_red=0.01,  # M
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
                           current=0.05,                     # A
                           )

# putting it all together
# all_results = protocol.run(cell_model=cell,
#                            degradation=deg,
#                            duration=1000,   # cycle time to simulate (s)
#                            )
#
# print(all_results.cycle_capacity[:5])

# should print out: [4.667500000000133, 9.357499999999945, 9.379499999999972, 9.379499999999972, 9.379499999999972]

########################################################################
# define the battery design parameters
cell = ZeroDModel(cls_volume=0.005,     # L
                  ncls_volume=0.05,     # L
                  cls_start_c_ox=0.01,   # M
                  cls_start_c_red=0.01,  # M
                  ncls_start_c_ox=0.01,  # M
                  ncls_start_c_red=0.01,  # M
                  init_ocv=1.2,         # V
                  resistance=1.0,       # ohms
                  k_0_cls=1e-3,         # cm/s
                  k_0_ncls=1e-3,        # cm/s
                  n_cls=1,              # electrons
                  n_ncls=1,             # electrons
                  )

# define degradation mechanisms
deg = ChemicalDegradation(rate_order=1,
                          rate_constant=1e-6,  # 1/s
                          species='red',
                          )

# define cycling protocol
protocol = ConstantCurrentConstantVoltage(voltage_limit_charge=1.5,voltage_limit_discharge=1.0,
                                          current_cutoff_charge=0.005,current_cutoff_discharge=-0.005,
                                          current=0.2)

# putting it all together
all_results = protocol.run(cell_model=cell,
                           degradation=deg,
                           duration=100,   # cycle time to simulate (s)
                           )
print(all_results.cycle_capacity[:5])
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot(all_results.times, all_results.current_profile)
ax.plot(all_results.times, all_results.cell_v_profile)
plt.show()

# should print out: [4.822421276429994, 9.634124272850148, 9.634583124882912, 9.63407360168263, 9.634085221274113]

#########################################

# define the battery design parameters
cell = ZeroDModel(cls_volume=0.005,     # L
                  ncls_volume=0.05,     # L
                  cls_start_c_ox=0.01,   # M
                  cls_start_c_red=0.01,  # M
                  ncls_start_c_ox=0.01,  # M
                  ncls_start_c_red=0.01,  # M
                  init_ocv=0.0,         # V
                  resistance=1.0,       # ohms
                  k_0_cls=1e-3,         # cm/s
                  k_0_ncls=1e-3,        # cm/s
                  n_cls=1,              # electrons
                  n_ncls=1,             # electrons
                  )

# define degradation mechanisms
deg = AutoReduction(rate_constant=5e-5)

# define cycling protocol
protocol = ConstantCurrent(voltage_cutoff_charge=0.2,       # V
                           voltage_cutoff_discharge=-0.2,   # V
                           current=0.05,                     # A
                           )

# putting it all together
# all_results = protocol.run(cell_model=cell,
#                            degradation=deg,
#                            duration=1000,   # cycle time to simulate (s)
#                            )
#
# print(all_results.cycle_capacity[:5])

# should print out: [4.656500000000139, 9.400999999999998, 9.33649999999992, 9.423500000000026, 9.335999999999919]
