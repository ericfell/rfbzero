import time

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

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
k_species = 5e-3

# for crossover
membrane_thickness = 183 / 10000  # cm, nafion 117
membrane_c = area / membrane_thickness
p_ox = 1.0e-6  # cm^2/s
p_red = 1.0e-6  # cm^2/s
##############################

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
protocol = ConstantCurrent(voltage_limit_charge=0.2,      # V
                           voltage_limit_discharge=-0.2,  # V
                           current=0.2,                   # A
                           )

# putting it all together
all_results = protocol.run(cell_model=cell,
                           degradation=deg,
                           duration=1000,   # cycle time to simulate (s)
                           )
