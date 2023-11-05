

# RFBzero

`rfbzero` is a Python package for zero dimensional simulation of electrochemical cycling in redox flow batteries (RFBs). 

This package includes modules to describe the initial flow cell setup, chemical and electrochemical properties of the redox-active electrolytes being cycled, cell cycling protocol selection, and optional inputs for capacity degradation mechanisms and active species crossover.

## Installation

`rfbzero` can be installed from GitHub:

```bash
git clone git@github.com:ericfell/rfbzero.git
```

##  Package Structure

- `redox_flow_cell.py`: Configures flow cell and electrolyte parameters
- `experiment.py`: Specifies electrochemical cycling protocol
- `degradation.py`: Optional degradation mechanism functions
- `crossover.py`: Optional crossover mechanism


## Examples

To simulate a full cell with constant current (CC) cycling:

```python
from redox_flow_cell import ZeroDModel
from experiment import ConstantCurrent


# define full cell and electrolyte parameters
cell = ZeroDModel(
    cls_volume=0.005,       # liters
    ncls_volume=0.010,      # liters
    cls_start_c_ox=0.01,    # molar
    cls_start_c_red=0.01,   # molar
    ncls_start_c_ox=0.01,   # molar
    ncls_start_c_red=0.01,  # molar
    init_ocv=1.0,           # volts
    resistance=0.5,         # ohms
    k_0_cls=1e-3,           # cm/s
    k_0_ncls=1e-3,          # cm/s
)

# define cycling protocol
protocol = ConstantCurrent(
    voltage_limit_charge=1.5,      # volts
    voltage_limit_discharge=0.5,   # volts
    current=0.1,                    # amps
)

# simulate cell via protocol for 1000 seconds
results = protocol.run(cell_model=cell, duration=1000)
```

or a symmetric cell with constant current constant voltage (CCCV) cycling, with an electrolyte undergoing a first order chemical degradation in the reduced form:

```python
from redox_flow_cell import ZeroDModel
from degradation import ChemicalDegradation
from experiment import ConstantCurrentConstantVoltage


# define symmetric cell and electrolyte parameters
cell = ZeroDModel(
    cls_volume=0.005,       # liters
    ncls_volume=0.010,      # liters
    cls_start_c_ox=0.01,    # molar
    cls_start_c_red=0.01,   # molar
    ncls_start_c_ox=0.01,   # molar
    ncls_start_c_red=0.01,  # molar
    init_ocv=0.0,           # volts
    resistance=0.5,         # ohms
    k_0_cls=1e-3,           # cm/s
    k_0_ncls=1e-3,          # cm/s
)

# define cycling protocol
protocol = ConstantCurrentConstantVoltage(
    voltage_limit_charge=0.2,               # volts
    voltage_limit_discharge=-0.2,           # volts
    current_limit_charge=0.005,            # amps
    current_limit_discharge=-0.005,        # amps
    current=0.1,                            # amps
)

# define first order chemical degradation mechanism
chem_deg = ChemicalDegradation(rate_order=1, rate_constant=1e-5, species='red')

# simulate cell via protocol for 1000 seconds
results = protocol.run(cell_model=cell, degradation=chem_deg, duration=1000)
```



## License
[MIT](https://choosealicense.com/licenses/mit/) 
