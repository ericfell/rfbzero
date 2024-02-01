[![Documentation Status](https://readthedocs.org/projects/rfbzero/badge/?version=latest)](https://rfbzero.readthedocs.io/en/latest/?badge=latest)  [![codecov](https://codecov.io/github/ericfell/rfbzero/graph/badge.svg)](https://codecov.io/github/ericfell/rfbzero)  [![Checked with mypy](http://www.mypy-lang.org/static/mypy_badge.svg)](http://mypy-lang.org/)  [![linting: pylint](https://img.shields.io/badge/linting-pylint-yellowgreen)](https://github.com/pylint-dev/pylint)



RFBzero
------------

`rfbzero` is a Python package for zero dimensional simulation of electrochemical cycling in redox flow batteries (RFBs). 

This package includes modules to describe the initial flow cell setup, chemical and electrochemical properties of the redox-active electrolytes being cycled, cell cycling protocol selection, and optional inputs for capacity degradation mechanisms and active species crossover.

## Installation

`rfbzero` can be installed from [PyPI](https://pypi.org/project/rfbzero/) using pip.

```bash
pip install rfbzero
```

See [Getting started with `rfbzero.py`](https://rfbzero.readthedocs.io/en/latest/getting-started.html) for instructions on simulating RFBs.

### Dependencies

`rfbzero.py` requires:

- Python (>=3.10)
- SciPy

### Examples and Documentation

Several simulated RFB examples can be found in a notebook in `docs/source/examples` (requires Jupyter notebook/lab). The documentation can be found at [rfbzero.readthedocs](https://rfbzero.readthedocs.io/en/latest/index.html).

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
    volume_cls=0.005,   # liters
    volume_ncls=0.010,  # liters
    c_ox_cls=0.01,      # molar
    c_red_cls=0.01,     # molar
    c_ox_ncls=0.01,     # molar
    c_red_ncls=0.01,    # molar
    ocv_50_soc=1.0,     # volts
    resistance=0.5,     # ohms
    k_0_cls=1e-3,       # cm/s
    k_0_ncls=1e-3,      # cm/s
)

# define cycling protocol
protocol = ConstantCurrent(
    voltage_limit_charge=1.5,      # volts
    voltage_limit_discharge=0.5,   # volts
    current=0.1,                   # amps
)

# simulate cell via protocol for 1000 seconds
results = protocol.run(cell_model=cell, duration=1000)
```

or a symmetric cell with constant current constant voltage (CCCV) cycling, with an electrolyte undergoing a first order chemical degradation in the reduced form:

```python
from redox_flow_cell import ZeroDModel
from degradation import ChemicalDegradationReduced
from experiment import ConstantCurrentConstantVoltage


# define symmetric cell and electrolyte parameters
cell = ZeroDModel(
    volume_cls=0.005,   # liters
    volume_ncls=0.010,  # liters
    c_ox_cls=0.01,      # molar
    c_red_cls=0.01,     # molar
    c_ox_ncls=0.01,     # molar
    c_red_ncls=0.01,    # molar
    ocv_50_soc=0.0,     # volts
    resistance=0.5,     # ohms
    k_0_cls=1e-3,       # cm/s
    k_0_ncls=1e-3,      # cm/s
)

# define cycling protocol
protocol = ConstantCurrentConstantVoltage(
    voltage_limit_charge=0.2,         # volts
    voltage_limit_discharge=-0.2,     # volts
    current_cutoff_charge=0.005,      # amps
    current_cutoff_discharge=-0.005,  # amps
    current=0.1,                      # amps
)

# define first order chemical degradation mechanism
chem_deg = ChemicalDegradationReduced(rate_order=1, rate_constant=1e-5)

# simulate cell via protocol for 1000 seconds
results = protocol.run(cell_model=cell, degradation=chem_deg, duration=1000)
```



## License
[MIT](https://choosealicense.com/licenses/mit/) 
