import pytest
#import numpy as np

from rfbzero.experiment import ConstantCurrent, ConstantCurrentConstantVoltage, ConstantVoltage
from rfbzero.redox_flow_cell import ZeroDModel
from rfbzero.degradation import ChemicalDegradation, AutoOxidation, AutoReduction, MultiDegradationMechanism#, Dimerization
from rfbzero.crossover import Crossover


class TestCyclingProtocolResults:


    @pytest.mark.skip(reason="not implemented yet")
    def test_state_of_charge(self):
        r_cls = 0.2
        o_cls = 0.8
        r_ncls = 0.3
        o_ncls = 0.7

        soc_c, soc_n = CyclingProtocolResults.state_of_charge(o_cls, r_cls, o_ncls, r_ncls)

        assert soc_c == 20.0
        assert soc_n == 30.0

    @pytest.mark.skip(reason="not implemented yet")
    def test_structure_data(self):
        raise NotImplementedError


class TestConstantCurrent:

    @pytest.mark.skip(reason="not implemented yet")
    def test_cc(self):
        raise NotImplementedError


class TestConstantVoltage:

    def test_cv1(self):
        # define the battery design parameters
        cell = ZeroDModel(cls_volume=0.005,  # L
                          ncls_volume=0.05,  # L
                          cls_start_c_ox=0.01,  # M
                          cls_start_c_red=0.01,  # M
                          ncls_start_c_ox=0.01,  # M
                          ncls_start_c_red=0.01,  # M
                          init_ocv=0.0,  # V
                          resistance=1.0,  # ohms
                          k_0_cls=1e-3,  # cm/s
                          k_0_ncls=1e-3,  # cm/s
                          n_cls=1,  # electrons
                          n_ncls=1,  # electrons
                          )

        # define degradation mechanisms
        deg = ChemicalDegradation(rate_order=1,
                                  rate_constant=1e-5,  # 1/s
                                  species='red',
                                  )

        protocol = ConstantCurrentConstantVoltage(
            voltage_limit_charge=0.2,  # volts
            voltage_limit_discharge=-0.2,  # volts
            current_cutoff_charge=0.005,  # amps
            current_cutoff_discharge=-0.005,  # amps
            current=1.5,  # amps
        )

        # putting it all together
        all_results = protocol.run(cell_model=cell,
                                   degradation=deg,
                                   duration=1000,  # cycle time to simulate (s)
                                   )
        expected = [4.809806770737516, 9.616034709809167, 9.61843069453628, 9.611898460659235, 9.611890636955804]

        assert all_results.half_cycle_capacity[:5] == expected


class TestConstantCurrentConstantVoltage:

    @pytest.mark.skip(reason="not implemented yet")
    def test_cccv(self):
        raise NotImplementedError




