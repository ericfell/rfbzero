import pytest
import numpy as np

from rfbzero.experiment import CyclingProtocolResults, CyclingProtocol, ConstantCurrent, ConstantCurrentConstantVoltage
#from experiment import CyclingProtocolResults, CyclingProtocol, ConstantCurrent, ConstantCurrentConstantVoltage
from rfbzero.redox_flow_cell import ZeroDModel
#from redox_flow_cell import ZeroDModel

class TestCyclingProtocolResults:

    @pytest.mark.skip(reason="not implemented yet")
    def test_results_init(self):
        raise NotImplementedError

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


class TestCyclingProtocol:

    @pytest.mark.skip(reason="not implemented yet")
    def test_abstract_class_init(self):
        raise NotImplementedError

    def test_abstract_class_run(self):
        test_model = ZeroDModel(cls_volume=6, ncls_volume=0.01, cls_start_c_ox=0.01, cls_start_c_red=0.01,
                                ncls_start_c_ox=0.01, ncls_start_c_red=0.01, init_ocv=1.0, resistance=1, k_0_cls=1e-3,
                                k_0_ncls=1e-3)
        with pytest.raises(NotImplementedError):
            abs_class = CyclingProtocol(current=1, charge_first=True)
            abs_class.run(duration=100, cell_model=test_model)

    @pytest.mark.skip(reason="not implemented yet")
    def test_validate_model(self):
        raise NotImplementedError

    def test_current_direction(self):
        x = CyclingProtocol(current=1, charge_first=True)
        y = CyclingProtocol(current=1, charge_first=False)
        assert x.current_direction() == 1
        assert y.current_direction() == -1


class TestConstantCurrent:

    @pytest.mark.skip(reason="not implemented yet")
    def test_cc_init(self):
        raise NotImplementedError

    @pytest.mark.skip(reason="not implemented yet")
    def test_cc_run(self):  # BIG ONE
        raise NotImplementedError


class TestConstantCurrentConstantVoltage:

    @pytest.mark.skip(reason="not implemented yet")
    def test_cccv_init(self):
        raise NotImplementedError

    @pytest.mark.skip(reason="not implemented yet")
    def test_cccv_run(self):  # BIG ONE
        raise NotImplementedError

    @pytest.mark.skip(reason="not implemented yet")
    def test_get_min_current(self):
        raise NotImplementedError



