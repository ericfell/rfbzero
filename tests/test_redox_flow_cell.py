import pytest
import numpy as np

#from rfbzero.rfbzero.redox_flow_cell import ZeroDModel
from redox_flow_cell import ZeroDModel # uncomment here when running test harness


class TestClassRedoxFlowCell:

    def test_class_init(self):
        with pytest.raises(ValueError):
            ZeroDModel(cls_volume=6, ncls_volume=0.01, cls_start_c_ox=0.01, cls_start_c_red=0.01, ncls_start_c_ox=0.01,
                       ncls_start_c_red=0.01, init_ocv=1.0, resistance=1, k_0_cls=1e-3, k_0_ncls=1e-3)

    def test_exchange_current(self):
        cell = ZeroDModel(cls_volume=0.005, ncls_volume=0.01, cls_start_c_ox=0.01, cls_start_c_red=0.01,
                          ncls_start_c_ox=0.01, ncls_start_c_red=0.01, init_ocv=1.0,resistance=1,k_0_cls=1e-3,
                          k_0_ncls=1e-3, n_ncls=2)
        i_0_cls, i_0_ncls = cell._exchange_current()
        assert np.isclose(i_0_cls, 0.12543093175)
        assert np.isclose(i_0_ncls, 0.25086186351)

    def test_limiting_current(self):
        limiting_c = 0.2
        cell = ZeroDModel(cls_volume=0.005, ncls_volume=0.01, cls_start_c_ox=0.01, cls_start_c_red=0.01,
                          ncls_start_c_ox=0.01, ncls_start_c_red=0.01, init_ocv=1.0, resistance=1, k_0_cls=1e-3,
                          k_0_ncls=1e-3, n_ncls=2)

        i_lim = cell._limiting_current(limiting_c)
        assert np.isclose(i_lim, 77.1882656)










    def test_cell_voltage(self):
        ocv = 1.0
        loss = 0.2

        cell_v1 = ZeroDModel.cell_voltage(ocv, loss, True)
        cell_v2 = ZeroDModel.cell_voltage(ocv, loss, False)

        assert cell_v1 == 1.2
        assert cell_v2 == 0.8





