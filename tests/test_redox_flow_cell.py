import pytest
import numpy as np

from rfbzero.redox_flow_cell import ZeroDModel
#from redox_flow_cell import ZeroDModel # uncomment here when running test harness


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

    def test_limiting_concentration(self):
        cell = ZeroDModel(cls_volume=0.005, ncls_volume=0.01, cls_start_c_ox=0.01, cls_start_c_red=0.02,
                          ncls_start_c_ox=0.02, ncls_start_c_red=0.01, init_ocv=1.0, resistance=1, k_0_cls=1e-3,
                          k_0_ncls=1e-3, n_ncls=2)

        i_lim_cls, i_lim_ncls = cell.limiting_concentration(True)
        assert np.isclose(i_lim_cls, 3.85941328)
        assert np.isclose(i_lim_ncls, 7.71882656)

        i_lim_cls, i_lim_ncls = cell.limiting_concentration(False)
        assert np.isclose(i_lim_cls, 7.71882656)
        assert np.isclose(i_lim_ncls, 15.43765312)

    def test_activation_overpotential(self):
        cell = ZeroDModel(cls_volume=0.005, ncls_volume=0.01, cls_start_c_ox=0.01, cls_start_c_red=0.01,
                          ncls_start_c_ox=0.01, ncls_start_c_red=0.01, init_ocv=1.0, resistance=1, k_0_cls=1e-3,
                          k_0_ncls=1e-3, n_ncls=2)
        current = 1
        i_0_cls = 0.01
        i_0_ncls = 0.01
        n_activation = cell._activation_overpotential(current, i_0_cls, i_0_ncls)
        assert np.isclose(n_activation, 0.177392243)

    def test_negative_concentrations(self):
        raise NotImplementedError

    def test_mass_transport_overpotential(self):
        raise NotImplementedError

    def test_total_overpotential(self):
        raise NotImplementedError

    def test_open_circuit_voltage(self):
        raise NotImplementedError

    def test_cell_voltage(self):
        raise NotImplementedError

    def test_coulomb_counter(self):
        raise NotImplementedError









    def test_cell_voltage(self):
        ocv = 1.0
        loss = 0.2

        cell_v1 = ZeroDModel.cell_voltage(ocv, loss, True)
        cell_v2 = ZeroDModel.cell_voltage(ocv, loss, False)

        assert cell_v1 == 1.2
        assert cell_v2 == 0.8





