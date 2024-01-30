import pytest
import numpy as np

from rfbzero.redox_flow_cell import ZeroDModel


class TestClassRedoxFlowCell:

    @pytest.mark.parametrize("v_cls,v_ncls,ox_cls,red_cls,ox_ncls,red_ncls,ocv,res,k_c,k_n,a_c,a_n,n_c,n_n",
                             [(6, 0.01, 0.01, 0.01, 0.01, 0.01, 1, 1, 1e-3, 1e-3, 0.5, 0.5, 1, 1),
                              (0.001, 0.1, -0.01, 0.01, 0.01, 0.01, 1, 1, 1e-3, 1e-3, 0.5, 0.5, 1, 1),
                              (0.001, 0.1, 0.1, 0.001, -0.01, 0.01, 1, 1, 1e-3, 1e-3, 0.5, 0.5, 1, 1),
                              (0.001, 0, 0.1, 0.001, -0.01, 0.01, 1, 1, 1e-3, 1e-3, 0.5, 0.5, 1, 1),
                              (0.001, 0, 0.1, 0.001, -0.01, 0.01, -1, 0, 1e-3, 1e-3, 0.5, 0.5, 1, 1),
                              (0.001, 0.1, 0.1, 0.001, 0.01, 0.01, 1, 0.1, 1e-3, 1e-3, -0.5, 0.5, 1, 1),
                              (0.001, 0.1, 0.1, 0.001, 0.01, 0.01, 1, 0.1, 1e-3, 1e-3, 0.5, -0.5, 1, 1),
                              (0.001, 0.1, 0.1, 0.001, 0.01, 0.01, 1, 0.1, 1e-3, 1e-3, 0.5, 1.5, 1, 1),
                              (0.001, 0.1, 0.1, 0.001, 0.01, 0.01, 1, 0.1, 1e-3, 1e-3, 1.5, 0.5, 1, 1),
                              (0.001, 0.1, 0.1, 0.001, 0.01, 0.01, 1, 0.1, 1e-3, 1e-3, 0.5, 0.5, 1.5, 1),
                              (0.001, 0.1, 0.1, 0.001, 0.01, 0.01, 1, 0.1, 1e-3, 1e-3, 0.5, 0.5, 1, -1),
                              (0.001, 0.1, 0.1, 0.001, 0.01, 0.01, 0, 0.1, 1e-3, 1e-3, 0.5, 0.5, 1, 2),
                              (0.001, 0.1, 0.1, 0.001, 0.01, 0.01, 1, 0.1, 1e-3, 1e-3, 0.5, 0.5, 1, 0),
                              (0.001, 0.1, 0.1, 0.001, 0.01, 0.01, 0, 0.1, 1e-3, 1e-3, 0.5, 0.5, 2, 3),
                              (0.001, 0.1, 0.1, 0.001, 0.01, 0.01, 1, 0.1, 1e-3, 1e-3, 0.5, 0.5, 2.5, 3),
                              (0.001, 0.1, 0.1, 0.001, 0.01, 0.01, 1, 0.1, 1e-3, 1e-3, 0.5, 0.5, 1.5, 1.5),
                              (0.2, 0.01, 15, 12, 0.01, 0.01, 1.5, 0.1, 1e-3, 1e-3, 0.5, 0.5, 5, 1),
                              (0.2, 0.01, 0.01, 0.01, 0.01, 0.01, 0.0, 0.1, 1e-3, 1e-3, 0.5, 0.5, 1, 1),
                              (0.002, 0.01, 0.01, 0.01, 0.01, 0.01, 1.0, 0.1, 1e-3, 1e-3, 0.5, 0.5, -2, 1),
                              (0.002, 0.01, 0.01, 0.01, 0.01, 0.01, 1.0, 0.1, 1e-3, 1e-3, 0.5, 0.5, 2, -11),
                              ])
    def test_class_init(self, v_cls, v_ncls, ox_cls, red_cls, ox_ncls, red_ncls, ocv, res, k_c, k_n, a_c, a_n, n_c, n_n):
        with pytest.raises(ValueError):
            ZeroDModel(volume_cls=v_cls,
                       volume_ncls=v_ncls,
                       c_ox_cls=ox_cls,
                       c_red_cls=red_cls,
                       c_ox_ncls=ox_ncls,
                       c_red_ncls=red_ncls,
                       ocv_50_soc=ocv,
                       resistance=res,
                       k_0_cls=k_c,
                       k_0_ncls=k_n,
                       alpha_cls=a_c,
                       alpha_ncls=a_n,
                       num_electrons_cls=n_c,
                       num_electrons_ncls=n_n,
                       )

    @pytest.mark.parametrize("v_cls,v_ncls,ox_cls,red_cls,ox_ncls,red_ncls,ocv,res,k_c,k_n,time_i",
                             [(0.01, 0.05, 0.01, 0.01, 0.01, 0.01, 1, 1, 1e-3, 1e-3, 1.0),
                              (0.01, 0.05, 0.01, 0.01, 0.01, 0.01, 1, 1, 1e-3, 1e-3, 5.0),
                              (0.01, 0.05, 0.01, 0.01, 0.01, 0.01, 1, 1, 1e-3, 1e-3, 11.2),
                              ])
    def test_time_step_warning(self,v_cls,v_ncls,ox_cls,red_cls,ox_ncls,red_ncls,ocv,res,k_c,k_n,time_i,capsys):
        ZeroDModel(volume_cls=v_cls,
                   volume_ncls=v_ncls,
                   c_ox_cls=ox_cls,
                   c_red_cls=red_cls,
                   c_ox_ncls=ox_ncls,
                   c_red_ncls=red_ncls,
                   ocv_50_soc=ocv,
                   resistance=res,
                   k_0_cls=k_c,
                   k_0_ncls=k_n,
                   time_step=time_i
                   )

        warn_out = "WARNING: 'time_step' >= 1 second will result in very coarse data.\
                  \nzero-D model approaches theory as time step decreases."
        captured = capsys.readouterr()
        assert captured.out.strip() == warn_out

    def test_exchange_current(self):
        cell = ZeroDModel(volume_cls=0.005, volume_ncls=0.01, c_ox_cls=0.01, c_red_cls=0.01,
                          c_ox_ncls=0.01, c_red_ncls=0.01, ocv_50_soc=1.0, resistance=1, k_0_cls=1e-3,
                          k_0_ncls=1e-3, num_electrons_ncls=2)
        i_0_cls, i_0_ncls = cell._exchange_current()
        assert np.isclose(i_0_cls, 0.12543093175)
        assert np.isclose(i_0_ncls, 0.25086186351)

    def test_limiting_current(self):
        limiting_c = 0.2
        cell = ZeroDModel(volume_cls=0.005, volume_ncls=0.01, c_ox_cls=0.01, c_red_cls=0.01,
                          c_ox_ncls=0.01, c_red_ncls=0.01, ocv_50_soc=1.0, resistance=1, k_0_cls=1e-3,
                          k_0_ncls=1e-3, num_electrons_ncls=2)

        i_lim = cell._limiting_current(limiting_c)
        assert np.isclose(i_lim, 77.1882656)

    def test_limiting_concentration(self):
        cell = ZeroDModel(volume_cls=0.005, volume_ncls=0.01, c_ox_cls=0.01, c_red_cls=0.02,
                          c_ox_ncls=0.02, c_red_ncls=0.01, ocv_50_soc=1.0, resistance=1, k_0_cls=1e-3,
                          k_0_ncls=1e-3, num_electrons_ncls=2)

        i_lim_cls, i_lim_ncls = cell._limiting_concentration(True)
        assert np.isclose(i_lim_cls, 3.85941328)
        assert np.isclose(i_lim_ncls, 7.71882656)

        i_lim_cls, i_lim_ncls = cell._limiting_concentration(False)
        assert np.isclose(i_lim_cls, 7.71882656)
        assert np.isclose(i_lim_ncls, 15.43765312)

    def test_activation_overpotential(self):
        cell = ZeroDModel(volume_cls=0.005, volume_ncls=0.01, c_ox_cls=0.01, c_red_cls=0.01,
                          c_ox_ncls=0.01, c_red_ncls=0.01, ocv_50_soc=1.0, resistance=1, k_0_cls=1e-3,
                          k_0_ncls=1e-3, num_electrons_ncls=2)
        current = 1
        i_0_cls = 0.01
        i_0_ncls = 0.01
        n_activation = cell._activation_overpotential(current, i_0_cls, i_0_ncls)
        assert np.isclose(n_activation, 0.177392243)

    def test_cell_voltage(self):
        ocv = 1.0
        loss = 0.2

        cell_v1 = ZeroDModel._cell_voltage(ocv, loss, True)
        cell_v2 = ZeroDModel._cell_voltage(ocv, loss, False)

        assert cell_v1 == 1.2
        assert cell_v2 == 0.8





