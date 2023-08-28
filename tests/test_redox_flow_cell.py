import pytest

from rfbzero.redox_flow_cell import ZeroDModel


def test_zerodmodel():
    """
    Tests for ZeroDmodel class methods

    """

    """
    def test_current_direction():
        x = ZeroDModel.current_direction(True)
        y = ZeroDModel.current_direction(False)
        assert x == 1
        assert y == -1

    ######
    test_current_direction()
    """
    def test_cell_voltage():
        ocv = 1.0
        loss = 0.2

        cell_v1 = ZeroDModel.cell_voltage(ocv, loss, True)
        cell_v2 = ZeroDModel.cell_voltage(ocv, loss, False)

        assert cell_v1 == 1.2
        assert cell_v2 == 0.8

    def test_state_of_charge():
        r_cls = 0.2
        o_cls = 0.8
        r_ncls = 0.3
        o_ncls = 0.7

        soc_c, soc_n = ZeroDModel.state_of_charge(o_cls, r_cls, o_ncls, r_ncls)

        assert soc_c == 20.0
        assert soc_n == 30.0

    #
    test_cell_voltage()
    test_state_of_charge()
    print('All test passed!')

######
test_zerodmodel()
