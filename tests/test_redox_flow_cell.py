import pytest

from rfbzero.rfbzero.zeroD_model_1e_vs_1e import ZeroDModel

def test_ZeroDModel():
    """
    Tests for ZeroDmodel class methods

    """

    def test_current_direction():
        x = ZeroDModel.current_direction(True)
        y = ZeroDModel.current_direction(False)
        assert x == 1
        assert y == -1


    ######
    test_current_direction()
    print('All test passed!')


######
test_ZeroDModel()
