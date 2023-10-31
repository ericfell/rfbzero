import pytest
import numpy as np

from rfbzero.crossover import Crossover
#from crossover import Crossover # uncomment here when running test harness


class TestClassCrossover:

    def test_init(self):
        with pytest.raises(ValueError):
            test_cross = Crossover(membrane_constant=0,
                                   permeability_ox=1,
                                   permeability_red=1)

        with pytest.raises(ValueError):
            test_cross = Crossover(membrane_constant=1,
                                   permeability_ox=-1,
                                   permeability_red=1)

        with pytest.raises(ValueError):
            test_cross = Crossover(membrane_constant=1,
                                   permeability_ox=0,
                                   permeability_red=0)

    def test_crossover(self):
        test_cross = Crossover(membrane_constant=1,
                               permeability_ox=2,
                               permeability_red=1)

        expected = (0.502, 0.9995, 1.499, 0.50025, -0.002, 0.0005)
        assert expected == test_cross.crossover(c_ox_cls=0.5, c_red_cls=1, c_ox_ncls=1.5, c_red_ncls=0.5,
                                                               vol_cls=1, vol_ncls=2, timestep=1)

