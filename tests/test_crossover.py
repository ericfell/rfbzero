import pytest
import numpy as np

from rfbzero.crossover import Crossover


class TestClassCrossover:

    def test_init(self):
        with pytest.raises(ValueError):
            test_cross = Crossover(membrane_thickness=0,
                                   permeability_ox=1,
                                   permeability_red=1)

        with pytest.raises(ValueError):
            test_cross = Crossover(membrane_thickness=1,
                                   permeability_ox=-1,
                                   permeability_red=1)

        with pytest.raises(ValueError):
            test_cross = Crossover(membrane_thickness=1,
                                   permeability_ox=0,
                                   permeability_red=0)

    def test_crossover(self):
        test_cross = Crossover(membrane_thickness=10000,
                               permeability_ox=2,
                               permeability_red=1)

        expected = (0.51, 0.9975, 1.495, 0.50125, -0.01, 0.0025)
        assert expected == test_cross.crossover(
            geometric_area=5.0,
            c_ox_cls=0.5,
            c_red_cls=1,
            c_ox_ncls=1.5,
            c_red_ncls=0.5,
            cls_volume=1,
            ncls_volume=2,
            timestep=1,
        )

