import pytest
import numpy as np

from rfbzero.rfbzero.degradation import (DegradationMechanism, ChemicalDegradation, AutoOxidation, AutoReduction,
                                         MultiDegradationMechanism)
#from degradation import (DegradationMechanism, ChemicalDegradation, AutoOxidation, AutoReduction,
#                                         MultiDegradationMechanism) # uncomment here when running test harness


class TestDegradationMechanism:

    def test_abstract_class_init(self):
        with pytest.raises(NotImplementedError):
            abs_class = DegradationMechanism()
            abs_class.degrade(c_ox=1, c_red=1, timestep=1)


class TestChemicalDegradation:

    def test_chem_deg_init(self):
        raise NotImplementedError

    def test_chem_deg_degrade(self):
        raise NotImplementedError


class TestAutoOxidation:

    def test_auto_ox_init(self):
        raise NotImplementedError

    def test_auto_ox_degrade(self):
        raise NotImplementedError


class TestAutoReduction:

    def test_auto_red_init(self):
        raise NotImplementedError

    def test_auto_red_degrade(self):
        raise NotImplementedError


class TestMultiDegradationMechanism:

    def test_multi_deg_init(self):
        raise NotImplementedError

    def test_multi_deg_degrade(self):
        raise NotImplementedError
