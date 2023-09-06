import pytest
#import numpy as np

#from rfbzero.rfbzero.degradation import (DegradationMechanism, ChemicalDegradation, AutoOxidation, AutoReduction, MultiDegradationMechanism)
from degradation import (DegradationMechanism, ChemicalDegradation, AutoOxidation, AutoReduction, MultiDegradationMechanism) # uncomment here when running test harness


class TestDegradationMechanism:

    def test_abstract_class_init(self):
        with pytest.raises(TypeError):
            DegradationMechanism()


class TestChemicalDegradation:
    @pytest.mark.parametrize("order,constant,spec",
                             [(-1, 0.1, 'red'),
                              ('one', 0.1, 'red'),
                              (1, -0.1, 'red'),
                              (1, 0.1, 'blue')])
    def test_chem_deg_init(self, order, constant, spec):
        with pytest.raises(ValueError):
            ChemicalDegradation(rate_order=order,
                                rate_constant=constant,
                                species=spec)

    def test_chem_deg_degrade(self):
        test_chemdeg = ChemicalDegradation(rate_order=1,
                                           rate_constant=0.1,
                                           species='red')
        c_o, c_r = test_chemdeg.degrade(c_ox=1, c_red=0.5, timestep=0.1)
        assert c_o == 1
        assert c_r == 0.495


class TestAutoOxidation:

    @pytest.mark.parametrize("constant", [0, -1, -0.01])
    def test_auto_ox_init(self, constant):
        with pytest.raises(ValueError):
            AutoOxidation(rate_constant=constant)

    def test_auto_ox_degrade(self):
        test_autoox = AutoOxidation(rate_constant=0.1)
        c_o, c_r = test_autoox.degrade(c_ox=1, c_red=0.5, timestep=0.1)
        assert c_o == 1.005
        assert c_r == 0.495


class TestAutoReduction:

    @pytest.mark.parametrize("constant", [0, -1, -0.01])
    def test_auto_red_init(self, constant):
        with pytest.raises(ValueError):
            AutoReduction(rate_constant=constant)

    def test_auto_red_degrade(self):
        test_autored = AutoReduction(rate_constant=0.1)
        c_o, c_r = test_autored.degrade(c_ox=1, c_red=0.5, timestep=0.1)
        assert c_o == 0.99
        assert c_r == 0.51


class TestMultiDegradationMechanism:

    #@pytest.mark.parametrize("mech", [(AutoOxidation, AutoReduction)])
    def test_multi_deg_init(self):#, mech):
        #with pytest.raises(ValueError):
        #    MultiDegradationMechanism(mechanisms_list=mech)
        raise NotImplementedError
    def test_multi_deg_degrade(self):
        raise NotImplementedError
