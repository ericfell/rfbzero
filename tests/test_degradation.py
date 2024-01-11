import pytest
#import numpy as np

from rfbzero.degradation import (DegradationMechanism, ChemicalDegradation, AutoOxidation, AutoReduction,
                                 MultiDegradationMechanism, AutoReductionO2Release, Dimerization)


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

        test_chemdeg2 = ChemicalDegradation(rate_order=1,
                                            rate_constant=0.1,
                                            species='ox')
        c_ox, c_red = test_chemdeg2.degrade(c_ox=1, c_red=0.5, timestep=0.1)
        assert c_ox == 0.99
        assert c_red == 0.5


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


class TestAutoReductionO2Release:

    @pytest.mark.parametrize("constant,factor", [(-1, -1), (0, -1), (-0.01, 11)])
    def test_auto_red_o2_init(self, constant, factor):
        with pytest.raises(ValueError):
            AutoReductionO2Release(rate_constant=constant, release_factor=factor)

    def test_auto_red_o2_degrade(self):
        test_autoredo2 = AutoReductionO2Release(rate_constant=0.1, release_factor=0.01)
        c_o, c_r = test_autoredo2.degrade(c_ox=1, c_red=0.5, timestep=0.1)
        assert c_o == 0.99
        assert c_r == 0.51
        assert test_autoredo2.rate_constant == 0.0999


class TestMultiDegradationMechanism:

    @pytest.mark.parametrize("mech", [(1, True)])
    def test_multi_deg_init(self, mech):
        with pytest.raises(ValueError):
            MultiDegradationMechanism(mechanisms=mech)
        #raise NotImplementedError
    def test_multi_deg_degrade(self):
        mech_1 = AutoReduction(rate_constant=0.1)
        mech_2 = AutoOxidation(rate_constant=0.1)
        mech_list = [mech_1, mech_2]
        test_multi = MultiDegradationMechanism(mechanisms=mech_list)
        c_o, c_r = test_multi.degrade(c_ox=1, c_red=0.5, timestep=0.1)
        assert c_o == 0.9951
        assert c_r == 0.5049


class TestDimerization:

    @pytest.mark.parametrize("k_f,k_b,c_dim", [(-1,-1,-1), (1,-1,1), (4,11,-1)])
    def test_dimerization_init(self, k_f, k_b, c_dim):
        with pytest.raises(ValueError):
            Dimerization(forward_rate_constant=k_f, backward_rate_constant=k_b, c_dimer=c_dim)

    def test_dimerization_degrade(self):
        test_dimerize = Dimerization(forward_rate_constant=2, backward_rate_constant=1, c_dimer=0.5)
        c_o, c_r = test_dimerize.degrade(c_ox=1, c_red=0.5, timestep=1)
        assert c_o == 0.5
        assert c_r == 0.0
