import pytest
import numpy as np

from rfbzero.degradation import (DegradationMechanism, ChemicalDegradationOxidized, ChemicalDegradationReduced,
                                 AutoOxidation, AutoReduction, Dimerization, MultiDegradationMechanism)


class TestDegradationMechanism:

    def test_abstract_class_init(self):
        with pytest.raises(TypeError):
            DegradationMechanism()


class TestChemicalDegradation:
    @pytest.mark.parametrize("order,constant",
                             [(-1, 0.1),
                              (-1, -0.1),
                              (1, -0.1,)])
    def test_chem_deg_init(self, order, constant):
        with pytest.raises(ValueError):
            ChemicalDegradationReduced(rate_order=order, rate_constant=constant)

        with pytest.raises(ValueError):
            ChemicalDegradationOxidized(rate_order=order, rate_constant=constant)

    def test_chem_deg_degrade(self):
        test_chemdeg = ChemicalDegradationReduced(rate_order=1, rate_constant=0.1)
        delta_c_o, delta_c_r = test_chemdeg.degrade(c_ox=1, c_red=0.5, time_step=0.1)
        assert np.isclose(delta_c_o, 0.0)
        assert np.isclose(delta_c_r, -0.005)

        test_chemdeg2 = ChemicalDegradationOxidized(rate_order=1, rate_constant=0.1)
        delta_c_ox, delta_c_red = test_chemdeg2.degrade(c_ox=1, c_red=0.5, time_step=0.1)
        assert np.isclose(delta_c_ox, -0.01)
        assert np.isclose(delta_c_red, 0.0)


class TestAutoOxidation:

    @pytest.mark.parametrize("constant", [0, -1, -0.01])
    def test_auto_ox_init(self, constant):
        with pytest.raises(ValueError):
            AutoOxidation(rate_constant=constant)

    @pytest.mark.parametrize("c_oxid, stoich", [(0, 1), (-1, -1), (2, -0.01)])
    def test_auto_ox_stoich(self, c_oxid, stoich):
        with pytest.raises(ValueError):
            AutoOxidation(rate_constant=0.1, c_oxidant=c_oxid, oxidant_stoich=stoich)

    def test_auto_ox_degrade(self):
        test_autoox = AutoOxidation(rate_constant=0.1)
        delta_c_o, delta_c_r = test_autoox.degrade(c_ox=1, c_red=0.5, time_step=0.1)
        assert np.isclose(delta_c_o, 0.005)
        assert np.isclose(delta_c_r, -0.005)


class TestAutoReduction:

    @pytest.mark.parametrize("constant", [0, -1, -0.01])
    def test_auto_red_init(self, constant):
        with pytest.raises(ValueError):
            AutoReduction(rate_constant=constant)

    @pytest.mark.parametrize("c_reduct, stoich", [(0, 5), (-11, -1), (32, -0.61)])
    def test_auto_red_stoich(self, c_reduct, stoich):
        with pytest.raises(ValueError):
            AutoReduction(rate_constant=0.1, c_reductant=c_reduct, reductant_stoich=stoich)

    def test_auto_red_degrade(self):
        test_autored = AutoReduction(rate_constant=0.1)
        delta_c_o, delta_c_r = test_autored.degrade(c_ox=1, c_red=0.5, time_step=0.1)
        assert np.isclose(delta_c_o, -0.01)
        assert np.isclose(delta_c_r, 0.01)


class TestDimerization:

    @pytest.mark.parametrize("k_f,k_b,c_dim", [(-1,-1,-1), (1,-1,1), (4,11,-1)])
    def test_dimerization_init(self, k_f, k_b, c_dim):
        with pytest.raises(ValueError):
            Dimerization(forward_rate_constant=k_f, backward_rate_constant=k_b, c_dimer=c_dim)

    def test_dimerization_degrade(self):
        test_dimerize = Dimerization(forward_rate_constant=2, backward_rate_constant=1, c_dimer=0.5)
        delta_c_o, delta_c_r = test_dimerize.degrade(c_ox=1, c_red=0.5, time_step=1)
        assert np.isclose(delta_c_o, -0.5)
        assert np.isclose(delta_c_r, -0.5)


class TestMultiDegradationMechanism:

    @pytest.mark.parametrize("mech", [(1, True)])
    def test_multi_deg_init(self, mech):
        with pytest.raises(ValueError):
            MultiDegradationMechanism(mechanisms=mech)

    def test_multi_deg_degrade(self):
        mech_1 = AutoReduction(rate_constant=0.1)
        mech_2 = AutoOxidation(rate_constant=0.1)
        mech_list = [mech_1, mech_2]
        test_multi = MultiDegradationMechanism(mechanisms=mech_list)
        delta_c_o, delta_c_r = test_multi.degrade(c_ox=1, c_red=0.5, time_step=0.1)
        assert np.isclose(delta_c_o, -0.005)
        assert np.isclose(delta_c_r, 0.005)
