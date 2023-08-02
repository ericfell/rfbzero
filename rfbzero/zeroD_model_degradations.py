
"""
DEGRADATION FUNCTIONS CALLED BY ZeroDmodel CLASS
"""

def degradation_mechanism(conc_ox, conc_red, timestep, *args, **kwargs):
    if not args:
        return conc_ox, conc_red

    for func, params in zip(args, kwargs.values()):
        conc_ox, conc_red = func(conc_ox, conc_red, timestep, *params)
    return conc_ox, conc_red
##############################################################################
##############################################################################


def chemical_degradation(conc_ox_t: float, conc_red_t: float, timestep: float, species: str = 'red',
                         rate_order: int = 1, rate: float = 0.0) -> tuple[float, float]:
    # nth order degradation of n[] --> n[]?
    assert species == 'red' or 'ox', "Options are 'red' or 'ox' "
    
    if rate == 0.0:
        return conc_ox_t, conc_red_t

    if species == 'red':
        concentration_red = conc_red_t - (timestep * rate * (conc_red_t**rate_order))
        return conc_ox_t, concentration_red
    else:
        concentration_ox = conc_ox_t - (timestep * rate * (conc_ox_t**rate_order))
        return concentration_ox, conc_red_t


def auto_oxidation(conc_ox_t: float, conc_red_t: float, timestep: float,
                   rate: float = 0.0) -> tuple[float, float]:
    # assumes first order process: red --> ox
    delta_conc = timestep * rate * conc_red_t

    concentration_red = conc_red_t - delta_conc
    concentration_ox = conc_ox_t + delta_conc
    return concentration_ox, concentration_red


def auto_reduction(conc_ox_t: float, conc_red_t: float, timestep: float,
                   rate: float = 0.0) -> tuple[float, float]:
    # assumes first order process: ox --> red
    delta_conc = timestep * rate * conc_ox_t

    concentration_red = conc_red_t + delta_conc
    concentration_ox = conc_ox_t - delta_conc
    return concentration_ox, concentration_red


"""
def auto_red_test(conc_ox_t, conc_red_t, timestep, extra_ratio=1, rate=0):
    # assumes first order process: ox --> red
    delta_conc = timestep*rate*conc_ox_t*extra_ratio
    concentration_red = conc_red_t + delta_conc
    concentration_ox = conc_ox_t - delta_conc
    return concentration_ox, concentration_red
"""

# def oxygen_oxidation(conc_o2, conc_ox_t, conc_red_t, timestep, rate=0):
    

# def potential_dependent(potential, conc_ox_t, conc_red_t, timestep, rate=0):
    
