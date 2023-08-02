
"""
CROSSOVER FUNCTIONS CALLED BY ZeroDmodel CLASS
"""

# membrane_const is geo Area divide by membrane thickness
# could change that so just thickness is input, geo area already in model

def crossover(conc_ox_t_CLS: float, conc_red_t_CLS: float,
              conc_ox_t_NCLS: float, conc_red_t_NCLS: float,
              timestep: float, vol_CLS: float, vol_NCLS: float,
              membrane_const: float = 1.0, P_ox: float = 0.0,
              P_red: float = 0.0) -> tuple[float, float, float, float, float, float]:
    # skip all this if no crossover is desired
    if (P_ox == 0.0) and (P_red == 0.0):
        return (conc_ox_t_CLS, conc_red_t_CLS,
                conc_ox_t_NCLS, conc_red_t_NCLS, 0.0, 0.0)

    # sets volume term based on side that is source of crossing species
    if conc_red_t_CLS < conc_red_t_NCLS:
        volume_red = vol_NCLS * 1000  # converts liters to cm^3
    else:
        volume_red = vol_CLS * 1000  # converts liters to cm^3

    if conc_ox_t_CLS < conc_ox_t_NCLS:
        volume_ox = vol_NCLS * 1000  # converts liters to cm^3
    else:
        volume_ox = vol_CLS * 1000  # converts liters to cm^3

    # amount added/subtracted from concentrations, units of deltas are mol/L 
    delta_red = (timestep * (P_red * membrane_const / volume_red)
                 * (conc_red_t_CLS - conc_red_t_NCLS))
    delta_ox = (timestep * (P_ox * membrane_const / volume_ox)
                * (conc_ox_t_CLS - conc_ox_t_NCLS))

    # update concentrations based on permeation
    concentration_red_CLS = conc_red_t_CLS - delta_red
    concentration_red_NCLS = conc_red_t_NCLS + delta_red

    concentration_ox_CLS = conc_ox_t_CLS - delta_ox
    concentration_ox_NCLS = conc_ox_t_NCLS + delta_ox

    return (concentration_ox_CLS, concentration_red_CLS,
            concentration_ox_NCLS, concentration_red_NCLS,
            delta_ox, delta_red)

