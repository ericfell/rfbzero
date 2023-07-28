# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 11:46:42 2023

@author: erico
"""

# CROSSOVER FUNCTIONS CALLED BY ZERO-D MAIN CLASS


def crossover_mechanism(conc_ox_CLS, conc_red_CLS, conc_ox_NCLS, conc_red_NCLS,
                        timestep, *args, **kwargs):
    if not args:
        return conc_ox_CLS, conc_red_CLS, conc_ox_NCLS, conc_red_NCLS, 0.0, 0.0
    #
    for func, params in zip(args, kwargs.values()):
        (conc_ox_CLS, conc_red_CLS,
         conc_ox_NCLS, conc_red_NCLS,
         delta_ox, delta_red) = func(conc_ox_CLS, conc_red_CLS,
                                     conc_ox_NCLS, conc_red_NCLS,
                                     timestep, *params)
    return (conc_ox_CLS, conc_red_CLS, 
            conc_ox_NCLS, conc_red_NCLS, delta_ox, delta_red)

######################

# membrane_const is geo Area divide by membrane thickness
def crossover(conc_ox_t_CLS, conc_red_t_CLS, 
                   conc_ox_t_NCLS, conc_red_t_NCLS,
                   timestep, membrane_const=1, P_ox=0.0, P_red=0.0,
                   vol_CLS=1, vol_NCLS=1):

    # skip all this if no crossover is desired
    if (P_ox == 0.0) and (P_red == 0.0):
        return (conc_ox_t_CLS, conc_red_t_CLS,
                conc_ox_t_NCLS, conc_red_t_NCLS, 0.0, 0.0)
    
    # set the volume term based on side that is source of crossing species
    if conc_red_t_CLS < conc_red_t_NCLS:
        volume_red = vol_NCLS*1000 # converts liters to cm^3
    else:
        volume_red = vol_CLS*1000 # converts liters to cm^3
    
    if conc_ox_t_CLS < conc_ox_t_NCLS:
        volume_ox = vol_NCLS*1000 # converts liters to cm^3
    else:
        volume_ox = vol_CLS*1000 # converts liters to cm^3
    
    
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
    
  