"""
Function used in zero dimensional model

"""


def experiment_time(duration, time_step):
    steps = int((duration / time_step)) + 1
    times = [0.0 + x * time_step for x in range(steps)]
    return times


#def coulomb_counter(self, current: float, volume_CLS: float, volume_NCLS: float, conc_0_ox_CLS: float,
#                    conc_0_red_CLS: float, conc_0_ox_NCLS: float, conc_0_red_NCLS: float, CLS: bool) -> tuple[float, float]: