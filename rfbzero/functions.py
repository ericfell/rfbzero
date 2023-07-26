"""
Function used in zero dimensional model

"""


def experiment_time(duration, time_step):
    steps = int((duration / time_step)) + 1
    times = [0.0 + x * time_step for x in range(steps)]
    return times