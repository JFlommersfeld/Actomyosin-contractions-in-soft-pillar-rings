import numpy as np


def discretize_curve(time_points, displacements, step_size):
    discretized_time_points = []
    discretized_displacements = []
    idx = 0
    for time_point, disp in zip(time_points, displacements):
        if time_point >= idx * step_size:
            discretized_time_points.append(time_point)
            discretized_displacements.append(disp)
            idx += 1
    return discretized_time_points, discretized_displacements


def velos_numerical(times, displacements):
    num_of_time_points = len(times)
    velos = [None] * num_of_time_points

    dt = (times[1] - times[0]) * 60
    velos[0] = (displacements[1] - displacements[0]) / dt

    for idx in range(1, num_of_time_points - 1):
        velos[idx] = (displacements[idx + 1] - displacements[idx - 1]) / (2. * dt)
    velos[-1] = (displacements[-1] - displacements[-2]) / dt
    return np.asarray(velos)
