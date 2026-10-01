import numpy as np


def power(x, A, theta):
    return A * x**theta


def power_log(x, A, theta, beta):
    return A * x**theta * np.log(x)**beta


def two_power(x, A, theta, B, phi):
    return A * x**theta + B * x**phi


def five_sixths(x, A):
    return A * x**(5 / 6)


def five_sixths_plus_two_thirds(x, A, B):
    return A * x**(5 / 6) + B * x**(2 / 3)


# Add new conjectural asymptotic functions below.
#
# Example:
#
# def isaac_model(x, A, B):
#     return A * x**(2/3) * np.log(x) + B * x**0.5
