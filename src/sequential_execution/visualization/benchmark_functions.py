import numpy as np

def spherical_function(x, y=0):
    return x**2 + y**2

def rastrigin_function(x):
    return 10 + x**2 - 10 * np.cos(2 * np.pi * x)
