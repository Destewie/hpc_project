import numpy as np

def min_spherical_function(x, y=0):
    return x**2 + y**2

def max_spherical_function(x, y=0):
    return -(x**2 + y**2)

def rastrigin_function(*args):
    if len(args) == 1:
        x = args[0]
        return 10 + x**2 - 10 * np.cos(2 * np.pi * x)
    elif len(args) == 2:
        x, y = args
        return 20 + x**2 + y**2 - 10 * (np.cos(2 * np.pi * x) + np.cos(2 * np.pi * y))
    else:
        raise TypeError("Expected 1 or 2 arguments.")
