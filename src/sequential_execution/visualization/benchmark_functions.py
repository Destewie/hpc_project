import numpy as np

# TODO: Generalize in n dimensions
def min_spherical_function(x, y=0):
    return x**2 + y**2

# TODO: Generalize in n dimensions
def max_spherical_function(x, y=0):
    return -(x**2 + y**2)

# TODO: Generalize in n dimensions
def min_rastrigin_function(*args):
    if len(args) == 1:
        x = args[0]
        return 10 + x**2 - 10 * np.cos(2 * np.pi * x)
    elif len(args) == 2:
        x, y = args
        return 20 + x**2 + y**2 - 10 * (np.cos(2 * np.pi * x) + np.cos(2 * np.pi * y))
    else:
        raise TypeError("Expected 1 or 2 arguments.")

# TODO: Generalize in n dimensions
def min_rosenbrock_function(x, y):
    return 100 * (y - x**2)**2 + (1 - x)**2

# TODO: Generalize in n dimensions
def min_ackley_function(x, y):
    return -20 * np.exp(-0.2 * np.sqrt(0.5 * (x**2 + y**2))) - np.exp(0.5 * (np.cos(2 * np.pi * x) + np.cos(2 * np.pi * y))) + np.e + 20

# def min_rastrigin_function(x, A=10):
#     n = len(x)  # Number of dimensions
#     return A * n + sum(xi**2 - A * np.cos(2 * np.pi * xi) for xi in x)