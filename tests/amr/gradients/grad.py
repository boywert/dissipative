import numpy as np

theta = np.deg2rad(57)

def analytic_rho(x, y):
    return np.cos(theta)*x + np.sin(theta)*y

def analytic_grad_rho(x, y):
    return (
            np.cos(theta)*np.ones_like(x),
            np.sin(theta)*np.ones_like(x),
            np.zeros_like(x) )
