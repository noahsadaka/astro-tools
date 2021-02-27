import numpy as np
import math
from math import sqrt as sqrt
from math import cos as cos
from math import sin as sin
from math import acos as acos
from math import asin as asin
from math import tan as tan
from math import atan as atan
pi = math.pi


def char_quants(mu1, mu2, a):
    """
    Calculate the characteristic quantities of a primary system 
    in the CR3BP

    Inputs:
        mu1: GM of the largest primary body
        mu2: GM of the other primary body
        a: distance between them (typically the SMA)
    Outputs:
        l_star: characteristic length
        m_star: characteristic mass
        t_star: characteristic time
        mu: mass ratio
    """
    G =  6.674e-11/(1000**3)
    mu_star = (mu1+mu2)
    m_star = mu_star/G
    l_star = a
    N = math.sqrt(mu_star/l_star**3)
    t_star = 1/N
    mu = mu1/mu_star
    return [l_star, m_star, t_star, mu]

def lagrange_points(mu, tol=1e-8):
    """
    Returns the x,y coordinates of the Lagrange points in a CR3BP
        system

    Inputs:
        mu: mass ratio of the system
        tol: tolerance of the newton-raphson solver solving each 
            point's location
    Outputs:
        Li: [x,y] coordinates of each point
    """
    L1 = L1_calc(mu, tol)
    L2 = L2_calc(mu, tol)
    L3 = L3_calc(mu, tol)
    L4, L5 = L45_calc(mu)
    return L1, L2, L3, L4, L5

def L1_calc(mu, tol=1e-8):
    """
    Returns the x,y coordinates of the L1 in a CR3BP system

    Inputs:
        mu: mass ratio of the system
        tol: tolerance of the newton-raphson solver solving the 
            point's location
    Outputs:
        Li: [x,y] coordinates of the L1 point
    """
    xn=mu
    curtol = 1
    while curtol > tol:
        x_n1 = xn - ((1-mu)/(1-xn)**2 - mu/xn**2 - 1 + mu + xn)/(2*(1-mu)/(1-xn)**3 + 2*mu/xn**3 + 1)
        curtol = abs(x_n1-xn)
        xn = x_n1
    return [1-mu-xn, 0]

def L2_calc(mu, tol=1e-8):
    """
    Returns the x,y coordinates of the L2 in a CR3BP system

    Inputs:
        mu: mass ratio of the system
        tol: tolerance of the newton-raphson solver solving the 
            point's location
    Outputs:
        Li: [x,y] coordinates of the L2 point
    """
    xn=mu
    curtol = 1
    while curtol > tol:
        x_n1 = xn - ((1-mu)/(xn+1)**2 + mu/xn**2 - 1 + mu - xn)/(-2*(1-mu)/(xn+1)**3 - 2*mu/xn**3 - 1)
        curtol = abs(x_n1-xn)
        xn = x_n1
    return [1-mu+xn, 0]

def L3_calc(mu, tol=1e-8):
    """
    Returns the x,y coordinates of the L3 in a CR3BP system

    Inputs:
        mu: mass ratio of the system
        tol: tolerance of the newton-raphson solver solving the 
            point's location
    Outputs:
        Li: [x,y] coordinates of the L3 point
    """
    curtol = 1
    xn = mu
    while curtol > rol:
        x_n1 = xn - (-(1-mu)/(xn)**2 - (mu)/(xn+1)**2 + xn + mu)/(2*(1-mu)/(xn)**3 + 2*mu/(xn+1)**3+1)
        curtol = abs(x_n1-xn)
        xn = x_n1
    return [-mu-xn, 0]

def L45_calc(mu):
    """
    Returns the x,y coordinates of the L4 and L5 in a CR3BP system

    Inputs:
        mu: mass ratio of the system
    Outputs:
        Li: [x,y] coordinates of the L4 and L5 points
    """
    L4 = [0.5-mu,math.sin(math.radians(60))]
    L5 = [0.5-mu,-math.sin(math.radians(60))]
    return L4, L5

def U_star(d, r, x, y, mu):
    """
    Calculates the pseudopotential energy

    Inputs:
        d: distance from the primary body
        r: distance from the secondary body
        x: x-position in space
        y: y-position in space
        mu: mass ratio

    Outputs:
        U: pseudopotential energy
    """
    U = (1-mu)/(d) + mu/r + 0.5*(x**2+y**2)
    return U
