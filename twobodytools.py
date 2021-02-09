import numpy as np
import matplotlib.pyplot as plt
import math
pi = math.pi
from math import sqrt as sqrt
from math import cos as cos
from math import sin as sin
from math import acos as acos
from math import asin as asin
from math import tan as tan
from math import atan as atan

def BP2_position_vector_rotating(theta, e, p, theta_r, d_omega=0, time=math.nan, n=math.nan):
    """
    Return a position vector in relative x-y coordinates.

        theta: true anomaly [rad] (float)
        e: eccentricity (float)
        p: semi-latus rectum [km] (float)
        theta_r: true anomaly of the body fixed in the relative frame
        d_omega: AoP [rad] (float)
        time: time past periapsis [s] (float)
        n: mean motion [km/s] (float)

        r_x_r, r_y_r: radius positions in relative x and y coordinates [km] (float)
    """
    if not math.isnan(time):
        if math.isnan(n):
            print('ERROR: Must provide a mean motion if giving a time')
            exit()
        theta = BP2_theta_from_time(time, e, n)
    r_mag = p/(1+e*math.cos(theta))
    r_e = r_mag * math.cos(theta+d_omega) 
    r_p = r_mag * math.sin(theta+d_omega) 
    r_x_r = cos(theta_r)*r_e + sin(theta_r)*r_p 
    r_y_r = -sin(theta_r)*r_e + cos(theta_r)*r_p 

    return r_x_r, r_y_r


def BP2_position_vector_inertial(theta, e, p, d_omega=0, time=math.nan, n=0):
    """
    Return a position vector in e and p components.

        theta: true anomaly [rad] (float)
        e: eccentricity (float)
        p: semi-latus rectum [km] (float)
        d_omega: AoP [rad] (float)
        time: time past periapsis [s] (float)
        n: mean motion [km/s] (float)

        r_e, r_p: radius positions in e and p coordinates [km] (float)
    """
    if not math.isnan(time):
        theta = BP2_theta_from_time(time, e, n)
    r_mag = p/(1+e*math.cos(theta)) # position vector magnitude [km]
    r_e = r_mag * math.cos(theta+d_omega) # position vector e-component [km]
    r_p = r_mag * math.sin(theta+d_omega) # position vector p-component [km]
    return r_e, r_p

def BP2_theta_from_time(time, e, n):
    """
    Get a true anomaly from time past periapsis

        time: time past periapsis [s] (float)
        e: eccentricity (float)
        n: mean motion [km/s] (float)

        theta: true anomaly [rad] (float)
    """

    M = n*time
    theta = keplereqn(e, M)
    return theta

def BP2_time_from_theta(theta, e, n):
    """
    Get the time past periapsis when at a certain true anomaly

    theta: true anomaly [rad] (float)
    e: eccentricity (float)
    n: mean motion [km/s] (float)
    """

    E = 2*atan(tan(theta/2) * sqrt((1-e)/(1+e)))
    time = 1/n * (E - e*sin(E))
    if time < 0:
        time = 2*pi/n - abs(time)
    return time

def keplereqn(e, M):
    En = M
    d_E = 100
    while d_E > 1e-6:
        En1 = En - (En - e*sin(En) - M)/(1 - e*cos(En))
        d_E = abs(En1-En)
        En = En1
    theta = 2*atan(sqrt((1+e)/(1-e)) * tan(En/2))
    return theta


