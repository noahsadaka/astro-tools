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


def generate_inertial_data(ecc, a, theta_start=0, theta_end=2*pi, delta_omega=0, IP=math.nan, mu=math.nan):
    """
    Get a vector for the plot of a body in an inertial frame
    Inputs:
        ecc: eccentricity of the resonant body
        a: SMA of the orbit
        theta_start: starting TA for the reference body
        theta_end: ending TA of the reference body
        delta_omega: difference in AOP of the resonant body wrt the reference
                    body

    Outputs:
        r: array of the resonant body in the inertial frame
    """
    if not math.isnan(IP):
        if math.isnan(mu):
            print('Error, must give a mu if using IP to get a')
        else:
            n = 2*pi/IP
            a = (mu/n**2)**(1/3)
    p = a*(1-ecc**2)
    n_pts = 1000
    r = []
    for i, theta in enumerate(np.linspace(theta_start, theta_end, n_pts)):
        r.append(BP2_position_vector_inertial(theta, ecc, p, d_omega=delta_omega))
    r = np.asarray(r)
    return r


def generate_rotating_data(resonance, ecc, mu, a_ref, theta_start,
                           theta_end, delta_omega=0):
    """
    Get a vector for the plot of a body in the relative 2 body problem
    rotating in the Sun-Ref frame, where Ref is some reference body

    Inputs:
        resonance: Resonance ratio wrt the ref body
        ecc: eccentricity of the resonant body
        mu: gravitational parameter of the central body
        a_ref: SMA of the reference body orbit
        theta_start: starting TA for the reference body
        theta_end: ending TA of the reference body
        delta_omega: difference in AOP of the resonant body wrt the reference
                    body

    Outputs:
        rA_r: array of the resonant body in the Sun-Ref frame
    """
    # Calculate orbital params
    IP_ref = 2*pi/sqrt(mu/a_ref**3)
    IP_A = IP_ref * resonance
    n_ref = 2*pi/IP_ref
    n_A = 2*pi/IP_A
    a_A = (mu/n_A**2)**(1/3)
    e_ref = 0
    p_A = a_A*(1-ecc**2)

    # Define points on orbit and set up problem
    n_pts = 2500
    rA_r = []

    # Generate points
    for ind, theta in enumerate(np.linspace(theta_start, theta_end, n_pts)):
        time = BP2_time_from_theta(theta, e_ref, n_ref)
        time = update_time_along_orbit(time, theta, IP_ref)
        rA_r.append(BP2_position_vector_rotating(theta, ecc, p_A, theta,
                    time=time, n=n_A, d_omega=delta_omega))
    rA_r = np.asarray(rA_r)

    return rA_r


def update_time_along_orbit(time, theta, IP):
    num_periods = math.floor(theta/(2*math.pi))
    for i in range(1, num_periods+1):
        if theta > i*2*pi:
            time += IP
    return time


def BP2_position_vector_rotating(theta, e, p, theta_r, d_omega=0,
                                 time=math.nan, n=math.nan):
    """
    Return a position vector in relative x-y coordinates.

        theta: true anomaly [rad] (float)
        e: eccentricity (float)
        p: semi-latus rectum [km] (float)
        theta_r: true anomaly of the body fixed in the relative frame
        d_omega: AoP [rad] (float)
        time: time past periapsis [s] (float)
        n: mean motion [km/s] (float)

        r_x_r, r_y_r: radius positions in relative x and y coordinates [km]
                      (float)
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
    r_mag = p/(1+e*math.cos(theta))  # position vector magnitude [km]
    r_e = r_mag * math.cos(theta+d_omega)  # position vector e-component [km]
    r_p = r_mag * math.sin(theta+d_omega)  # position vector p-component [km]
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
