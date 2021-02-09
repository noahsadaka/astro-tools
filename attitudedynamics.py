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


def DCM_from_quaternion(e):
    """
    Converts Euler parameters to the row vector direction cosine C

    Inputs
        e: 4-vector
    Outputs
        C: 3x3 matrix with direction cosine elements
    """
    check_quaternion(e)

    C11 = 1 - 2*e[1]**2 - 2*e[2]**2
    C12 = 2*(e[0]*e[1] - e[2]*e[3])
    C13 = 2*(e[0]*e[2] + e[1]*e[3])
    C21 = 2*(e[0]*e[1] + e[2]*e[3])
    C22 = 1 - 2*e[2]**2 - 2*e[0]**2
    C23 = 2*(e[1]*e[2] - e[0]*e[3])
    C31 = 2*(e[2]*e[0] - e[1]*e[3])
    C32 = 2*(e[1]*e[2] + e[0]*e[3])
    C33 = 1 - 2*e[0]**2 - 2*e[1]**2

    DCM = np.array([[C11, C12, C13], [C21, C22, C23], [C31, C32, C33]])
    check_DCM(DCM)
    return DCM


def quaternion_from_DCM(DCM):
    """
    Converts a DCM into Euler Parameters

    Inputs:
        DCM: a 3x3 direction cosine matrix
    Outputs
        q: a 4-vector quaternion
    """

    check_DCM(DCM)
    E_4 = 0.5*math.sqrt(1 + DCM[0, 0] + DCM[1, 1] + DCM[2, 2])
    E_1 = (DCM[2, 1]-DCM[1, 2])/(4*E_4)
    E_2 = (DCM[0, 2]-DCM[2, 0])/(4*E_4)
    E_3 = (DCM[1, 0]-DCM[0, 1])/(4*E_4)

    q = np.array([E_1, E_2, E_3, E_4])
    check_quaternion(q)
    return q


def quaternion_from_lambda_theta(lam, theta):
    """
    Converts Euler axis/Euler Angle to Euler parameters

    Inputs
        lam: lambda 3-vector
        theta: an angle [rad]
    outputs
        q: 4-vector of the quaternion
    """

    check_lambda(lam)
    q = np.append(lam*sin(theta/2), cos(theta/2))
    check_quaternion(q)
    return q


def lambda_theta_from_quaternion(e):
    """
    Converts Euler parameters to the Euler axis/angle

    Inputs
        e: 4-vector
    Outputs
        theta: an angle [rad]
        lam: a 3-vector
    """
    check_quaternion(e)
    theta = 2*acos(e[3])
    norm = sqrt(e[0]**2 + e[1]**2 + e[2]**2)
    lam = np.array([e[0], e[1], e[2]])/norm
    check_lambda(lam)
    return lam, theta


def check_lambda(lam):
    """
    Checks that lambda is a unit vector
    Inputs
        lam: a 3-vector
    """

    tol = 1 - abs(lam[0]**2 + lam[1]**2 + lam[2]**2)
    if tol > 1e-6:
        print(f'Lambda vector is not a unit vector {lam}')


def check_DCM(DCM):
    """
    Checks that the DCM satisfies the orthogonality conditions
    Inputs
        DCM: 3x3 matrix
    """
    badflag = 0
    for i in range(3):
        tol = 1-abs(sqrt(DCM[i, 0]**2 + DCM[i, 1]**2 + DCM[i, 2]**2))
        if tol > 1e-6:
            badflag = 1
    if badflag == 1:
        print(f'DCM does not satisfy orthogonality conditions {DCM}')


def check_quaternion(e):
    """
    Checks that the quaternion satisfies the constraint equation
    Inputs
        e: 4-vector
    """
    summation = e[0]**2 + e[1]**2 + e[2]**2 + e[3]**2
    tol = 1 - abs(summation)
    if tol > 1e-6:
        print(f'Quaternion constraint not met for e = {e}')
