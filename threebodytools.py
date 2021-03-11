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
    mu = mu2/mu_star
    return [l_star, m_star, t_star, mu]

def lagrange_points(mu, tol=1e-8):
    """
    Returns the x,y coordinates of the Lagrange points in a CR3BP
        system as well as their JC

    Inputs:
        mu: mass ratio of the system
        tol: tolerance of the newton-raphson solver solving each 
            point's location
    Outputs:
        Li: ([x,y], JC) coordinates of each point and the JC
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
    coord = [1-mu-xn, 0]
    JC = coord[0]**2 + coord[1]**2 + 2*(1-mu)/np.linalg.norm(get_d(coord[0], coord[1], mu)) + 2*mu/np.linalg.norm(get_r(coord[0], coord[1], mu))
    return (coord, JC)

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
    coord = [1-mu+xn, 0]
    JC = coord[0]**2 + coord[1]**2 + 2*(1-mu)/np.linalg.norm(get_d(coord[0], coord[1], mu)) + 2*mu/np.linalg.norm(get_r(coord[0], coord[1], mu))
    return (coord, JC)

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
    while curtol > tol:
        x_n1 = xn - (-(1-mu)/(xn)**2 - (mu)/(xn+1)**2 + xn + mu)/(2*(1-mu)/(xn)**3 + 2*mu/(xn+1)**3+1)
        curtol = abs(x_n1-xn)
        xn = x_n1
    coord = [-mu-xn, 0]
    JC = coord[0]**2 + coord[1]**2 + 2*(1-mu)/np.linalg.norm(get_d(coord[0], coord[1], mu)) + 2*mu/np.linalg.norm(get_r(coord[0], coord[1], mu))
    return (coord, JC)

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

    JC = L4[0]**2 + L4[1]**2 + 2*(1-mu)/np.linalg.norm(get_d(L4[0], L4[1], mu)) + 2*mu/np.linalg.norm(get_r(L4[0], L4[1], mu))
    
    return (L4, JC), (L5, JC)

def U_star(x, y, mu):
    """
    Calculates the pseudopotential energy

    Inputs:
        x: x-position in space
        y: y-position in space
        mu: mass ratio

    Outputs:
        U: pseudopotential energy
    """
    d = np.linalg.norm(get_d(x, y, mu))
    r = np.linalg.norm(get_r(x, y, mu))
    U = (1-mu)/(d) + mu/r + 0.5*(x**2+y**2)
    return U

def get_r(x, y, mu, z=0):
    """
    Returns the r vector

    Inputs:
        x: x-position in nd units
        y: y-position in nd units
        z: z-position in nd units (default: 0)
        mu: mass ratio

    Outputs
        r: vector
    """
    return np.array([x-1+mu, y, z])


def get_d(x, y, mu, z=0):
    """
    Returns the d vector

    Inputs:
        x: x-position in nd units
        y: y-position in nd units
        z: z-position in nd units (default: 0)
        mu: mass ratio

    Outputs
        d: vector
    """
    return np.array([x+mu, y, z])

# ZVC STUFF

def ZVC_iterator(JC, mu):
    Li = lagrange_points(mu)
    if JC > Li[0][1]:
         x_vals, y_vals, slope_val = ZVC_L1(Li, JC, mu)
    elif JC <= Li[0][1] and JC > Li[1][1]:
         x_vals, y_vals, slope_val = ZVC_L1_L2(Li, JC, mu)
    elif JC <= Li[1][1] and JC > Li[2][1]:
         x_vals, y_vals, slope_val = ZVC_L2_L3(Li, JC, mu)
    elif JC <=Li[2][1] and JC > Li[3][1]:
        x_vals, y_vals, slope_val = ZVC_L3_L45(Li, JC, mu)
    elif JC <= Li[3][1]:
        print('smaller than L4 and L5, out of plane')
        x_vals = np.array([])
        y_vals = np.array([])
        slope_val = np.array([])
    else:
        print("I wasn't expecting this!")
    return x_vals, y_vals, slope_val


def sloper(slope, delta_min = 0.00001, delta_max=0.001):
    # note: delta_min was originally 0.0001 and that worked for everything up to L3
    slope_limit=15 # originally 3
    #print(abs(slope))
    if abs(slope) > slope_limit:
        return delta_min
    else:
        return delta_max - abs(slope)*(delta_max-delta_min)/slope_limit
        
def x_step(x, guess, JC, mu):
    yn = guess
    tol = 1
    it=0
    while tol >1e-12:
        yn_1 = yn - (x**2 + yn**2 + 2*(1-mu)/math.sqrt((x+mu)**2 + yn**2) + 2*mu/math.sqrt((x-1+mu)**2 + yn**2) - JC)\
           /(2*yn - 2*yn*(1-mu)/((x+mu)**2 + yn**2)**(3/2) - 2*yn*mu/((x-1+mu)**2 + yn**2)**(3/2))
        tol = abs(yn_1-yn)
        yn=yn_1
        it+=1
        if it>40:
            return np.nan
    return yn

def y_step(y, guess, JC, mu):
    xn = guess
    tol = 1
    it=0
    while tol >1e-12:
        xn_1 = xn - (xn**2 + y**2 + 2*(1-mu)/math.sqrt((xn+mu)**2 + y**2) + 2*mu/math.sqrt((xn-1+mu)**2 + y**2) - JC)\
            /(2*xn - 2*(xn+mu)*(1-mu)/((xn+mu)**2 + y**2)**(3/2) - 2*(xn-1+mu)*mu/((xn-1+mu)**2 + y**2)**(3/2))
        tol = abs(xn_1-xn)
        xn=xn_1
        it+=1
        if it>40:
            return np.nan
    return xn   


def ZVC_L1(Li, JC, mu):
    slope_val = []
    delta = 0.00000001
    guess = 0.01
    
    # right loop, left
    x_vals = np.array([1-mu])
    y_vals = np.array([abs(x_step(x_vals[0], guess, JC, mu))])
    x_vals = np.append(x_vals, 1 - mu - delta)
    y_vals = np.append(y_vals, abs(x_step(x_vals[-1], guess, JC, mu)))
    x_v, y_v, s_v = curve_walker(mu, JC, x_vals, y_vals, slope_val)
    x_vals = np.append(x_vals, x_v)
    y_vals = np.append(y_vals, y_v)
    slope_val.extend(s_v)
    
    # right loop, right
    delta = 0.00000001
    x_vals = np.append(x_vals, 1-mu + delta)
    y_vals =  np.append(y_vals, abs(x_step(x_vals[-1], guess, JC, mu)))
    x_v, y_v, s_v = curve_walker(mu, JC, x_vals, y_vals, slope_val, sign=-1, slope_lim=200)
    x_vals = np.append(x_vals, x_v)
    y_vals = np.append(y_vals, y_v)
    slope_val.extend(s_v)

    # Left loop, right
    delta = 0.001
    guess = 0.1
    x_vals = np.append(x_vals, mu + delta)
    y_vals =  np.append(y_vals, abs(x_step(x_vals[-1], guess, JC, mu)))
    x_vals = np.append(x_vals, mu + delta+delta)
    y_vals =  np.append(y_vals, abs(x_step(x_vals[-1], guess, JC, mu)))
    x_v, y_v, s_v = curve_walker(mu, JC, x_vals, y_vals, slope_val, sign=-1, slope_lim=20)
    x_vals = np.append(x_vals, x_v)
    y_vals = np.append(y_vals, y_v)
    slope_val.extend(s_v)
 
    # Left loop, left
    delta = 0.001
    guess = 0.1
    x_vals = np.append(x_vals, mu - delta)
    y_vals =  np.append(y_vals, abs(x_step(x_vals[-1], guess, JC, mu)))
    x_vals = np.append(x_vals, mu - 2*delta)
    y_vals =  np.append(y_vals, abs(x_step(x_vals[-1], guess, JC, mu)))
    x_v, y_v, s_v = curve_walker(mu, JC, x_vals, y_vals, slope_val, slope_lim=20)
    x_vals = np.append(x_vals, x_v)
    y_vals = np.append(y_vals, y_v)
    slope_val.extend(s_v)

    # Outer loop, left pick left of L3 as starting point
    nanswitch=False
    delta = 0.001
    guess = 1
    x_vals = np.append(x_vals, Li[2][0][0])
    y_vals =  np.append(y_vals, abs(x_step(x_vals[-1], guess, JC, mu)))
    x_vals = np.append(x_vals, Li[2][0][0] - delta)
    y_vals =  np.append(y_vals, abs(x_step(x_vals[-1], guess, JC, mu)))
    x_v, y_v, s_v = curve_walker(mu, JC, x_vals, y_vals, slope_val)
    x_vals = np.append(x_vals, x_v)
    y_vals = np.append(y_vals, y_v)
    slope_val.extend(s_v)
 
    # outer loop, right
    delta = 0.001
    guess = 1
    x_vals = np.append(x_vals, Li[2][0][0])
    y_vals =  np.append(y_vals, abs(x_step(x_vals[-1], guess, JC, mu)))
    x_vals = np.append(x_vals, Li[2][0][0] + delta)
    y_vals =  np.append(y_vals, abs(x_step(x_vals[-1], guess, JC, mu)))
    x_v, y_v, s_v = curve_walker(mu, JC, x_vals, y_vals, slope_val, sign=-1, slope_lim=20)
    x_vals = np.append(x_vals, x_v)
    y_vals = np.append(y_vals, y_v)
    slope_val.extend(s_v)
    return x_vals, y_vals, slope_val


def curve_walker(mu, JC, x_vals, y_vals, slope_val, sign=1, slope_lim=50):
    nanswitch = 0
    cusp = False
    while nanswitch==False and y_vals[-1] >= 0:
        slope = (y_vals[-2]-y_vals[-1])/(x_vals[-2]-x_vals[-1])
        slope_val.append(slope)
        delta = sloper(slope)
        if abs(slope)<slope_lim and cusp == False:
            x_n = x_vals[-1]-sign*delta
            y_val = x_step(x_n, y_vals[-1], JC, mu)
            if np.isnan(y_val) == False:
                x_vals = np.append(x_vals, x_n)
                y_vals = np.append(y_vals, y_val)
            else:
                nanswitch = True
        else:
            cusp = True
            y_n = y_vals[-1]-delta
            x_val = y_step(y_n, x_vals[-1], JC, mu)
            if np.isnan(x_val) == False:
                x_vals = np.append(x_vals, x_val)
                y_vals = np.append(y_vals, y_n)
            else:
                nanswitch = True
    return x_vals, y_vals, slope_val

def ZVC_L1_L2(Li, JC, mu):
    
    # Inner loop, left
    slope_val = []
    delta = 0.0001
    guess = 0.01
    x_vals = np.array([Li[0][0][0]])
    y_vals = np.array([abs(x_step(x_vals[0], guess, JC, mu))])
    x_vals = np.append(x_vals, Li[0][0][0] - delta)
    y_vals = np.append(y_vals, abs(x_step(x_vals[-1], guess, JC, mu)))
    x_v, y_v, s_v = curve_walker(mu, JC, x_vals, y_vals, slope_val)
    x_vals = np.append(x_vals, x_v)
    y_vals = np.append(y_vals, y_v)
    slope_val.extend(s_v)
    
    # Inner loop, right
    delta = 0.0001
    x_vals = np.append(x_vals, Li[0][0][0] + delta)
    y_vals =  np.append(y_vals, abs(x_step(x_vals[-1], guess, JC, mu)))
    x_v, y_v, s_v = curve_walker(mu, JC, x_vals, y_vals, slope_val, sign=-1, slope_lim=10)
    x_vals = np.append(x_vals, x_v)
    y_vals = np.append(y_vals, y_v)
    slope_val.extend(s_v)
    
    # outer loop, left
    delta = 0.00001
    x_vals = np.append(x_vals, Li[2][0][0]-delta)
    y_vals = np.append(y_vals, abs(x_step(x_vals[-1], 0.5, JC, mu)))
    x_v, y_v, s_v = curve_walker(mu, JC, x_vals, y_vals, slope_val, slope_lim=200)
    x_vals = np.append(x_vals, x_v)
    y_vals = np.append(y_vals, y_v)
    slope_val.extend(s_v)
   
    ## outer loop, right
    delta = 0.001
    x_vals = np.append(x_vals, Li[2][0][0]+delta)
    y_vals = np.append(y_vals, abs(x_step(x_vals[-1], 0.5, JC, mu)))
    x_v, y_v, s_v = curve_walker(mu, JC, x_vals, y_vals, slope_val, sign=-1, slope_lim=20)
    x_vals = np.append(x_vals, x_v)
    y_vals = np.append(y_vals, y_v)
    slope_val.extend(s_v)

    return x_vals, y_vals, slope_val
 

def ZVC_L2_L3(Li, JC, mu):
    slope_val = []
    delta = 0.00000001
    guess = 0.01
    # moving right and under from L3
    x_vals = np.array([Li[2][0][0]])
    y_vals = np.array([abs(x_step(x_vals[0], guess, JC, mu))])
    x_vals = np.append(x_vals, Li[2][0][0] + delta)
    y_vals = np.append(y_vals, abs(x_step(x_vals[-1], guess, JC, mu)))
    nanswitch = 0
    cusp = False
    while nanswitch==False and y_vals[-1] >= 0:
        slope = (y_vals[-2]-y_vals[-1])/(x_vals[-2]-x_vals[-1])
        slope_val.append(slope)
        delta = sloper(slope)
        if abs(slope)<50 and cusp == False: # originally 30
            x_n = x_vals[-1]+delta
            y_val = x_step(x_n, y_vals[-1], JC, mu)
            if np.isnan(y_val) == False:
                x_vals = np.append(x_vals, x_n)
                y_vals = np.append(y_vals, y_val)
            else:
                nanswitch = True
        else:
            cusp=True
            if abs(slope) > 5:
                y_n = y_vals[-1]-delta
                x_val = y_step(y_n, x_vals[-1], JC, mu)
                if np.isnan(x_val) == False:
                    x_vals = np.append(x_vals, x_val)
                    y_vals = np.append(y_vals, y_n)
                else:
                    nanswitch = True
            else:
                x_n = x_vals[-1]-delta
                y_val = x_step(x_n, y_vals[-1], JC, mu)
                if np.isnan(y_val) == False:
                    x_vals = np.append(x_vals, x_n)
                    y_vals = np.append(y_vals, y_val)
                else:
                    nanswitch = True
    
    # Moving left from L3
    delta = 0.00000001
    x_vals = np.append(x_vals, Li[2][0][0] - delta)
    y_vals =  np.append(y_vals, abs(x_step(x_vals[-1], guess, JC, mu)))
    x_v, y_v, s_v = curve_walker(mu, JC, x_vals, y_vals, slope_val)
    x_vals = np.append(x_vals, x_v)
    y_vals = np.append(y_vals, y_v)
    slope_val.extend(s_v)
    return x_vals, y_vals, slope_val


def ZVC_L3_L45(Li, JC, mu):
    slope_val = []
    delta = 0.00000001
    guess = Li[3][0][1]+0.000000001
    # upper left
    x_vals = np.array([Li[3][0][0]])
    y_vals = np.array([abs(x_step(x_vals[0], guess, JC, mu))])
    x_vals = np.append(x_vals, Li[3][0][0] - delta)
    y_vals = np.append(y_vals, abs(x_step(x_vals[-1], guess+0.001, JC, mu)))
    x_v, y_v, s_v = curve_walker(mu, JC, x_vals, y_vals, slope_val)
    x_vals = np.append(x_vals, x_v)
    y_vals = np.append(y_vals, y_v)
    slope_val.extend(s_v)
    
    # upper right and underneath
    delta = 0.00000001
    x_vals = np.append(x_vals, Li[3][0][0] + delta)
    y_vals =  np.append(y_vals, abs(x_step(x_vals[-1], guess, JC, mu)))
    nanswitch = 0
    cusp = False
    while nanswitch==False and y_vals[-1] >= 0:
        slope = (y_vals[-2]-y_vals[-1])/(x_vals[-2]-x_vals[-1])
        slope_val.append(slope)
        delta = sloper(slope)
        if abs(slope)<2 and cusp == False:
            x_n = x_vals[-1]+delta
            y_val = x_step(x_n, y_vals[-1], JC, mu)
            if np.isnan(y_val) == False:
                x_vals = np.append(x_vals, x_n)
                y_vals = np.append(y_vals, y_val)
            else:
                nanswitch = True
        else:
            cusp=True
            if abs(slope) >= 2:
                y_n = y_vals[-1]-delta
                x_val = y_step(y_n, x_vals[-1], JC, mu)
                if np.isnan(x_val) == False:
                    x_vals = np.append(x_vals, x_val)
                    y_vals = np.append(y_vals, y_n)
                else:
                    nanswitch = True
            else:
                x_n = x_vals[-1]-delta
                y_val = x_step(x_n, y_vals[-1], JC, mu)
                if np.isnan(y_val) == False:
                    x_vals = np.append(x_vals, x_n)
                    y_vals = np.append(y_vals, y_val)
                else:
                    nanswitch = True
    return x_vals, y_vals, slope_val   
