# -*- coding: utf-8 -*-
"""
Created on 12/02/2021

Analytical solution for Richards' equatuion 1D


@author: Niccolo` Tubini, Riccardo Rigon
@license: this work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License

Reference: Srivastava, R., & Yeh, T. C. J. (1991). Analytical solutions for one‐dimensional, transient infiltration toward the water table in homogeneous and layered soils. Water Resources Research, 27(5), 753-762.

"""

import numpy as np

def solution_homogeneous_soil(L_star, z_star, Ks, alpha, theta_s, theta_r, qA_star, qB_star, psi0, times):
    """ Compute the analytical solution for the case of homogeneous soil
    
    [1] Srivastava, R., & Yeh, T. C. J. (1991). Analytical solutions for one‐dimensional, transient   infiltration toward the water table in homogeneous and layered soils. Water Resources Research, 27(5), 753-762.
    

    Parameters
    ----------
    L_star: float 
            dimensional depth to water table 
    z_star: numpy.array (float) 
            dimensional coordinate z 
    KS: float
            saturated hydraulic conductivity
    alpha: double
            is a soil pore-size distribution parameter representing the rate of reduction in hydraulic conductivity or moisture
    theta_s: float
            adimensional water content at saturation
    theta_r: float
            adimensional residual water content
    qA_star: float 
            dimensional initial water flux at soil surface 
    qB_star: float 
            dimensional  water flux at soil surface for time greater than 0 
    psi0: float 
            prescribed water suction at the water table    
    KMAX: int
            number of solution points
    times: array float
            times at which compute the solution
        
    Returns
    -------
    psi_0: numpy.array float
            initial condition for water suction
    theta_0: numpy.array flaot
            initial condition for water content
    psi: dict
            solution for water suction
            dict of{int: numpy.array}
    theta: dict
            solution for water content
            dict of{int: numpy.array}
    q: 2darray float
            solution outflow at the water table
    """
        
    # compute adimensional quantities Eqs. (4)
    z = z_star*alpha
    L = L_star*alpha
    qA = qA_star/Ks
    qB = qB_star/Ks
    t = alpha*Ks*times/(theta_s-theta_r)
    
    # Eq. (6a)
    K_0 = qA - (qA-np.exp(alpha*psi0))*np.exp(-z)
    
    psi_0 = np.log(K_0)/alpha
    theta_0 = theta_r + (theta_s-theta_r)*K_0
    
    # Eq. (10)
    f = lambda l: np.tan(l*L) + 2*l 
    
    # compute positive roots of the characteristic equation, Eq. (10)
    # Eq. (10) is solved with the bisection algorithm
    intervals = np.arange(-np.pi/2/(L), 200*np.pi/2/(L), np.pi/(L))
    lambda_n = np.zeros(len(intervals)-1)
    for i in range(0,len(intervals)-1):
        lambda_n[i] = bisection(f, intervals[i]+ 1e-9, intervals[i+1]- 1e-9, 100)
    
    # compute adimensional solution Eq. (12a)
    psi = {}
    theta = {}
    for i in range(0,len(t)):
        summation_n = 0.0
        for nn in range(1,len(lambda_n)):
            summation_n = summation_n + np.sin(lambda_n[nn]*z)*np.sin(lambda_n[nn]*L)*np.exp(-lambda_n[nn]**2*t[i]) / (1 + (L/2) + 2*L*lambda_n[nn]**2)

        K = qB - (qB - np.exp(alpha*psi0))*np.exp(-z) - 4*(qB - qA)*np.exp((L - z)/2)*np.exp(-t[i]/4)*summation_n
        
        psi[times[i]] = np.log(K)/alpha
        theta[times[i]] = theta_r + (theta_s-theta_r)*K
 

    
    # compute outflow at the water table (z=0) Eq. (12b)
    tmp_times_star = np.linspace(times[0], times[len(times)-1], 101)
    tmp_t = alpha*Ks*tmp_times_star/(theta_s-theta_r)

    q = np.zeros((2,len(tmp_t)))
    for i in range(0, len(tmp_t)):
        summation_n = 0.0
        for nn in range(1,len(lambda_n)):
            summation_n = summation_n + ( lambda_n[nn]*np.sin(lambda_n[nn]*L)*np.exp(-lambda_n[nn]**2*tmp_t[i]) )/( 1 + L/2 + 2*lambda_n[nn]**2*L )
        
        q[1,i] = Ks*qB - 4*Ks*(qB-qA)*np.exp(L/2)*np.exp(-tmp_t[i]/4)*summation_n
        
    q[0,:] = tmp_times_star
                 
    
    
    return (psi_0, theta_0, psi, theta, q)




def solution_layered_soil(L1_star, L2_star, z_star, Ks1, Ks2, alpha1, alpha2, theta_s1, theta_r1, theta_s2, theta_r2, qA1_star, qB1_star, qA2_star, qB2_star, psi0, times):
    """ Compute the analytical solution for the case of layered soil
    
    [1] Srivastava, R., & Yeh, T. C. J. (1991). Analytical solutions for one‐dimensional, transient   infiltration toward the water table in homogeneous and layered soils. Water Resources Research, 27(5), 753-762.
    

    Parameters
    ----------
    L1_star: float 
            dimensional thickness of the lower layer
    L2_star: float 
            dimensional thickness of the upper layer
    z_star: numpy.array
            dimensional coordinate
    Ks1: float
            saturated hydraulic conductivity of the lower layer
    Ks2: float
            saturated hydraulic conductivity of the upper layer
    alpha1: double
            is a soil pore-size distribution parameter representing the rate of reduction in hydraulic conductivity or moisture for the lower layer
    alpha2: double
            is a soil pore-size distribution parameter representing the rate of reduction in hydraulic conductivity or moisture for the upper layer
    theta_s1: float
            adimensional water content at saturation of the lower layer
    theta_r1: float
            adimensional residual water content of the lower layer
    theta_s2: float
            adimensional water content at saturation of the upper layer
    theta_r2: float
            adimensional residual water content of the upper layer
    qA1_star: float 
            dimensional initial water flux for the lower layer
    qB1_star: float 
            dimensional  water flux for time greater than 0  for the lower layer
    qA2_star: float 
            dimensional initial water flux for the upper layer
    qB2_star: float 
            dimensional  water flux for time greater than 0  for the upper layer
    psi0: float 
            prescribed water suction at the water table    
    times: array float
            times at which compute the solution
        
    Returns
    -------
    psi_0: numpy.array float
            initial condition for water suction
    theta_0: numpy.array flaot
            initial condition for water content
    psi: dict
            solution for water suction
            dict of{int: numpy.array}
    theta: dict
            solution for water content
            dict of{int: numpy.array}
    q: 2darray float
            solution outflow at the water table
    """
    
#     z_star = np.linspace(-L1_star,L2_star,KMAX)
#     z_plot = z_star+L1_star
    
    z1_star = z_star[z_star<0]
    z2_star = z_star[z_star>=0]
    
    # compute adimensional quantities Eqs. (4)
    L1 = L1_star*alpha1 # cm
    L2 = L2_star*alpha2
    qA1 = qA1_star/Ks1
    qA2 = qA2_star/Ks2
    qB1 = qB1_star/Ks1
    qB2 = qB2_star/Ks2
    z1 = z1_star*alpha1
    z2 = z2_star*alpha2
    t = alpha1*Ks1*times/(theta_s1-theta_r1)
    beta = (alpha1*Ks1*(theta_s2-theta_r2))/(alpha2*Ks2*(theta_s1-theta_r1))
    
    # Eq. (15a)
    K1_0 = qA1-(qA1-np.exp(alpha1*0))*np.exp(-L1-z1)
    psi1_0 = 1/alpha1*np.log(K1_0)
    theta1_0 = theta_r1 + (theta_s1-theta_r1)*K1_0
    
    # Eq. (16a)
    K2_0 =  qA2 - (qA2 - (qA1 - (qA1-np.exp(alpha1*0))*np.exp(-L1))**(alpha2/alpha1))*np.exp(-z2)
    psi2_0 = 1/alpha2*np.log(K2_0)
    theta2_0 = theta_r2 + (theta_s2-theta_r2)*K2_0
    
    psi_0 = np.concatenate([psi1_0, psi2_0])
    theta_0 = np.concatenate([theta1_0, theta2_0])
    # Case A
    # Eq. (21b)
    mu = lambda l: (beta*l**2 + (beta-1)/4)**(1/2)
    
    # Eq.(21a)
    f = lambda l: ( np.sin(l*L1) + 2*l*np.cos(l*L1) )*( np.sin(mu(l)*L2) + 2*mu(l)*np.cos(mu(l)*L2) ) - ( 1+4*(mu(l))**2 )*Ks2/Ks1*np.sin(l*L1)*np.sin(mu(l)*L2)
    
    # Compute zeros of Eq.(21a)
    l = np.linspace(0,100,50000)
    Alambda_n = []
    for i in range(0,len(l)-1):
        if(f(l[i])*f(l[i+1])<0):
            Alambda_n.append(bisection(f,l[i],l[i+1],100))
    
    # Case B
    Blambda_n = []
    if(beta>1):
        # Eq. (24b)
        mu = lambda l: ((beta-1)/(4) - beta*l**2)**(1/2)
        # Eq. (24a)
        f = lambda l: (np.sinh(l*L1) + 2*l*np.cosh(l*L1))*(np.sin(mu(l)*L2) + 2*mu(l)*np.cos(mu(l)*L2)) - (1+4*(mu(l)**2))*Ks2/Ks1*np.sinh(l*L1)*np.sin(mu(l)*L2)
        
        # Compute zeros of Eq. (24a)
        l = np.linspace(0,10,5000)
        for i in range(0,len(l)-1):
            if(f(l[i])*f(l[i+1])<0):
                Blambda_n.append(bisection(f,l[i],l[i+1],100))
                
    # Case C
    Cmu_n = []
    if(beta<1):
        print('here')
        # Eq. (27b)
        l = lambda mu: ((1-beta)/(4*beta) - (mu**2)/beta)**(1/2)
        # Eq. (27a)
        f = lambda mu: ( np.sin(l(mu)*L1) + 2*l(mu)*np.cos(l(mu)*L1) ) * ( np.sinh(mu*L2) + 2*mu*np.cosh(mu*L2) ) - (1-4*mu**2)*Ks2/Ks1*np.sin(l(mu)*L1)*np.sinh(mu*L2)
    
        # Compute zeros of Eq. (27a)
        mu = np.linspace(0,np.sqrt(1-beta)/2,10000)
        for i in range(0,len(mu)-1):
            if(f(mu[i])*f(mu[i+1])<0):
                Cmu_n.append(bisection(f,mu[i],mu[i+1],100))
    
    
    # Compute solution Eq. (30a) Eq. (30b)
    psi = {}
    theta = {}
    psi1 = {}
    theta1 = {}
    psi2 = {}
    theta2 = {}
    K1 = np.zeros(len(z1))
    K2 = np.zeros(len(z2))
    for i in range(0,len(t)):
    
        if(beta>1):
            [RA1, RA2] = compute_RA(L1, L2, Ks1, Ks2, alpha1, alpha2, beta, z1, z2, t[i], Alambda_n)
            [RBRC1, RBRC2] = compute_RB(L1, L2, Ks1, Ks2, alpha1, alpha2, beta, z1, z2, t[i], Blambda_n)

        if(beta<1):
            [RA1, RA2] = compute_RA(L1, L2, Ks1, Ks2, alpha1, alpha2, beta, z1, z2, t[i], Alambda_n)
            [RBRC1, RBRC2] = compute_RC(L1, L2, Ks1, Ks2, alpha1, alpha2, beta, z1, z2, t[i], Cmu_n)

#         print(RA1)
#         print(RA2)
#         print(RBRC1)
#         print(RBRC2)
        K2 = qB2 - (qB2 - qB1 + (qB1-np.exp(alpha2*psi0))*np.exp(-L1))*np.exp(-z2) - 4*(qB1-qA1)*np.exp((L2-z2)/2)*(RA2+RBRC2)    
        
        K1 = qB1 - (qB1-np.exp(alpha1*psi0))*np.exp(-(L1+z1)) - 4*(qB1-qA1)*np.exp((L2-z1)/2)*(RA1+RBRC1)
        
        psi1[times[i]] = np.log(K1)/alpha1
        theta1[times[i]] = theta_r1 + (theta_s1-theta_r1)*K1 
        
        psi2[times[i]] = np.log(K2)/alpha2
        theta2[times[i]] = theta_r2 + (theta_s2-theta_r2)*K2 
        
#         tmp = np.concatenate([np.array(np.log(K1)/alpha1),np.array(np.log(K2)/alpha2)])
                             
        psi[times[i]] = np.concatenate([np.log(K1)/alpha1, np.log(K2)/alpha2])
                            
        theta[times[i]] = np.concatenate([theta_r1+(theta_s1-theta_r1)*K1, theta_r2+(theta_s2-theta_r2)*K2])
        
        if i==0:
            KK =K1

    return (psi_0, theta_0, psi, theta)


def compute_RA(L1, L2, Ks1, Ks2, alpha1, alpha2, beta, z1, z2, t, lambda_n):
    """ Compute RA1 and RA2
    
    [1] Srivastava, R., & Yeh, T. C. J. (1991). Analytical solutions for one‐dimensional, transient   infiltration toward the water table in homogeneous and layered soils. Water Resources Research, 27(5), 753-762.
    

    Parameters
    ----------
    L1: float 
            adimensional thickness of the lower layer
    L2: float 
            sdimensional thickness of the upper layer
    Ks1: float
            saturated hydraulic conductivity of the lower layer
    Ks2: float
            saturated hydraulic conductivity of the upper layer
    alpha1: double
            is a soil pore-size distribution parameter representing the rate of reduction in hydraulic conductivity or moisture for the lower layer
    alpha2: double
            is a soil pore-size distribution parameter representing the rate of reduction in hydraulic conductivity or moisture for the upper layer
    beta: float
            
    z1: numpy.array
            adimensional z coordinate of the lower layer
    z2: numpy.array
            adimensional z coordinate of the upper layer
    t: float
            adimensional time
    lambda_n: list
            solution of Eq. (21a)
        
    Returns
    -------
    RA1: float
            sum of residues Eq. (22a)
    RA2: flaot
            sum of residues Eq. (22b)
    """
    
    # Eq. (21b)
    mu = lambda l: (beta*l**2 + (beta-1)/4)**(1/2)
    
    # Eq. (23b)
    A_ss = lambda l: 4*beta*l*Ks2/Ks1 + l*(L1 + beta*L2)
    
    # Eq. (23c)
    A_sc = lambda l: Ks2/Ks1*(beta*l*L2)/(2*mu(l))*(1+4*(mu(l))**2) + 2*l*mu(l)*L1 - (beta*l)/(2*mu(l))*(2+L2)

    # Eq. (23d)
    A_cs = lambda l: Ks2/Ks1*L1/2*(1+4*(mu(l))**2) + 2*beta*l**2*L2 - 1 - L1/2
    
    # Eq. (23e)
    A_cc = lambda l: -mu(l)*(2+L1) - beta*l**2/mu(l)*(2+L2)

    # Eq. (23.a)
    D_n = lambda l: (1/4+l**2)/(mu(l)*l)*( A_ss(l)*np.sin(l*L1)*np.sin(mu(l)*L2) + A_sc(l)*np.sin(l*L1)*np.cos(mu(l)*L2) + A_cs(l)*np.cos(l*L1)*np.sin(mu(l)*L2) + A_cc(l)*np.cos(l*L1)*np.cos(mu(l)*L2) )

    # Eq. (22a)
    RA1 = 0.0
    for nn in range(0,len(lambda_n)):
        D_nn = D_n(lambda_n[nn])
        RA1 = RA1 + np.sin(lambda_n[nn]*(L1+z1))*np.exp(-(1/4+lambda_n[nn]**2)*t)/D_nn
    
    # Eq. (22b)
    RA2 = 0.0
    for nn in range(0,len(lambda_n)):
        mu_nn = mu(lambda_n[nn])
        D_nn = D_n(lambda_n[nn])
        RA2 = RA2 + ( np.sin(lambda_n[nn]*L1)*( np.sin(mu_nn*(L2-z2)) + 2*mu_nn*np.cos(mu_nn*(L2-z2)) )*np.exp(-(1/4+lambda_n[nn]**2)*t) ) * ( D_nn*( np.sin(mu_nn*L2) + 2*mu_nn*np.cos(mu_nn*L2) ) )**(-1)
        
    
    return (RA1, RA2)


def compute_RB(L1, L2, Ks1, Ks2, alpha1, alpha2, beta, z1, z2, t, lambda_n):
    """ Compute RB1 and RB2
    
    [1] Srivastava, R., & Yeh, T. C. J. (1991). Analytical solutions for one‐dimensional, transient   infiltration toward the water table in homogeneous and layered soils. Water Resources Research, 27(5), 753-762.
    

    Parameters
    ----------
    L1: float 
            adimensional thickness of the lower layer
    L2: float 
            sdimensional thickness of the upper layer
    Ks1: float
            saturated hydraulic conductivity of the lower layer
    Ks2: float
            saturated hydraulic conductivity of the upper layer
    alpha1: double
            is a soil pore-size distribution parameter representing the rate of reduction in hydraulic conductivity or moisture for the lower layer
    alpha2: double
            is a soil pore-size distribution parameter representing the rate of reduction in hydraulic conductivity or moisture for the upper layer
    beta: float
            
    z1: numpy.array
            adimensional z coordinate of the lower layer
    z2: numpy.array
            adimensional z coordinate of the upper layer
    t: float
            adimensional time
    lambda_n: list
            solution of Eq. (24a)

    Returns
    -------
    RB1: float
            sum of residues Eq. (25a)
    RB2: flaot
            sum of residues Eq. (25b)
    """

    # Eq. (24b)
    mu = lambda l: ((beta-1)/(4) - beta*l**2)**(1/2)
    
    # Eq. (26b)
    B_ss = lambda l: 4*beta*l*Ks2/Ks1 + l*(L1+beta*L2)
    
    # Eq. (26c)
    B_sc = lambda l: Ks2/Ks1*beta*l*L2/(2*mu(l))*(1+4*(mu(l))**2) + 2*l*mu(l)*L1 - beta*l/(2*mu(l))*(2+L2)

    # Eq. (26d)
    B_cs = lambda l: -Ks2/Ks1*L1/2*(1+4*(mu(l))**2) + 2*beta*l**2*L2 + 1 + L1/2
    
    # Eq. (26e)
    B_cc = lambda l: mu(l)*(2+L1) - beta*l**2/mu(l)*(2+L2)

    # Eq. (26a)
    D_n = lambda l: (1/4-l**2)/(mu(l)*l) * ( B_ss(l)*np.sinh(l*L1) * np.sin(mu(l)*L2) + B_sc(l)*np.sinh(l*L1)*np.cos(mu(l)*L2) + B_cs(l)*np.cosh(l*L1)*np.sin(mu(l)*L2) + B_cc(l)*np.cosh(l*L1)*np.cos(mu(l)*L2) )


    # Eq. (25a)
    RB1 = 0.0
    for nn in range(0,len(lambda_n)):
#         mu_nn = mu(lambda_n[nn])
        D_nn = D_n(lambda_n[nn])
        RB1 = RB1 + np.sinh(lambda_n[nn]*(L1+z1))*np.exp(-(1/4-lambda_n[nn]**2)*t)/D_nn
       
    # Eq. (25b)
    RB2 = 0.0
    for nn in range(0,len(lambda_n)):
        mu_nn = mu(lambda_n[nn])
        D_nn = D_n(lambda_n[nn])
        RB2 = RB2 + ( np.sinh(lambda_n[nn]*L1)*( np.sin(mu(lambda_n[nn])*(L2-z2)) + 2*mu(lambda_n[nn])*np.cos(mu(lambda_n[nn])*(L2-z2)) )*np.exp(-(1/4-lambda_n[nn]**2)*t) )  * ( D_nn*( np.sin(mu(lambda_n[nn])*L2) + 2*mu(lambda_n[nn])*np.cos(mu(lambda_n[nn])*L2) ) )**(-1)
                 
    return (RB1, RB2)


def compute_RC(L1, L2, Ks1, Ks2, alpha1, alpha2, beta, z1, z2, t, mu_n):
    """ Compute RC1 and RC2
    
    [1] Srivastava, R., & Yeh, T. C. J. (1991). Analytical solutions for one‐dimensional, transient   infiltration toward the water table in homogeneous and layered soils. Water Resources Research, 27(5), 753-762.
    

    Parameters
    ----------
    L1: float 
            adimensional thickness of the lower layer
    L2: float 
            sdimensional thickness of the upper layer
    Ks1: float
            saturated hydraulic conductivity of the lower layer
    Ks2: float
            saturated hydraulic conductivity of the upper layer
    alpha1: double
            is a soil pore-size distribution parameter representing the rate of reduction in hydraulic conductivity or moisture for the lower layer
    alpha2: double
            is a soil pore-size distribution parameter representing the rate of reduction in hydraulic conductivity or moisture for the upper layer
    beta: float
            
    z1: numpy.array
            adimensional z coordinate of the lower layer
    z2: numpy.array
            adimensional z coordinate of the upper layer
    t: float
            adimensional time
    mu_n: list
            solution of Eq. (27a)
        
    Returns
    -------
    RC1: float
            sum of residues Eq. (28a)
    RC2: flaot
            sum of residues Eq. (28b)
    """
    
    # Eq. (27b)
    l = lambda mu: ((1-beta)/(4*beta) - (mu**2)/beta)**(1/2)
    
    # Eq. (29b)
    C_ss = lambda mu: 4*beta*l(mu)*Ks2/Ks1 + l(mu)*(L1+beta*L2)

    # Eq. (29c)
    C_sc = lambda mu: -Ks2/Ks1*(beta*l(mu)*L2)/(2*mu)*(1-4*mu**2) + 2*l(mu)*mu*L1 + (beta*l(mu))/(2*mu)*(2+L2)

    # Eq. (29d)
    C_cs = lambda mu: Ks2/Ks1*L1/2*(1-4*mu**2) + 2*beta*(l(mu))**2*L2 - 1 - L2/2
    
    # Eq. (29e)
    C_cc = lambda mu: -mu*(2+L1) + beta*(l(mu))**2/mu*(2+L2)

    # Eq. (29a)
    D_n = lambda mu: (1/4+(l(mu))**2)/(mu*l(mu)) * ( C_ss(mu)*np.sin(l(mu)*L1)*np.sinh(mu*L2) + C_sc(mu)*np.sin(l(mu)*L1)*np.cosh(mu*L2) + C_cs(mu)*np.cos(l(mu)*L1)*np.sinh(mu*L2) + C_cc(mu)*np.cos(l(mu)*L1)*np.cosh(mu*L2) )

    # Eq. (28a)
    RC1 = 0.0
    for nn in range(0,len(mu_n)):
        l_nn = l(mu_n[nn])
        D_nn = D_n(mu_n[nn])
        RC1 = RC1 + np.sin(l_nn*(L1+z1))*np.exp(-(1/4+l_nn**2)*t)/D_nn
    
    # Eq. (28b)
    RC2 = 0.0
    for nn in range(0,len(mu_n)):
        l_nn = l(mu_n[nn])
        D_nn = D_n(mu_n[nn])
        RC2 = RC2 + ( np.sin(l_nn*L1)*( np.sinh(mu_n[nn]*(L2-z2)) + 2*mu_n[nn]*np.cosh(mu_n[nn]*(L2-z2)) )*np.exp(-(1/4+l_nn**2)*t) )  * ( D_nn*( np.sinh(mu_n[nn]*L2) + 2*mu_n[nn]*np.cosh(mu_n[nn]*L2) ) )**(-1)
      
    
    return (RC1, RC2)



def bisection(f,a,b,N):
    '''Approximate solution of f(x)=0 on interval [a,b] by bisection method.

    Parameters
    ----------
    f : function
        The function for which we are trying to approximate a solution f(x)=0.
    a,b : numbers
        The interval in which to search for a solution. The function returns
        None if f(a)*f(b) >= 0 since a solution is not guaranteed.
    N : (positive) integer
        The number of iterations to implement.

    Returns
    -------
    x_N : number
        The midpoint of the Nth interval computed by the bisection method. The
        initial interval [a_0,b_0] is given by [a,b]. If f(m_n) == 0 for some
        midpoint m_n = (a_n + b_n)/2, then the function returns this solution.
        If all signs of values f(a_n), f(b_n) and f(m_n) are the same at any
        iteration, the bisection method fails and return None.

    Examples
    --------
    >>> f = lambda x: x**2 - x - 1
    >>> bisection(f,1,2,25)
    1.618033990263939
    >>> f = lambda x: (2*x - 1)*(x - 3)
    >>> bisection(f,0,1,10)
    0.5
    
    https://www.math.ubc.ca/~pwalls/math-python/roots-optimization/bisection/
    '''
    if f(a)*f(b) >= 0:
        print("Bisection method fails.")
        return None
    a_n = a
    b_n = b
    for n in range(1,N+1):
        m_n = (a_n + b_n)/2
        f_m_n = f(m_n)
        if f(a_n)*f_m_n < 0:
            a_n = a_n
            b_n = m_n
        elif f(b_n)*f_m_n < 0:
            a_n = m_n
            b_n = b_n
        elif f_m_n == 0:
#             print("Found exact solution.")
            return m_n
        else:
            print("Bisection method fails.")
            return None
    return (a_n + b_n)/2