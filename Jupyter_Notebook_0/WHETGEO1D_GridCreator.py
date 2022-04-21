# -*- coding: utf-8 -*-
"""
Created on 10/20/2020

This is used to create the grid for WHETGEO 1D model.

@author: Niccolo` Tubini, Concetta D'Amato, Riccardo Rigon
@license: creative commons 4.0
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy.interpolate import interp1d
from netCDF4 import Dataset

def grid1D(data_grid, dz_min, b, dz_max, grid_type, **kwargs):
    '''
    This function creates the geometry of 1D grid for a finite volume numerical.     
    
    :param data_grid: pandas dataframe containg the grid_input_file.csv
    :type data_grid: pandas dataframe

     :param dz_min: thickness of the first layer
    :type dz_min: float

    :param b: growth rate, range [0, 1[.
    :type b: float

    :param grid_type: type of the grid. The grid can be created using picewise mesh spacing, classical, or an exponential function, exponential
    
    :param \**kwargs:
    See below
    
    :Keyword Arguments
        * *param shallow_water* (boolean) add a node for surface water

    return:
    
    KMAX: number of control volumes
    type NMAX: int

    eta: vertical coordinate of control volumes centroids positive upward with origin set at soil surface.
    type eta: array

    eta_dual: vertical coordinate of control volumes interfaces positive upward with origin set at soil surface.
    type eta_dual: array

    space_delta: is the distance between two adjacent control volumes. 
                This quantity is used to compute gradients.    
    type space_delta: array
    
    z: vertical coordinate of control volumes centroids positive upward with origin set at soil column bottom. This is the spatial coordinate used to to write the equation
    type z: array

    z_dual: vertical coordinate of control volumes interfaces positive upward with origin set at soil column bottom. 
    type zDual: array
        
    '''
    shallow_water = kwargs.get('shallow_water', False)
    
    if(grid_type=='classical'):
        [KMAX, eta, eta_dual, space_delta, z, z_dual, control_volume] = build_grid(data_grid)
    elif(grid_type=='exponential'):
        [KMAX, eta, eta_dual, space_delta, z, z_dual, control_volume] = build_grid_exponential(data_grid, dz_min, b)
    elif(grid_type=='mixed'):
        [KMAX, eta, eta_dual, space_delta, z, z_dual, control_volume] = build_grid_mixed(data_grid, dz_min, b, dz_max)    
    else:
        ('Check grid_type')
        
    if(shallow_water==True):
        KMAX = KMAX+1
        eta = np.append(eta, [0])
        z = np.append(z, - data_grid['eta'][np.size(data_grid['eta'])-1])
        control_volume = np.append(control_volume, [1])
        space_delta = np.append(space_delta, [-9999.0])
        
    return [KMAX, eta, eta_dual, space_delta, z, z_dual, control_volume]


def build_grid(data_grid):
    '''
    This function creates the geometry of 1D grid for a finite volume numerical. The discretizion is 
	with constant mesh spacing in each layers.
    scheme.
    
    
    :param data_grid: pandas dataframe containg the grid_input_file.csv
	:type data_grid: pandas dataframe
	

    return:
    
	KMAX: number of control volumes
	type KMAX: int
	
	VECTOR_LENGTH: length of the array. This value is equal to NMAX is size_factor 1, suggested whenever there is no need to regrid.
	type VECTOR_LENGTH: int
	
    eta: vertical coordinate of control volumes centroids positive upward with origin set at soil surface.
	type eta: array
	
    eta_dual: vertical coordinate of control volumes interfaces positive upward with origin set at soil surface.
	type eta_dual: array

	space_delta: is the distance between two adjacent control volumes. 
                This quantity is used to compute gradients.    
	type space_delta: array
    
    z: vertical coordinate of control volumes centroids positive upward with origin set at soil column bottom. This is the spatial coordinate used to to write the equation
	type z: array
    
    z_dual: vertical coordinate of control volumes interfaces positive upward with origin set at soil column bottom. 
	type zDual: array		
	
    control_volume: dimension of control volume 
	type control_volume: array		
        
    '''
    # get the number of control volumes
    KMAX = int(data_grid['K'].sum())-int(len(data_grid[data_grid.Type == 'M']))
    
    # list containing centroids coordinates
    tmp_eta = []
    # list containing control volumes interface coordinates
    tmp_eta_dual = []
    
    # array containing centroids coordinates measured along eta
    eta = np.zeros(KMAX,dtype=float)
    # array containing control volumes interface coordinates measured along eta
    eta_dual = np.zeros(KMAX+1,dtype=float)
    # array containing centroids coordinates measured along z
    z = np.zeros(KMAX,dtype=float)
    # array containing control volumes interface coordinates measured along z
    z_dual = np.zeros(KMAX+1,dtype=float)
    # array containing distances between centroids (used to compute gradient)
    space_delta = np.zeros(KMAX+1,dtype=float)
    # array containing control volume size
    control_volume = np.zeros(KMAX,dtype=float)



    
    for i in range(np.size(data_grid.index)-1,0,-1):
		
        if data_grid['Type'][i]=='L' and data_grid['Type'][i-1]=='L':
			
            deta = ( data_grid['eta'][i]-data_grid['eta'][i-1])/data_grid['K'][i-1]
            tmp_eta=np.append(tmp_eta, np.linspace(data_grid['eta'][i]-deta/2,data_grid['eta'][i-1]+deta/2,num=data_grid['K'][i-1],endpoint=True) )
            tmp_eta_dual=np.append(tmp_eta_dual, np.linspace(data_grid['eta'][i],data_grid['eta'][i-1],num=data_grid['K'][i-1]+1,endpoint=True) )
			
        elif data_grid['Type'][i]=='L' and data_grid['Type'][i-1]=='M':
			
            deta = ( data_grid['eta'][i]-data_grid['eta'][i-1])/data_grid['K'][i-1]
            tmp_eta=np.append(tmp_eta, np.linspace(data_grid['eta'][i]-deta/2,data_grid['eta'][i-1],num=data_grid['K'][i-1],endpoint=True) )
            tmp_eta_dual=np.append(tmp_eta_dual, np.linspace(data_grid['eta'][i],data_grid['eta'][i-1]+deta/2,num=data_grid['K'][i-1],endpoint=True) )
			
        elif data_grid['Type'][i]=='M' and data_grid['Type'][i-1]=='L':
			
            deta = ( data_grid['eta'][i]-data_grid['eta'][i-1])/data_grid['K'][i-1]
            tmp_eta=np.append(tmp_eta, np.linspace(data_grid['eta'][i],data_grid['eta'][i-1]+deta/2,num=data_grid['K'][i-1],endpoint=True) )
            tmp_eta_dual=np.append(tmp_eta_dual, np.linspace(data_grid['eta'][i]-deta/2,data_grid['eta'][i-1],num=data_grid['K'][i-1],endpoint=True) )
			
        else:
            print("ERROR!!")  
        
    # to eliminate doubles
    tmp_eta=[ii for n,ii in enumerate(tmp_eta) if ii not in tmp_eta[:n]]
    tmp_eta_dual=[ii for n,ii in enumerate(tmp_eta_dual) if ii not in tmp_eta_dual[:n]]
    
   
    # move from list to array
    for i in range(0,len(tmp_eta)):
        eta[i] = tmp_eta[i]
        z[i] = tmp_eta[i] - data_grid['eta'][np.size(data_grid['eta'])-1]
        

    for i in range(0,len(tmp_eta_dual)):

        eta_dual[i] = tmp_eta_dual[i]
        z_dual[i] = tmp_eta_dual[i] - data_grid['eta'][np.size(data_grid['eta'])-1]

        if i==0:

            space_delta[i] = np.abs(eta_dual[i]-eta[i])

        elif i==np.size(eta_dual)-1:

            space_delta[i] = np.abs(eta_dual[i]-eta[i-1])

        else:

            space_delta[i] = np.abs(eta[i-1]-eta[i]) 
           
    for i in range(0,len(eta_dual)-1):
        control_volume[i] = np.abs(eta_dual[i]-eta_dual[i+1])		
   
    return [KMAX, eta, eta_dual, space_delta, z, z_dual, control_volume]


def build_grid_exponential(data_grid, dz_min, b):
    '''
    This function creates the geometry of 1D grid for a finite volume numerical. The discretizion is 
    with exponential function (Gubler S. et al. 2013, doi:10.5194/gmd-6-1319-2013).
    
    
    :param data_grid: pandas dataframe containg the grid_input_file.csv
    :type data_grid: pandas dataframe

     :param dz_min: thickness of the first layer
    :type dz_min: float

    :param b: growth rate, range [0, 1[.
    :type b: float

   
    
    return:
    
    KMAX: number of control volumes
    type NMAX: int

    eta: vertical coordinate of control volumes centroids positive upward with origin set at soil surface.
    type eta: array

    eta_dual: vertical coordinate of control volumes interfaces positive upward with origin set at soil surface.
    type eta_dual: array

    space_delta: is the distance between two adjacent control volumes. 
                This quantity is used to compute gradients.    
    type space_delta: array
    
    z: vertical coordinate of control volumes centroids positive upward with origin set at soil column bottom. This is the spatial coordinate used to to write the equation
    type z: array

    z_dual: vertical coordinate of control volumes interfaces positive upward with origin set at soil column bottom. 
    type zDual: array		
	
    control_volume: dimension of control volume 
	type control_volume: array	
        
    '''
  
    # list containing layer thickness
    tmp_dz = []

    z_max = -data_grid['eta'][np.size(data_grid['eta'])-1]
    dz_sum = 0
    k = 0
    while (z_max-dz_sum)>1E-12:

        tmp_dz.append(dz_min*(1+b)**k)
        dz_sum = dz_sum + dz_min*(1+b)**k

        if(dz_sum>z_max):

            tmp_dz.pop()
            dz_sum = dz_sum-dz_min*(1+b)**k
            tmp_dz.append(z_max-dz_sum)
            dz_sum = dz_sum + z_max-dz_sum

        k = k+1
        
    # get the number og control volumes
    KMAX = len(tmp_dz)
    dz = np.zeros(KMAX,dtype=float)
    
    for i in range(0,KMAX):

        dz[i] = tmp_dz[i]
    
    # array containing centroids coordinates measured along eta
    eta = np.zeros(KMAX,dtype=float)
    # array containing control volumes interface coordinates measured along eta
    eta_dual = np.zeros(KMAX+1,dtype=float)
    # array containing centroids coordinates measured along z
    z = np.zeros(KMAX,dtype=float)
    # array containing control volumes interface coordinates measured along z
    z_dual = np.zeros(KMAX+1,dtype=float)
    # array containing distances between centroids (used to compute gradient)
    space_delta = np.zeros(KMAX+1,dtype=float)
    # array containing control volume size
    control_volume = np.zeros(KMAX,dtype=float)
	
	
    tmp = 0
    for i in range(0,KMAX):
        z[i] = dz[KMAX-1-i]/2+tmp
        z_dual[i] = tmp
        tmp = tmp+dz[KMAX-1-i]
        eta[i] = -z_max + z[i]
        eta_dual[i] = -z_max + z_dual[i]

    z_dual[KMAX] = z_max
    eta_dual[KMAX] = 0.0       


    for i in range(0,KMAX+1):

        if i==0:

            space_delta[i] = np.abs(eta_dual[i]-eta[i])

        elif i==np.size(eta_dual)-1:

            space_delta[i] = np.abs(eta_dual[i]-eta[i-1])

        else:

            space_delta[i] = np.abs(eta[i-1]-eta[i]) 

    for i in range(0,len(eta_dual)-1):
        control_volume[i] = np.abs(eta_dual[i]-eta_dual[i+1])	
	

    return [KMAX, eta, eta_dual, space_delta, z, z_dual, control_volume]


def build_grid_mixed(data_grid, dz_min, b, dz_max):
    '''
    This function creates the geometry of 1D grid for a finite volume numerical. The discretizion is based
    on the exponential function (Gubler S. et al. 2013, doi:10.5194/gmd-6-1319-2013).
    Where the dz is larger that dz_max, the grid is switched to a constant spaced grid.
    
    
    :param data_grid: pandas dataframe containg the grid_input_file.csv
    :type data_grid: pandas dataframe

    :param dz_min: thickness of the first layer
    :type dz_min: float

    :param b: growth rate, range [0, 1[.
    :type b: float
    
    :param dz_max: maximum thickness of a layer
    :type dz_max: float


   
    return:
    
    KMAX: number of control volumes
    type NMAX: int

    eta: vertical coordinate of control volumes centroids positive upward with origin set at soil surface.
    type eta: array

    eta_dual: vertical coordinate of control volumes interfaces positive upward with origin set at soil surface.
    type eta_dual: array

    space_delta: is the distance between two adjacent control volumes. 
                This quantity is used to compute gradients.    
    type space_delta: array
    
    z: vertical coordinate of control volumes centroids positive upward with origin set at soil column bottom. This is the spatial coordinate used to to write the equation
    type z: array

    z_dual: vertical coordinate of control volumes interfaces positive upward with origin set at soil column bottom. 
    type zDual: array		
	
    control_volume: dimension of control volume 
	type control_volume: array	
        
    '''
    
    tmp_dz = []
    
    for layer in range(0,np.size(data_grid['eta'])-1):
        
        z_0 = -data_grid['eta'][layer]
        z_1 = -data_grid['eta'][layer+1]
        z_max = z_1-z_0

        dz_sum = 0
        k = 0
        tmp_layer_dz = []
        
        while (z_max-dz_sum)>1E-12:
            tmp_layer_dz.append(dz_min*(1+b)**k)
            dz_sum = dz_sum + dz_min*(1+b)**k

            if(dz_sum>z_max):
        
                tmp_layer_dz.pop()
                dz_sum = dz_sum-dz_min*(1+b)**k
                tmp_layer_dz.append(z_max-dz_sum)
                dz_sum = dz_sum + z_max-dz_sum

            k = k+1
        tmp_layer_dz_flipped = list(reversed(tmp_layer_dz))
                
        combined_tmp_layer_dz= np.minimum(tmp_layer_dz, tmp_layer_dz_flipped)
        
        combined_layer_dz = [x if (x<dz_max) else np.nan for x in combined_tmp_layer_dz]

        if np.isnan(combined_layer_dz).any()==False:
            compare = [1 if i/j==1 else 0 for i, j in zip(combined_layer_dz,tmp_layer_dz)]
            index = len(compare)-compare[::-1].index(1)-1
        else:
            index = combined_layer_dz[0:].index(np.nan)
        constant_layer_dz = np.ones(int (np.ceil( (z_max-np.nansum(combined_layer_dz))/dz_max))) * (z_max-np.nansum(combined_layer_dz))/np.ceil( (z_max-np.nansum(combined_layer_dz))/dz_max)

        tmp = []
        if index == 0:
            tmp.extend(constant_layer_dz)
        else:
            tmp.extend(combined_layer_dz[0:index])
            tmp.extend(list(constant_layer_dz))
            tmp.extend(combined_layer_dz[index+combined_layer_dz.count(np.nan):len(combined_layer_dz)] ) 

        tmp_dz.extend(tmp)

    # get the number og control volumes
    KMAX = len(tmp_dz)
    dz = np.zeros(KMAX,dtype=float)
    
    for i in range(0,KMAX):
        dz[i] = tmp_dz[i]
    
    # array containing centroids coordinates measured along eta
    eta = np.zeros(KMAX,dtype=float)
    # array containing control volumes interface coordinates measured along eta
    eta_dual = np.zeros(KMAX+1,dtype=float)
    # array containing centroids coordinates measured along z
    z = np.zeros(KMAX,dtype=float)
    # array containing control volumes interface coordinates measured along z
    z_dual = np.zeros(KMAX+1,dtype=float)
    # array containing distances between centroids (used to compute gradient)
    space_delta = np.zeros(KMAX+1,dtype=float)
    # array containing control volume size
    control_volume = np.zeros(KMAX,dtype=float)
	
    tmp = 0
    for i in range(0,KMAX):
        z[i] = dz[KMAX-1-i]/2+tmp
        z_dual[i] = tmp
        tmp = tmp+dz[KMAX-1-i]
        eta[i] = data_grid['eta'][np.size(data_grid['eta'])-1] + z[i]
        eta_dual[i] = data_grid['eta'][np.size(data_grid['eta'])-1] + z_dual[i]


    z_dual[KMAX] = -data_grid['eta'][np.size(data_grid['eta'])-1]#z_max
    eta_dual[KMAX] = 0.0       


    for i in range(0,KMAX+1):

        if i==0:

            space_delta[i] = np.abs(eta_dual[i]-eta[i])

        elif i==np.size(eta_dual)-1:

            space_delta[i] = np.abs(eta_dual[i]-eta[i-1])

        else:

            space_delta[i] = np.abs(eta[i-1]-eta[i]) 

    for i in range(0,len(eta_dual)-1):
        control_volume[i] = np.abs(eta_dual[i]-eta_dual[i+1])	
	
    return [KMAX, eta, eta_dual, space_delta, z, z_dual, control_volume]


def set_initial_condition(data, eta, psi_interp_model, T_interp_model, **kwargs):
    '''
    This function define the problem initial condition for water suction and temperature. The initial condition
    is interpolated starting from some pairs (eta,Psi0,T0) contained in a .csv file.
    The interpolation is performed using the class scipy.interpolate.interp1d
    
    
    :param data: pandas dataframe containg tuple of (eta, Psi0, T0)
    :type data_grid: pandas dataframe
    
    :param eta: vertical coordinate of control volume centroids. It is positive upward with 
        origin set at soil surface.
    :type eta: list
    
    :param interp_model: specifies the kind of interpolation as a string. 
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html#scipy.interpolate.interp1d
    :type ic_type: str

    :param psi_interp_model: specifies the kind of interpolation for water suction. 
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html#scipy.interpolate.interp1d
    :type psi_interp_model: str
    
    :param T_interp_model: specifies the kind of interpolation for temperature. 
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html#scipy.interpolate.interp1d
    :type psi_interp_model: str
    
    :param \**kwargs:
    See below
    
    :Keyword Arguments
        * *bounds_error* look at the documentation of scipy.interpolate.interp1d
        * *fill_value* look at the documentation of scipy.interpolate.interp1d
        * *shallow_water* (boolean) value to consider the shallow water
        * *param water_ponding_0* (numpy.float64) water ponding depth at the beginnig of the simulation   
        * *param T_water_ponding_0* (numpy.float64) temperature of water ponding depth at the beginnig of the simulation  
    
    return:
    
    psi_0: initial condition for water suction
    type psi_0: array

    T_0: initial condition for temperature
    type T_0: array
    '''
    
    bound_error = kwargs.get('bounds_error',False)
    fill_value =  kwargs.get('fill_value',np.nan)
    
    shallow_water = kwargs.get('shallow_water',False)
    water_ponding_0 = kwargs.get('water_ponding_0',np.nan)
    T_water_ponding_0 = kwargs.get('T_water_ponding_0',np.nan)
    
    eta_points = data['eta']
    psi_0_points = data['Psi0']
    T_0_points = data['T0']
    
    f = interp1d(eta_points, psi_0_points, kind=psi_interp_model, assume_sorted=False, bounds_error=bound_error, fill_value=fill_value)
    psi_0 = f(eta)
    
    f = interp1d(eta_points, T_0_points, kind=T_interp_model, assume_sorted=False, bounds_error=bound_error, fill_value=fill_value)
    T_0 = f(eta)
    
    if(shallow_water==True):
        psi_0[len(psi_0)-1] = water_ponding_0        
        T_0[len(T_0)-1] = T_water_ponding_0

    return [psi_0, T_0]




def set_parameters_richards(data_grid, data_parameter, data_dictionary, KMAX, eta):
    '''
    This function associate to each control volume a label that identifies 
    the rheology model, the set of parameters describing the soil type, and the max/min cell size 
    for regridding.
    
    :param data_grid: pandas dataframe containg the grid_input_file.csv
    :type data: pandas dataframe
    
    :param data_parameter: pandas dataframe containg the parameter_input_file.csv
    :type data: pandas dataframe
    
    :param data_dictionary: pandas dataframe containg a dictionary for the SWRC parameters
    :type data: pandas dataframe
    
    :param KMAX: number of control volumes.
    :type KMAX: int
    
    :param eta: vertical coordinate of control volume centroids. It is positive upward with 
        origin set at soil surface.
    :type eta: list
    
    return:
    
    equation_state_ID: vector containing the ID identifing the state equation of the k-th element
    type: array
    
    parameters_ID: vector containing the ID identifing the set of parameters of the k-th element
    type: array
        	
    theta_s: vector containing the adimensional water content at saturation
    type:array
    
    theta_r: vector containing the residual adimensional water content
    type:array
    
    par_1: vector containing the parameter 1 of the SWRC model
    type:array
    
    par_2: vector containing the parameter 2 of the SWRC model
    type:array
    
    par_3: vector containing the parameter 3 of the SWRC model
    type:array
    
    par_4: vector containing the parameter 4 of the SWRC model
    type:array
    
    par_5: vector containing the parameter 5 of the SWRC model
    type:array
    
    alpha_ss: vector containing the acquitard compressibility 
    type:array
    
    beta_ss: vector containing the water compressibility
    type:array
    
    ks: vector containing the saturated hydraulic conductivity
    type:array
    '''
    equation_state_ID = np.zeros(KMAX, dtype=float)
    parameter_ID = np.zeros(KMAX, dtype=float)
    
    theta_s = np.zeros(data_parameter.shape[0]+1, dtype=float)
    theta_r = np.zeros(data_parameter.shape[0]+1, dtype=float)
    par_1 = np.zeros(data_parameter.shape[0]+1, dtype=float)
    par_2 = np.zeros(data_parameter.shape[0]+1, dtype=float)
    par_3 = np.zeros(data_parameter.shape[0]+1, dtype=float)
    par_4 = np.zeros(data_parameter.shape[0]+1, dtype=float)
    par_5 = np.zeros(data_parameter.shape[0]+1, dtype=float)
    alpha_ss = np.zeros(data_parameter.shape[0]+1, dtype=float)
    beta_ss = np.zeros(data_parameter.shape[0]+1, dtype=float)
    ks = np.zeros(data_parameter.shape[0]+1, dtype=float)
    
    coord_layer = []
    tmp_equation_state_ID = []
    tmp_parameter_ID = []
    
    for i in data_grid.index:

        if data_grid['Type'][i] == 'L':
            
            coord_layer.append(data_grid['eta'][i])
            tmp_equation_state_ID.append(data_grid['equationStateID'][i])
            tmp_parameter_ID.append(data_grid['parameterID'][i])
           
     
    for i in range(np.size(coord_layer)-1,0,-1):
        
        for j in range(0,KMAX):
            
            if(eta[j]>coord_layer[i] and eta[j]<coord_layer[i-1] ):
                
                equation_state_ID[j] = tmp_equation_state_ID[i-1]
                parameter_ID[j] = tmp_parameter_ID[i-1]
        
    parameter_ID[KMAX-1] = parameter_ID[KMAX-2]   
    
    data_parameter.rename(columns=dict(data_dictionary.values), inplace=True)    

    columns = ['thetaS', 'thetaR', 'par1', 'par2', 'par3', 'par4', 'par5', 'alphaSpecificStorage', 'betaSpecificStorage', 'Ks']
    tmp_df = pd.DataFrame(columns=columns)
    tmp_df = pd.concat([tmp_df, data_parameter],sort=False)
    tmp_df.replace(np.nan, -9999.0, inplace=True)
    
    theta_s[1:] = tmp_df['thetaS'].to_numpy()
    theta_r[1:] = tmp_df['thetaR'].to_numpy()
    par_1[1:] = tmp_df['par1'].to_numpy()
    par_2[1:] = tmp_df['par2'].to_numpy()
    par_3[1:] = tmp_df['par3'].to_numpy()
    par_4[1:] = tmp_df['par4'].to_numpy()
    par_5[1:] = tmp_df['par5'].to_numpy()
    alpha_ss[1:] = tmp_df['alphaSpecificStorage'].to_numpy()
    beta_ss[1:] = tmp_df['betaSpecificStorage'].to_numpy()
    ks[1:] = tmp_df['Ks'].to_numpy()
    
    return [equation_state_ID, parameter_ID, theta_s, theta_r, par_1, par_2, par_3, par_4, par_5,
           alpha_ss, beta_ss, ks]


def set_parameters_richards_lysimeter(data_grid, data_parameter, data_dictionary, KMAX, eta):
    '''
    This function associate to each control volume a label that identifies 
    the rheology model, the set of parameters describing the soil type, and the max/min cell size 
    for regridding.
    
    :param data_grid: pandas dataframe containg the grid_input_file.csv
    :type data: pandas dataframe
    
    :param data_parameter: pandas dataframe containg the parameter_input_file.csv
    :type data: pandas dataframe
    
    :param data_dictionary: pandas dataframe containg a dictionary for the SWRC parameters
    :type data: pandas dataframe
    
    :param KMAX: number of control volumes.
    :type KMAX: int
    
    :param eta: vertical coordinate of control volume centroids. It is positive upward with 
        origin set at soil surface.
    :type eta: list
    
    return:
    
    equation_state_ID: vector containing the ID identifing the state equation of to model the k-th element
    type: array
    
    parameters_ID: vector containing the ID identifing the set of parameters of to model the k-th element
    type: array
        	
    theta_s: vector containing the adimensional water content at saturation
    type:array
    
    theta_r: vector containing the residual adimensional water content
    type:array
    
    theta_wp: vector containing the adimensional water content at the wilting point
    type:array
    
    theta_fc: vector containing the adimensional water content of the field capacity
    type:array
    
    par_1: vector containing the parameter 1 of the SWRC model
    type:array
    
    par_2: vector containing the parameter 2 of the SWRC model
    type:array
    
    par_3: vector containing the parameter 3 of the SWRC model
    type:array
    
    par_4: vector containing the parameter 4 of the SWRC model
    type:array
    
    par_5: vector containing the parameter 5 of the SWRC model
    type:array
    
    alpha_ss: vector containing the acquitard compressibility 
    type:array
    
    beta_ss: vector containing the water compressibility
    type:array
    
    ks: vector containing the saturated hydraulic conductivity
    type:array
    '''
    equation_state_ID = np.zeros(KMAX, dtype=float)
    parameter_ID = np.zeros(KMAX, dtype=float)
    
    theta_s = np.zeros(data_parameter.shape[0]+1, dtype=float)
    theta_r = np.zeros(data_parameter.shape[0]+1, dtype=float)
    theta_wp = np.zeros(data_parameter.shape[0]+1, dtype=float)
    theta_fc = np.zeros(data_parameter.shape[0]+1, dtype=float)
    par_1 = np.zeros(data_parameter.shape[0]+1, dtype=float)
    par_2 = np.zeros(data_parameter.shape[0]+1, dtype=float)
    par_3 = np.zeros(data_parameter.shape[0]+1, dtype=float)
    par_4 = np.zeros(data_parameter.shape[0]+1, dtype=float)
    par_5 = np.zeros(data_parameter.shape[0]+1, dtype=float)
    alpha_ss = np.zeros(data_parameter.shape[0]+1, dtype=float)
    beta_ss = np.zeros(data_parameter.shape[0]+1, dtype=float)
    ks = np.zeros(data_parameter.shape[0]+1, dtype=float)
    
    coord_layer = []
    tmp_equation_state_ID = []
    tmp_parameter_ID = []
    tmp_regrid_ID = []

    for i in data_grid.index:

        if data_grid['Type'][i] == 'L':
            
            coord_layer.append(data_grid['eta'][i])
            tmp_equation_state_ID.append(data_grid['equationStateID'][i])
            tmp_parameter_ID.append(data_grid['parameterID'][i])
           
     
    for i in range(np.size(coord_layer)-1,0,-1):
        
        for j in range(0,KMAX):
            
            if(eta[j]>coord_layer[i] and eta[j]<coord_layer[i-1] ):
                
                equation_state_ID[j] = tmp_equation_state_ID[i-1]
                parameter_ID[j] = tmp_parameter_ID[i-1]
     
    parameter_ID[KMAX-1] = parameter_ID[KMAX-2]
        
    data_parameter.rename(columns=dict(data_dictionary.values), inplace=True)    
    
    columns = ['thetaS', 'thetaR', 'thetaWP', 'thetaFC', 'par1', 'par2', 'par3', 'par4', 'par5', 'alphaSpecificStorage', 'betaSpecificStorage', 'Ks']
    tmp_df = pd.DataFrame(columns=columns)
    tmp_df = pd.concat([tmp_df, data_parameter],sort=False)
    tmp_df.replace(np.nan, -9999.0, inplace=True)
    
    theta_s[1:] = tmp_df['thetaS'].to_numpy()
    theta_r[1:] = tmp_df['thetaR'].to_numpy()
    theta_wp[1:] = tmp_df['thetaWP'].to_numpy()
    theta_fc[1:] = tmp_df['thetaFC'].to_numpy()
    par_1[1:] = tmp_df['par1'].to_numpy()
    par_2[1:] = tmp_df['par2'].to_numpy()
    par_3[1:] = tmp_df['par3'].to_numpy()
    par_4[1:] = tmp_df['par4'].to_numpy()
    par_5[1:] = tmp_df['par5'].to_numpy()
    alpha_ss[1:] = tmp_df['alphaSpecificStorage'].to_numpy()
    beta_ss[1:] = tmp_df['betaSpecificStorage'].to_numpy()
    ks[1:] = tmp_df['Ks'].to_numpy()
    
    return [equation_state_ID, parameter_ID, theta_s, theta_r, theta_wp, theta_fc, par_1, par_2, par_3, par_4, par_5,
           alpha_ss, beta_ss, ks]


def set_parameters_heat_diffusion(data_grid, data_parameter, data_dictionary, KMAX, eta):
    '''
    This function associate to each control volume a label that identifies 
    the rheology model, the set of parameters describing the soil type, and the max/min cell size 
    for regridding.
    
    :param data_grid: pandas dataframe containg the grid_input_file.csv
    :type data: pandas dataframe
    
    :param data_parameter: pandas dataframe containg the parameter_input_file.csv
    :type data: pandas dataframe
    
    :param data_dictionary: pandas dataframe containg a dictionary for the SWRC parameters
    :type data: pandas dataframe
    
    :param KMAX: number of control volumes.
    :type KMAX: int
    
    :param eta: vertical coordinate of control volume centroids. It is positive upward with 
        origin set at soil surface.
    :type eta: list
    
    return:
    
    equation_state_ID: vector containing the ID identifing the state equation of to model the k-th element
    type: array
    
    parameters_ID: vector containing the ID identifing the set of parameters of to model the k-th element
    type: array
        	
    soil_particles_density:
    type:array
    
    thermal_conductivity_soil_particles:
    type:array
    
    specific_thermal_capacity_soil_particles:
    type:array
    
    theta_s: vector containing the adimensional water content at saturation
    type:array
    
    theta_r: vector containing the residual adimensional water content
    type:array
    
    melting_temperature:
    type:array
    
    par_1: vector containing the parameter 1 of the SWRC model
    type:array
    
    par_2: vector containing the parameter 2 of the SWRC model
    type:array
    
    par_3: vector containing the parameter 3 of the SWRC model
    type:array
    
    par_4: vector containing the parameter 4 of the SWRC model
    type:array
    
    par_5: vector containing the parameter 5 of the SWRC model
    type:array
    
    alpha_ss: vector containing the acquitard compressibility 
    type:array
    
    beta_ss: vector containing the water compressibility
    type:array
    
    ks: vector containing the saturated hydraulic conductivity
    type:array
    '''
    equation_state_ID = np.zeros(KMAX, dtype=float)
    parameter_ID = np.zeros(KMAX, dtype=float)
    
    soil_particles_density = np.zeros(data_parameter.shape[0]+1, dtype=float)
    thermal_conductivity_soil_particles = np.zeros(data_parameter.shape[0]+1, dtype=float)
    specific_heat_capacity_soil_particles = np.zeros(data_parameter.shape[0]+1, dtype=float)
    theta_s = np.zeros(data_parameter.shape[0]+1, dtype=float)
    theta_r = np.zeros(data_parameter.shape[0]+1, dtype=float)
    melting_temperature = np.zeros(data_parameter.shape[0]+1, dtype=float)
    par_1 = np.zeros(data_parameter.shape[0]+1, dtype=float)
    par_2 = np.zeros(data_parameter.shape[0]+1, dtype=float)
    par_3 = np.zeros(data_parameter.shape[0]+1, dtype=float)
    par_4 = np.zeros(data_parameter.shape[0]+1, dtype=float)
    par_5 = np.zeros(data_parameter.shape[0]+1, dtype=float)
    alpha_ss = np.zeros(data_parameter.shape[0]+1, dtype=float)
    beta_ss = np.zeros(data_parameter.shape[0]+1, dtype=float)
    ks = np.zeros(data_parameter.shape[0]+1, dtype=float)
    
    coord_layer = []
    tmp_equation_state_ID = []
    tmp_parameter_ID = []

    for i in data_grid.index:

        if data_grid['Type'][i] == 'L':
            
            coord_layer.append(data_grid['eta'][i])
            tmp_equation_state_ID.append(data_grid['equationStateID'][i])
            tmp_parameter_ID.append(data_grid['parameterID'][i])
           
     
    for i in range(np.size(coord_layer)-1,0,-1):
        
        for j in range(0,KMAX):
            
            if(eta[j]>coord_layer[i] and eta[j]<coord_layer[i-1] ):
                
                equation_state_ID[j] = tmp_equation_state_ID[i-1]
                parameter_ID[j] = tmp_parameter_ID[i-1]
        
    data_parameter.rename(columns=dict(data_dictionary.values), inplace=True)   
    
    columns = ['spDensity', 'spConductivity', 'spSpecificHeatCapacity', 'thetaS', 'thetaR', 'par1', 'par2', 'par3', 'par4', 'par5', 'alphaSpecificStorage', 'betaSpecificStorage', 'Ks']
    tmp_df = pd.DataFrame(columns=columns)
    tmp_df = pd.concat([tmp_df, data_parameter],sort=False)
    tmp_df.replace(np.nan, -9999.0, inplace=True)
    
    soil_particles_density[1:] = tmp_df['spDensity'].to_numpy()
    thermal_conductivity_soil_particles[1:] = tmp_df['spConductivity'].to_numpy()
    specific_heat_capacity_soil_particles[1:] = tmp_df['spSpecificHeatCapacity'].to_numpy()
    theta_s[1:] = tmp_df['thetaS'].to_numpy()
    theta_r[1:] = tmp_df['thetaR'].to_numpy()
    par_1[1:] = tmp_df['par1'].to_numpy()
    par_2[1:] = tmp_df['par2'].to_numpy()
    par_3[1:] = tmp_df['par3'].to_numpy()
    par_4[1:] = tmp_df['par4'].to_numpy()
    par_5[1:] = tmp_df['par5'].to_numpy()
    alpha_ss[1:] = tmp_df['alphaSpecificStorage'].to_numpy()
    beta_ss[1:] = tmp_df['betaSpecificStorage'].to_numpy()
    ks[1:] = tmp_df['Ks'].to_numpy()

    return [equation_state_ID, parameter_ID, soil_particles_density, thermal_conductivity_soil_particles, specific_heat_capacity_soil_particles,
         theta_s, theta_r, melting_temperature, par_1, par_2, par_3, par_4, par_5, alpha_ss, beta_ss, ks]


def set_parameters_richards_heat_advection_diffusion(data_grid, data_parameter, data_dictionary, KMAX, eta):
    '''
    This function associate to each control volume a label that identifies 
    the rheology model, the set of parameters describing the soil type, and the max/min cell size 
    for regridding.
    
    :param data_grid: pandas dataframe containg the grid_input_file.csv
    :type data: pandas dataframe
    
    :param data_parameter: pandas dataframe containg the parameter_input_file.csv
    :type data: pandas dataframe
    
    :param data_dictionary: pandas dataframe containg a dictionary for the SWRC parameters
    :type data: pandas dataframe
    
    :param KMAX: number of control volumes.
    :type KMAX: int
    
    :param eta: vertical coordinate of control volume centroids. It is positive upward with 
        origin set at soil surface.
    :type eta: list
    
    return:
    
    equation_state_ID: vector containing the ID identifing the state equation of the k-th element
    type: array
    
    parameters_ID: vector containing the ID identifing the set of parameters of the k-th element
    type: array
    
    soil_particles_density:
    type:array
    
    thermal_conductivity_soil_particles:
    type:array
    
    specific_thermal_capacity_soil_particles:
    type:array
        	
    theta_s: vector containing the adimensional water content at saturation
    type:array
    
    theta_r: vector containing the residual adimensional water content
    type:array
    
    par_1: vector containing the parameter 1 of the SWRC model
    type:array
    
    par_2: vector containing the parameter 2 of the SWRC model
    type:array
    
    par_3: vector containing the parameter 3 of the SWRC model
    type:array
    
    par_4: vector containing the parameter 4 of the SWRC model
    type:array
    
    par_5: vector containing the parameter 5 of the SWRC model
    type:array
    
    alpha_ss: vector containing the acquitard compressibility 
    type:array
    
    beta_ss: vector containing the water compressibility
    type:array
    
    ks: vector containing the saturated hydraulic conductivity
    type:array
    '''
    equation_state_ID = np.zeros(KMAX, dtype=float)
    parameter_ID = np.zeros(KMAX, dtype=float)
    
    soil_particles_density = np.zeros(data_parameter.shape[0]+1, dtype=float)
    thermal_conductivity_soil_particles = np.zeros(data_parameter.shape[0]+1, dtype=float)
    specific_heat_capacity_soil_particles = np.zeros(data_parameter.shape[0]+1, dtype=float)
    theta_s = np.zeros(data_parameter.shape[0]+1, dtype=float)
    theta_r = np.zeros(data_parameter.shape[0]+1, dtype=float)
    par_1 = np.zeros(data_parameter.shape[0]+1, dtype=float)
    par_2 = np.zeros(data_parameter.shape[0]+1, dtype=float)
    par_3 = np.zeros(data_parameter.shape[0]+1, dtype=float)
    par_4 = np.zeros(data_parameter.shape[0]+1, dtype=float)
    par_5 = np.zeros(data_parameter.shape[0]+1, dtype=float)
    alpha_ss = np.zeros(data_parameter.shape[0]+1, dtype=float)
    beta_ss = np.zeros(data_parameter.shape[0]+1, dtype=float)
    ks = np.zeros(data_parameter.shape[0]+1, dtype=float)
    
    coord_layer = []
    tmp_equation_state_ID = []
    tmp_parameter_ID = []
    
    for i in data_grid.index:

        if data_grid['Type'][i] == 'L':
            
            coord_layer.append(data_grid['eta'][i])
            tmp_equation_state_ID.append(data_grid['equationStateID'][i])
            tmp_parameter_ID.append(data_grid['parameterID'][i])
           
     
    for i in range(np.size(coord_layer)-1,0,-1):
        
        for j in range(0,KMAX):
            
            if(eta[j]>coord_layer[i] and eta[j]<coord_layer[i-1] ):
                
                equation_state_ID[j] = tmp_equation_state_ID[i-1]
                parameter_ID[j] = tmp_parameter_ID[i-1]
        
    parameter_ID[KMAX-1] = parameter_ID[KMAX-2]   
    
    data_parameter.rename(columns=dict(data_dictionary.values), inplace=True)    

    columns = ['spDensity', 'spConductivity', 'spSpecificHeatCapacity', 'thetaS', 'thetaR', 'par1', 'par2', 'par3', 'par4', 'par5', 'alphaSpecificStorage', 'betaSpecificStorage', 'Ks']

    tmp_df = pd.DataFrame(columns=columns)
    tmp_df = pd.concat([tmp_df, data_parameter],sort=False)
    tmp_df.replace(np.nan, -9999.0, inplace=True)
    
    soil_particles_density[1:] = tmp_df['spDensity'].to_numpy()
    thermal_conductivity_soil_particles[1:] = tmp_df['spConductivity'].to_numpy()
    specific_heat_capacity_soil_particles[1:] = tmp_df['spSpecificHeatCapacity'].to_numpy()
    theta_s[1:] = tmp_df['thetaS'].to_numpy()
    theta_r[1:] = tmp_df['thetaR'].to_numpy()
    par_1[1:] = tmp_df['par1'].to_numpy()
    par_2[1:] = tmp_df['par2'].to_numpy()
    par_3[1:] = tmp_df['par3'].to_numpy()
    par_4[1:] = tmp_df['par4'].to_numpy()
    par_5[1:] = tmp_df['par5'].to_numpy()
    alpha_ss[1:] = tmp_df['alphaSpecificStorage'].to_numpy()
    beta_ss[1:] = tmp_df['betaSpecificStorage'].to_numpy()
    ks[1:] = tmp_df['Ks'].to_numpy()
    
    return [equation_state_ID, parameter_ID, soil_particles_density, thermal_conductivity_soil_particles, specific_heat_capacity_soil_particles, theta_s, theta_r, par_1, par_2, par_3, par_4, par_5,
           alpha_ss, beta_ss, ks]


def calibration_point_index(data_grid, eta):
    '''
    This function identifies the index of the calibration points.
    
    :param data_grid: pandas dataframe containg the grid_input_file.csv
    :type data: pandas dataframe

    :param eta: vertical coordinate of control volume centroids. It is positive upward with 
        origin set at soil surface.
    :type eta: list
    
    return:
    
    control_volume_index: vector containing the index of calibration (measurament) points)
    type: array
 
    '''
    
    measurament_eta = []
    for i in np.where(data_grid['Type']=='M'):
        measurament_eta = np.append(measurament_eta, data_grid['eta'][i])
    
    control_volume_index = []
    for i in measurament_eta:
        control_volume_index = np.append(control_volume_index, np.where(eta==i))

    return control_volume_index


