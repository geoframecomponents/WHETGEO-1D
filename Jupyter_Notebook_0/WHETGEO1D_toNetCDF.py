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


def write_grid_netCDF_richards(eta, eta_dual, z, z_dual, space_delta, control_volume, control_volume_index, psi_0, T_0, equation_state_ID, parameter_ID, KMAX, theta_s, theta_r, par_1, par_2, par_3, par_4, par_5, alpha_ss, beta_ss, ks,
					  output_file_name, output_title, output_institution, output_summary, output_date,
					  grid_input_file_name, parameter_input_file_name):
    '''
    Save all grid data in a NetCDF file for Richards simulation
    
    :param eta: vertical coordinate of control volume centroids. It is positive upward with.
        origin set at soil surface.
    :type eta: array
    
    :param eta_dual: vertical coordinate of control volume interface. It is positive upward with.
        origin set at soil surface.
    :type eta_dual: array
    
    :param z: vertical coordinate of control volume centroids. It is positive upward with.
        origin set at soil column bottom.
    :type z: array
    
    :param z_dual: vertical coordinate of control volume interfaces. It is positive upward with.
        origin set at soil column bottom.
    :type z_dual: array
    
    :param space_delta: is the distance between two adjacent control volumes.
        This quantity is used to compute gradients
    :type space_delta: array
    
    :param control_volume: dimension of each control volume.
    :type control_volume: array.
    
    :param control_volume_index: index of the calibration points. 
    :type control_volume_index: array.
    
    :param psi_0: water suction initial condition.
    :type ic: array.
    
    :param T_0: temperature initial condition.
    :type ic: array.
    
    :param equation_state_ID: containing a label for each control volume defining the equation state of the k-th element.
    :type equation_state_ID: array.
    
    :param parameter_ID: containing a label for each control volume defining the parameter set to be used for the k-th element.
    :type parameter_ID: array.
    	
    :param KMAX: number of control volumes.
    :type KMAX: int.
    	
    :param theta_s: vector containing the adimensional water content at saturation
    :type theta_s: array.
    
    :param theta_r: vector containing the residual adimensional water content
    :type theta_r: array.
    
    :param par_1: vector containing the parameter 1 of the SWRC model
    :type par_1: array.
    
    :param par_2: vector containing the parameter 2 of the SWRC model
    :type par_2: array.
    
    :param par_3: vector containing the parameter 3 of the SWRC model
    :type par_3: array.
    
    :param par_4: vector containing the parameter 4 of the SWRC model
    :type par_4: array.
    
    :param par_5: vector containing the parameter 5 of the SWRC model
    :type par_5: array.
    
    :param alpha_ss: vector containing the acquitard compressibility
    :type alpha_ss: array.
    
    :param beta_ss: vector containing the water compressibility
    :type beta_ss: array.
    
    :param ks: vector containing the saturated hydraulic conductivity
    :type ks: array.

    :param output_file_name: 
    :type output_file_name: str
    
    :param output_title: 
    :type output_title: str
    
    :param output_institution: 
    :type param output_institution: str
    
    :param output_summary: 
    :type output_summary: str
    
    :param output_date: 
    :type output_date: str
    
    :param input_file_name: 
    :type input_file_name: str

    '''

    # the output array to write will be nx x ny
    dim = np.size(eta);
    dim1 = np.size(eta_dual);
    dim2 = np.size(space_delta);
    dim_parameter = np.size(par_1) 
    dim_control_volume_index = max(np.size(control_volume_index),1)
    
    # open a new netCDF file for writing.
    ncfile = Dataset(output_file_name,'w') 
    
    # Create global attributes
    ncfile.title = output_title + '\\n' + 'grid input file' + grid_input_file_name + 'parameter input file' + parameter_input_file_name
    ncfile.institution =  output_institution
    ncfile.summary = output_summary
    #ncfile.acknowledgment = ""
    ncfile.date_created = output_date
    
    # create the z dimensions.
    ncfile.createDimension('z',dim)
    ncfile.createDimension('z_dual',dim1)
    ncfile.createDimension('space_delta',dim2)
    ncfile.createDimension('parameter',dim_parameter)
    ncfile.createDimension('scalar',1)
    ncfile.createDimension('control_volume_index',dim_control_volume_index)
    
    # create the variable
    # first argument is name of variable, second is datatype, third is
    # a tuple with the names of dimensions.
    data_KMAX = ncfile.createVariable('KMAX','i4',('scalar'))
    data_KMAX.unit = '-'
    
    data_eta = ncfile.createVariable('eta','f8',('z'))
    data_eta.unit = 'm'
    data_eta.long_name = '\u03b7 coordinate of volume centroids: zero is at soil surface and and positive upward'
    
    data_eta_dual = ncfile.createVariable('etaDual','f8',('z_dual'))
    data_eta_dual.unit = 'm'
    data_eta_dual.long_name = '\u03b7 coordinate of volume interfaces: zero is at soil surface and and positive upward. '
    
    data_z = ncfile.createVariable('z','f8',('z'))
    data_z.unit = 'm'
    data_z.long_name = 'z coordinate  of volume centroids: zero is at the bottom of the column and and positive upward'
    
    data_z_dual = ncfile.createVariable('zDual','f8',('z_dual'))
    data_z_dual.unit = 'm'
    data_z_dual.long_name = 'z coordinate of volume interfaces: zero is at soil surface and and positive upward.'
    
    data_psi_0 = ncfile.createVariable('psi0','f8',('z'))
    data_psi_0.units = 'm'
    data_psi_0.long_name = 'Water suction initial condition'
    
    data_T_0 = ncfile.createVariable('T0','f8',('z'))
    data_T_0.units = 'K'
    data_T_0.long_name = 'Temperature initial condition'
    	
    data_space_delta = ncfile.createVariable('spaceDelta','f8',('space_delta'))
    data_space_delta.unit = 'm'
    data_space_delta.long_name = 'Distance between consecutive controids, is used to compute gradients'
    
    data_control_volume = ncfile.createVariable('controlVolume','f8',('z'))
    data_control_volume.unit = 'm'
    data_control_volume.long_name = 'Control volume size'
    
    data_control_volume_index = ncfile.createVariable('controlVolumeIndex','i4',('control_volume_index'))
    data_control_volume_index.unit = ''
    data_control_volume_index.long_name = 'Index of control volumes for calibration'
    
    data_equation_state_ID = ncfile.createVariable('equationStateID','i4',('z'))
    data_equation_state_ID.units = '-'
    data_equation_state_ID.long_name = 'label describing the equation state of the k-th element'
    
    data_parameter_ID = ncfile.createVariable('parameterID','i4',('z'))
    data_parameter_ID.units = '-'
    data_parameter_ID.long_name = 'label identifying the set of parameters'

    data_theta_s = ncfile.createVariable('thetaS','f8',('parameter'))
    data_theta_s.units = '-'
    data_theta_s.long_name = 'adimensional water content at saturation'
    
    data_theta_r = ncfile.createVariable('thetaR','f8',('parameter'))
    data_theta_r.units = '-'
    data_theta_r.long_name = 'residual adimensional water content'

    data_par_1 = ncfile.createVariable('par1SWRC','f8',('parameter'))
    data_par_1.units = '-'
    data_par_1.long_name = 'SWRC parameter'
    
    data_par_2 = ncfile.createVariable('par2SWRC','f8',('parameter'))
    data_par_2.units = '-'
    data_par_2.long_name = 'SWRC parameter'
    
    data_par_3 = ncfile.createVariable('par3SWRC','f8',('parameter'))
    data_par_3.units = '-'
    data_par_3.long_name = 'SWRC parameter'
    
    data_par_4 = ncfile.createVariable('par4SWRC','f8',('parameter'))
    data_par_4.units = '-'
    data_par_4.long_name = 'SWRC parameter'
    
    data_par_5 = ncfile.createVariable('par5SWRC','f8',('parameter'))
    data_par_5.units = '-'
    data_par_5.long_name = 'SWRC parameter'
    
    data_alpha_ss = ncfile.createVariable('alphaSpecificStorage','f8',('parameter'))
    data_alpha_ss.units = '1/Pa'
    data_alpha_ss.long_name = 'acquitard compressibility'
    
    data_beta_ss = ncfile.createVariable('betaSpecificStorage','f8',('parameter'))
    data_beta_ss.units = '1/Pa'
    data_beta_ss.long_name = 'water compressibility'
    
    data_ks = ncfile.createVariable('ks','f8',('parameter'))
    data_ks.units = 'm/s'
    data_ks.long_name = 'saturated hydraulic conductivity'
    
    ## write data to variable.

    data_KMAX[0] = KMAX

    for i in range(0,dim):
        data_eta[i] = eta[i]
        data_z[i] = z[i]
        data_control_volume[i] = control_volume[i]
        data_psi_0[i] = psi_0[i]
        data_T_0[i] = T_0[i]
        data_equation_state_ID[i] = equation_state_ID[i]
        data_parameter_ID[i] = parameter_ID[i]
    
    for i in range(0,dim1):
        data_eta_dual[i] = eta_dual[i]
        data_z_dual[i] = z_dual[i]
        data_space_delta[i] = space_delta[i]
        
    for i in range(0,dim2):
        data_space_delta[i] = space_delta[i]
        
    for i in range(0,dim_parameter):
        data_theta_s[i] = theta_s[i]
        data_theta_r[i] = theta_r[i]
        data_par_1[i] = par_1[i]
        data_par_2[i] = par_2[i]
        data_par_3[i] = par_3[i]
        data_par_4[i] = par_4[i]
        data_par_5[i] = par_5[i]
        data_alpha_ss[i] = alpha_ss[i]
        data_beta_ss[i] = beta_ss[i]
        data_ks[i] = ks[i]
        
#     if not control_volume_index:
    if control_volume_index is None:
        data_control_volume_index = -9999
    else:
        for i in range(0,len(control_volume_index)):
            data_control_volume_index[i] = control_volume_index[i]
    
    ## close the file.
    ncfile.close()
    print ('\n\n***SUCCESS writing!  '+ output_file_name)

    
    return


def write_grid_netCDF_richards_lysimeter(eta, eta_dual, z, z_dual, space_delta, control_volume, control_volume_index, psi_0, T_0, equation_state_ID, parameter_ID, KMAX, theta_s,              
					  theta_r, theta_wp, theta_fc, par_1, par_2, par_3, par_4, par_5, alpha_ss, beta_ss, ks,
					  output_file_name, output_title, output_institution, output_summary, output_date,
					  grid_input_file_name, parameter_input_file_name):
    '''
    Save all grid data in a NetCDF file for Richards with lysimeter simulation
    
    :param eta: vertical coordinate of control volume centroids. It is positive upward with.
        origin set at soil surface.
    :type eta: array
    
    :param eta_dual: vertical coordinate of control volume interface. It is positive upward with.
        origin set at soil surface.
    :type eta_dual: array
    
    :param z: vertical coordinate of control volume centroids. It is positive upward with.
        origin set at soil column bottom.
    :type z: array
    
    :param z_dual: vertical coordinate of control volume interfaces. It is positive upward with.
        origin set at soil column bottom.
    :type z_dual: array
    
    :param space_delta: is the distance between two adjacent control volumes.
        This quantity is used to compute gradients
    :type space_delta: array
    
    :param control_volume: dimension of each control volume.
    :type control_volume: array.
    
    :param control_volume_index: index of the calibration points. 
    :type control_volume_index: array.
    
    :param psi_0: water suction initial condition.
    :type ic: array.
    
    :param T_0: temperature initial condition.
    :type ic: array.
    
    :param equation_state_ID: containing a label for each control volume defining the type of the equation state of the k-th element.
    :type equation_state_ID: array.
    
    :param parameter_ID: containing a label for each control volume defining the parameter set to be used.
    :type parameter_ID: array.
    	
    :param KMAX: number of control volumes.
    :type KMAX: int.
    	
    :param theta_s: vector containing the adimensional water content at saturation
    :type theta_s: array.
    
    :param theta_r: vector containing the residual adimensional water content
    :type theta_r: array.
    
    :param theta_wp: vector containing the adimensional water content at the wilting point
    :type theta_wp: array.
    
    :param theta_fc: vector containing the adimensional water content of the field capacity
    :type theta_fc: array.
    
    :param par_1: vector containing the parameter 1 of the SWRC model
    :type par_1: array.
    
    :param par_2: vector containing the parameter 2 of the SWRC model
    :type par_2: array.
    
    :param par_3: vector containing the parameter 3 of the SWRC model
    :type par_3: array.
    
    :param par_4: vector containing the parameter 4 of the SWRC model
    :type par_4: array.
    
    :param par_5: vector containing the parameter 5 of the SWRC model
    :type par_5: array.
    
    :param alpha_ss: vector containing the acquitard compressibility
    :type alpha_ss: array.
    
    :param beta_ss: vector containing the water compressibility
    :type beta_ss: array.
    
    :param ks: vector containing the saturated hydraulic conductivity
    :type ks: array.

    :param output_file_name: 
    :type output_file_name: str
    
    :param output_title: 
    :type output_title: str
    
    :param output_institution: 
    :type param output_institution: str
    
    :param output_summary: 
    :type output_summary: str
    
    :param output_date: 
    :type output_date: str
    
    :param input_file_name: 
    :type input_file_name: str

    '''

    # the output array to write will be nx x ny
    dim = np.size(eta);
    dim1 = np.size(eta_dual);
    dim2 = np.size(space_delta);
    dim_parameter = np.size(par_1) 
    dim_control_volume_index = min(np.size(control_volume_index),1)
    
    # open a new netCDF file for writing.
    ncfile = Dataset(output_file_name,'w') 
    
    # Create global attributes
    ncfile.title = output_title + '\\n' + 'grid input file' + grid_input_file_name + 'parameter input file' + parameter_input_file_name
    ncfile.institution =  output_institution
    ncfile.summary = output_summary
    #ncfile.acknowledgment = ""
    ncfile.date_created = output_date
    
    # create the z dimensions.
    ncfile.createDimension('z',dim)
    ncfile.createDimension('z_dual',dim1)
    ncfile.createDimension('space_delta',dim2)
    ncfile.createDimension('parameter',dim_parameter)
    ncfile.createDimension('scalar',1)
    ncfile.createDimension('control_volume_index',dim_control_volume_index)
    
    # create the variable
    # first argument is name of variable, second is datatype, third is
    # a tuple with the names of dimensions.
    data_KMAX = ncfile.createVariable('KMAX','i4',('scalar'))
    data_KMAX.unit = '-'
    
    data_eta = ncfile.createVariable('eta','f8',('z'))
    data_eta.unit = 'm'
    data_eta.long_name = '\u03b7 coordinate of volume centroids: zero is at soil surface and and positive upward'
    
    data_eta_dual = ncfile.createVariable('etaDual','f8',('z_dual'))
    data_eta_dual.unit = 'm'
    data_eta_dual.long_name = '\u03b7 coordinate of volume interfaces: zero is at soil surface and and positive upward. '
    
    data_z = ncfile.createVariable('z','f8',('z'))
    data_z.unit = 'm'
    data_z.long_name = 'z coordinate  of volume centroids: zero is at the bottom of the column and and positive upward'
    
    data_z_dual = ncfile.createVariable('zDual','f8',('z_dual'))
    data_z_dual.unit = 'm'
    data_z_dual.long_name = 'z coordinate of volume interfaces: zero is at soil surface and and positive upward.'
    
    data_psi_0 = ncfile.createVariable('psi0','f8',('z'))
    data_psi_0.units = 'm'
    data_psi_0.long_name = 'Water suction initial condition'
    
    data_T_0 = ncfile.createVariable('T0','f8',('z'))
    data_T_0.units = 'K'
    data_T_0.long_name = 'Temperature initial condition'
    	
    data_space_delta = ncfile.createVariable('spaceDelta','f8',('space_delta'))
    data_space_delta.unit = 'm'
    data_space_delta.long_name = 'Distance between consecutive controids, is used to compute gradients'
    
    data_control_volume = ncfile.createVariable('controlVolume','f8',('z'))
    data_control_volume.unit = 'm'
    data_control_volume.long_name = 'Control volume size'
    
    data_control_volume_index = ncfile.createVariable('controlVolumeIndex','i4',('control_volume_index'))
    data_control_volume_index.unit = ''
    data_control_volume_index.long_name = 'Index of control volumes for calibration'
    
    data_equation_state_ID = ncfile.createVariable('equationStateID','i4',('z'))
    data_equation_state_ID.units = '-'
    data_equation_state_ID.long_name = 'label describing the equation state of the k-th element'
    
    data_parameter_ID = ncfile.createVariable('parameterID','i4',('z'))
    data_parameter_ID.units = '-'
    data_parameter_ID.long_name = 'label identifying the set of parameters'

    data_theta_s = ncfile.createVariable('thetaS','f8',('parameter'))
    data_theta_s.units = '-'
    data_theta_s.long_name = 'adimensional water content at saturation'
    
    data_theta_r = ncfile.createVariable('thetaR','f8',('parameter'))
    data_theta_r.units = '-'
    data_theta_r.long_name = 'residual adimensional water content'
    
    data_theta_wp = ncfile.createVariable('thetaWP','f8',('parameter'))
    data_theta_wp.units = '-'
    data_theta_wp.long_name = 'adimensional water content at wilting point'
    
    data_theta_fc = ncfile.createVariable('thetaFC','f8',('parameter'))
    data_theta_fc.units = '-'
    data_theta_fc.long_name = 'adimensional water content at field capacity'

    data_par_1 = ncfile.createVariable('par1SWRC','f8',('parameter'))
    data_par_1.units = '-'
    data_par_1.long_name = 'SWRC parameter'
    
    data_par_2 = ncfile.createVariable('par2SWRC','f8',('parameter'))
    data_par_2.units = '-'
    data_par_2.long_name = 'SWRC parameter'
    
    data_par_3 = ncfile.createVariable('par3SWRC','f8',('parameter'))
    data_par_3.units = '-'
    data_par_3.long_name = 'SWRC parameter'
    
    data_par_4 = ncfile.createVariable('par4SWRC','f8',('parameter'))
    data_par_4.units = '-'
    data_par_4.long_name = 'SWRC parameter'
    
    data_par_5 = ncfile.createVariable('par5SWRC','f8',('parameter'))
    data_par_5.units = '-'
    data_par_5.long_name = 'SWRC parameter'
    
    data_alpha_ss = ncfile.createVariable('alphaSpecificStorage','f8',('parameter'))
    data_alpha_ss.units = '1/Pa'
    data_alpha_ss.long_name = 'acquitard compressibility'
    
    data_beta_ss = ncfile.createVariable('betaSpecificStorage','f8',('parameter'))
    data_beta_ss.units = '1/Pa'
    data_beta_ss.long_name = 'water compressibility'
    
    data_ks = ncfile.createVariable('ks','f8',('parameter'))
    data_ks.units = 'm/s'
    data_ks.long_name = 'saturated hydraulic conductivity'
    
    ## write data to variable.

    data_KMAX[0] = KMAX

    for i in range(0,dim):
        data_eta[i] = eta[i]
        data_z[i] = z[i]
        data_control_volume[i] = control_volume[i]
        data_psi_0[i] = psi_0[i]
        data_T_0[i] = T_0[i]
        data_equation_state_ID[i] = equation_state_ID[i]
        data_parameter_ID[i] = parameter_ID[i]
    
    for i in range(0,dim1):
        data_eta_dual[i] = eta_dual[i]
        data_z_dual[i] = z_dual[i]
        data_space_delta[i] = space_delta[i]
        
    for i in range(0,dim2):
        data_space_delta[i] = space_delta[i]
        
    for i in range(0,dim_parameter):
        data_theta_s[i] = theta_s[i]
        data_theta_r[i] = theta_r[i]
        data_theta_wp[i] = theta_wp[i]
        data_theta_fc[i] = theta_fc[i]
        data_par_1[i] = par_1[i]
        data_par_2[i] = par_2[i]
        data_par_3[i] = par_3[i]
        data_par_4[i] = par_4[i]
        data_par_5[i] = par_5[i]
        data_alpha_ss[i] = alpha_ss[i]
        data_beta_ss[i] = beta_ss[i]
        data_ks[i] = ks[i]
        
    if not control_volume_index:
        data_control_volume_index[0] = -9999
    else:
        for i in range(0,len(control_volume_index)):
            data_control_volume_index[i] = control_volume_index[i]
    
    ## close the file.
    ncfile.close()
    print ('\n\n***SUCCESS writing!  '+ output_file_name)

    
    return

def write_grid_netCDF_heat_diffusion(eta, eta_dual, z, z_dual, space_delta, control_volume, control_volume_index, psi_0, T_0, equation_state_ID, parameter_ID, KMAX, soil_particles_density,              
    				  thermal_conductivity_soil_particles, 
    				  specific_heat_capacity_soil_particles, theta_s, theta_r, melting_temperature, par_1, par_2, par_3, par_4, par_5, alpha_ss, beta_ss, ks,
    				  output_file_name, output_title, output_institution, output_summary, output_date,
    				  grid_input_file_name, parameter_input_file_name):
    '''
    Save all grid data in a NetCDF file
    
    :param eta: vertical coordinate of control volume centroids. It is positive upward with.
        origin set at soil surface.
    :type eta: array
    
    :param eta_dual: vertical coordinate of control volume interface. It is positive upward with.
        origin set at soil surface.
    :type eta_dual: array
    
    :param z: vertical coordinate of control volume centroids. It is positive upward with.
        origin set at soil column bottom.
    :type z: array
    
    :param z_dual: vertical coordinate of control volume interfaces. It is positive upward with.
        origin set at soil column bottom.
    :type z_dual: array
    
    :param space_delta: is the distance between two adjacent control volumes.
        This quantity is used to compute gradients
    :type space_delta: array
    
    :param control_volume: dimension of each control volume.
    :type control_volume: array.
    
    :param control_volume_index: index of the calibration points. 
    :type control_volume_index: array.
    
    :param psi_0: water suction initial condition.
    :type ic: array.
    
    :param T_0: temperature initial condition.
    :type ic: array.
    
    :param equation_state_ID: containing a label for each control volume defining the type of the equation state of the k-th element.
    :type equation_state_ID: array.
    
    :param parameter_ID: containing a label for each control volume defining the parameter set to be used.
    :type parameter_ID: array.
    	
    :param KMAX: number of control volumes.
    :type KMAX: int.
    	
    :param soil_particles_density: array containing the soil particles density.
    :type soil_particles_density: array.
    
    :param thermal_conductivity_soil_particles: array containing the thermal conductivity of soil particles.
    :type thermal_conductivity_soil_particles: array.
    
    :param specific_heat_capacity_soil_particles: array containing the specific heat capacity of soil particles.
    :type specific_heat_capacity_soil_particles: array.
    
    :param theta_s: array containing the values of the water content at saturation.
    :type theta_s: array.
    
    :param theta_r: array containing the values of the residual water content.
    :type theta_r: array.

    :param melting_temperature: array containing the melting temperature.
    :type melting_temperature: array.
    
    :param par1: array containing the values of the SFCC parameter.
    :type par1: array.
    
    :param par2: array containing the values of the SFCC parameter.
    :type par2: array.
    
    :param par3: array containing the values of the SFCC parameter.
    :type par3: array.
    
    :param par4: array containing the values of the SFCC parameter.
    :type par4: array.
    
    :param par5: array containing the values of the SFCC parameter.
    :type par5: array.
    
    :param alpha_ss: vector containing the acquitard compressibility
    :type alpha_ss: array.
    
    :param beta_ss: vector containing the water compressibility
    :type beta_ss: array.
    
    :param ks: vector containing the saturated hydraulic conductivity
    :type ks: array.

    :param output_file_name: 
    :type output_file_name: str
    
    :param output_title: 
    :type output_title: str
    
    :param output_institution: 
    :type param output_institution: str
    
    :param output_summary: 
    :type output_summary: str
    
    :param output_date: 
    :type output_date: str
    
    :param input_file_name: 
    :type input_file_name: str
    
    '''
    
    # the output array to write will be nx x ny
    dim = np.size(eta)
    dim1 = np.size(eta_dual)
    dim2 = np.size(space_delta)
    dim_parameter = np.size(par_1) 
    dim_control_volume_index = min(np.size(control_volume_index),1)
    print(dim1)
    
    
    # open a new netCDF file for writing.
    ncfile = Dataset(output_file_name,'w') 
    
    # Create global attributes
    ncfile.title = output_title + '\\n' + 'grid input file' + grid_input_file_name + 'parameter input file' + parameter_input_file_name
    ncfile.institution =  output_institution
    ncfile.summary = output_summary
    #ncfile.acknowledgment = ""
    ncfile.date_created = output_date
      
    # create the z dimensions.
    ncfile.createDimension('z',dim)
    ncfile.createDimension('z_dual',dim1)
    ncfile.createDimension('space_delta',dim2)
    ncfile.createDimension('parameter',dim_parameter)
    ncfile.createDimension('scalar',1)
    ncfile.createDimension('control_volume_index',dim_control_volume_index)
    
    # create the variable
    # first argument is name of variable, second is datatype, third is
    # a tuple with the names of dimensions.
    data_KMAX = ncfile.createVariable('KMAX','i4',('scalar'))
    data_KMAX.unit = '-'
    
    data_eta = ncfile.createVariable('eta','f8',('z'))
    data_eta.unit = 'm'
    data_eta.long_name = '\u03b7 coordinate of volume centroids: zero is at soil surface and and positive upward'
    
    data_eta_dual = ncfile.createVariable('etaDual','f8',('z_dual'))
    data_eta_dual.unit = 'm'
    data_eta_dual.long_name = '\u03b7 coordinate of volume interfaces: zero is at soil surface and and positive upward. '
    
    data_z = ncfile.createVariable('z','f8',('z'))
    data_z.unit = 'm'
    data_z.long_name = 'z coordinate  of volume centroids: zero is at the bottom of the column and and positive upward'
    
    data_z_dual = ncfile.createVariable('zDual','f8',('z_dual'))
    data_z_dual.unit = 'm'
    data_z_dual.long_name = 'z coordinate of volume interfaces: zero is at soil surface and and positive upward.'

    data_psi_0 = ncfile.createVariable('psi0','f8',('z'))
    data_psi_0.units = 'm'
    data_psi_0.long_name = 'Water suction initial condition'
    
    data_T_0 = ncfile.createVariable('T0','f8',('z'))
    data_T_0.units = 'K'
    data_T_0.long_name = 'Temperature initial condition'
    	
    data_space_delta = ncfile.createVariable('spaceDelta','f8',('z_dual'))
    data_space_delta.unit = 'm'
    data_space_delta.long_name = 'Distance between consecutive controids, is used to compute gradients'
    
    data_control_volume = ncfile.createVariable('controlVolume','f8',('z'))
    data_control_volume.unit = 'm'
    data_control_volume.long_name = 'Control volume size'
    
    data_control_volume_index = ncfile.createVariable('controlVolumeIndex','i4',('control_volume_index'))
    data_control_volume_index.unit = ''
    data_control_volume_index.long_name = 'Index of control volumes for calibration'
    
    data_equation_state_ID = ncfile.createVariable('equationStateID','i4',('z'))
    data_equation_state_ID.units = '-'
    data_equation_state_ID.long_name = 'label describing the equation state of the k-th element'
    
    data_parameter_ID = ncfile.createVariable('parameterID','i4',('z'))
    data_parameter_ID.units = '-'
    data_parameter_ID.long_name = 'label identifying the set of parameters'

    data_soil_particles_density = ncfile.createVariable('soilParticlesDensity','f8',('parameter'))
    data_soil_particles_density.units = 'kg/m3'
    data_soil_particles_density.long_name = 'density of soil particles'
    
    data_thermal_conductivity_soil_particles = ncfile.createVariable('thermalConductivitySoilParticles','f8',('parameter'))
    data_thermal_conductivity_soil_particles.units = 'W/m2'
    data_thermal_conductivity_soil_particles.long_name = 'thermal conductivity of soil particles'
    
    data_specific_heat_capacity_soil_particles = ncfile.createVariable('specificThermalCapacitySoilParticles','f8',('parameter'))
    data_specific_heat_capacity_soil_particles.units = 'J/kg m3'
    data_specific_heat_capacity_soil_particles.long_name = 'specific thermal capacity of soil particles'
    
    data_theta_s = ncfile.createVariable('thetaS','f8',('parameter'))
    data_theta_s.units = '-'
    data_theta_s.long_name = 'adimensional water content at saturation'
    
    data_theta_r = ncfile.createVariable('thetaR','f8',('parameter'))
    data_theta_r.units = '-'
    data_theta_r.long_name = 'adimensional residual water content'
    
    data_melting_temperature = ncfile.createVariable('meltingTemperature','f8',('parameter'))
    data_melting_temperature.units = 'K'
    data_melting_temperature.long_name = 'melting temperature of soil water'
    
    data_par_1 = ncfile.createVariable('par1SWRC','f8',('parameter'))
    data_par_1.units = '-'
    data_par_1.long_name = 'SWRC parameter'
    
    data_par_2 = ncfile.createVariable('par2SWRC','f8',('parameter'))
    data_par_2.units = '-'
    data_par_2.long_name = 'SWRC parameter'
    
    data_par_3 = ncfile.createVariable('par3SWRC','f8',('parameter'))
    data_par_3.units = '-'
    data_par_3.long_name = 'SWRC parameter'
    
    data_par_4 = ncfile.createVariable('par4SWRC','f8',('parameter'))
    data_par_4.units = '-'
    data_par_4.long_name = 'SWRC parameter'
    
    data_par_5 = ncfile.createVariable('par5SWRC','f8',('parameter'))
    data_par_5.units = '-'
    data_par_5.long_name = 'SWRC parameter'
    
    data_alpha_ss = ncfile.createVariable('alphaSpecificStorage','f8',('parameter'))
    data_alpha_ss.units = '1/Pa'
    data_alpha_ss.long_name = 'acquitard compressibility'
    
    data_beta_ss = ncfile.createVariable('betaSpecificStorage','f8',('parameter'))
    data_beta_ss.units = '1/Pa'
    data_beta_ss.long_name = 'water compressibility'
    
    data_ks = ncfile.createVariable('ks','f8',('parameter'))
    data_ks.units = 'm/s'
    data_ks.long_name = 'saturated hydraulic conductivity'
        
    
    ## write data to variable.

    data_KMAX[0] = KMAX

    for i in range(0,dim):
        data_eta[i] = eta[i]
        data_z[i] = z[i]
        data_control_volume[i] = control_volume[i]
        data_psi_0[i] = psi_0[i]
        data_T_0[i] = T_0[i]
        data_equation_state_ID[i] = equation_state_ID[i]
        data_parameter_ID[i] = parameter_ID[i]
    	
    for i in range(0,dim1):
        data_eta_dual[i] = eta_dual[i]
        data_z_dual[i] = z_dual[i]
        data_space_delta[i] = space_delta[i]
    	
    for i in range(0,dim_parameter):
        data_soil_particles_density[i] = soil_particles_density[i]
        data_thermal_conductivity_soil_particles[i] = thermal_conductivity_soil_particles[i]
        data_specific_heat_capacity_soil_particles[i] = specific_heat_capacity_soil_particles[i]
        data_theta_s[i] = theta_s[i]
        data_theta_r[i] = theta_r[i]
        data_melting_temperature[i] = melting_temperature[i]
        data_par_1[i] = par_1[i]
        data_par_2[i] = par_2[i]
        data_par_3[i] = par_3[i]
        data_par_4[i] = par_4[i]
        data_par_5[i] = par_5[i]
        data_alpha_ss[i] = alpha_ss[i]
        data_beta_ss[i] = beta_ss[i]
        data_ks[i] = ks[i]
        
    print(par_1)
    ## close the file.
    ncfile.close()
    print ('\n\n***SUCCESS writing!  '+ output_file_name)

	
    return


def write_grid_netCDF_richards_heat_advection_diffusion(eta, eta_dual, z, z_dual, space_delta, control_volume, control_volume_index, psi_0, T_0, equation_state_ID, parameter_ID, KMAX, soil_particles_density,
thermal_conductivity_soil_particles, 
specific_heat_capacity_soil_particles,theta_s, theta_r, par_1, par_2, par_3, par_4, par_5, alpha_ss, beta_ss, ks,
					  output_file_name, output_title, output_institution, output_summary, output_date,
					  grid_input_file_name, parameter_input_file_name):
    '''
    Save all grid data in a NetCDF file for Richards simulation
    
    :param eta: vertical coordinate of control volume centroids. It is positive upward with.
        origin set at soil surface.
    :type eta: array
    
    :param eta_dual: vertical coordinate of control volume interface. It is positive upward with.
        origin set at soil surface.
    :type eta_dual: array
    
    :param z: vertical coordinate of control volume centroids. It is positive upward with.
        origin set at soil column bottom.
    :type z: array
    
    :param z_dual: vertical coordinate of control volume interfaces. It is positive upward with.
        origin set at soil column bottom.
    :type z_dual: array
    
    :param space_delta: is the distance between two adjacent control volumes.
        This quantity is used to compute gradients
    :type space_delta: array
    
    :param control_volume: dimension of each control volume.
    :type control_volume: array.
    
    :param control_volume_index: index of the calibration points. 
    :type control_volume_index: array.
    
    :param psi_0: water suction initial condition.
    :type ic: array.
    
    :param T_0: temperature initial condition.
    :type ic: array.
    
    :param equation_state_ID: containing a label for each control volume defining the equation state of the k-th element.
    :type equation_state_ID: array.
    
    :param parameter_ID: containing a label for each control volume defining the parameter set to be used for the k-th element.
    :type parameter_ID: array.
    	
    :param KMAX: number of control volumes.
    :type KMAX: int.
    
    :param soil_particles_density: array containing the soil particles density.
    :type soil_particles_density: array.
    
    :param thermal_conductivity_soil_particles: array containing the thermal conductivity of soil particles.
    :type thermal_conductivity_soil_particles: array.
    
    :param specific_heat_capacity_soil_particles: array containing the specific heat capacity of soil particles.
    :type specific_heat_capacity_soil_particles: array.
    	
    :param theta_s: vector containing the adimensional water content at saturation
    :type theta_s: array.
    
    :param theta_r: vector containing the residual adimensional water content
    :type theta_r: array.
    
    :param par_1: vector containing the parameter 1 of the SWRC model
    :type par_1: array.
    
    :param par_2: vector containing the parameter 2 of the SWRC model
    :type par_2: array.
    
    :param par_3: vector containing the parameter 3 of the SWRC model
    :type par_3: array.
    
    :param par_4: vector containing the parameter 4 of the SWRC model
    :type par_4: array.
    
    :param par_5: vector containing the parameter 5 of the SWRC model
    :type par_5: array.
    
    :param alpha_ss: vector containing the acquitard compressibility
    :type alpha_ss: array.
    
    :param beta_ss: vector containing the water compressibility
    :type beta_ss: array.
    
    :param ks: vector containing the saturated hydraulic conductivity
    :type ks: array.

    :param output_file_name: 
    :type output_file_name: str
    
    :param output_title: 
    :type output_title: str
    
    :param output_institution: 
    :type param output_institution: str
    
    :param output_summary: 
    :type output_summary: str
    
    :param output_date: 
    :type output_date: str
    
    :param input_file_name: 
    :type input_file_name: str

    '''

    # the output array to write will be nx x ny
    dim = np.size(eta);
    dim1 = np.size(eta_dual);
    dim2 = np.size(space_delta);
    dim_parameter = np.size(par_1) 
    dim_control_volume_index = max(np.size(control_volume_index),1)
    
    # open a new netCDF file for writing.
    ncfile = Dataset(output_file_name,'w') 
    
    # Create global attributes
    ncfile.title = output_title + '\\n' + 'grid input file' + grid_input_file_name + 'parameter input file' + parameter_input_file_name
    ncfile.institution =  output_institution
    ncfile.summary = output_summary
    #ncfile.acknowledgment = ""
    ncfile.date_created = output_date
    
    # create the z dimensions.
    ncfile.createDimension('z',dim)
    ncfile.createDimension('z_dual',dim1)
    ncfile.createDimension('space_delta',dim2)
    ncfile.createDimension('parameter',dim_parameter)
    ncfile.createDimension('scalar',1)
    ncfile.createDimension('control_volume_index',dim_control_volume_index)
    
    # create the variable
    # first argument is name of variable, second is datatype, third is
    # a tuple with the names of dimensions.
    data_KMAX = ncfile.createVariable('KMAX','i4',('scalar'))
    data_KMAX.unit = '-'
    
    data_eta = ncfile.createVariable('eta','f8',('z'))
    data_eta.unit = 'm'
    data_eta.long_name = '\u03b7 coordinate of volume centroids: zero is at soil surface and and positive upward'
    
    data_eta_dual = ncfile.createVariable('etaDual','f8',('z_dual'))
    data_eta_dual.unit = 'm'
    data_eta_dual.long_name = '\u03b7 coordinate of volume interfaces: zero is at soil surface and and positive upward. '
    
    data_z = ncfile.createVariable('z','f8',('z'))
    data_z.unit = 'm'
    data_z.long_name = 'z coordinate  of volume centroids: zero is at the bottom of the column and and positive upward'
    
    data_z_dual = ncfile.createVariable('zDual','f8',('z_dual'))
    data_z_dual.unit = 'm'
    data_z_dual.long_name = 'z coordinate of volume interfaces: zero is at soil surface and and positive upward.'
    
    data_psi_0 = ncfile.createVariable('psi0','f8',('z'))
    data_psi_0.units = 'm'
    data_psi_0.long_name = 'Water suction initial condition'
    
    data_T_0 = ncfile.createVariable('T0','f8',('z'))
    data_T_0.units = 'K'
    data_T_0.long_name = 'Temperature initial condition'
    	
    data_space_delta = ncfile.createVariable('spaceDelta','f8',('space_delta'))
    data_space_delta.unit = 'm'
    data_space_delta.long_name = 'Distance between consecutive controids, is used to compute gradients'
    
    data_control_volume = ncfile.createVariable('controlVolume','f8',('z'))
    data_control_volume.unit = 'm'
    data_control_volume.long_name = 'Control volume size'
    
    data_control_volume_index = ncfile.createVariable('controlVolumeIndex','i4',('control_volume_index'))
    data_control_volume_index.unit = ''
    data_control_volume_index.long_name = 'Index of control volumes for calibration'
    
    data_equation_state_ID = ncfile.createVariable('equationStateID','i4',('z'))
    data_equation_state_ID.units = '-'
    data_equation_state_ID.long_name = 'label describing the equation state of the k-th element'
    
    data_parameter_ID = ncfile.createVariable('parameterID','i4',('z'))
    data_parameter_ID.units = '-'
    data_parameter_ID.long_name = 'label identifying the set of parameters'

    data_soil_particles_density = ncfile.createVariable('soilParticlesDensity','f8',('parameter'))
    data_soil_particles_density.units = 'kg/m3'
    data_soil_particles_density.long_name = 'density of soil particles'
    
    data_thermal_conductivity_soil_particles = ncfile.createVariable('thermalConductivitySoilParticles','f8',('parameter'))
    data_thermal_conductivity_soil_particles.units = 'W/m2'
    data_thermal_conductivity_soil_particles.long_name = 'thermal conductivity of soil particles'
    
    data_specific_heat_capacity_soil_particles = ncfile.createVariable('specificThermalCapacitySoilParticles','f8',('parameter'))
    data_specific_heat_capacity_soil_particles.units = 'J/kg m3'
    data_specific_heat_capacity_soil_particles.long_name = 'specific thermal capacity of soil particles'
   
    data_theta_s = ncfile.createVariable('thetaS','f8',('parameter'))
    data_theta_s.units = '-'
    data_theta_s.long_name = 'adimensional water content at saturation'
    
    data_theta_r = ncfile.createVariable('thetaR','f8',('parameter'))
    data_theta_r.units = '-'
    data_theta_r.long_name = 'residual adimensional water content'
    
    data_melting_temperature = ncfile.createVariable('meltingTemperature','f8',('parameter'))
    data_melting_temperature.units = 'K'
    data_melting_temperature.long_name = 'melting temperature of soil water'

    data_par_1 = ncfile.createVariable('par1SWRC','f8',('parameter'))
    data_par_1.units = '-'
    data_par_1.long_name = 'SWRC parameter'
    
    data_par_2 = ncfile.createVariable('par2SWRC','f8',('parameter'))
    data_par_2.units = '-'
    data_par_2.long_name = 'SWRC parameter'
    
    data_par_3 = ncfile.createVariable('par3SWRC','f8',('parameter'))
    data_par_3.units = '-'
    data_par_3.long_name = 'SWRC parameter'
    
    data_par_4 = ncfile.createVariable('par4SWRC','f8',('parameter'))
    data_par_4.units = '-'
    data_par_4.long_name = 'SWRC parameter'
    
    data_par_5 = ncfile.createVariable('par5SWRC','f8',('parameter'))
    data_par_5.units = '-'
    data_par_5.long_name = 'SWRC parameter'
    
    data_alpha_ss = ncfile.createVariable('alphaSpecificStorage','f8',('parameter'))
    data_alpha_ss.units = '1/Pa'
    data_alpha_ss.long_name = 'acquitard compressibility'
    
    data_beta_ss = ncfile.createVariable('betaSpecificStorage','f8',('parameter'))
    data_beta_ss.units = '1/Pa'
    data_beta_ss.long_name = 'water compressibility'
    
    data_ks = ncfile.createVariable('ks','f8',('parameter'))
    data_ks.units = 'm/s'
    data_ks.long_name = 'saturated hydraulic conductivity'
    
    ## write data to variable.

    data_KMAX[0] = KMAX

    for i in range(0,dim):
        data_eta[i] = eta[i]
        data_z[i] = z[i]
        data_control_volume[i] = control_volume[i]
        data_psi_0[i] = psi_0[i]
        data_T_0[i] = T_0[i]
        data_equation_state_ID[i] = equation_state_ID[i]
        data_parameter_ID[i] = parameter_ID[i]
    
    for i in range(0,dim1):
        data_eta_dual[i] = eta_dual[i]
        data_z_dual[i] = z_dual[i]
        data_space_delta[i] = space_delta[i]
        
    for i in range(0,dim2):
        data_space_delta[i] = space_delta[i]
        
    for i in range(0,dim_parameter):
        data_soil_particles_density[i] = soil_particles_density[i]
        data_thermal_conductivity_soil_particles[i] = thermal_conductivity_soil_particles[i]
        data_specific_heat_capacity_soil_particles[i] = specific_heat_capacity_soil_particles[i]
        data_theta_s[i] = theta_s[i]
        data_theta_r[i] = theta_r[i]
        data_melting_temperature[i] = -9999.0
        data_par_1[i] = par_1[i]
        data_par_2[i] = par_2[i]
        data_par_3[i] = par_3[i]
        data_par_4[i] = par_4[i]
        data_par_5[i] = par_5[i]
        data_alpha_ss[i] = alpha_ss[i]
        data_beta_ss[i] = beta_ss[i]
        data_ks[i] = ks[i]
        
#     if not control_volume_index:
    if control_volume_index is None:
        data_control_volume_index = -9999
    else:
        for i in range(0,len(control_volume_index)):
            data_control_volume_index[i] = control_volume_index[i]
    
    ## close the file.
    ncfile.close()
    print ('\n\n***SUCCESS writing!  '+ output_file_name)

    
    return