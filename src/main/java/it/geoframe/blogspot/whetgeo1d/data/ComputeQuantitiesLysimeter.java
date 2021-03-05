/*
 * GNU GPL v3 License
 *
 * Copyright 2019 Concetta D'Amato
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package it.geoframe.blogspot.whetgeo1d.data;

import oms3.annotations.Author;
import oms3.annotations.Bibliography;
import oms3.annotations.Description;
import oms3.annotations.Documentation;
import oms3.annotations.Keywords;
import oms3.annotations.License;

@Description("This class compute the quantities necessary to run the lysimeter.")
@Documentation("")
@Author(name = "Concetta D'Amato, Niccolo' Tubini and Riccardo Rigon", contact = "")
@Keywords("Richards equation, numerical solver, finite volume ")
@Bibliography("")
//@Label()
//@Name()
//@Status()
@License("General Public License Version 3 (GPLv3)")
public class ComputeQuantitiesLysimeter {
	
	
	private ProblemQuantities variables;
	private Geometry geometry;
	
	@Description("Values of the adimensional water content at the wilting point")
	private double[] thetaWP;
	
	@Description("Values of the adimensional water content at the field capacity")
	private double[] thetaFC;

	
	public ComputeQuantitiesLysimeter(double[] thetaWP, double[] thetaFC) {
		
		variables = ProblemQuantities.getInstance();
		geometry = Geometry.getInstance();
		
		this.thetaWP = thetaWP;
		this.thetaFC = thetaFC;
	
	}
	
	/*
	 * FIXME: check the extreme of the for loops in order to deal with Neumann BC and 
	 * Dirichlet BC.
	 */
	public void computeEvapoTranspirations(int KMAX, double tTimeStep, double timeDelta, double[] stressedETs) {
			
		for(int element = 0; element < KMAX-2; element++) {
			variables.ETs[element] = stressedETs[element]/1000/tTimeStep*timeDelta;
		}
		
	}
	
	/*
	 * FIXME: check the extreme of the for loops in order to deal with Neumann BC and 
	 * Dirichlet BC.
	 */
	public void checkEvapoTranspirations(int KMAX) {
		
		variables.sumETs = 0.0;
		
		for(int element = 0; element < KMAX-1; element++) {
			
			if(variables.volumes[element] > thetaWP[variables.parameterID[element]]*geometry.controlVolume[element]) {
				if(variables.ETs[element] > (variables.volumes[element] - thetaWP[variables.parameterID[element]]*geometry.controlVolume[element])) {
					variables.ETs[element] = variables.volumes[element] - thetaWP[variables.parameterID[element]]*geometry.controlVolume[element];
					/*
					 * FIXME: add a warning if this happens?
					 */
					System.out.println("\tETs volume larger than water volume available for evapotranspiration, ETs reduced.");
				}
			} else {
				if(variables.ETs[element] != 0.0) {
					System.out.println("\tBroker error, g!=0");
					
					variables.ETs[element] = 0.0;
					System.out.println("\tETs set to 0");
				} else {
					System.out.println("\tETs is already 0");
				}
			}
			
			variables.sumETs = variables.sumETs + variables.ETs[element];
		}
		
		
	}
	
}
