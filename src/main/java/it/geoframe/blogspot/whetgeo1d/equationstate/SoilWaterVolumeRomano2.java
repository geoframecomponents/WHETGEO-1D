/*
 * GNU GPL v3 License
 *
 * Copyright 2019 Niccolo` Tubini
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

/**
 * 
 */
package it.geoframe.blogspot.whetgeo1d.equationstate;

import it.geoframe.blogspot.whetgeo1d.data.Geometry;
import it.geoframe.blogspot.whetgeo1d.data.ProblemQuantities;
import rheology.Rheology;
import stateequation.*;
import rootfinding.Bisection;

/**
 * @author Niccolo` Tubini
 *
 */
public class SoilWaterVolumeRomano2 extends StateEquation {

	private Rheology soilModel;
	protected Geometry geometry;
	protected ProblemQuantities variables;
	private Bisection bisection;

	public SoilWaterVolumeRomano2(Rheology rheology) {
		this.soilModel = rheology;
		this.bisection = new Bisection(this);
		this.geometry = Geometry.getInstance();
		this.variables = ProblemQuantities.getInstance();
	}



	@Override
	public double stateEquation(double x, double y, int id, int element) {
		
		return soilModel.f(x, y, id)*geometry.controlVolume[element];

	}


	@Override
	public double dStateEquation(double x, double y, int id, int element) {

		return soilModel.df(x, y, id)*geometry.controlVolume[element];

	}


	@Override
	public double ddStateEquation(double x, double y, int id, int element) {

		return soilModel.ddf(x, y, id)*geometry.controlVolume[element];

	}


	@Override
	public double p(double x, double y, int id, int element) {

		if(x<=variables.xStar1[element]) {
			
			return dStateEquation(x, y, id, element);  
			
		} else if(variables.xStar1[element]<x && x<variables.xStar3[element]) {
			
			return this.dStateEquation(variables.xStar1[element], y, id, element);
			
		} else if (variables.xStar3[element]<=x && x<=variables.xStar2[element]) {
			
			return this.dStateEquation(variables.xStar1[element], y, id, element) + this.dStateEquation(x, y, id, element) ;
			
		} else if (variables.xStar2[element]<x && x<0){
			
			return this.dStateEquation(variables.xStar1[element], y, id, element) + this.dStateEquation(variables.xStar2[element], y, id, element);
			
		} else {
			
			return this.dStateEquation(variables.xStar1[element], y, id, element) + this.dStateEquation(variables.xStar2[element], y, id, element) + 
					+ this.dStateEquation(x, y, id, element);
			
		}

	}
	

	@Override
	public double pIntegral(double x, double y, int id, int element) {

		if(x <=variables.xStar1[element]) {
			
			return this.stateEquation(x, y, id, element); 
			
		} else if(variables.xStar1[element]<x && x<=variables.xStar3[element]) {
			
			return this.stateEquation(variables.xStar1[element], y, id, element) 
					+ this.dStateEquation(variables.xStar1[element], y, id, element)*(x - variables.xStar1[element]);
			
		} else if (variables.xStar3[element]<x && x<=variables.xStar2[element]) {
			
			return this.stateEquation(variables.xStar1[element], y, id, element) 
					+ this.dStateEquation(variables.xStar1[element], y, id, element)*(x-variables.xStar1[element])
					+ this.stateEquation(x, y, id, element) - this.stateEquation(variables.xStar3[element], y, id, element);
			
		} else if (variables.xStar2[element]<x && x<0){
			
			return this.stateEquation(variables.xStar1[element], y, id, element)
					+ this.dStateEquation(variables.xStar1[element], y, id, element)*(variables.xStar2[element]-variables.xStar1[element])
					+ this.stateEquation(variables.xStar2[element], y, id, element) - this.stateEquation(variables.xStar3[element], y, id, element)
					+ this.dStateEquation(variables.xStar2[element], y, id, element)*(x-variables.xStar2[element]);
			
		} else {
			
			return this.stateEquation(variables.xStar1[element], y, id, element) 
					+ this.dStateEquation(variables.xStar1[element], y, id, element)*(variables.xStar2[element]-variables.xStar1[element]) 
					+ this.stateEquation(variables.xStar2[element], y, id, element) - this.stateEquation(variables.xStar3[element], y, id, element) 
					+ this.dStateEquation(variables.xStar2[element], y, id, element)*(0-variables.xStar2[element])
					+ this.dStateEquation(0, y, id, element)*(x-0);

		}

	}
    


	@Override
	public void computeXStar(double y, int id, int element) {
		
		double x1 = super.parameters.par4[id]*Math.exp(-Math.pow(super.parameters.par2[id],2));
		variables.xStar1[element] = bisection.findZero(x1*1.1, x1*0.9, y, id, element);
		double x2 = super.parameters.par5[id]*Math.exp(-Math.pow(super.parameters.par3[id],2));
		variables.xStar2[element] = bisection.findZero(x2*1.1, x2*0.9, y, id, element);
		variables.xStar3[element] = bisection.findZero(variables.xStar1[element]*0.9, variables.xStar2[element]*1.1, y, id, element);
	}
	
	@Override
	public double initialGuess(double x, int id, int element) {
		
		return Math.min(variables.waterSuctions[element], variables.xStar1[element]);
		
	}



}
