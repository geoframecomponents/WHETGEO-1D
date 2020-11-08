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
import it.geoframe.blogspot.closureequation.closureequation.ClosureEquation;
import it.geoframe.blogspot.closureequation.equationstate.EquationState;
import it.geoframe.blogspot.numerical.rootfinding.Bisection;

/**
 * @author Niccolo` Tubini
 *
 */
public class SoilWaterVolumeRomano extends EquationState {

	private Geometry geometry;
	private ProblemQuantities variables;
	private Bisection bisection;

	public SoilWaterVolumeRomano(ClosureEquation closureEquation) {
		super(closureEquation);
		this.bisection = new Bisection(this);
		this.geometry = Geometry.getInstance();
		this.variables = ProblemQuantities.getInstance();
	}



	@Override
	public double equationState(double x, double y, int id, int element) {
		
		return super.closureEquation.f(x, y, id)*geometry.controlVolume[element];

	}


	@Override
	public double dEquationState(double x, double y, int id, int element) {

		return super.closureEquation.df(x, y, id)*geometry.controlVolume[element];

	}


	@Override
	public double ddEquationState(double x, double y, int id, int element) {

		return super.closureEquation.ddf(x, y, id)*geometry.controlVolume[element];

	}


	@Override
	public double p(double x, double y, int id, int element) {

		if(x<=variables.xStar1[element]) {
			
			return dEquationState(x, y, id, element);  
			
		} else if(variables.xStar1[element]<x && x<variables.xStar3[element]) {
			
			return this.dEquationState(variables.xStar1[element], y, id, element);
			
		} else if (variables.xStar3[element]<=x && x<=variables.xStar2[element]) {
			
			return this.dEquationState(x, y, id, element) + (this.dEquationState(variables.xStar1[element], y, id, element)-this.dEquationState(variables.xStar3[element], y, id, element));
			
		} else if (variables.xStar2[element]<x && x<0){
			
			return this.dEquationState(variables.xStar2[element], y, id, element) + (this.dEquationState(variables.xStar1[element], y, id, element)-this.dEquationState(variables.xStar3[element], y, id, element));
			
		} else {
			
			return this.dEquationState(variables.xStar2[element], y, id, element) + (this.dEquationState(variables.xStar1[element], y, id, element)-this.dEquationState(variables.xStar3[element], y, id, element))
					+ this.dEquationState(x, y, id, element);
			
		}

	}
	

	@Override
	public double pIntegral(double x, double y, int id, int element) {

		if(x <=variables.xStar1[element]) {
			
			return this.equationState(x, y, id, element); 
			
		} else if(variables.xStar1[element]<x && x<=variables.xStar3[element]) {
			
			return this.equationState(variables.xStar1[element], y, id, element) 
					+ this.dEquationState(variables.xStar1[element], y, id, element)*(x - variables.xStar1[element]);
			
		} else if (variables.xStar3[element]<x && x<=variables.xStar2[element]) {
			
			return this.equationState(variables.xStar1[element], y, id, element) 
					+ this.dEquationState(variables.xStar1[element], y, id, element)*(x-variables.xStar1[element])
					- this.dEquationState(variables.xStar3[element], y, id, element)*(x-variables.xStar3[element])
					+ this.equationState(x, y, id, element) - this.equationState(variables.xStar3[element], y, id, element);
			
		} else if (variables.xStar2[element]<x && x<0){
			
			return this.equationState(variables.xStar1[element], y, id, element)
					+ this.dEquationState(variables.xStar1[element], y, id, element)*(variables.xStar2[element]-variables.xStar1[element])
					- this.dEquationState(variables.xStar3[element], y, id, element)*(variables.xStar2[element]-variables.xStar3[element])
					+ this.equationState(variables.xStar2[element], y, id, element) - this.equationState(variables.xStar3[element], y, id, element)
					+ this.dEquationState(variables.xStar2[element], y, id, element)*(x-variables.xStar2[element]);
			
		} else {
			
			return this.equationState(variables.xStar1[element], y, id, element) 
					+ this.dEquationState(variables.xStar1[element], y, id, element)*(variables.xStar2[element]-variables.xStar1[element]) 
					- this.dEquationState(variables.xStar3[element], y, id, element)*(variables.xStar2[element]-variables.xStar3[element])
					+ this.equationState(variables.xStar2[element], y, id, element) - this.equationState(variables.xStar3[element], y, id, element) 
					+ this.dEquationState(variables.xStar2[element], y, id, element)*(0-variables.xStar2[element])
					+ this.dEquationState(0, y, id, element)*(x-0);

		}

	}
    


	@Override
	public void computeXStar(double y, int id, int element) {
		
		double x1 = super.closureEquation.parameters.par4[id]*Math.exp(-Math.pow(super.closureEquation.parameters.par2[id],2));
		variables.xStar1[element] = bisection.findZero(x1*1.1, x1*0.9, y, id, element);
		double x2 = super.closureEquation.parameters.par5[id]*Math.exp(-Math.pow(super.closureEquation.parameters.par3[id],2));
		variables.xStar2[element] = bisection.findZero(x2*1.2, x2*0.8, y, id, element);
		variables.xStar3[element] = bisection.findZero(variables.xStar1[element]*0.9, variables.xStar2[element]*1.1, y, id, element);
	}
	
	
	@Override
	public double initialGuess(double x, int id, int element) {
		
		return Math.min(variables.waterSuctions[element], variables.xStar1[element]);
		
	}



}
