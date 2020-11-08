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


/**
 * @author Niccolo` Tubini
 *
 */
public class SoilWaterVolumeBrooksCorey extends EquationState {

	private Geometry geometry;
	private ProblemQuantities variables;

	public SoilWaterVolumeBrooksCorey(ClosureEquation closureEquation) {
		super(closureEquation);
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
		} else {
			return dEquationState(variables.xStar1[element], y, id, element);
		}

	}


	@Override
	public double pIntegral(double x, double y, int id, int element) {

		if(x<=variables.xStar1[element]) {
			return equationState(x, y, id, element);  
		} else {
			return equationState(variables.xStar1[element], y, id, element) + dEquationState(variables.xStar1[element], y, id, element)*(x-variables.xStar1[element]);
		}

	}



	@Override
	public void computeXStar(double y, int id, int element) {
		
		variables.xStar1[element] = super.closureEquation.parameters.par2[id];
		variables.xStar2[element] = -9999.0;
		variables.xStar3[element] = -9999.0;
		
	}
	
	
	@Override
	public double initialGuess(double x, int id, int element) {
		
		return Math.min(variables.waterSuctions[element], variables.xStar1[element]);
		
	}



}
