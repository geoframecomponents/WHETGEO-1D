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

import it.geoframe.blogspot.closureequation.ClosureEquation;
import it.geoframe.blogspot.equationstate.EquationState;
import it.geoframe.blogspot.whetgeo1d.data.ProblemQuantities;

/**
 * @author Niccolo` Tubini
 *
 */
public class WaterDepth extends EquationState {

	protected ProblemQuantities variables;

	public WaterDepth(ClosureEquation closureEquation) {
		super(closureEquation);
		this.variables = ProblemQuantities.getInstance();
	}



	@Override
	public double stateEquation(double x, double y, int id, int element) {
		
		if(x>0) {
			return x;
		} else {
			return 0.0;
		}
		
	}


	@Override
	public double dStateEquation(double x, double y, int id, int element) {

		if(x>0) {
			return 1.0;
		} else {
			return 0.0;
		}
		
	}


	@Override
	public double ddStateEquation(double x, double y,  int id, int element) {

		return 0.0;
	}


	@Override
	public double p(double x, double y, int id, int element) {

		return dStateEquation(x, y, id, element);

	}


	@Override
	public double pIntegral(double x, double y, int id, int element) {

		return stateEquation(x, y, id, element);

	}



	@Override
	public void computeXStar(double y, int id, int element) {
		
		variables.xStar1[element] = 0.0;
		variables.xStar2[element] = -9999.0;
		variables.xStar3[element] = -9999.0;
	}

	@Override
	public double initialGuess(double x, int id, int element) {
		
		return Math.max(variables.waterSuctions[element], 0.1);
		
	}


}