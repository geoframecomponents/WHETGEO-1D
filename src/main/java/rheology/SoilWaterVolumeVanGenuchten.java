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
package rheology;

import data.Geometry;
import data.ProblemQuantities;
import rheology.Rheology;
import stateequation.*;

/**
 * @author Niccolo` Tubini
 *
 */
public class SoilWaterVolumeVanGenuchten extends StateEquation {

	private Geometry geometry;
	private ProblemQuantities variables;


	public SoilWaterVolumeVanGenuchten(Rheology rheology) {
		super(rheology);
		this.geometry = Geometry.getInstance();
		this.variables = ProblemQuantities.getInstance();
	}



	@Override
	public double stateEquation(double x, double y, int id, int element) {
		
		return rheology.f(x, y, id)*geometry.controlVolume[element];

	}


	@Override
	public double dStateEquation(double x, double y, int id, int element) {

		return rheology.df(x, y, id)*geometry.controlVolume[element];

	}


	@Override
	public double ddStateEquation(double x, double y, int id, int element) {

		return rheology.ddf(x, y, id)*geometry.controlVolume[element];

	}


	@Override
	public double p(double x, double y, int id, int element) {

		if(x<=variables.xStar1[id]) {
			return dStateEquation(x, y, id, element);  
		} else {
			return dStateEquation(variables.xStar1[id], y, id, element);
		}

	}


	@Override
	public double pIntegral(double x, double y, int id, int element) {

		if(x<=variables.xStar1[id]) {
			return stateEquation(x, y, id, element);  
		} else {
			return stateEquation(variables.xStar1[id], y, id, element) + dStateEquation(variables.xStar1[id], y, id, element)*(x-variables.xStar1[id]);
		}

	}



	@Override
	public void computeXStar(double y, int id, int element) {
		
		variables.xStar1[element] = -1/super.rheology.parameters.par2[id]*Math.pow((super.rheology.parameters.par1[id]-1)/(super.rheology.parameters.par1[id]), 1/super.rheology.parameters.par1[id]);
		variables.xStar2[element] = -9999.0;
		variables.xStar3[element] = -9999.0;
		
	}

	@Override
	public double initialGuess(double x, int id, int element) {
		
		return Math.min(variables.waterSuctions[element], variables.xStar1[element]);
		
	}


}
