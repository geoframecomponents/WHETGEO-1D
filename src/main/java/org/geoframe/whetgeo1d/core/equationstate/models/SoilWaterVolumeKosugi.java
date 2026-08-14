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
package org.geoframe.whetgeo1d.core.equationstate.models;

import org.geoframe.closureequation.closureequation.ClosureEquation;
import org.geoframe.closureequation.closureequation.Parameters;
import org.geoframe.closureequation.equationstate.EquationState;
import org.geoframe.whetgeo1d.utils.GFGeometry;
import org.geoframe.whetgeo1d.utils.ProblemQuantities;

/**
 * @author Niccolo` Tubini
 *
 */
public class SoilWaterVolumeKosugi extends EquationState {

	private GFGeometry geometry;
	private ProblemQuantities variables;

	public SoilWaterVolumeKosugi(ClosureEquation closureEquation, GFGeometry geometry, ProblemQuantities variables) {
		super(closureEquation);
		this.geometry = geometry;
		this.variables = variables;
	}

	@Override
	public double equationState(double x, double y, int id, int element) {

		return super.closureEquation.f(x, y, id) * geometry.controlVolume[element];

	}

	@Override
	public double dEquationState(double x, double y, int id, int element) {

		return super.closureEquation.df(x, y, id) * geometry.controlVolume[element];

	}

	@Override
	public double ddEquationState(double x, double y, int id, int element) {

		return super.closureEquation.ddf(x, y, id) * geometry.controlVolume[element];

	}

	@Override
	public double p(double x, double y, int id, int element) {

		if (x <= variables.waterSuctionStar1[element]) {
			return dEquationState(x, y, id, element);
		} else if (x > variables.waterSuctionStar1[element] && x <= 0) {
			return dEquationState(variables.waterSuctionStar1[element], y, id, element);
		} else {
			return dEquationState(variables.waterSuctionStar1[element], y, id, element)
					+ dEquationState(x, y, id, element);
		}

	}

	@Override
	public double pIntegral(double x, double y, int id, int element) {

		if (x <= variables.waterSuctionStar1[element]) {
			return equationState(x, y, id, element);
		} else if (x > variables.waterSuctionStar1[element] && x <= 0) {
			return equationState(variables.waterSuctionStar1[element], y, id, element)
					+ dEquationState(variables.waterSuctionStar1[element], y, id, element)
							* (x - variables.waterSuctionStar1[element]);
		} else {
			return equationState(variables.waterSuctionStar1[element], y, id, element)
					+ dEquationState(variables.waterSuctionStar1[element], y, id, element)
							* (x - variables.waterSuctionStar1[element])
					+ dEquationState(x, y, id, element) * (x - 0);
		}

	}

	@Override
	public void computeXStar(double y, int id, int element) {

		Parameters parameters = closureEquation.getParameters();
		variables.waterSuctionStar1[element] = parameters.par1[id]
				/ Math.exp(Math.pow(parameters.par2[id], 2));
		variables.waterSuctionStar2[element] = -9999.0;
		variables.waterSuctionStar3[element] = -9999.0;

	}

	@Override
	public double initialGuess(double x, int id, int element) {

//		return Math.min(variables.waterSuctions[element], variables.waterSuctionStar1[element]);
		return Math.min(x, variables.waterSuctionStar1[element]);

	}

}
