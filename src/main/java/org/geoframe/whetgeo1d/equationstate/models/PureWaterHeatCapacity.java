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
package org.geoframe.whetgeo1d.equationstate.models;

import org.geoframe.closureequation.closureequation.ClosureEquation;
import org.geoframe.closureequation.closureequation.Parameters;
import org.geoframe.closureequation.equationstate.EquationState;
import org.geoframe.whetgeo1d.utils.GFGeometry;
import org.geoframe.whetgeo1d.utils.ProblemQuantities;

/**
 * @author Niccolo` Tubini
 *
 */
public class PureWaterHeatCapacity extends EquationState {

	private ProblemQuantities variables;

	public PureWaterHeatCapacity(ClosureEquation closureEquation, GFGeometry geometry, ProblemQuantities variables) {
		super(closureEquation);
		this.variables = variables;
	}

	@Override
	public double equationState(double x, double y, int id, int element) {

		Parameters parameters = super.closureEquation.getParameters();
		return parameters.waterDensity
				* parameters.specificThermalCapacityWater * closureEquation.f(y, id);

	}

	@Override
	public double dEquationState(double x, double y, int id, int element) {

		return -9999.0;

	}

	@Override
	public double ddEquationState(double x, double y, int id, int element) {

		return -9999.0;

	}

	@Override
	public double p(double x, double y, int id, int element) {

		return -9999.0;

	}

	@Override
	public double pIntegral(double x, double y, int id, int element) {

		return -9999.0;

	}

	@Override
	public void computeXStar(double y, int id, int element) {

		variables.temperatureStar1[element] = -9999.0;
		variables.temperatureStar2[element] = -9999.0;
		variables.temperatureStar3[element] = -9999.0;

	}

	@Override
	public double initialGuess(double x, int id, int element) {

		return -9999.0;

	}

}
