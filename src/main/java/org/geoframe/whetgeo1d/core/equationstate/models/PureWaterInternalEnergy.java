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
public class PureWaterInternalEnergy extends EquationState {

	private GFGeometry geometry;
	private ProblemQuantities variables;
	private double epsilon;

	public PureWaterInternalEnergy(ClosureEquation closureEquation, GFGeometry geometry, ProblemQuantities variables) {
		super(closureEquation);
		this.geometry = geometry;
		this.variables = variables;
	}

	@Override
	public double equationState(double x, double y, int id, int element) {

		Parameters parameters = closureEquation.getParameters();
		epsilon = parameters.referenceTemperatureInternalEnergy
				- parameters.meltingTemperature[id];

		if (x <= parameters.meltingTemperature[id]) {
			return (parameters.iceDensity
					* parameters.specificThermalCapacityIce
					* (x - parameters.referenceTemperatureInternalEnergy))
					* geometry.controlVolume[element];
		} else if (x > 273.15) {
			return (parameters.waterDensity
					* parameters.specificThermalCapacityWater
					* (x - parameters.referenceTemperatureInternalEnergy)
					+ parameters.waterDensity * parameters.latentHeatFusion)
					* geometry.controlVolume[element];
		} else {
			return parameters.iceDensity
					* parameters.specificThermalCapacityIce
					* (parameters.meltingTemperature[id]
							- parameters.referenceTemperatureInternalEnergy)
					* geometry.controlVolume[element]
					+ (parameters.waterDensity * parameters.latentHeatFusion
							* geometry.controlVolume[element]
							- parameters.iceDensity
									* parameters.specificThermalCapacityIce
									* (parameters.meltingTemperature[id]
											- parameters.referenceTemperatureInternalEnergy)
									* geometry.controlVolume[element])
							/ epsilon * (x - (parameters.meltingTemperature[id]));

		}

	}

	@Override
	public double dEquationState(double x, double y, int id, int element) {

		Parameters parameters = closureEquation.getParameters();
		epsilon = parameters.referenceTemperatureInternalEnergy
				- parameters.meltingTemperature[id];

		if (x <= parameters.meltingTemperature[id]) {
			return parameters.iceDensity
					* parameters.specificThermalCapacityIce * geometry.controlVolume[element];
		} else if (x > 273.15) {
			return parameters.waterDensity
					* parameters.specificThermalCapacityWater * geometry.controlVolume[element];
		} else {
			return (parameters.waterDensity * parameters.latentHeatFusion
					* geometry.controlVolume[element]
					- parameters.iceDensity
							* parameters.specificThermalCapacityIce
							* (parameters.meltingTemperature[id]
									- parameters.referenceTemperatureInternalEnergy)
							* geometry.controlVolume[element])
					/ epsilon;

		}

	}

	@Override
	public double ddEquationState(double x, double y, int id, int element) {

		return 0.0;

	}

	@Override
	public double p(double x, double y, int id, int element) {

		Parameters parameters = closureEquation.getParameters();
		epsilon = parameters.referenceTemperatureInternalEnergy
				- parameters.meltingTemperature[id];

		if (x <= variables.temperatureStar1[element]) {
			return dEquationState(x, y, id, element);
		} else {
			return (parameters.waterDensity * parameters.latentHeatFusion
					* geometry.controlVolume[element]
					- parameters.iceDensity
							* parameters.specificThermalCapacityIce
							* (parameters.meltingTemperature[id]
									- parameters.referenceTemperatureInternalEnergy)
							* geometry.controlVolume[element])
					/ epsilon;
		}

	}

	@Override
	public double pIntegral(double x, double y, int id, int element) {

		Parameters parameters = closureEquation.getParameters();
		epsilon = parameters.referenceTemperatureInternalEnergy
				- parameters.meltingTemperature[id];

		if (x <= variables.temperatureStar1[element]) {
			return equationState(x, y, id, element);
		} else {
			return parameters.iceDensity
					* parameters.specificThermalCapacityIce
					* (parameters.meltingTemperature[id]
							- parameters.referenceTemperatureInternalEnergy)
					* geometry.controlVolume[element]
					+ (parameters.waterDensity * parameters.latentHeatFusion
							* geometry.controlVolume[element]
							- parameters.iceDensity
									* parameters.specificThermalCapacityIce
									* (parameters.meltingTemperature[id]
											- parameters.referenceTemperatureInternalEnergy)
									* geometry.controlVolume[element])
							/ epsilon * (x - parameters.meltingTemperature[id]);
		}

	}

	@Override
	public void computeXStar(double y, int id, int element) {

		variables.temperatureStar1[element] = closureEquation.getParameters().meltingTemperature[id];
		variables.temperatureStar2[element] = -9999.0;
		variables.temperatureStar3[element] = -9999.0;

	}

	@Override
	public double initialGuess(double x, int id, int element) {

		return Math.min(variables.temperatures[element], variables.temperatureStar1[element]);

	}

}
