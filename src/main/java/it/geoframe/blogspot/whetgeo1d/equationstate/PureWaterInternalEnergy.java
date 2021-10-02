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

import it.geoframe.blogspot.closureequation.closureequation.ClosureEquation;
import it.geoframe.blogspot.closureequation.equationstate.EquationState;
import it.geoframe.blogspot.whetgeo1d.data.Geometry;
import it.geoframe.blogspot.whetgeo1d.data.ProblemQuantities;

/**
 * @author Niccolo` Tubini
 *
 */
public class PureWaterInternalEnergy extends EquationState {

	private Geometry geometry;
	private ProblemQuantities variables;
	private double epsilon;

	public PureWaterInternalEnergy(ClosureEquation closureEquation) {
		super(closureEquation);
		this.geometry = Geometry.getInstance();
		this.variables = ProblemQuantities.getInstance();
	}



	@Override
	public double equationState(double x, double y, int id, int element) {
		
		epsilon = super.closureEquation.parameters.referenceTemperatureInternalEnergy-super.closureEquation.parameters.meltingTemperature[id];

		if(x<=super.closureEquation.parameters.meltingTemperature[id]) {
			return  ( super.closureEquation.parameters.iceDensity*super.closureEquation.parameters.specificThermalCapacityIce*(x-super.closureEquation.parameters.referenceTemperatureInternalEnergy) ) * geometry.controlVolume[element];
		} else if(x>273.15) {
			return  ( super.closureEquation.parameters.waterDensity*super.closureEquation.parameters.specificThermalCapacityWater*(x-super.closureEquation.parameters.referenceTemperatureInternalEnergy)
						+ super.closureEquation.parameters.waterDensity*super.closureEquation.parameters.latentHeatFusion ) * geometry.controlVolume[element];
		} else {
			return  super.closureEquation.parameters.iceDensity*super.closureEquation.parameters.specificThermalCapacityIce*(super.closureEquation.parameters.meltingTemperature[id]-super.closureEquation.parameters.referenceTemperatureInternalEnergy)*geometry.controlVolume[element] +
					(super.closureEquation.parameters.waterDensity*super.closureEquation.parameters.latentHeatFusion*geometry.controlVolume[element] - super.closureEquation.parameters.iceDensity*super.closureEquation.parameters.specificThermalCapacityIce*(super.closureEquation.parameters.meltingTemperature[id]-super.closureEquation.parameters.referenceTemperatureInternalEnergy)*geometry.controlVolume[element] )/epsilon*(x-(super.closureEquation.parameters.meltingTemperature[id]));
		
		}
		
	}


	@Override
	public double dEquationState(double x, double y, int id, int element) {

		epsilon = super.closureEquation.parameters.referenceTemperatureInternalEnergy-super.closureEquation.parameters.meltingTemperature[id];

		if(x<=super.closureEquation.parameters.meltingTemperature[id]) {
			return  super.closureEquation.parameters.iceDensity*super.closureEquation.parameters.specificThermalCapacityIce * geometry.controlVolume[element];
		} else if(x>273.15) {
			return  super.closureEquation.parameters.waterDensity*super.closureEquation.parameters.specificThermalCapacityWater * geometry.controlVolume[element];
		} else {
			return ( super.closureEquation.parameters.waterDensity*super.closureEquation.parameters.latentHeatFusion*geometry.controlVolume[element] 
					- 	super.closureEquation.parameters.iceDensity*super.closureEquation.parameters.specificThermalCapacityIce*(super.closureEquation.parameters.meltingTemperature[id]-super.closureEquation.parameters.referenceTemperatureInternalEnergy)*geometry.controlVolume[element]  ) /epsilon;
							
		}

	}


	@Override
	public double ddEquationState(double x, double y, int id, int element) {

		return 0.0;

	}


	@Override
	public double p(double x, double y, int id, int element) {

		epsilon = super.closureEquation.parameters.referenceTemperatureInternalEnergy-super.closureEquation.parameters.meltingTemperature[id];
		
		if(x<=variables.temperatureStar1[element]) {
			return dEquationState(x, y, id, element);
		} else {
			return ( super.closureEquation.parameters.waterDensity*super.closureEquation.parameters.latentHeatFusion*geometry.controlVolume[element] 
					- 	super.closureEquation.parameters.iceDensity*super.closureEquation.parameters.specificThermalCapacityIce*(super.closureEquation.parameters.meltingTemperature[id]-super.closureEquation.parameters.referenceTemperatureInternalEnergy)*geometry.controlVolume[element]  ) /epsilon;
		}

	}


	@Override
	public double pIntegral(double x, double y, int id, int element) {

		epsilon = super.closureEquation.parameters.referenceTemperatureInternalEnergy-super.closureEquation.parameters.meltingTemperature[id];

		if(x<=variables.temperatureStar1[element]) {
			return equationState(x, y, id, element);
		} else {
			return super.closureEquation.parameters.iceDensity*super.closureEquation.parameters.specificThermalCapacityIce*(super.closureEquation.parameters.meltingTemperature[id]-super.closureEquation.parameters.referenceTemperatureInternalEnergy)*geometry.controlVolume[element] +
					(super.closureEquation.parameters.waterDensity*super.closureEquation.parameters.latentHeatFusion*geometry.controlVolume[element] - super.closureEquation.parameters.iceDensity*super.closureEquation.parameters.specificThermalCapacityIce*(super.closureEquation.parameters.meltingTemperature[id]-super.closureEquation.parameters.referenceTemperatureInternalEnergy)*geometry.controlVolume[element] )/epsilon*(x-super.closureEquation.parameters.meltingTemperature[id]);
		}

	}



	@Override
	public void computeXStar(double y, int id, int element) {
		
		variables.temperatureStar1[element] = super.closureEquation.parameters.meltingTemperature[id];
		variables.temperatureStar2[element] = -9999.0;
		variables.temperatureStar3[element] = -9999.0;
		
	}

	@Override
	public double initialGuess(double x, int id, int element) {
		
		return Math.min(variables.temperatures[element], variables.temperatureStar1[element]);
		
	}


}
