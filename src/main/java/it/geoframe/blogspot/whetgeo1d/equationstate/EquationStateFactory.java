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

package it.geoframe.blogspot.whetgeo1d.equationstate;

import it.geoframe.blogspot.closureequation.closureequation.ClosureEquation;
import it.geoframe.blogspot.closureequation.equationstate.EquationState;

public class EquationStateFactory {
	
	
	public EquationState create(String model, ClosureEquation closureEquation) {
		
		EquationState myModel = null;
		
		if(model.equalsIgnoreCase("Van Genuchten") || model.equalsIgnoreCase("VanGenuchten") || model.equalsIgnoreCase("VG")) {
			myModel = new SoilWaterVolumeVanGenuchten(closureEquation);
		} else if(model.equalsIgnoreCase("Brooks Corey") || model.equalsIgnoreCase("BrooksCorey") || model.equalsIgnoreCase("BC")) {
			myModel = new SoilWaterVolumeBrooksCorey(closureEquation);
		} else if(model.equalsIgnoreCase("Kosugi") ) {
			myModel = new SoilWaterVolumeKosugi(closureEquation);
		} else if(model.equalsIgnoreCase("Romano") ) {
			myModel = new SoilWaterVolumeRomano(closureEquation);
		} else if(model.equalsIgnoreCase("Gardner") ) {
			myModel = new SoilWaterVolumeGardner(closureEquation);
		} else if(model.equalsIgnoreCase("Water Depth") || model.equalsIgnoreCase("WaterDepth") ) {
			myModel = new WaterDepth(closureEquation);
		} else if(model.equalsIgnoreCase("Soil Internal Energy") || model.equalsIgnoreCase("SoilInternalEnergy")) {
			myModel = new SoilInternalEnergy(closureEquation);
		} else if(model.equalsIgnoreCase("Freezing Soil Internal Energy") || model.equalsIgnoreCase("FreezingSoilInternalEnergy")) {
			myModel = new FreezingSoilInternalEnergy(closureEquation);
		} else if(model.equalsIgnoreCase("Soil heat capacity") || model.equalsIgnoreCase("SoilHeatCapacity")) {
			myModel = new SoilHeatCapacity(closureEquation);
		} else if(model.equalsIgnoreCase("Water heat capacity") || model.equalsIgnoreCase("WaterHeatCapacity")) {
			myModel = new PureWaterHeatCapacity(closureEquation);
		} else {
			System.out.println("\n\n\tERROR: please check the stateEquationModel");
		}
		return myModel;
	}

}
