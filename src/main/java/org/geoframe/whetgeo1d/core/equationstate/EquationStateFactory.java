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

package org.geoframe.whetgeo1d.core.equationstate;

import org.geoframe.closureequation.closureequation.ClosureEquation;
import org.geoframe.closureequation.equationstate.EquationState;
import org.geoframe.whetgeo1d.core.equationstate.models.FreezingSoilInternalEnergy;
import org.geoframe.whetgeo1d.core.equationstate.models.PureWaterHeatCapacity;
import org.geoframe.whetgeo1d.core.equationstate.models.PureWaterInternalEnergy;
import org.geoframe.whetgeo1d.core.equationstate.models.SoilHeatCapacity;
import org.geoframe.whetgeo1d.core.equationstate.models.SoilInternalEnergy;
import org.geoframe.whetgeo1d.core.equationstate.models.SoilWaterVolumeBrooksCorey;
import org.geoframe.whetgeo1d.core.equationstate.models.SoilWaterVolumeGardner;
import org.geoframe.whetgeo1d.core.equationstate.models.SoilWaterVolumeKosugi;
import org.geoframe.whetgeo1d.core.equationstate.models.SoilWaterVolumeRomano;
import org.geoframe.whetgeo1d.core.equationstate.models.SoilWaterVolumeVanGenuchten;
import org.geoframe.whetgeo1d.core.equationstate.models.WaterDepth;
import org.geoframe.whetgeo1d.utils.GFGeometry;
import org.geoframe.whetgeo1d.utils.ProblemQuantities;

public class EquationStateFactory {
	public enum StateEquationModel {
		VAN_GENUCHTEN, BROOKS_COREY, KOSUGI, ROMANO, GARDNER, WATER_DEPTH, SOIL_INTERNAL_ENERGY,
		FREEZING_SOIL_INTERNAL_ENERGY, SOIL_HEAT_CAPACITY, WATER_HEAT_CAPACITY, WATER_INTERNAL_ENERGY
	}

	public EquationState create(StateEquationModel model, ClosureEquation closureEquation, GFGeometry geometry,
			ProblemQuantities problemQuantities) {

		switch (model) {
		case VAN_GENUCHTEN:
			return new SoilWaterVolumeVanGenuchten(closureEquation, geometry, problemQuantities);

		case BROOKS_COREY:
			return new SoilWaterVolumeBrooksCorey(closureEquation, geometry, problemQuantities);

		case KOSUGI:
			return new SoilWaterVolumeKosugi(closureEquation, geometry, problemQuantities);

		case ROMANO:
			return new SoilWaterVolumeRomano(closureEquation, geometry, problemQuantities);

		case GARDNER:
			return new SoilWaterVolumeGardner(closureEquation, geometry, problemQuantities);

		case WATER_DEPTH:
			return new WaterDepth(closureEquation, problemQuantities);

		case SOIL_INTERNAL_ENERGY:
			return new SoilInternalEnergy(closureEquation, geometry, problemQuantities);

		case FREEZING_SOIL_INTERNAL_ENERGY:
			return new FreezingSoilInternalEnergy(closureEquation, geometry, problemQuantities);

		case SOIL_HEAT_CAPACITY:
			return new SoilHeatCapacity(closureEquation, geometry, problemQuantities);

		case WATER_HEAT_CAPACITY:
			return new PureWaterHeatCapacity(closureEquation, geometry, problemQuantities);

		case WATER_INTERNAL_ENERGY:
			return new PureWaterInternalEnergy(closureEquation, geometry, problemQuantities);

		default:
			throw new IllegalArgumentException("Unsupported StateEquationModel: " + model);
		}
	}

}
