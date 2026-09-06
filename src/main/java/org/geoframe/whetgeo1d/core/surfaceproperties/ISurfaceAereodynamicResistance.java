/*
 * GNU GPL v3 License
 *
 * Copyright 2019  Niccolo` Tubini
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
package org.geoframe.whetgeo1d.core.surfaceproperties;

import org.geoframe.whetgeo1d.core.surfaceproperties.aerodynamicresistance.SurfaceAereodynamicResistanceNeutralCondition;
import org.geoframe.whetgeo1d.core.utils.EnumUtils;

/**
 * Evaluate the aerodynamic resistance
 * 
 * @author Niccolo' Tubini
 *
 */
public interface ISurfaceAereodynamicResistance {
	enum SurfaceAereodynamicResistanceType {
		NeutralCondition;
		
		public static SurfaceAereodynamicResistanceType fromString(String value) {
	        return EnumUtils.fromString(SurfaceAereodynamicResistanceType.class, value);
	    }
	}

	double constantVonKarman = 0.4;

	double evaluate(double referenceHeight, double surfaceRoughness, double windVelocity, double zeroPlaneDisplacement);

	public static ISurfaceAereodynamicResistance create(SurfaceAereodynamicResistanceType type) {

		ISurfaceAereodynamicResistance surfaceAereodynamic = null;
		if (type == SurfaceAereodynamicResistanceType.NeutralCondition) {
			surfaceAereodynamic = new SurfaceAereodynamicResistanceNeutralCondition();
		} else {
			throw new IllegalArgumentException("Check the surface aereodynamic resistance model.");
		}

		return surfaceAereodynamic;

	}
}
