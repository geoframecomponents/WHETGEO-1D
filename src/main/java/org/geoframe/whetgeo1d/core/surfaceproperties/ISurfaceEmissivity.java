/*
 * GNU GPL v3 License
 *
 * Copyright 2020  Niccolo` Tubini
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

import org.geoframe.whetgeo1d.core.surfaceproperties.emissivity.SurfaceEmissivityConstant;
import org.geoframe.whetgeo1d.core.surfaceproperties.emissivity.SurfaceEmissivityVanBavelHillel;
import org.geoframe.whetgeo1d.core.utils.EnumUtils;

/**
 * Interface for the evaluation of the surface emissivity.
 * 
 * @author Niccolo' Tubini
 *
 */
public interface ISurfaceEmissivity {
	enum SurfaceEmissivityModel {
		CONSTANT, VAN_BAVEL_HILLEL;
		
		public static SurfaceEmissivityModel fromString(String value) {
	        return EnumUtils.fromString(SurfaceEmissivityModel.class, value);
	    }
	}

	/**
	 * 
	 * TODO this interface should probably be remove as the two implmentations make
	 * use, one of the first parameter ignoring the second... and the other one
	 * ignoring the first.
	 */
	double evaluate(double waterContent, double par1);

	/**
	 * Factory method to create a surface emissivity object based on the selected model.
	 * 
	 * @param model the surface emissivity model to be used for the evaluation of the surface emissivity.
	 * @return an instance of ISurfaceEmissivity based on the selected model.
	 */
	public static ISurfaceEmissivity create(SurfaceEmissivityModel model) {
		ISurfaceEmissivity surfaceEmissivity = null;
		if (model == SurfaceEmissivityModel.CONSTANT) {
			surfaceEmissivity = new SurfaceEmissivityConstant();
		} else if (model == SurfaceEmissivityModel.VAN_BAVEL_HILLEL) {
			surfaceEmissivity = new SurfaceEmissivityVanBavelHillel();
		} else {
			throw new IllegalArgumentException("Please check the surface emissivity model.");
		}

		return surfaceEmissivity;
	}

}
