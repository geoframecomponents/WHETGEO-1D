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
package org.geoframe.whetgeo1d.surfaceproperties;

import org.geoframe.whetgeo1d.surfaceproperties.vaporresistance.SurfaceWaterVaporResistanceFeddes;
import org.geoframe.whetgeo1d.utils.EnumUtils;

/**
 * Interface to evaluate the surface water vapor resistance.
 * 
 * @author Niccolo' Tubini
 */
public interface ISurfaceWaterVaporResistance {
	enum SurfaceWaterVaporResistanceModel {
		FEDDES;
		
		public static SurfaceWaterVaporResistanceModel fromString(String value) {
	        return EnumUtils.fromString(SurfaceWaterVaporResistanceModel.class, value);
	    }
	}
	
    /**
     * TODO javadoc to be checked.
     * 
     * Evaluates the vapor resistance reduction factor.
     *
     * @param h1 upper (wet) suction threshold
     * @param h2 upper bound of optimal suction range
     * @param h3 lower bound of optimal suction range
     * @param h4 lower (dry) suction threshold
     * @param waterSuction current soil water suction (pressure head)
     * @return a dimensionless reduction factor in the range {@code [0, 1]}
     */
	double evaluate(double h1, double h2, double h3, double h4, double waterSuction);
	
	/**
	 * Factory method to create a surface water vapor resistance object based on the selected model.
	 * 
	 * @param model the surface water vapor resistance model to be used for the evaluation of the surface water vapor resistance.
	 * @return an instance of ISurfaceWaterVaporResistance based on the selected model.
	 */
	public static ISurfaceWaterVaporResistance create(SurfaceWaterVaporResistanceModel model) {
		ISurfaceWaterVaporResistance surfaceWaterVaporResistance = null;
		if (model == SurfaceWaterVaporResistanceModel.FEDDES) {
			surfaceWaterVaporResistance = new SurfaceWaterVaporResistanceFeddes();
		} else {
			throw new IllegalArgumentException("Please check the surface water vapor resistance model.");
		}
		return surfaceWaterVaporResistance;
	}
	
}
