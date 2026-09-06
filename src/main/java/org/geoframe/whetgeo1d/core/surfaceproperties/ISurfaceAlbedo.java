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

import org.geoframe.whetgeo1d.core.surfaceproperties.albedo.SurfaceAlbedoConstant;
import org.geoframe.whetgeo1d.core.surfaceproperties.albedo.SurfaceAlbedoVanBavelHillel;
import org.geoframe.whetgeo1d.core.utils.EnumUtils;

/**
 * Interface for the evaluation of the surface albedo.
 * 
 * @author Niccolo' Tubini
 *
 */
public interface ISurfaceAlbedo {
	enum SurfaceAlbedoModel {
		CONSTANT, VAN_BAVEL_HILLEL;

		public static SurfaceAlbedoModel fromString(String value) {
			return EnumUtils.fromString(SurfaceAlbedoModel.class, value);
		}
	}

	double evaluate(double waterContent, double par1);

	public static ISurfaceAlbedo create(SurfaceAlbedoModel type) {

		ISurfaceAlbedo surfaceAlbedo = null;
		if (type == SurfaceAlbedoModel.CONSTANT) {
			surfaceAlbedo = new SurfaceAlbedoConstant();
		} else if (type == SurfaceAlbedoModel.VAN_BAVEL_HILLEL) {
			surfaceAlbedo = new SurfaceAlbedoVanBavelHillel();
		} else {
			throw new IllegalArgumentException("Check the surface albedo model.");
		}

		return surfaceAlbedo;

	}

}
