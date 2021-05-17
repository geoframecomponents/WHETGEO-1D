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
package it.geoframe.blogspot.whetgeo1d.surfaceproperties;

/**
 * A simple design factory to create a SurfaceAlbedo objects.
 */

public class SurfaceAlbedoFactory {
	
	
	public SurfaceAlbedo create (String type) {

		SurfaceAlbedo surfaceAlbedo = null;
		if(type.equalsIgnoreCase("Constant") ){
			surfaceAlbedo = new SurfaceAlbedoConstant();
		} else if(type.equalsIgnoreCase("van Bavel Hillel") || type.equalsIgnoreCase("vanBavelHillel")){
			surfaceAlbedo = new SurfaceAlbedoVanBavelHillel();
		} else {
			System.out.println("\n\n\tERROR: please check the surface albedo model.");
		}

		return surfaceAlbedo;

		}
}
