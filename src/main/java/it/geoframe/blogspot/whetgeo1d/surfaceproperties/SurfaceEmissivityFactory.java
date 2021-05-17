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
 * A simple design factory to create a SurfaceEmissivity objects.
 */

public class SurfaceEmissivityFactory {
	
	
	public SurfaceEmissivity create (String type) {

		SurfaceEmissivity surfaceEmissivity = null;
		if(type.equalsIgnoreCase("Constant") ){
			surfaceEmissivity = new SurfaceEmissivityConstant();
		} else if(type.equalsIgnoreCase("van Bavel Hillel") || type.equalsIgnoreCase("vanBavelHillel")){
			surfaceEmissivity = new SurfaceEmissivityVanBavelHillel();
		} else {
			System.out.println("\n\n\tERROR: please check the surface emissivity model.");
		}

		return surfaceEmissivity;

		}
}
