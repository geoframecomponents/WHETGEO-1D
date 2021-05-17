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
 * van Bavel, C.H.M., and D.I. Hillel. 1976. Calculating potential andactual evaporation from a bare soil surface by simulation 
 * of concurrent flow of water and heat. Agric. For. Meteorol. 17:453–476.
 * 
 * @author Niccolo' Tubini
 *
 */
public class SurfaceEmissivityVanBavelHillel extends SurfaceEmissivity {

	
	public double evaluate(double waterContent, double par1) {
		return Math.min(0.9 + 0.18*waterContent, 1.0);
	}
	
}
