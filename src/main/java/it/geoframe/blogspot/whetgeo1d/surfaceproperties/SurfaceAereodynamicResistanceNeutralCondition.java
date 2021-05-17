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
package it.geoframe.blogspot.whetgeo1d.surfaceproperties;

/**
 * Evaluate the aerodynamic resistance for neutral condition
 * 
 * Liu, S., Lu, L., Mao, D., & Jia, L. (2007). Evaluating parameterizations of aerodynamic resistance to heat transfer 
 * using field measurements. Hydrology and earth system sciences, 11(2), 769-783.
 * 
 * @author Niccolo' Tubini
 *
 */
public class SurfaceAereodynamicResistanceNeutralCondition extends SurfaceAereodynamicResistance{

	
	public double evaluate(double referenceHeight, double surfaceRoughness, double windVelocity, double zeroPlaneDisplacement) {
		
		if(windVelocity==0) {
			return -9999;
		} else {
			return 1/(Math.pow(constantVonKarman,2)*windVelocity) * Math.pow(Math.log((referenceHeight-zeroPlaneDisplacement)/surfaceRoughness), 2);
		}
	}

	
}
