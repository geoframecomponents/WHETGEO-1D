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
public class SurfaceWaterVaporResistanceFeddes extends SurfaceWaterVaporResistance{

	
	public double evaluate(double h1, double h2, double h3, double h4, double waterSuction) {
		
		if(waterSuction>h1 || waterSuction<h4) {
			return 0;
		} else if (waterSuction>=h4 && waterSuction<h3){
			return (waterSuction-h4)/(h3-h4);
		} else if (waterSuction>h2 && waterSuction<=h1) {
			return (waterSuction-h2)/(h1-h2);
		} else {
			return 1.0;
		}
	}

	
}
