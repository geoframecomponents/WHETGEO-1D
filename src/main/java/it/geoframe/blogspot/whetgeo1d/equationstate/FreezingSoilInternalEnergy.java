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

/**
 * 
 */
package it.geoframe.blogspot.whetgeo1d.equationstate;

import it.geoframe.blogspot.closureequation.closureequation.ClosureEquation;
import it.geoframe.blogspot.closureequation.equationstate.EquationState;
import it.geoframe.blogspot.numerical.rootfinding.Bisection;
import it.geoframe.blogspot.whetgeo1d.data.Geometry;
import it.geoframe.blogspot.whetgeo1d.data.ProblemQuantities;

/**
 * @author Niccolo` Tubini
 *
 */
public class FreezingSoilInternalEnergy extends EquationState {

	private Geometry geometry;
	private ProblemQuantities variables;
	private Bisection bisection;


	private double f;
	private double df;
	private double ddf;
	private double meltingTemperature;


	public FreezingSoilInternalEnergy(ClosureEquation closureEquation) {
		super(closureEquation);
		this.geometry = Geometry.getInstance();
		this.variables = ProblemQuantities.getInstance();
		this.bisection = new Bisection(this);
	}



	@Override
	public double equationState(double x, double y, int id, int element) {
		
		f = closureEquation.f(y, x, id);
		return (  ( super.closureEquation.parameters.specificThermalCapacitySoilParticles[id]*super.closureEquation.parameters.soilParticlesDensity[id]*(1.0-super.closureEquation.parameters.thetaS[id]) + 
				super.closureEquation.parameters.waterDensity*super.closureEquation.parameters.specificThermalCapacityWater*f +
				super.closureEquation.parameters.iceDensity*super.closureEquation.parameters.specificThermalCapacityIce*(closureEquation.f(y, 273.2, id)-f) )*(x-super.closureEquation.parameters.referenceTemperatureInternalEnergy)
				+ super.closureEquation.parameters.waterDensity*super.closureEquation.parameters.latentHeatFusion*f  ) * geometry.controlVolume[element];

	}


	@Override
	public double dEquationState(double x, double y, int id, int element) {

		f = closureEquation.f(y, x, id);
		df = closureEquation.df(y, x, id);
		return ( super.closureEquation.parameters.specificThermalCapacitySoilParticles[id]*super.closureEquation.parameters.soilParticlesDensity[id]*(1.0-super.closureEquation.parameters.thetaS[id]) + 
				super.closureEquation.parameters.waterDensity*super.closureEquation.parameters.specificThermalCapacityWater*f +
				super.closureEquation.parameters.iceDensity*super.closureEquation.parameters.specificThermalCapacityIce*(closureEquation.f(y, 273.2, id)-f)  +
				+ (super.closureEquation.parameters.waterDensity*super.closureEquation.parameters.specificThermalCapacityWater - super.closureEquation.parameters.iceDensity*super.closureEquation.parameters.specificThermalCapacityIce)*df*(x-super.closureEquation.parameters.referenceTemperatureInternalEnergy)
				+ super.closureEquation.parameters.waterDensity*super.closureEquation.parameters.latentHeatFusion*df ) * geometry.controlVolume[element];

	}


	@Override
	public double ddEquationState(double x, double y, int id, int element) {

		f = closureEquation.f(y, x, id);
		df = closureEquation.df(y, x, id);
		ddf = closureEquation.ddf(y, x, id);
				
		return ( 2*(super.closureEquation.parameters.waterDensity*super.closureEquation.parameters.specificThermalCapacityWater - super.closureEquation.parameters.iceDensity*super.closureEquation.parameters.specificThermalCapacityIce)*df +
				(super.closureEquation.parameters.waterDensity*super.closureEquation.parameters.specificThermalCapacityWater - super.closureEquation.parameters.iceDensity*super.closureEquation.parameters.specificThermalCapacityIce)*ddf*(x-super.closureEquation.parameters.referenceTemperatureInternalEnergy) +
				 super.closureEquation.parameters.waterDensity*super.closureEquation.parameters.latentHeatFusion*ddf  ) * geometry.controlVolume[element];

	}


	@Override
	public double p(double x, double y, int id, int element) {

		if(x<=variables.temperatureStar1[element]) {
//			System.out.println(x+"\t"+variables.temperatureStar1[element]+"\t"+ dEquationState(x, y, id, element));
			return dEquationState(x, y, id, element);  
		} else {
//			System.out.println(x+"\t"+variables.temperatureStar1[element]+"\t"+ dEquationState(variables.temperatureStar1[element], y, id, element));
			return dEquationState(variables.temperatureStar1[element], y, id, element);
		}  

	}


	@Override
	public double pIntegral(double x, double y, int id, int element) {

		if(x<=variables.temperatureStar1[element]) {
			return equationState(x, y, id, element);  
		} else {
			return equationState(variables.temperatureStar1[element], y, id, element) + dEquationState(variables.temperatureStar1[element], y, id, element)*(x-variables.temperatureStar1[element]);
		}
	}



	@Override
	public void computeXStar(double y, int id, int element) {
		meltingTemperature = super.closureEquation.parameters.referenceTemperatureInternalEnergy + 9.81*super.closureEquation.parameters.referenceTemperatureInternalEnergy/super.closureEquation.parameters.latentHeatFusion*y;
//		variables.temperatureStar1[element] =  bisection.findZero(273.15-10, meltingTemperature-1e-4, y, id, element);
		variables.temperatureStar1[element] =  findZero(273.15-10, meltingTemperature-1e-9, y, id, element);
//		System.out.println(y + "\t" + variables.temperatures[element] +"\t"+ closureEquation.f(y, 273, id));
		variables.temperatureStar2[element] = -9999.0;
		variables.temperatureStar3[element] = -9999.0;
		
	}

	@Override
	public double initialGuess(double x, int id, int element) {
		
		return Math.min(variables.temperatures[element], variables.temperatureStar1[element]);
		
	}
	
	private double findZero(double a, double b, double y, int id, int element) {
		int counter;
		double fa;
		double fb;
		double c;
		double fc;
		double tolerance = 1e-12;
		
		
		counter = 1;

		fa = ddEquationState(a, y, id, element);
		fb = ddEquationState(b, y, id, element);
		c = (a+b)/2;
		fc = ddEquationState(c, y, id, element);

		if(fc == 0.0) {
			return c;
		} else {
			while(Math.abs(a-b) > tolerance) {


				if(fa*fb > 0.0) {
//					System.out.println("\tBISECTION: error finding zero.");
					return meltingTemperature;
				}
//				System.out.println("a "+a+" b "+b+" fa*fb "+fa*fb);
				c = (a+b)/2;
				fc = ddEquationState(c, y, id, element);

				if(fc == 0.0) {
					return c;
				} else {
					if(fa*fc < 0.0) {
						b = c;
						fb = fc;
					} else {
						a = c;
						fa = fc;
					}
				}

				if(counter>150) {
					System.out.println("\tBISECTION: reached 250 iteration |a-b| = "+Math.abs(a-b));
					return c;
				}
				counter ++;
			}
		}
//		System.out.println("element: " +element+ "counter: "+counter+", intervall: "+Math.abs(a-b));

		return c;

	}


}
