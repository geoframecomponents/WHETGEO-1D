/*
 * GNU GPL v3 License
 *
 * Copyright 2021 Niccolo` Tubini
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

package it.geoframe.blogspot.whetgeo1d.data;


public class ProblemQuantities {
	
	private static ProblemQuantities uniqueInstance;
	
	public static ProblemQuantities getInstance() {
		return uniqueInstance;
	}
	
	public static ProblemQuantities getInstance(double[] icWaterSuction, double[] icTemperature, int[] equationStateID, int[] parameterID) {
		if (uniqueInstance == null) {
			uniqueInstance = new ProblemQuantities(icWaterSuction, icTemperature, equationStateID, parameterID);
		}
		return uniqueInstance;
	}
	
	public static ProblemQuantities getInstance(double[] icWaterSuction, double[] icTemperature, double[] icConcentration, int[] equationStateID, int[] parameterID) {
		if (uniqueInstance == null) {
			uniqueInstance = new ProblemQuantities(icWaterSuction, icTemperature, icConcentration, equationStateID, parameterID);
		}
		return uniqueInstance;
	}
	
	
	public double[] waterSuctions;
	public double[] temperatures;
	
	public double[] thetas;
    public double[] thetasNew;
    public double[] saturationDegree;
    public double[] iceContent;
	public double[] kappas;
	public double[] kappasInterface;
	public double[] volumes;
	public double[] volumesNew;
	public double[] darcyVelocities;
	public double[] darcyVelocitiesCapillary;
	public double[] darcyVelocitiesGravity;
	public double[] poreVelocities;
	public double[] celerities;        // Rasmussen et al. 2000
	public double[] kinematicRatio;  // Rasmussen et al. 2000
	public double[] waterSuctionStar1;
	public double[] waterSuctionStar2;
	public double[] waterSuctionStar3;
	
	public double[] internalEnergys;
	public double[] internalEnergysNew;
	public double[] heatCapacitys;
	public double[] heatCapacitysNew;
	public double[] waterCapacityTransported;
	public double[] lambdas;
	public double[] lambdasInterface;
	public double[] conductionHeatFluxs;
	public double[] advectionHeatFluxs;
	public double[] heatFluxs;
	public double[] heatSourcesSinksTerm;
	public double[] temperatureStar1;
	public double[] temperatureStar2;
	public double[] temperatureStar3;
	
	public int[] equationStateID;
	public int[] parameterID;
	
	public double runOff;
	public double kappaBottom;
	public double waterVolume;
	public double waterVolumeNew;
	public double errorVolume;
	public double volumeLost;
	public double richardsTopBCValue;
	public double richardsBottomBCValue;
	
	public double internalEnergy;
	public double internalEnergyNew;
	public double errorInternalEnergy;
	public double internalEnergyTopBCValue;
	public double shortWaveIn;
	public double longWaveIn;
	public double shortWaveOut;
	public double longWaveOut;
	public double airT;
	public double windVelocity;
	public double potentialLatentHeatFlux;
	public double actualLatentHeatFlux;
	public double sensibleHeatFlux;
	public double referenceHeight;
	public double surfaceRoughness;
	public double surfaceZeroHeightDisplacement;
	public double h1;
	public double h2;
	public double h3;
	public double h4;
	public double sensibleHeatCoefficient;
	public double internalEnergyBottomBCValue;
	public double heatFluxTop;
	public double heatFluxBottom;
	public double surfaceAlbedo;
	public double surfaceEmissivity;
	public double linearizedLongWaveOut;
	public final double constantStefanBoltzmann = 5.67*Math.pow(10, -8);
	public final double airDensity = 1.225;
	public final double airHeatCapacity = 1000;
	

	public double sumETs;
	public double[] ETs;
	
	public double[] concentrations;
	public double soluteTopBCValue;
	public double soluteBottomBCValue;
	public double waterVolumeConcentration;
	public double waterVolumeConcentrationNew;
	public double[] waterVolumeConcentrations;
	public double[] waterVolumeConcentrationsNew;
	
	public double[] dispersionCoefficients;
	public double[] dispersionFactors;
	public double[] thetasInterface;
	
	public double averageSoluteConcentration;
	public double averageWaterVolumeSoluteConcentration;
	public double[] dispersionSoluteFluxes;
	public double[] advectionSoluteFluxes;
	public double[] soluteFluxes;
	public double errorWaterVolumeConcentration;
	public double[] soluteQuantitiesTransported;
	public double[] soluteSourcesSinksTerm;
	public double sumSoluteSourceSinkTerm;
	public double[] timeVariationWaterVolumesConcentration;
	public double[] tortuosityFactors;
	public double[] tortuosityFactorsInterface;
	
	
	
	private ProblemQuantities(double[] icWaterSuction, double[] icTemperature, int[] equationStateID, int[] parameterID) {
		
		waterSuctions = icWaterSuction.clone();
		temperatures = icTemperature.clone();
		thetas = new double[icWaterSuction.length];
		thetasNew = new double[icWaterSuction.length];
		saturationDegree = new double[icWaterSuction.length];
		iceContent = new double[icWaterSuction.length];
		kappas = new double[icWaterSuction.length];
		kappasInterface = new double[icWaterSuction.length+1];
		volumes = new double[icWaterSuction.length];
		volumesNew = new double[icWaterSuction.length];
		darcyVelocities = new double[icWaterSuction.length+1];
		darcyVelocitiesCapillary = new double[icWaterSuction.length+1];
		darcyVelocitiesGravity = new double[icWaterSuction.length+1];
		poreVelocities = new double[icWaterSuction.length+1];
		celerities = new double[icWaterSuction.length+1];
		kinematicRatio = new double[icWaterSuction.length+1];
		ETs = new double[icWaterSuction.length];
		waterSuctionStar1 = new double[icWaterSuction.length];
		waterSuctionStar2 = new double[icWaterSuction.length];
		waterSuctionStar3 = new double[icWaterSuction.length];
		
		internalEnergys = new double[icWaterSuction.length];
		internalEnergysNew = new double[icWaterSuction.length];
		heatCapacitys = new double[icWaterSuction.length];
		heatCapacitysNew = new double[icWaterSuction.length];
		waterCapacityTransported = new double[icWaterSuction.length];
		lambdas = new double[icWaterSuction.length];
		lambdasInterface = new double[icWaterSuction.length+1];
		conductionHeatFluxs = new double[icWaterSuction.length+1];
		advectionHeatFluxs = new double[icWaterSuction.length+1];
		heatFluxs = new double[icWaterSuction.length+1];
		heatSourcesSinksTerm = new double[icWaterSuction.length];
		temperatureStar1 = new double[icWaterSuction.length];
		temperatureStar2 = new double[icWaterSuction.length];
		temperatureStar3 = new double[icWaterSuction.length];
		
		
		this.equationStateID = equationStateID.clone();
		this.parameterID = parameterID.clone();
		
		
	}
	private ProblemQuantities(double[] icWaterSuction, double[] icTemperature, double[] icConcentration, int[] equationStateID, int[] parameterID) {
		waterSuctions = icWaterSuction.clone();
		temperatures = icTemperature.clone();
		
		thetas = new double[icWaterSuction.length];
		thetasNew = new double[icWaterSuction.length];
		saturationDegree = new double[icWaterSuction.length];
		iceContent = new double[icWaterSuction.length];
		kappas = new double[icWaterSuction.length];
		kappasInterface = new double[icWaterSuction.length+1];
		volumes = new double[icWaterSuction.length];
		volumesNew = new double[icWaterSuction.length];
		darcyVelocities = new double[icWaterSuction.length+1];
		darcyVelocitiesCapillary = new double[icWaterSuction.length+1];
		darcyVelocitiesGravity = new double[icWaterSuction.length+1];
		poreVelocities = new double[icWaterSuction.length+1];
		celerities = new double[icWaterSuction.length+1];
		kinematicRatio = new double[icWaterSuction.length+1];
		ETs = new double[icWaterSuction.length];
		waterSuctionStar1 = new double[icWaterSuction.length];
		waterSuctionStar2 = new double[icWaterSuction.length];
		waterSuctionStar3 = new double[icWaterSuction.length];
		
		internalEnergys = new double[icWaterSuction.length];
		internalEnergysNew = new double[icWaterSuction.length];
		heatCapacitys = new double[icWaterSuction.length];
		heatCapacitysNew = new double[icWaterSuction.length];
		waterCapacityTransported = new double[icWaterSuction.length];
		lambdas = new double[icWaterSuction.length];
		lambdasInterface = new double[icWaterSuction.length+1];
		conductionHeatFluxs = new double[icWaterSuction.length+1];
		advectionHeatFluxs = new double[icWaterSuction.length+1];
		heatFluxs = new double[icWaterSuction.length+1];
		heatSourcesSinksTerm = new double[icWaterSuction.length];
		temperatureStar1 = new double[icWaterSuction.length];
		temperatureStar2 = new double[icWaterSuction.length];
		temperatureStar3 = new double[icWaterSuction.length];
		
		concentrations = icConcentration.clone();
		
		waterVolumeConcentrations = new double[icWaterSuction.length];
		waterVolumeConcentrationsNew = new double[icWaterSuction.length];
		dispersionCoefficients = new double[icWaterSuction.length+1];
		dispersionFactors = new double[icWaterSuction.length+1];
		thetasInterface = new double[icWaterSuction.length+1];
		dispersionSoluteFluxes = new double[icWaterSuction.length+1];
		advectionSoluteFluxes = new double[icWaterSuction.length+1];
		soluteFluxes = new double[icWaterSuction.length+1];
		soluteQuantitiesTransported = new double[icWaterSuction.length];
		soluteSourcesSinksTerm = new double[icWaterSuction.length];
		tortuosityFactors = new double[icWaterSuction.length];
		tortuosityFactorsInterface = new double[icWaterSuction.length+1];
		timeVariationWaterVolumesConcentration = new double[icWaterSuction.length];
		

			
		this.equationStateID = equationStateID.clone();
		this.parameterID = parameterID.clone();
	}

}




