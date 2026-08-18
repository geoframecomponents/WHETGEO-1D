/*
  * GNU GPL v3 License
 *
 * Copyright 2020 Niccolo` Tubini
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

package org.geoframe.solvers.toreview;

import java.io.File;
import java.nio.file.Files;
import java.util.HashMap;

import org.geoframe.whetgeo.WGTestCase;
import org.geoframe.whetgeo1d.core.boundaryconditions.IBoundaryCondition;
import org.geoframe.whetgeo1d.solvers.HeatDiffusionSolverWithSurfaceEnergyBalance1D;
import org.hortonmachine.gears.io.geoframe.HeatDiffusionBufferWithSurfaceEnergyBudget1D;
import org.hortonmachine.gears.io.geoframe.ReadNetCDFHeatDiffusionGrid1D;
import org.hortonmachine.gears.io.geoframe.ReadNetCDFHeatDiffusionOutput1D;
import org.hortonmachine.gears.io.geoframe.WriteNetCDFHeatDiffusionWithSurfaceEnergyBudget1DDouble;
import org.hortonmachine.gears.io.timedependent.OmsTimeSeriesIteratorReader;

/**
 * Test the {@link TestHeatDiffusionSurfaceEnergyBalance} module.
 * 
 * 
 * @author Niccolo' Tubini
 */
public class TestHeatDiffusionSurfaceEnergyBalance  extends WGTestCase {

	public void testHeatDiffusionSurfaceEnergyBalance() throws Exception {
		String startDate = "2003-01-01 00:00";
		String endDate = "2007-01-01 00:00";
		int timeStepMinutes = 60;
		String fId = "ID";
		
		String pathAirT = getRes("/input/TimeSeries/airT_T0135.csv");
		String pathWindVelocity = getRes("/input/TimeSeries/windVelocity_T0135.csv");
		String pathSW = getRes("/input/TimeSeries/TotalSolarRadiation_T0135.csv");
		String pathLW = getRes("/input/TimeSeries/LWDownwelling_T0135.csv");
		String pathLE = getRes("/input/TimeSeries/LatentHeat_PT_T0135.csv");
		String pathBottomBC = getRes("/input/TimeSeries/noFlux_T0135.csv");
		String pathSaveDates = getRes("/input/TimeSeries/saveDates_T0135.csv"); 
		
		String pathGrid =  getRes("/input/Grid_NetCDF/Heat_diffusion.nc");
		
		File tempFile = Files.createTempFile("Sim_heat_diffusion", ".nc").toFile();
		String pathOutput = tempFile.getAbsolutePath();
		
		var bottomBC = IBoundaryCondition.DiffusionBoundaryConditionType.BOTTOM_NEUMANN;

		String outputDescription = "\n"
				+ "Pure heat diffusion driven by the surface energy budget. Soil is saturated.";
		
		int writeFrequency = 10000;

		OmsTimeSeriesIteratorReader airTReader = getTimeseriesReader(pathAirT, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader windVelocityReader = getTimeseriesReader(pathWindVelocity, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader swReader = getTimeseriesReader(pathSW, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader lwReader = getTimeseriesReader(pathLW, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader leReader = getTimeseriesReader(pathLE, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader bottomBCReader = getTimeseriesReader(pathBottomBC, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader saveDatesReader = getTimeseriesReader(pathSaveDates, fId, startDate, endDate, timeStepMinutes);

		HeatDiffusionBufferWithSurfaceEnergyBudget1D buffer = new HeatDiffusionBufferWithSurfaceEnergyBudget1D();
		WriteNetCDFHeatDiffusionWithSurfaceEnergyBudget1DDouble writeNetCDF = new WriteNetCDFHeatDiffusionWithSurfaceEnergyBudget1DDouble();
		ReadNetCDFHeatDiffusionGrid1D readNetCDF = new ReadNetCDFHeatDiffusionGrid1D();
		HeatDiffusionSolverWithSurfaceEnergyBalance1D solver = new HeatDiffusionSolverWithSurfaceEnergyBalance1D();
		
		
		readNetCDF.gridFilename = pathGrid;
		
		readNetCDF.read();
		
		
		solver.z = readNetCDF.z;
		solver.spaceDeltaZ = readNetCDF.spaceDelta;
		solver.psiIC = readNetCDF.psi;
		solver.temperature = readNetCDF.temperatureIC;
		solver.controlVolume = readNetCDF.controlVolume;
		solver.soilParticlesDensity = readNetCDF.soilParticlesDensity;
		solver.thermalConductivitySoilParticles = readNetCDF.soilParticlesThermalConductivity;
		solver.specificThermalCapacitySoilParticles = readNetCDF.soilParticlesSpecificHeatCapacity;
		solver.meltingTemperature = readNetCDF.meltingTemperature;
		solver.ks = readNetCDF.Ks;
		solver.thetaS = readNetCDF.thetaS;
		solver.thetaR = readNetCDF.thetaR;
		solver.par1SWRC = readNetCDF.par1SWRC;
		solver.par2SWRC = readNetCDF.par2SWRC;
		solver.par3SWRC = readNetCDF.par3SWRC;
		solver.par4SWRC = readNetCDF.par4SWRC;
		solver.par5SWRC = readNetCDF.par5SWRC;
		solver.alphaSpecificStorage = readNetCDF.alphaSS;
		solver.betaSpecificStorage = readNetCDF.betaSS;
		solver.surfaceAlbedo = 0.1;
		solver.surfaceEmissivity = 0.95;
		solver.referenceHeight = 10;
		solver.surfaceRoughness = 0.01;
		solver.surfaceZeroHeightDisplacement = 0.0;
		solver.inEquationStateID = readNetCDF.equationStateID;
		solver.inParameterID = readNetCDF.parameterID;
		solver.typeClosureEquation = new String[] {"Van Genuchten"};
		solver.typeEquationState = new String[] {"SoilInternalEnergy"};
		solver.typeThermalConductivity = new String[] {"Cosenza"};
		solver.interfaceThermalConductivityModel = "max";
		solver.bottomBCType = bottomBC;
		solver.surfaceAlbedoType = "Constant";
		solver.surfaceEmissivityType = "Constant";
		solver.surfaceAereodynamicResistanceType = "NeutralCondition";
		solver.surfaceWaterVaporResistanceType = "Feddes";
		solver.h1 = 0.1;
		solver.h2 = -5.0;
		solver.h3 = -15;
		solver.h4 = -50;

		solver.stationID = 135;
		solver.delta = 0;
		solver.tTimeStep = 3600;
		solver.timeDelta = 3600;
		solver.newtonTolerance = Math.pow(10,-5);
		solver.nestedNewton = 1;
		solver.picardIteration = 1;

		buffer.writeFrequency = writeFrequency;
		
		writeNetCDF.fileName = pathOutput;
		writeNetCDF.briefDescritpion = outputDescription;
		writeNetCDF.pathGrid = pathGrid;
		writeNetCDF.pathBottomBC = pathBottomBC; 
		writeNetCDF.bottomBC = bottomBC.name();
		writeNetCDF.swrcModel = "VG";
		writeNetCDF.soilThermalConductivityModel = "Cosenza";
		writeNetCDF.interfaceConductivityModel = "max";
		writeNetCDF.writeFrequency = writeFrequency;
		writeNetCDF.spatialCoordinate = readNetCDF.eta;
		writeNetCDF.dualSpatialCoordinate = readNetCDF.etaDual;	
		writeNetCDF.controlVolume = readNetCDF.controlVolume;
		writeNetCDF.psi = readNetCDF.psi;
		writeNetCDF.temperatureIC = readNetCDF.temperatureIC;
		writeNetCDF.timeUnits = "Minutes since 01/01/1970 00:00:00 UTC";
		writeNetCDF.timeZone = "UTC"; 
		writeNetCDF.fileSizeMax = 1000;
		
		while( swReader.doProcess  ) {
		
			
			swReader.nextRecord();	
			HashMap<Integer, double[]> bCValueMap = swReader.outData;
			solver.inShortWave= bCValueMap;
			
			lwReader.nextRecord();	
			bCValueMap = lwReader.outData;
			solver.inLongWave= bCValueMap;
			
			leReader.nextRecord();	
			bCValueMap = leReader.outData;
			solver.inPotentialLatentHeatFlux= bCValueMap;
			
			airTReader.nextRecord();	
			bCValueMap = airTReader.outData;
			solver.inAirT= bCValueMap;
			
			windVelocityReader.nextRecord();	
			bCValueMap = windVelocityReader.outData;
			solver.inWindVelocity= bCValueMap;

			bottomBCReader.nextRecord();
			bCValueMap = bottomBCReader.outData;
			solver.inBottomBC = bCValueMap;

			saveDatesReader.nextRecord();
			bCValueMap = saveDatesReader.outData;
			solver.inSaveDate = bCValueMap;
			
			solver.inCurrentDate = swReader.tCurrent;
			solver.solve();

			
			buffer.inputDate = solver.inCurrentDate;
			buffer.doProcessBuffer = solver.doProcessBuffer;
			buffer.inputVariable = solver.outputToBuffer;
			
			buffer.solve();
			

			writeNetCDF.variables = buffer.myVariable;
			writeNetCDF.doProcess = swReader.doProcess;
			writeNetCDF.writeNetCDF();


		}

		swReader.close();
		lwReader.close();
		bottomBCReader.close();
				
		/*
		 * ASSERT 
		 */
//		System.out.println("Assert");
//		ReadNetCDFHeatDiffusionOutput1D readTestData = new ReadNetCDFHeatDiffusionOutput1D();
//		readTestData.gridFilename = "resources/Output/Check_heat_diffusion_0000.nc";
//		readTestData.read();
//		
//		ReadNetCDFHeatDiffusionOutput1D readSimData = new ReadNetCDFHeatDiffusionOutput1D();
//		readSimData.gridFilename = pathOutput.replace(".nc","_0000.nc");
//		readSimData.read();
//
//		for(int k=0; k<readSimData.temperature[(readSimData.temperature.length)-1].length; k++) {
//			if(Math.abs(readSimData.temperature[(readSimData.temperature.length)-1][k]-readTestData.temperature[(readTestData.temperature.length)-1][k])>Math.pow(10,-11)) {
//				System.out.println("\n\n\t\tERROR: temperature mismatch");
//			}
//		}

	}
}
