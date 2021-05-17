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

package heat;

import java.net.URISyntaxException;
import java.util.*;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorReader;

import it.geoframe.blogspot.buffer.buffertowriter.HeatDiffusionFreezingThawingBufferWithSurfaceEnergyBudget1D;
import it.geoframe.blogspot.whetgeo1d.heatsolver.HeatDiffusionFreezingThawingSolverWithSurfaceEnergyBalance1DMain;
import it.geoframe.blogspot.whetgeo1d.heatsolver.HeatDiffusionFreezingThawingSolverWithSurfaceEnergyBalance1DMain2;
import it.geoframe.blogspot.whetgeo1d.heatsolver.HeatDiffusionSolverWithSurfaceEnergyBalance1DMain;
import it.geoframe.blogspot.whetgeo1d.heatsolver.SFCC;
import it.geoframe.blogpsot.netcdf.monodimensionalproblemtimedependent.ReadNetCDFHeatDiffusionGrid1D;
import it.geoframe.blogpsot.netcdf.monodimensionalproblemtimedependent.WriteNetCDFHeatDiffusionFreezingThawingWithSurfaceEnergyBudget1DDouble;

import org.junit.Test;

/**
 * Test the {@link TestSFCC} module.
 * 
 * 
 * @author Niccolo' Tubini
 */
public class TestSFCC {

	@Test
	public void Test() throws Exception {


		String startDate = "2003-01-01 00:00";
		String endDate = "2003-01-01 00:00";
		int timeStepMinutes = 60;
		String fId = "ID";
		String oms ="C:/Users/Niccolo/OMS/OMS_WHETGEO1D_surface_energy_balance";
		String pathAirT = oms+"/data/TimeSeries/Roncafort_y20032021/airT_T0135.csv";
		String pathWindVelocity = oms+"/data/TimeSeries/Roncafort_y20032021/windVelocity_T0135.csv";
		String pathSW = oms+"/data/TimeSeries/Roncafort_y20032021/TotalSolarRadiation_T0135.csv";
		String pathLW = oms+"/output/Roncafort/LWDownwelling_T0135.csv";
		String pathLE = oms+"/output/Roncafort/LatentHeat_PT_T0135.csv";
		String pathBottomBC = oms+"/data/TimeSeries/Roncafort_y20032021/tmpNoFlux_T0135.csv";
		String pathSaveDates = oms+"/data/TimeSeries/Roncafort_y20032021/tmpsavedates_T0135.csv"; 
		String pathGrid =  oms+"/data/Grid_NetCDF/New_diffusion_soil_class_wetwet.nc";
		String pathOutput = "resources/output/provaFTSurfaceEnergyBalance_FT_LE.nc";
		
//		String topBC = "Top dirichlet";
		String bottomBC = "Bottom Neumann";

		String outputDescription = "\n"
				+ "Initial condition constant temperature\n		"
				+ "DeltaT: 3600s\n		"
				+ "Picard iteration: 1\n		";
		
		int writeFrequency = 1000000;

		OmsTimeSeriesIteratorReader airTReader = getTimeseriesReader(pathAirT, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader windVelocityReader = getTimeseriesReader(pathWindVelocity, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader swReader = getTimeseriesReader(pathSW, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader lwReader = getTimeseriesReader(pathLW, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader leReader = getTimeseriesReader(pathLE, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader bottomBCReader = getTimeseriesReader(pathBottomBC, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader saveDatesReader = getTimeseriesReader(pathSaveDates, fId, startDate, endDate, timeStepMinutes);

		HeatDiffusionFreezingThawingBufferWithSurfaceEnergyBudget1D buffer = new HeatDiffusionFreezingThawingBufferWithSurfaceEnergyBudget1D();
		WriteNetCDFHeatDiffusionFreezingThawingWithSurfaceEnergyBudget1DDouble writeNetCDF = new WriteNetCDFHeatDiffusionFreezingThawingWithSurfaceEnergyBudget1DDouble();
		ReadNetCDFHeatDiffusionGrid1D readNetCDF = new ReadNetCDFHeatDiffusionGrid1D();
		
		SFCC solver = new SFCC();
		
		
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
		solver.beta0 = -766.45;
		solver.referenceTemperatureSWRC = 278.15;
		solver.typeClosureEquation = new String[] {"VanGenuchtenDallAmico"};
		solver.typeEquationState = new String[] {"FreezingSoilInternalEnergy"};
		solver.typeThermalConductivity = new String[] {"Cosenza"};
		solver.interfaceThermalConductivityModel = "max";
		solver.bottomBCType = bottomBC;
		solver.surfaceAlbedoType = "Constant";
		solver.surfaceEmissivityType = "Constant";
		solver.surfaceAereodynamicResistanceType = "Neutral";
		solver.surfaceWaterVaporResistanceType = "Feddes";
		solver.h1 = -0.5;
		solver.h2 = -5;
		solver.h3 = -15.0;
		solver.h4 = -50.0;
		
		solver.stationID = 135;
		solver.delta = 0;
		solver.tTimeStep = 3600;
		solver.timeDelta = 3600;
		solver.newtonTolerance = 0.003337000000000;//1e-11/(1000*333700);
		solver.nestedNewton = 1;
		solver.picardIteration = 1;

		buffer.writeFrequency = writeFrequency;
		
		writeNetCDF.fileName = pathOutput;
		writeNetCDF.briefDescritpion = outputDescription;
		writeNetCDF.pathGrid = pathGrid;
		writeNetCDF.pathBottomBC = pathBottomBC; 
//		writeNetCDF.pathTopBC = pathTopBC; 
		writeNetCDF.bottomBC = bottomBC;
//		writeNetCDF.topBC = topBC;
		writeNetCDF.swrcModel = "VG";
		writeNetCDF.soilThermalConductivityModel = "Cosenza";
		writeNetCDF.interfaceConductivityModel = "max";
		writeNetCDF.writeFrequency = writeFrequency;
		writeNetCDF.spatialCoordinate = readNetCDF.eta;
		writeNetCDF.dualSpatialCoordinate = readNetCDF.etaDual;	
		writeNetCDF.controlVolume = readNetCDF.controlVolume;
		writeNetCDF.psi = readNetCDF.psi;
		writeNetCDF.temperatureIC = readNetCDF.temperatureIC;
		writeNetCDF.timeUnits = "Minutes since 01/01/1970 01:00:00 UTC";
		writeNetCDF.timeZone = "UTC"; 
		writeNetCDF.fileSizeMax = 10000;
		
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
			System.out.println(swReader.tCurrent);
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
//		ReadNetCDFRichardsOutput1D readTestData = new ReadNetCDFRichardsOutput1D();
//		readTestData.richardsOutputFilename = "resources/Output/Check_VG.nc";
//		readTestData.read();
//		
//		ReadNetCDFRichardsOutput1D readSimData = new ReadNetCDFRichardsOutput1D();
//		readSimData.richardsOutputFilename = pathOutput.replace(".nc","_0000.nc");
//		readSimData.read();
//
//		for(int k=0; k<readSimData.psi[(readSimData.psi.length)-1].length; k++) {
//			if(Math.abs(readSimData.psi[(readSimData.psi.length)-1][k]-readTestData.psi[(readTestData.psi.length)-1][k])>Math.pow(10,-11)) {
//				System.out.println("\n\n\t\tERROR: psi mismatch");
//			}
//		}

	}

	private OmsTimeSeriesIteratorReader getTimeseriesReader( String inPath, String id, String startDate, String endDate,
			int timeStepMinutes ) throws URISyntaxException {
		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file = inPath;
		reader.idfield = "ID";
		reader.tStart = startDate;
		reader.tTimestep = timeStepMinutes;
		reader.tEnd = endDate;
		reader.fileNovalue = "-9999";
		reader.initProcess();
		return reader;
	}
}
