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

package org.geoframe.solvers;

import java.io.File;
import java.nio.file.Files;
import java.util.HashMap;

import org.geoframe.whetgeo.WGTestCase;
import org.geoframe.whetgeo1d.core.boundaryconditions.IBoundaryCondition;
import org.geoframe.whetgeo1d.core.boundaryconditions.IBoundaryCondition.DiffusionBoundaryConditionType;
import org.geoframe.whetgeo1d.solvers.HeatDiffusionSolver1D;
import org.hortonmachine.gears.io.geoframe.HeatDiffusionBuffer1D;
import org.hortonmachine.gears.io.geoframe.ReadNetCDFHeatDiffusionGrid1D;
import org.hortonmachine.gears.io.geoframe.WriteNetCDFHeatDiffusion1DDouble;
import org.hortonmachine.gears.io.timedependent.OmsTimeSeriesIteratorReader;

/**
 * Test the {@link TestHeatDiffusion} module.
 * 
 * 
 * @author Niccolo' Tubini
 */
public class TestHeatDiffusion extends WGTestCase{

	public void testHeatDiffusion() throws Exception {


		String startDate = "2015-01-15 00:00";
		String endDate = "2015-12-15 00:00";
		int timeStepMinutes = 60;
		String fId = "ID";
				
		String pathTopBC = getRes("/input/TimeSeries/tempTop.csv");
		String pathBottomBC = getRes("/input/TimeSeries/tempBottom.csv");
//		String pathSaveDates = getRes("/input/TimeSeries/save.csv"); 
		String pathGrid = "/home/hydrologis/TMP/UNITN/whetgeo1d/Heat_diffusion.nc";// getRes("/input/Grid_NetCDF/Heat_diffusion.nc");
		
		File tempFile = Files.createTempFile("Sim_heat_diffusion", ".nc").toFile();
		String pathOutput = tempFile.getAbsolutePath();
		
		var topBC = IBoundaryCondition.DiffusionBoundaryConditionType.TOP_DIRICHLET;
		var bottomBC = IBoundaryCondition.DiffusionBoundaryConditionType.BOTTOM_DIRICHLET;

		String outputDescription = "\n"
				+ "Initial condition constant temperature\n		"
				+ "DeltaT: 3600s\n		"
				+ "Picard iteration: 1\n		";
		
		int writeFrequency = 100;
		
		OmsTimeSeriesIteratorReader topBCReader = getTimeseriesReader(pathTopBC, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader bottomBCReader = getTimeseriesReader(pathBottomBC, fId, startDate, endDate, timeStepMinutes);
//		OmsTimeSeriesIteratorReader saveDatesReader = getTimeseriesReader(pathSaveDates, fId, startDate, endDate, timeStepMinutes);

		HeatDiffusionBuffer1D buffer = new HeatDiffusionBuffer1D();
		WriteNetCDFHeatDiffusion1DDouble writeNetCDF = new WriteNetCDFHeatDiffusion1DDouble();
		ReadNetCDFHeatDiffusionGrid1D readNetCDF = new ReadNetCDFHeatDiffusionGrid1D();
		
		HeatDiffusionSolver1D solver = new HeatDiffusionSolver1D();
		
		
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
		solver.inEquationStateID = readNetCDF.equationStateID;
		solver.inParameterID = readNetCDF.parameterID;
		solver.beta0 = -766.45;
		solver.referenceTemperatureSWRC = 278.15;
		solver.typeClosureEquation = new String[] {"Van Genuchten"};
		solver.typeEquationState = new String[] {"SoilInternalEnergy"};
		solver.typeThermalConductivity = new String[] {"Cosenza"};
		solver.interfaceThermalConductivityModel = "max";
		solver.topBCType = topBC;
		solver.bottomBCType = bottomBC;
		solver.delta = 0;
		solver.tTimeStep = 3600;
		solver.timeDelta = 3600;
		solver.newtonTolerance = Math.pow(10,-6);
		solver.nestedNewton = 1;
		solver.picardIteration = 1;
		solver.stationID = 0;

		buffer.writeFrequency = writeFrequency;
		
		writeNetCDF.fileName = pathOutput;
		writeNetCDF.briefDescritpion = outputDescription;
		writeNetCDF.pathGrid = pathGrid;
		writeNetCDF.pathBottomBC = pathBottomBC; 
		writeNetCDF.pathTopBC = pathTopBC; 
		writeNetCDF.bottomBC = bottomBC.toString();
		writeNetCDF.topBC = topBC.toString();
		writeNetCDF.swrcModel = "VG";
		writeNetCDF.soilThermalConductivityModel = "Mualem VG no temperature";
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
		
		while( topBCReader.doProcess  ) {
		
			
			topBCReader.nextRecord();	
			HashMap<Integer, double[]> bCValueMap = topBCReader.outData;
			solver.inTopBC= bCValueMap;


			bottomBCReader.nextRecord();
			bCValueMap = bottomBCReader.outData;
			solver.inBottomBC = bCValueMap;

//			saveDatesReader.nextRecord();
//			bCValueMap = saveDatesReader.outData;
//			solver.inSaveDate = bCValueMap;
			
			solver.inCurrentDate = topBCReader.tCurrent;
			
			solver.solve();

			
			buffer.inputDate = solver.inCurrentDate;
			buffer.doProcessBuffer = true;// solver.doProcessBuffer;
			buffer.inputVariable = solver.outputToBuffer;
			
			buffer.solve();
			

			writeNetCDF.variables = buffer.myVariable;
			writeNetCDF.doProcess = topBCReader.doProcess;
			writeNetCDF.writeNetCDF();


		}

		topBCReader.close();
		bottomBCReader.close();
				

	}

	
}
