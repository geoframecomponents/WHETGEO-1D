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

import it.geoframe.blogspot.buffer.buffertowriter.HeatDiffusionBuffer1D;
import it.geoframe.blogspot.buffer.buffertowriter.HeatDiffusionFreezingThawingBufferWithSurfaceEnergyBudget1D;
import it.geoframe.blogspot.whetgeo1d.heatsolver.HeatDiffusionFreezingThawingSolver1DMain;
import it.geoframe.blogspot.whetgeo1d.heatsolver.HeatDiffusionSolver1DMain;
import it.geoframe.blogpsot.netcdf.monodimensionalproblemtimedependent.ReadNetCDFHeatDiffusionGrid1D;
import it.geoframe.blogpsot.netcdf.monodimensionalproblemtimedependent.WriteNetCDFHeatDiffusion1DDouble;
import it.geoframe.blogpsot.netcdf.monodimensionalproblemtimedependent.WriteNetCDFHeatDiffusionFreezingThawingWithSurfaceEnergyBudget1DDouble;

import org.junit.Test;

/**
 * Test the {@link TestHeatDiffusionFrozenSoilLake} module.
 * 
 * 
 * @author Niccolo' Tubini
 */
public class TestHeatDiffusionFrozenSoilLake {

	@Test
	public void Test() throws Exception {


		String startDate = "1950-01-01 00:00";
		String endDate = "2020-01-01 00:00";
		int timeStepMinutes = 60;
		String fId = "ID";
		String oms_path = "C:/Users/Niccolo/OMS/OMS_Project_WHETGEO1D_mio/";		
		String pathTopBC = oms_path+"data/TimeSeries/prova.csv";
		String pathBottomBC = oms_path+"data/TimeSeries/prova_bottom.csv";
		String pathSaveDates = oms_path+"data/TimeSeries/prova_save_dates.csv"; 
		String pathGrid =  oms_path+"data/Grid_NetCDF/Heat_conduction_frozensoil_lake.nc";
		String pathOutput = oms_path+"output/Sim_heat_diffusion_frozensoil_lake.nc";
		
		String topBC = "Top dirichlet";
		String bottomBC = "Bottom Neumann";

		String outputDescription = "\n"
				+ "Initial condition constant temperature\n		"
				+ "DeltaT: 3600s\n		"
				+ "Picard iteration: 1\n		";
		
		int writeFrequency = 3000;
		
		OmsTimeSeriesIteratorReader topBCReader = getTimeseriesReader(pathTopBC, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader bottomBCReader = getTimeseriesReader(pathBottomBC, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader saveDatesReader = getTimeseriesReader(pathSaveDates, fId, startDate, endDate, timeStepMinutes);

		HeatDiffusionFreezingThawingBufferWithSurfaceEnergyBudget1D buffer = new HeatDiffusionFreezingThawingBufferWithSurfaceEnergyBudget1D();
		WriteNetCDFHeatDiffusionFreezingThawingWithSurfaceEnergyBudget1DDouble writeNetCDF = new WriteNetCDFHeatDiffusionFreezingThawingWithSurfaceEnergyBudget1DDouble();
		ReadNetCDFHeatDiffusionGrid1D readNetCDF = new ReadNetCDFHeatDiffusionGrid1D();
		
		HeatDiffusionFreezingThawingSolver1DMain solver = new HeatDiffusionFreezingThawingSolver1DMain();
		
		
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
		solver.meltingTemperature = new double[] {-9999,273.149,-9999};//readNetCDF.meltingTemperature;
		solver.beta0 = -766.45;
		solver.referenceTemperatureSWRC = 278.15;
		solver.typeClosureEquation = new String[] {"None","VanGenuchtenDallAmico"};
		solver.typeEquationState = new String[] {"WaterInternalEnergy","FreezingSoilInternalEnergy"};
		solver.typeThermalConductivity = new String[] {"Water","Cosenza"};
		solver.interfaceThermalConductivityModel = "harmonic mean";
		solver.topBCType = topBC;
		solver.bottomBCType = bottomBC;
		solver.delta = 0;
		solver.tTimeStep = 3600;
		solver.timeDelta = 3600;
		solver.newtonTolerance = 3*Math.pow(10,-3);
		solver.nestedNewton = 1;
		solver.picardIteration = 1;
		solver.stationID = 0;

		buffer.writeFrequency = writeFrequency;
		
		writeNetCDF.fileName = pathOutput;
		writeNetCDF.briefDescritpion = outputDescription;
		writeNetCDF.pathGrid = pathGrid;
		writeNetCDF.pathBottomBC = pathBottomBC; 
		writeNetCDF.pathTopBC = pathTopBC; 
		writeNetCDF.bottomBC = bottomBC;
		writeNetCDF.topBC = topBC;
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

			saveDatesReader.nextRecord();
			bCValueMap = saveDatesReader.outData;
			solver.inSaveDate = bCValueMap;
			
			solver.inCurrentDate = topBCReader.tCurrent;
			
			solver.solve();

			
			buffer.inputDate = solver.inCurrentDate;
			buffer.doProcessBuffer = solver.doProcessBuffer;
			buffer.inputVariable = solver.outputToBuffer;
			
			buffer.solve();
			

			writeNetCDF.variables = buffer.myVariable;
			writeNetCDF.doProcess = topBCReader.doProcess;
			writeNetCDF.writeNetCDF();


		}

		topBCReader.close();
		bottomBCReader.close();
				

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
