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
import org.geoframe.whetgeo1d.core.boundaryconditions.IBoundaryCondition.RichardsBoundaryConditionType;
import org.geoframe.whetgeo1d.solvers.RichardsSolver1D;
import org.hortonmachine.gears.io.geoframe.ReadNetCDFRichardsGrid1D;
import org.hortonmachine.gears.io.geoframe.ReadNetCDFRichardsOutput1D;
import org.hortonmachine.gears.io.geoframe.RichardsBuffer1D;
import org.hortonmachine.gears.io.geoframe.WriteNetCDFRichards1DDouble;
import org.hortonmachine.gears.io.timedependent.OmsTimeSeriesIteratorReader;

/**
 * Test the {@link TestVanGenuchten} module.
 * 
 * This test consider an initial hydrostatic condition with Neumann boundary condition at the 
 * top and free drainage at the bottom. 
 * 
 * @author Niccolo' Tubini
 */
public class TestVanGenuchten extends WGTestCase {

	public void testVanGenuchten() throws Exception {


		String startDate = "2015-01-15 00:00";
		String endDate = "2015-12-15 00:00";
		int timeStepMinutes = 60;
		String fId = "ID";
				
		String pathTopBC = getRes("/input/TimeSeries/precip.csv");
		String pathBottomBC = getRes("/input/TimeSeries/bottom.csv");
		String pathSaveDates = getRes("/input/TimeSeries/save.csv"); 
		String pathGrid =  getRes("/input/Grid_NetCDF/RichardsCoupled_VG.nc");
		File tempFile = Files.createTempFile("Sim_RichardsCoupled_VG_", ".nc").toFile();
		String pathOutput = tempFile.getAbsolutePath();
		
		
		var topBC = RichardsBoundaryConditionType.TOP_COUPLED;
		var bottomBC = RichardsBoundaryConditionType.BOTTOM_FREE_DRAINAGE;

		String outputDescription = "\n"
				+ "Initial condition hydrostatic no ponding\n		"
				+ "DeltaT: 1800s\n		"
				+ "Picard iteration: 1\n		";
		
		int writeFrequency = 1000000;
		
		OmsTimeSeriesIteratorReader topBCReader = getTimeseriesReader(pathTopBC, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader bottomBCReader = getTimeseriesReader(pathBottomBC, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader saveDatesReader = getTimeseriesReader(pathSaveDates, fId, startDate, endDate, timeStepMinutes);

		RichardsBuffer1D buffer = new RichardsBuffer1D();
		WriteNetCDFRichards1DDouble writeNetCDF = new WriteNetCDFRichards1DDouble();
		ReadNetCDFRichardsGrid1D netcdfReader = new ReadNetCDFRichardsGrid1D();
		
		RichardsSolver1D R1DSolver = new RichardsSolver1D();
		
		
		netcdfReader.richardsGridFilename = pathGrid;
		
		netcdfReader.read();
		
		
		R1DSolver.z = netcdfReader.z;
		R1DSolver.spaceDeltaZ = netcdfReader.spaceDelta;
		R1DSolver.psiIC = netcdfReader.psiIC;
		R1DSolver.temperature = netcdfReader.temperature;
		R1DSolver.controlVolume = netcdfReader.controlVolume;
		R1DSolver.ks = netcdfReader.Ks;
		R1DSolver.thetaS = netcdfReader.thetaS;
		R1DSolver.thetaR = netcdfReader.thetaR;
		R1DSolver.par1SWRC = netcdfReader.par1SWRC;
		R1DSolver.par2SWRC = netcdfReader.par2SWRC;
		R1DSolver.par3SWRC = netcdfReader.par3SWRC;
		R1DSolver.par4SWRC = netcdfReader.par4SWRC;
		R1DSolver.par5SWRC = netcdfReader.par5SWRC;
		R1DSolver.alphaSpecificStorage = netcdfReader.alphaSS;
		R1DSolver.betaSpecificStorage = netcdfReader.betaSS;
		R1DSolver.inEquationStateID = netcdfReader.equationStateID;
		R1DSolver.inParameterID = netcdfReader.parameterID;
		

		R1DSolver.typeClosureEquation = new String[] {"Water Depth", "Van Genuchten"};
		R1DSolver.typeEquationState = new String[] {"Water Depth", "Van Genuchten"};
		
		R1DSolver.typeUHCModel = new String[] {"", "Mualem Van Genuchten"};
		R1DSolver.typeUHCTemperatureModel = "notemperature"; //"Ronan1998";
		R1DSolver.interfaceHydraulicConductivityModel = "max";
		
		R1DSolver.topBCType = topBC;
		R1DSolver.bottomBCType = bottomBC;
		
		R1DSolver.beta0 = -766.45;
		R1DSolver.referenceTemperatureSWRC = 278.15;
		
		R1DSolver.maxPonding = 0.0;
		
		R1DSolver.tTimeStep = 3600;
		R1DSolver.timeDelta = 1800;
		
		R1DSolver.newtonTolerance = 0.00000000001;//Math.pow(10,-10);
		R1DSolver.nestedNewton = 1;
		R1DSolver.delta = 0;
		
		R1DSolver.picardIteration = 1;

		
		
		buffer.writeFrequency = writeFrequency;
		
		
		
		writeNetCDF.fileName = pathOutput;
		writeNetCDF.briefDescritpion = outputDescription;
		writeNetCDF.pathGrid = pathGrid;
		writeNetCDF.pathBottomBC = pathBottomBC; 
		writeNetCDF.pathTopBC = pathTopBC; 
		writeNetCDF.bottomBC = bottomBC.name();
		writeNetCDF.topBC = topBC.name();
		writeNetCDF.swrcModel = "VG";
		writeNetCDF.soilHydraulicConductivityModel = "Mualem VG no temperature";
		writeNetCDF.interfaceConductivityModel = "max";
		writeNetCDF.writeFrequency = writeFrequency;
		writeNetCDF.spatialCoordinate = netcdfReader.eta;
		writeNetCDF.dualSpatialCoordinate = netcdfReader.etaDual;	
		writeNetCDF.controlVolume = netcdfReader.controlVolume;
		writeNetCDF.psiIC = netcdfReader.psiIC;
		writeNetCDF.temperature = netcdfReader.temperature;
		writeNetCDF.outVariables = new String[] {"darcyVelocity"};
		writeNetCDF.timeUnits = "Minutes since 01/01/1970 00:00:00 UTC";
		writeNetCDF.timeZone = "UTC"; 
		writeNetCDF.fileSizeMax = 10000;
		
		while( topBCReader.doProcess  ) {
			topBCReader.nextRecord();	
			HashMap<Integer, double[]> bCValueMap = topBCReader.outData;
			R1DSolver.inTopBC= bCValueMap;


			bottomBCReader.nextRecord();
			bCValueMap = bottomBCReader.outData;
			R1DSolver.inBottomBC = bCValueMap;

			saveDatesReader.nextRecord();
			bCValueMap = saveDatesReader.outData;
			R1DSolver.inSaveDate = bCValueMap;
			
			R1DSolver.inCurrentDate = topBCReader.tCurrent;
			R1DSolver.solve();

			
			buffer.inputDate = R1DSolver.inCurrentDate;
			buffer.doProcessBuffer = R1DSolver.doProcessBuffer;
			buffer.inputVariable = R1DSolver.outputToBuffer;
			
			buffer.solve();
			

			writeNetCDF.variables = buffer.myVariable;
			writeNetCDF.doProcess = topBCReader.doProcess;
			writeNetCDF.writeNetCDF();

		}

		topBCReader.close();
		bottomBCReader.close();
				
		/*
		 * ASSERT 
		 */
		ReadNetCDFRichardsOutput1D readTestData = new ReadNetCDFRichardsOutput1D();
		readTestData.richardsOutputFilename = getRes("/output/Check_RichardsCoupled_VG.nc");
		readTestData.read();
		
		ReadNetCDFRichardsOutput1D readSimData = new ReadNetCDFRichardsOutput1D();
		readSimData.richardsOutputFilename = pathOutput.replace(".nc","_0000.nc"); // check the proper file of the output series
		readSimData.read();

		for(int k=0; k<readSimData.psi[(readSimData.psi.length)-1].length; k++) {
			double psiSim = readSimData.psi[(readSimData.psi.length)-1][k];
			double psiTest = readTestData.psi[(readTestData.psi.length)-1][k];
			assertEquals("Error in output data check at k = " + k, psiTest, psiSim, 1e-7);
		}
		
	}

//	private OmsTimeSeriesIteratorReader getTimeseriesReader( String inPath, String id, String startDate, String endDate,
//			int timeStepMinutes ) throws URISyntaxException {
//		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
//		reader.file = inPath;
//		reader.idfield = "ID";
//		reader.tStart = startDate;
//		reader.tTimestep = timeStepMinutes;
//		reader.tEnd = endDate;
//		reader.fileNovalue = "-9999";
//		reader.initProcess();
//		return reader;
//	}
}
