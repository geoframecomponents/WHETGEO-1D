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

package org.geoframe.richards;

import java.nio.file.Files;
import java.util.HashMap;

import org.geoframe.whetgeo.WGTestCase;
import org.geoframe.whetgeo1d.boundaryconditions.IBoundaryCondition.RichardsBoundaryConditionType;
import org.geoframe.whetgeo1d.richardssolver.RichardsSolver1DMain;
import org.hortonmachine.gears.io.geoframe.ReadNetCDFRichardsGrid1D;
import org.hortonmachine.gears.io.geoframe.ReadNetCDFRichardsOutput1D;
import org.hortonmachine.gears.io.geoframe.RichardsBuffer1D;
import org.hortonmachine.gears.io.geoframe.WriteNetCDFRichards1DDouble;
import org.hortonmachine.gears.io.timedependent.OmsTimeSeriesIteratorReader;
/**
 * Test the {@link TestBrooksCorey} module.
 * 
 * This test consider an initial hydrostatic condition with Neumann boundary condition at the 
 * top and free drainage at the bottom. 
 * 
 * @author Niccolo' Tubini
 */
public class TestBrooksCorey extends WGTestCase {

	
	public void testBrooksCorey() throws Exception {


		String startDate = "2015-01-15 00:00";
		String endDate = "2015-12-15 00:00";
		int timeStepMinutes = 60;
		String fId = "ID";
				
		String pathTopBC = getRes("/input/TimeSeries/precip.csv");
		String pathBottomBC = getRes("/input/TimeSeries/bottom.csv");
		String pathSaveDates = getRes("/input/TimeSeries/save.csv"); 
		String pathGrid =  getRes("/input/Grid_NetCDF/RichardsCoupled_BC.nc");
		String pathOutput = Files.createTempFile("Sim_RichardsCoupled_BC_", ".nc").toString();
		
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
		ReadNetCDFRichardsGrid1D readNetCDF = new ReadNetCDFRichardsGrid1D();
		
		RichardsSolver1DMain r1DSolver = new RichardsSolver1DMain();
		
		readNetCDF.richardsGridFilename = pathGrid;
		
		readNetCDF.read();
		
		
		r1DSolver.z = readNetCDF.z;
		r1DSolver.spaceDeltaZ = readNetCDF.spaceDelta;
		r1DSolver.psiIC = readNetCDF.psiIC;
		r1DSolver.temperature = readNetCDF.temperature;
		r1DSolver.controlVolume = readNetCDF.controlVolume;
		r1DSolver.ks = readNetCDF.Ks;
		r1DSolver.thetaS = readNetCDF.thetaS;
		r1DSolver.thetaR = readNetCDF.thetaR;
		r1DSolver.par1SWRC = readNetCDF.par1SWRC;
		r1DSolver.par2SWRC = readNetCDF.par2SWRC;
		r1DSolver.par3SWRC = readNetCDF.par3SWRC;
		r1DSolver.par4SWRC = readNetCDF.par4SWRC;
		r1DSolver.par5SWRC = readNetCDF.par5SWRC;
		r1DSolver.alphaSpecificStorage = readNetCDF.alphaSS;
		r1DSolver.betaSpecificStorage = readNetCDF.betaSS;
		r1DSolver.inEquationStateID = readNetCDF.equationStateID;
		r1DSolver.inParameterID = readNetCDF.parameterID;
		r1DSolver.beta0 = -766.45;
		r1DSolver.referenceTemperatureSWRC = 278.15;
		r1DSolver.maxPonding = 0.0;
		r1DSolver.typeClosureEquation = new String[] {"Water Depth", "Brooks Corey"};
		r1DSolver.typeEquationState = new String[] {"Water Depth", "Brooks Corey"};
		r1DSolver.typeUHCModel = new String[] {"", "Mualem Brooks Corey"};
		r1DSolver.typeUHCTemperatureModel = "notemperature"; //"Ronan1998";
		r1DSolver.interfaceHydraulicConductivityModel = "max";
		r1DSolver.topBCType = topBC;
		r1DSolver.bottomBCType = bottomBC;
		r1DSolver.delta = 0;
		r1DSolver.tTimeStep = 3600;
		r1DSolver.timeDelta = 1800;
		r1DSolver.newtonTolerance = 0.00000000001;//Math.pow(10,-10);
		r1DSolver.nestedNewton = 1;
		r1DSolver.picardIteration = 1;

		buffer.writeFrequency = writeFrequency;
		
		writeNetCDF.fileName = pathOutput;
		writeNetCDF.briefDescritpion = outputDescription;
		writeNetCDF.pathGrid = pathGrid;
		writeNetCDF.pathBottomBC = pathBottomBC; 
		writeNetCDF.pathTopBC = pathTopBC; 
		writeNetCDF.bottomBC = bottomBC.name();
		writeNetCDF.topBC = topBC.name();
		writeNetCDF.swrcModel = "BC";
		writeNetCDF.soilHydraulicConductivityModel = "Mualem BC no temperature";
		writeNetCDF.interfaceConductivityModel = "max";
		writeNetCDF.writeFrequency = writeFrequency;
		writeNetCDF.spatialCoordinate = readNetCDF.eta;
		writeNetCDF.dualSpatialCoordinate = readNetCDF.etaDual;	
		writeNetCDF.controlVolume = readNetCDF.controlVolume;
		writeNetCDF.psiIC = readNetCDF.psiIC;
		writeNetCDF.temperature = readNetCDF.temperature;
		writeNetCDF.outVariables = new String[] {"darcy_velocity"};
		writeNetCDF.timeUnits = "Minutes since 01/01/1970 01:00:00 UTC";
		writeNetCDF.timeZone = "UTC"; 
		writeNetCDF.fileSizeMax = 10000;
		
		while( topBCReader.doProcess  ) {
		
			
			topBCReader.nextRecord();	
			HashMap<Integer, double[]> bCValueMap = topBCReader.outData;
			r1DSolver.inTopBC= bCValueMap;


			bottomBCReader.nextRecord();
			bCValueMap = bottomBCReader.outData;
			r1DSolver.inBottomBC = bCValueMap;

			saveDatesReader.nextRecord();
			bCValueMap = saveDatesReader.outData;
			r1DSolver.inSaveDate = bCValueMap;
			
			r1DSolver.inCurrentDate = topBCReader.tCurrent;
			
			r1DSolver.solve();

			
			buffer.inputDate = r1DSolver.inCurrentDate;
			buffer.doProcessBuffer = r1DSolver.doProcessBuffer;
			buffer.inputVariable = r1DSolver.outputToBuffer;
			
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
		System.out.println("Assert");
		ReadNetCDFRichardsOutput1D readTestData = new ReadNetCDFRichardsOutput1D();
		readTestData.richardsOutputFilename = getRes("/Output/Check_RichardsCoupled_BC.nc");
		readTestData.read();
		
		ReadNetCDFRichardsOutput1D readSimData = new ReadNetCDFRichardsOutput1D();
		readSimData.richardsOutputFilename = pathOutput;//.replace(".nc","_0000.nc");
		readSimData.read();

		for(int k=0; k<readSimData.psi[(readSimData.psi.length)-1].length; k++) {
			if(Math.abs(readSimData.psi[(readSimData.psi.length)-1][k]-readTestData.psi[(readTestData.psi.length)-1][k])>Math.pow(10,-11)) {
				System.out.println("\n\n\t\tERROR: psi mismatch");
			}
		}

	}


}
