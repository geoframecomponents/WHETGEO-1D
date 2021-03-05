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

package richards;

import java.net.URISyntaxException;
import java.util.*;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorWriter;

import it.geoframe.blogspot.buffer.buffertowriter.RichardsBuffer1D;
import it.geoframe.blogspot.whetgeo1d.richardssolver.RichardsSolver1DMain;
import it.geoframe.blogspot.buffer.buffercalibrationpoints.BufferCalibrationRichards1D;
import it.geoframe.blogspot.buffer.bufferparameter.BufferParameterRichards1D;
import it.geoframe.blogpsot.netcdf.monodimensionalproblemtimedependent.ReadNetCDFRichardsGrid1D;
import it.geoframe.blogpsot.netcdf.monodimensionalproblemtimedependent.WriteNetCDFRichards1DDouble;

import org.junit.Test;

/**
 * Test the {@link TestVanGenuchtenBufferParameter} module.
 * 
 * This test consider an initial hydrostatic condition with Neumann boundary condition at the 
 * top and free drainage at the bottom. 
 * 
 * Here the simulated values in the measurament point are saved in a .csv file.
 * 
 * @author Niccolo' Tubini
 */
public class TestVanGenuchtenBufferParameter {

	@Test
	public void Test() throws Exception {


		String startDate = "2015-01-15 00:00";
		String endDate = "2015-12-15 00:00";
		int timeStepMinutes = 60;
		String fId = "ID";
				
		String pathTopBC = "resources/input/TimeSeries/precip.csv";
		String pathBottomBC = "resources/input/TimeSeries/bottom.csv";
		String pathSaveDates = "resources/input/TimeSeries/save.csv"; 
		String pathGrid =  "C:/Users/Niccolo/OMS/OMS_Project_WHETGEO1D/data/Grid_NetCDF/xLUCA2.nc";
		String pathOutput = "resources/output/cancella.nc";
		
		String pathCalibrationPointPsi = "resources/output/cancella_psi.csv";
		String pathCalibrationPointTheta = "resources/output/cancella_theta.csv";
		
		String topBC = "Top Neumann";
		String bottomBC = "Bottom free drainage";

		String outputDescription = "\n"
				+ "Initial condition hydrostatic no ponding\n		"
				+ "DeltaT: 1800s\n		"
				+ "Picard iteration: 1\n		";
		
		int writeFrequency = 1000000;
		
		OmsTimeSeriesIteratorReader topBCReader = getTimeseriesReader(pathTopBC, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader bottomBCReader = getTimeseriesReader(pathBottomBC, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader saveDatesReader = getTimeseriesReader(pathSaveDates, fId, startDate, endDate, timeStepMinutes);
		
		OmsTimeSeriesIteratorWriter writerCalibrationPointsPsi = getTimeseriesWriter(pathCalibrationPointPsi, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorWriter writerCalibrationPointsTheta = getTimeseriesWriter(pathCalibrationPointTheta, fId, startDate, endDate, timeStepMinutes);

		RichardsBuffer1D buffer = new RichardsBuffer1D();
		WriteNetCDFRichards1DDouble writeNetCDF = new WriteNetCDFRichards1DDouble();
		ReadNetCDFRichardsGrid1D readNetCDF = new ReadNetCDFRichardsGrid1D();
	
		BufferCalibrationRichards1D bufferCalibration = new BufferCalibrationRichards1D();

		BufferParameterRichards1D bufferParameter = new BufferParameterRichards1D();
		
		RichardsSolver1DMain R1DSolver = new RichardsSolver1DMain();
		
		
		bufferParameter.thetaS1 = 0.43;
		bufferParameter.thetaR1 = 0.045;
		bufferParameter.par1SWRC1 = 2.68;
		bufferParameter.par2SWRC1 = 14.5;
		bufferParameter.par3SWRC1 = 0.0;
		bufferParameter.par4SWRC1 = 0.0;
		bufferParameter.par5SWRC1 = 0.0;
		bufferParameter.alphaSpecificStorage1 = 0.0;
		bufferParameter.betaSpecificStorage1 = 0.0;
		bufferParameter.ks1 = 8.25e-05;
		bufferParameter.solve();
		
		
		readNetCDF.richardsGridFilename = pathGrid;
		
		readNetCDF.read();
		
		
		R1DSolver.z = readNetCDF.z;
		R1DSolver.spaceDeltaZ = readNetCDF.spaceDelta;
		R1DSolver.psiIC = readNetCDF.psiIC;
		R1DSolver.temperature = readNetCDF.temperature;
		R1DSolver.controlVolume = readNetCDF.controlVolume;
//		R1DSolver.ks = readNetCDF.Ks;
//		R1DSolver.thetaS = readNetCDF.thetaS;
//		R1DSolver.thetaR = readNetCDF.thetaR;
//		R1DSolver.par1SWRC = readNetCDF.par1SWRC;
//		R1DSolver.par2SWRC = readNetCDF.par2SWRC;
//		R1DSolver.par3SWRC = readNetCDF.par3SWRC;
//		R1DSolver.par4SWRC = readNetCDF.par4SWRC;
//		R1DSolver.par5SWRC = readNetCDF.par5SWRC;
//		R1DSolver.alphaSpecificStorage = readNetCDF.alphaSS;
//		R1DSolver.betaSpecificStorage = readNetCDF.betaSS;
		
		R1DSolver.ks = bufferParameter.ks;
		R1DSolver.thetaS = bufferParameter.thetaS;
		R1DSolver.thetaR = bufferParameter.thetaR;
		R1DSolver.par1SWRC = bufferParameter.par1SWRC;
		R1DSolver.par2SWRC = bufferParameter.par2SWRC;
		R1DSolver.par3SWRC = bufferParameter.par3SWRC;
		R1DSolver.par4SWRC = bufferParameter.par4SWRC;
		R1DSolver.par5SWRC = bufferParameter.par5SWRC;
		R1DSolver.alphaSpecificStorage = bufferParameter.alphaSpecificStorage;
		R1DSolver.betaSpecificStorage = bufferParameter.betaSpecificStorage;
		
		R1DSolver.inEquationStateID = readNetCDF.equationStateID;
		R1DSolver.inParameterID = readNetCDF.parameterID;
		R1DSolver.beta0 = -766.45;
		R1DSolver.referenceTemperatureSWRC = 278.15;
		R1DSolver.maxPonding = 0.0;
		R1DSolver.typeClosureEquation = new String[] {"Van Genuchten"};
		R1DSolver.typeEquationState = new String[] {"Van Genuchten"};
		R1DSolver.typeUHCModel = new String[] {"Mualem Van Genuchten"};
		R1DSolver.typeUHCTemperatureModel = "notemperature"; //"Ronan1998";
		R1DSolver.interfaceHydraulicConductivityModel = "max";
		R1DSolver.topBCType = topBC;
		R1DSolver.bottomBCType = bottomBC;
		R1DSolver.delta = 0;
		R1DSolver.tTimeStep = 3600;
		R1DSolver.timeDelta = 1800;
		R1DSolver.newtonTolerance = 0.00000000001;//Math.pow(10,-10);
		R1DSolver.nestedNewton = 1;
		R1DSolver.picardIteration = 1;
		
		buffer.writeFrequency = writeFrequency;
		
		bufferCalibration.controlVolumeIndex = readNetCDF.controlVolumeIndex;
		
		writeNetCDF.fileName = pathOutput;
		writeNetCDF.briefDescritpion = outputDescription;
		writeNetCDF.pathGrid = pathGrid;
		writeNetCDF.pathBottomBC = pathBottomBC; 
		writeNetCDF.pathTopBC = pathTopBC; 
		writeNetCDF.bottomBC = bottomBC;
		writeNetCDF.topBC = topBC;
		writeNetCDF.swrcModel = "VG";
		writeNetCDF.soilHydraulicConductivityModel = "Mualem VG no temperature";
		writeNetCDF.interfaceConductivityModel = "max";
		writeNetCDF.writeFrequency = writeFrequency;
		writeNetCDF.spatialCoordinate = readNetCDF.eta;
		writeNetCDF.dualSpatialCoordinate = readNetCDF.etaDual;	
		writeNetCDF.controlVolume = readNetCDF.controlVolume;
		writeNetCDF.psiIC = readNetCDF.psiIC;
		writeNetCDF.temperature = readNetCDF.temperature;
		writeNetCDF.outVariables = new String[] {"darcy_velocity"};
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
			System.out.println(topBCReader.tCurrent);
			R1DSolver.solve();

			bufferCalibration.inputVariable = R1DSolver.outputToBuffer;
			
			bufferCalibration.solve();
			
			buffer.inputDate = R1DSolver.inCurrentDate;
			buffer.doProcessBuffer = R1DSolver.doProcessBuffer;
			buffer.inputVariable = R1DSolver.outputToBuffer;
			
			buffer.solve();
			

			writeNetCDF.variables = buffer.myVariable;
			writeNetCDF.doProcess = topBCReader.doProcess;
			writeNetCDF.writeNetCDF();
			
			writerCalibrationPointsPsi.inData = bufferCalibration.simulatedPsi;
			writerCalibrationPointsTheta.inData = bufferCalibration.simulatedTheta;
			
			writerCalibrationPointsPsi.writeNextLine();
			writerCalibrationPointsTheta.writeNextLine();

		}

		topBCReader.close();
		bottomBCReader.close();
		saveDatesReader.close();
		
		writerCalibrationPointsPsi.close();
		writerCalibrationPointsTheta.close();
				

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
	
	private OmsTimeSeriesIteratorWriter getTimeseriesWriter( String inPath, String id, String startDate, String endDate,
			int timeStepMinutes ) throws URISyntaxException {
		OmsTimeSeriesIteratorWriter writer = new OmsTimeSeriesIteratorWriter();
		writer.file = inPath;
		writer.tStart = startDate;
		writer.tTimestep = timeStepMinutes;
		writer.fileNovalue = "-9999";
		return writer;
	}
}
