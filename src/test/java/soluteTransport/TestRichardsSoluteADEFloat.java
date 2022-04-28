/*
  * GNU GPL v3 License
 *
 * Copyright 2022 Concetta D'Amato
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

package soluteTransport;

import static java.lang.Math.pow;

import java.net.URISyntaxException;
import java.util.*;
import org.hortonmachine.gears.io.timedependent.OmsTimeSeriesIteratorReader;

import it.geoframe.blogspot.buffer.buffertowriter.*;
import it.geoframe.blogspot.whetgeo1d.solutetransport.*;
import it.geoframe.blogpsot.netcdf.monodimensionalproblemtimedependent.*;


import org.junit.Test;

/**
 * Test the {@link TestRichardsSoluteADEFloat} module.
 * 
 * 
 * @author Concetta D'Amato and Niccolo' Tubini
 */
public class TestRichardsSoluteADEFloat {

	@Test
	public void Test() throws Exception {


		String startDate = "2003-01-01 01:00";
		String endDate = "2003-09-10 05:00";
		int timeStepMinutes = 60;
		String fId = "ID";
		String lab = "01Float";
				
		String pathSoluteTopBC = "resources/input/TimeSeries/ConcTopVar2_T0135.csv";
		String pathSoluteBottomBC = "resources/input/TimeSeries/ConcBottom0_T0135.csv";
		String pathRichardsTopBC = "resources/input/TimeSeries/PrecipVar2_T0135.csv";
		String pathRichardsBottomBC = "resources/input/TimeSeries/noFlux_T0135.csv";
		String pathSaveDates = "resources/input/TimeSeries/saveAll_T0135.csv"; 
		String pathGrid =  "data/Grid_NetCDF/Grid_WHETGEO_Richards_Solute_2704.nc";
		String pathOutput = "resources/output/Sim_RichardsADE_test2704_"+lab+".nc";
		
		//Solute boundary conditions
		String topSoluteBC = "Top dirichlet";
		String bottomSoluteBC = "Bottom dirichlet";
		
		//Richards boundary conditions
		String topRichardsBC = "Top Coupled";
		String bottomRichardsBC = "Bottom Free Drainage"; //"Bottom Free drainage"

		String outputDescription = "\n"
				+ "Richards' equation coupled with the solute advection-dispersion equation";
		
		int writeFrequency = 1000000;
		
		OmsTimeSeriesIteratorReader topSoluteBCReader = getTimeseriesReader(pathSoluteTopBC, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader bottomSoluteBCReader = getTimeseriesReader(pathSoluteBottomBC, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader topRichardsBCReader = getTimeseriesReader(pathRichardsTopBC, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader bottomRichardsBCReader = getTimeseriesReader(pathRichardsBottomBC, fId, startDate, endDate, timeStepMinutes);
		OmsTimeSeriesIteratorReader saveDatesReader = getTimeseriesReader(pathSaveDates, fId, startDate, endDate, timeStepMinutes);

		SoluteAdvectionDispersionBuffer1D buffer = new SoluteAdvectionDispersionBuffer1D();
		WriteNetCDFRichardsSoluteADE1DFloat writeNetCDF = new WriteNetCDFRichardsSoluteADE1DFloat();
		ReadNetCDFRichardsSoluteADEGrid1D readNetCDF = new ReadNetCDFRichardsSoluteADEGrid1D();
		
		RichardsConservativeSoluteADESolver1DMain solver = new RichardsConservativeSoluteADESolver1DMain();
		
		//double[] concentrationIC = {0,0,0,0,0,0};
		//solver.concentrationIC = concentrationIC;
		
		//solver.molecularDiffusion = 5.91667 * pow(10,-8); 
		//solver.longitudinalDispersivity = 0.05;
		
		//solver.molecularDiffusion = 0; 
		//solver.longitudinalDispersivity = 0;
		
		
		//solver.longitudinalDispersivity = 1;
		
		
		readNetCDF.richardsGridFilename = pathGrid;
		
		readNetCDF.read();
		
		
		solver.z = readNetCDF.z;
		solver.spaceDeltaZ = readNetCDF.spaceDelta;
		solver.psiIC = readNetCDF.psiIC;
		solver.temperatureIC = readNetCDF.temperature;
		solver.controlVolume = readNetCDF.controlVolume;
		solver.concentrationIC = readNetCDF.concentrationIC;
		//solver.soilParticlesDensity = readNetCDF.soilParticlesDensity;
		//solver.thermalConductivitySoilParticles = readNetCDF.soilParticlesThermalConductivity;
		//solver.specificThermalCapacitySoilParticles = readNetCDF.soilParticlesSpecificHeatCapacity;
		
		solver.molecularDiffusion = readNetCDF.molecularDiffusion;
		solver.longitudinalDispersivity = readNetCDF.longitudinalDispersivity;
		
		solver.ks = readNetCDF.Ks;
		solver.thetaS = readNetCDF.thetaS;
		solver.thetaR = readNetCDF.thetaR;
		solver.par1SWRC = readNetCDF.par1SWRC;
		solver.par2SWRC = readNetCDF.par2SWRC;
		solver.par3SWRC = readNetCDF.par3SWRC;
		solver.par4SWRC = readNetCDF.par4SWRC;
		solver.par5SWRC = readNetCDF.par5SWRC;
		//solver.meltingTemperature = readNetCDF.meltingTemperature;
		solver.alphaSpecificStorage = readNetCDF.alphaSS;
		solver.betaSpecificStorage = readNetCDF.betaSS;
		solver.inEquationStateID = readNetCDF.equationStateID;
		solver.inParameterID = readNetCDF.parameterID;
		solver.typeClosureEquation = new String[] {"Water depth", "Van Genuchten"};
		
		solver.typeRichardsEquationState = new String[] {"Water depth", "Van Genuchten"};
		solver.typeUHCModel = new String[] {"", "Mualem Van Genuchten"};
		solver.interfaceHydraulicConductivityModel = "max";
		solver.typeUHCTemperatureModel = "notemperature";
		
		solver.interfaceDispersionModel = "max";
		solver.maxPonding = 0;
		
		//solver.typeInternalEnergyEquationState = new String[] {"Water heat capacity", "SoilHeatCapacity"};
		//solver.typeThermalConductivity = new String[] {"Water", "Cosenza"};
		//solver.interfaceThermalConductivityModel = "max";
		solver.topRichardsBCType = topRichardsBC;
		solver.bottomRichardsBCType = bottomRichardsBC;
		solver.topSoluteBCType = topSoluteBC;
		solver.bottomSoluteBCType = bottomSoluteBC;
		solver.tTimeStep = 3600;
		solver.timeDelta = 3600;
		solver.newtonTolerance = Math.pow(10,-12);
		solver.nestedNewton = 1;
		solver.picardIteration = 1;
		solver.stationID = 135;

		buffer.writeFrequency = writeFrequency;
		
		writeNetCDF.fileName = pathOutput;
		writeNetCDF.briefDescritpion = outputDescription;
		writeNetCDF.pathGrid = pathGrid;
		writeNetCDF.pathSoluteBottomBC = ""; 
		writeNetCDF.pathSoluteTopBC = ""; 
		writeNetCDF.pathRichardsBottomBC = ""; 
		writeNetCDF.pathRichardsTopBC = ""; 
		writeNetCDF.bottomSoluteBC = "";
		writeNetCDF.topSoluteBC = "";
		writeNetCDF.bottomRichardsBC = "";
		writeNetCDF.topRichardsBC = "";
		writeNetCDF.swrcModel = "VG";
		writeNetCDF.soilHydraulicConductivityModel = "Mualem VG no temperature";
		writeNetCDF.interfaceHydraulicConductivityModel = "max";
		//writeNetCDF.soilThermalConductivityModel = "Cosenza";
		writeNetCDF.interfaceDispersionCoefficientModel = "max";
		writeNetCDF.writeFrequency = writeFrequency;
		writeNetCDF.spatialCoordinate = readNetCDF.eta;
		writeNetCDF.dualSpatialCoordinate = readNetCDF.etaDual;	
		writeNetCDF.controlVolume = readNetCDF.controlVolume;
		writeNetCDF.psi = readNetCDF.psiIC;
		writeNetCDF.concentrationIC = readNetCDF.concentrationIC;
		//writeNetCDF.rootIC = readNetCDF.rootIC;;
		writeNetCDF.timeUnits = "Minutes since 01/01/1970 00:00:00 UTC";
		writeNetCDF.timeZone = "UTC"; 
		writeNetCDF.fileSizeMax = 10000;
		
		while( topRichardsBCReader.doProcess  ) {
		
			
			topSoluteBCReader.nextRecord();	
			HashMap<Integer, double[]> bCValueMap = topSoluteBCReader.outData;
			solver.inSoluteTopBC= bCValueMap;


			bottomSoluteBCReader.nextRecord();
			bCValueMap = bottomSoluteBCReader.outData;
			solver.inSoluteBottomBC = bCValueMap;
			
			topRichardsBCReader.nextRecord();	
			bCValueMap = topRichardsBCReader.outData;
			solver.inRichardsTopBC= bCValueMap;


			bottomRichardsBCReader.nextRecord();
			bCValueMap = bottomRichardsBCReader.outData;
			solver.inRichardsBottomBC = bCValueMap;

			saveDatesReader.nextRecord();
			bCValueMap = saveDatesReader.outData;
			solver.inSaveDate = bCValueMap;
			
			solver.inCurrentDate = topRichardsBCReader.tCurrent;
			
			solver.solve();

			
			buffer.inputDate = solver.inCurrentDate;
			buffer.doProcessBuffer = solver.doProcessBuffer;
			buffer.inputVariable = solver.outputToBuffer;
			
			buffer.solve();
			

			writeNetCDF.variables = buffer.myVariable;
			writeNetCDF.doProcess = topRichardsBCReader.doProcess;
			writeNetCDF.writeNetCDF();


		}

		topRichardsBCReader.close();
		bottomRichardsBCReader.close();
		topSoluteBCReader.close();
		bottomSoluteBCReader.close();
						
		/*
		 * ASSERT 
		 */
		/*System.out.println("Assert");
		ReadNetCDFHeatAdvectionDiffusionOutput1D readTestData = new ReadNetCDFHeatAdvectionDiffusionOutput1D();
		readTestData.gridFilename = "resources/Output/Check_heat_advection_diffusion_0000.nc";
		readTestData.read();
		
		ReadNetCDFHeatAdvectionDiffusionOutput1D readSimData = new ReadNetCDFHeatAdvectionDiffusionOutput1D();
		readSimData.gridFilename = pathOutput.replace(".nc","_0000.nc");
		readSimData.read();

		for(int k=0; k<readSimData.temperature[(readSimData.temperature.length)-1].length; k++) {
			if(Math.abs(readSimData.temperature[(readSimData.temperature.length)-1][k]-readTestData.temperature[(readTestData.temperature.length)-1][k])>Math.pow(10,-11)) {
				System.out.println("\n\n\t\tERROR: temperature mismatch");
			}
		}
		
		for(int k=0; k<readSimData.psi[(readSimData.psi.length)-1].length; k++) {
			if(Math.abs(readSimData.psi[(readSimData.psi.length)-1][k]-readTestData.psi[(readTestData.psi.length)-1][k])>Math.pow(10,-11)) {
				System.out.println("\n\n\t\tERROR: psi mismatch");
			}
		}*/

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
