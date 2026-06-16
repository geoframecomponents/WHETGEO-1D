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

package org.geoframe.heatsolver;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Date;
import java.util.Map;

import org.geoframe.whetgeo.WGTestCase;
import org.geoframe.whetgeo1d.boundaryconditions.IBoundaryCondition.DiffusionBoundaryConditionType;
import org.geoframe.whetgeo1d.heatsolver.HeatDiffusionFreezingThawingSolverWithSurfaceEnergyBalance1DMain;
import org.hortonmachine.gears.io.geoframe.whetgeo.Whetgeo1DInputsHandler;
import org.hortonmachine.gears.io.geoframe.whetgeo.Whetgeo1DOutputsHandler;
import org.hortonmachine.gears.utils.time.ETimeUtilities;

/**
 * Test the heat diffusion freezing-thawing solver with surface energy balance
 * using GeoPackage inputs and outputs via {@link Whetgeo1DInputsHandler} /
 * {@link Whetgeo1DOutputsHandler}.
 *
 * @author Niccolo' Tubini, Andrea Antonello (https://g-ant.eu)
 */
public class TestHeatDiffusionFreezingThawingSurfaceEnergyBalanceGpkg extends WGTestCase {

	public void testHeatDiffusionFreezingThawingSurfaceEnergyBalance() throws Exception {

		String startDate = "2003-01-01 00:00";
		String endDate = "2007-01-01 00:00";

		String inputsPath = "/home/hydrologis/TMP/UNITN/whetgeo1d/TestHeatDiffusionFreezingThawingSurfaceEnergyBalance/HeatDiffusionSurfaceEnergyBalance.gpkg";
		String outputsPath ="/home/hydrologis/TMP/UNITN/whetgeo1d/TestHeatDiffusionFreezingThawingSurfaceEnergyBalance/HeatDiffusionFreezingThawingSurfaceEnergyBalance_output.gpkg";
		Files.deleteIfExists(Path.of(outputsPath));

		var bottomBC = DiffusionBoundaryConditionType.BOTTOM_NEUMANN;

		var inputsHandler = new Whetgeo1DInputsHandler(inputsPath);
		inputsHandler.read();

		HeatDiffusionFreezingThawingSolverWithSurfaceEnergyBalance1DMain solver = new HeatDiffusionFreezingThawingSolverWithSurfaceEnergyBalance1DMain();

		solver.z = inputsHandler.z;
		solver.spaceDeltaZ = inputsHandler.spaceDelta;
		solver.psiIC = inputsHandler.psi;
		solver.temperature = inputsHandler.temperatureIC;
		solver.controlVolume = inputsHandler.controlVolume;
		solver.soilParticlesDensity = inputsHandler.soilParticlesDensity;
		solver.thermalConductivitySoilParticles = inputsHandler.soilParticlesThermalConductivity;
		solver.specificThermalCapacitySoilParticles = inputsHandler.soilParticlesSpecificHeatCapacity;
		solver.meltingTemperature = inputsHandler.meltingTemperature;
		solver.ks = inputsHandler.Ks;
		solver.thetaS = inputsHandler.thetaS;
		solver.thetaR = inputsHandler.thetaR;
		solver.par1SWRC = inputsHandler.par1SWRC;
		solver.par2SWRC = inputsHandler.par2SWRC;
		solver.par3SWRC = inputsHandler.par3SWRC;
		solver.par4SWRC = inputsHandler.par4SWRC;
		solver.par5SWRC = inputsHandler.par5SWRC;
		solver.alphaSpecificStorage = inputsHandler.alphaSS;
		solver.betaSpecificStorage = inputsHandler.betaSS;
		solver.surfaceAlbedo = 0.1;
		solver.surfaceEmissivity = 0.95;
		solver.referenceHeight = 10;
		solver.surfaceRoughness = 0.01;
		solver.surfaceZeroHeightDisplacement = 0.0;
		solver.inEquationStateID = inputsHandler.equationStateID;
		solver.inParameterID = inputsHandler.parameterID;
		solver.typeClosureEquation = new String[] { "VanGenuchtenDallAmico" };
		solver.typeEquationState = new String[] { "FreezingSoilInternalEnergy" };
		solver.typeThermalConductivity = new String[] { "Cosenza" };
		solver.interfaceThermalConductivityModel = "max";
		solver.bottomBCType = bottomBC;
		solver.surfaceAlbedoType = "Constant";
		solver.surfaceEmissivityType = "Constant";
		solver.surfaceAereodynamicResistanceType = "NeutralCondition";
		solver.surfaceWaterVaporResistanceType = "Feddes";
		solver.h1 = 0.1;
		solver.h2 = -5;
		solver.h3 = -15.0;
		solver.h4 = -50.0;
		solver.stationID = 135;
		solver.delta = 0;
		solver.tTimeStep = 3600;
		solver.timeDelta = 3600;
		solver.newtonTolerance = 0.003337000000000;
		solver.nestedNewton = 1;
		solver.picardIteration = 1;

		try (var airTIter = inputsHandler.iterateTimeseries(Whetgeo1DInputsHandler.TABLE_TIMESERIES_AIR_T, startDate, endDate, 1000);
				var swIter = inputsHandler.iterateTimeseries(Whetgeo1DInputsHandler.TABLE_TIMESERIES_SW_RADIATION, startDate, endDate, 1000);
				var lwIter = inputsHandler.iterateTimeseries(Whetgeo1DInputsHandler.TABLE_TIMESERIES_LW_DOWNWELLING, startDate, endDate, 1000);
				var leIter = inputsHandler.iterateTimeseries(Whetgeo1DInputsHandler.TABLE_TIMESERIES_LATENT_HEAT, startDate, endDate, 1000);
				var windIter = inputsHandler.iterateTimeseries(Whetgeo1DInputsHandler.TABLE_TIMESERIES_WIND_VELOCITY, startDate, endDate, 1000);
				var bottomBCIter = inputsHandler.iterateTimeseries(Whetgeo1DInputsHandler.TABLE_TIMESERIES_NO_FLUX, startDate, endDate, 1000);
				var writer = new Whetgeo1DOutputsHandler(outputsPath, 500)) {

			writer.writeIntervalMinutes = 60 * 6; // write every 6 hours

			writer.eta = inputsHandler.eta;
			writer.etaDual = inputsHandler.etaDual;
			writer.controlVolume = inputsHandler.controlVolume;
			writer.psi = inputsHandler.psi;
			writer.temperatureIC = inputsHandler.temperatureIC;

			int iterCount = 0;
			while (airTIter.next()) {
				swIter.next();
				lwIter.next();
				leIter.next();
				windIter.next();
				bottomBCIter.next();

				long timestamp = airTIter.timestamp();
				solver.inAirT = Map.of(solver.stationID, airTIter.values());
				solver.inShortWave = Map.of(solver.stationID, swIter.values());
				solver.inLongWave = Map.of(solver.stationID, lwIter.values());
				solver.inPotentialLatentHeatFlux = Map.of(solver.stationID, leIter.values());
				solver.inWindVelocity = Map.of(solver.stationID, windIter.values());
				solver.inBottomBC = Map.of(solver.stationID, bottomBCIter.values());
				solver.inCurrentDate = ETimeUtilities.INSTANCE.TIME_FORMATTER_UTC.format(new Date(timestamp));
				iterCount++;
				if (iterCount % 100 == 0) {
					System.out.println(iterCount + ") Solving for timestamp: " + solver.inCurrentDate);
				}
				solver.solve();

				writer.timestamp = timestamp;
				writer.temperature = solver.outTemperature;
				writer.theta = solver.outTheta;
				writer.iceContent = solver.outIceContent;
				writer.internalEnergy = solver.outInternalEnergy;
				writer.heatFlux = solver.outConductionHeatFlux;
				writer.error = solver.outErrorInternalEnergy;
				writer.airT = solver.outAirT;
				writer.shortWaveOut = solver.outShortWaveOut;
				writer.shortWaveIn = solver.outShortWaveIn;
				writer.longWaveOut = solver.outLongWaveOut;
				writer.longWaveIn = solver.outLongWaveIn;
				writer.sensibleHeatFlux = solver.outSensibleHeatFlux;
				writer.latentHeatFlux = solver.outActualLatentHeatFlux;
				writer.heatFluxBottom = solver.outHeatFluxBottom;
				writer.write();
			}
		} // all iterators and writer closed here
	}
}
