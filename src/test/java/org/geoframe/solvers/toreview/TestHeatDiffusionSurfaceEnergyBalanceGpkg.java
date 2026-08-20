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

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Date;
import java.util.Map;

import org.geoframe.whetgeo.WGTestCase;
import org.geoframe.whetgeo1d.core.boundaryconditions.IBoundaryCondition.DiffusionBoundaryConditionType;
import org.geoframe.whetgeo1d.solvers.HeatDiffusionSolverWithSurfaceEnergyBalance1D;
import org.hortonmachine.dbs.compat.ADb;
import org.hortonmachine.dbs.compat.EDb;
import org.hortonmachine.dbs.compat.objects.QueryResult;
import org.geoframe.whetgeo1d.io.Whetgeo1DInputsHandler;
import org.geoframe.whetgeo1d.io.Whetgeo1DOutputsHandler;
import org.hortonmachine.gears.utils.time.ETimeUtilities;

/**
 * Test the heat diffusion solver with surface energy balance using GeoPackage
 * inputs and outputs via {@link Whetgeo1DInputsHandler} /
 * {@link Whetgeo1DOutputsHandler}.
 *
 * @author Niccolo' Tubini
 * @author Andrea Antonello (geopackage input/output adaption)
 */
public class TestHeatDiffusionSurfaceEnergyBalanceGpkg extends WGTestCase {

	public void testHeatDiffusionSurfaceEnergyBalance() throws Exception {

		String startDate = "2003-01-01 00:00";
		String endDate = "2004-01-01 00:00";

		String inputsPath = getRes("/input/gpkg/HeatDiffusionSurfaceEnergyBalance.gpkg");
		String outputsPath = getTmpPath("HeatDiffusionSurfaceEnergyBalance_output", "gpkg");
		Files.deleteIfExists(Path.of(outputsPath));

		var bottomBC = DiffusionBoundaryConditionType.BOTTOM_NEUMANN;

		var inputsHandler = new Whetgeo1DInputsHandler(inputsPath);
		inputsHandler.read();

		HeatDiffusionSolverWithSurfaceEnergyBalance1D solver = new HeatDiffusionSolverWithSurfaceEnergyBalance1D();

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
		solver.typeClosureEquation = new String[] { "Van Genuchten" };
		solver.typeEquationState = new String[] { "SoilInternalEnergy" };
		solver.typeThermalConductivity = new String[] { "Cosenza" };
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
		solver.newtonTolerance = Math.pow(10, -5);
		solver.nestedNewton = 1;
		solver.picardIteration = 1;

		try (var airTIter = inputsHandler.iterateTimeseries(Whetgeo1DInputsHandler.TABLE_TIMESERIES_AIR_T, startDate,
				endDate, 1000);
				var swIter = inputsHandler.iterateTimeseries(Whetgeo1DInputsHandler.TABLE_TIMESERIES_SW_RADIATION,
						startDate, endDate, 1000);
				var lwIter = inputsHandler.iterateTimeseries(Whetgeo1DInputsHandler.TABLE_TIMESERIES_LW_DOWNWELLING,
						startDate, endDate, 1000);
				var leIter = inputsHandler.iterateTimeseries(Whetgeo1DInputsHandler.TABLE_TIMESERIES_LATENT_HEAT,
						startDate, endDate, 1000);
				var windIter = inputsHandler.iterateTimeseries(Whetgeo1DInputsHandler.TABLE_TIMESERIES_WIND_VELOCITY,
						startDate, endDate, 1000);
				var bottomBCIter = inputsHandler.iterateTimeseries(Whetgeo1DInputsHandler.TABLE_TIMESERIES_NO_FLUX,
						startDate, endDate, 1000);
				var writer = new Whetgeo1DOutputsHandler(outputsPath, 50)) {

			writer.writeIntervalMinutes = 60 * 24; // write every day

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

			// check average temperature, internal_energy and theta for min and max eta values
			// through time in the output gpkg
			// TODO create a more meaningful testcase
			try (ADb db = EDb.GEOPACKAGE.getDb()) {
				db.open(outputsPath);
				String sql = """
						SELECT eta, avg(temperature), avg(internal_energy), avg(theta)
						FROM output_state
						WHERE eta = (SELECT min(eta) FROM output_state)
						   OR eta = (SELECT max(eta) FROM output_state)
						GROUP BY eta
						order by eta
						""";
				QueryResult result = db.getTableRecordsMapFromRawSql(sql, -1);
				assertEquals(-29.975, ((Number) result.data.get(0)[0]).doubleValue(), 0);
				assertEquals(-0.025, ((Number) result.data.get(1)[0]).doubleValue(), 0);
				assertEquals(285.1499998872765, ((Number) result.data.get(0)[1]).doubleValue(), 0.0001);
				assertEquals(282.9290082162099, ((Number) result.data.get(1)[1]).doubleValue(), 0.0001);
				assertEquals(8225163.982329215, ((Number) result.data.get(0)[2]).doubleValue(), 0.0001);
				assertEquals(7958993.834700684, ((Number) result.data.get(1)[2]).doubleValue(), 0.0001);
				assertEquals(0.38, ((Number) result.data.get(0)[3]).doubleValue(), 0.0001);
				assertEquals(0.38, ((Number) result.data.get(1)[3]).doubleValue(), 0.0001);
			}
		} // all iterators and writer closed here
	}
}
