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
import org.geoframe.whetgeo1d.boundaryconditions.IBoundaryCondition;
import org.geoframe.whetgeo1d.heatsolver.HeatDiffusionSolver1DMain;
import org.hortonmachine.dbs.compat.ADb;
import org.hortonmachine.dbs.compat.EDb;
import org.hortonmachine.dbs.compat.objects.QueryResult;
// import org.hortonmachine.dbs.compat.ADb;
// import org.hortonmachine.dbs.compat.EDb;
// import org.hortonmachine.gears.io.geoframe.HeatDiffusionBuffer1D;
// import org.hortonmachine.gears.io.geoframe.ReadNetCDFHeatDiffusionGrid1D;
// import org.hortonmachine.gears.io.geoframe.WriteNetCDFHeatDiffusion1DDouble;
// import org.hortonmachine.gears.io.geoframe.whetgeo.WhetgeoTemperatureIterator;
import org.hortonmachine.gears.io.geoframe.whetgeo.Whetgeo1DInputsHandler;
import org.hortonmachine.gears.io.geoframe.whetgeo.Whetgeo1DOutputsHandler;
// import org.hortonmachine.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import org.hortonmachine.gears.utils.time.ETimeUtilities;

/**
 * Test the HeatDiffusion module using GeoPackage inputs and outputs.
 *
 * @author Niccolo' Tubini
 */
public class TestHeatDiffusionGpkg extends WGTestCase {

	public void testHeatDiffusion() throws Exception {

		String startDate = "2013-12-15 01:00";
		String endDate = "2015-12-16 01:00";
		String inputsPath = getRes("/input/gpkg/HeatDiffusion.gpkg");
		String outputsPath = getTmpPath("HeatDiffusion_output", "gpkg");   
		Files.deleteIfExists(Path.of(outputsPath));

		var inputsHandler = new Whetgeo1DInputsHandler(inputsPath);
		inputsHandler.read();

		var topBC = IBoundaryCondition.DiffusionBoundaryConditionType.TOP_DIRICHLET;
		var bottomBC = IBoundaryCondition.DiffusionBoundaryConditionType.BOTTOM_NEUMANN;

		HeatDiffusionSolver1DMain solver = new HeatDiffusionSolver1DMain();

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
		solver.inEquationStateID = inputsHandler.equationStateID;
		solver.inParameterID = inputsHandler.parameterID;
		solver.beta0 = -766.45;
		solver.referenceTemperatureSWRC = 278.15;
		solver.typeClosureEquation = new String[] { "Van Genuchten" };
		solver.typeEquationState = new String[] { "SoilInternalEnergy" };
		solver.typeThermalConductivity = new String[] { "Cosenza" };
		solver.interfaceThermalConductivityModel = "max";
		solver.topBCType = topBC;
		solver.bottomBCType = bottomBC;
		solver.delta = 0;
		solver.tTimeStep = 3600;
		solver.timeDelta = 3600;
		solver.newtonTolerance = Math.pow(10, -6);
		solver.nestedNewton = 1;
		solver.picardIteration = 1;
		solver.stationID = 0;

		try (var topBCIterator = inputsHandler.iterateTimeseries("timeseries_temperature_bottom_interface", startDate,
				endDate, 1000);
				var bottomBCIterator = inputsHandler.iterateTimeseries("timeseries_temperature_top_interface",
						startDate, endDate, 1000);
				var writer = new Whetgeo1DOutputsHandler(outputsPath, 500)) {

			writer.writeIntervalMinutes = 60 * 24; // write every day

			writer.eta = inputsHandler.eta;
			writer.etaDual = inputsHandler.etaDual;
			writer.controlVolume = inputsHandler.controlVolume;
			writer.psi = inputsHandler.psi;
			writer.temperatureIC = inputsHandler.temperatureIC;

			int iterCount = 0;
			while (topBCIterator.next() && bottomBCIterator.next()) {

				long timestamp = topBCIterator.timestamp();
				solver.inTopBC = Map.of(solver.stationID, topBCIterator.values());

				timestamp = bottomBCIterator.timestamp();
				solver.inBottomBC = Map.of(solver.stationID, bottomBCIterator.values());

				solver.inCurrentDate = ETimeUtilities.INSTANCE.TIME_FORMATTER_UTC.format(new Date(timestamp));
				iterCount++;
				if (iterCount % 1000 == 0) {
					System.out.println(iterCount + ") Solving for timestamp: " + solver.inCurrentDate);
				}

				solver.solve();

				writer.timestamp = timestamp;
				writer.temperature = solver.outTemperature;
				writer.theta = solver.outTheta;
				writer.internalEnergy = solver.outInternalEnergy;
				writer.heatFlux = solver.outDiffusionHeatFlux;
				writer.error = solver.outErrorInternalEnergy;
				writer.topBC = solver.outHeatFluxTop;
				writer.bottomBC = solver.outHeatFluxBottom;
				writer.write();

			}
		} // iterators and writer all closed here; writer flushes remaining buffer
		
		
		// check average temperature and internal energy for min and max eta values
		// through time in the output gpkg
		// TODO create a more meaningful testcase
		ADb db = EDb.GEOPACKAGE.getDb();
		db.open(outputsPath);	
		String sql ="""
				SELECT eta, avg(temperature), avg(internal_energy)
				FROM output_state
				WHERE eta = (SELECT min(eta) FROM output_state)
				   OR eta = (SELECT max(eta) FROM output_state)
				GROUP BY eta order by eta
				""";
		QueryResult result = db.getTableRecordsMapFromRawSql(sql, -1);
		assertEquals(-29.975, ((Number)result.data.get(0)[0]).doubleValue(), 0);
		assertEquals(-0.025, ((Number)result.data.get(1)[0]).doubleValue(), 0);
		assertEquals(1068.7785001837876, ((Number)result.data.get(0)[1]).doubleValue(), 0.0001);
		assertEquals(285.1552027413406, ((Number)result.data.get(1)[1]).doubleValue(), 0.0001);
		assertEquals(1.3129997785765049E8, ((Number)result.data.get(0)[2]).doubleValue(), 0.0001);
		assertEquals(8226205.610392944, ((Number)result.data.get(1)[2]).doubleValue(), 0.0001);
		
	}
}
