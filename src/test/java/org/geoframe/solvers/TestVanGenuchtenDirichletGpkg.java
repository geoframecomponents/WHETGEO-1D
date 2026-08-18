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

import java.util.Date;
import java.util.HashMap;
import java.util.Map;

import org.geoframe.whetgeo.WGTestCase;
import org.geoframe.whetgeo1d.core.boundaryconditions.IBoundaryCondition.RichardsBoundaryConditionType;
import org.geoframe.whetgeo1d.solvers.RichardsSolver1D;
import org.hortonmachine.gears.io.geoframe.whetgeo.Whetgeo1DInputsHandler;
import org.hortonmachine.gears.io.geoframe.whetgeo.Whetgeo1DOutputsHandler;
import org.hortonmachine.gears.libs.monitor.LogProgressMonitor;
import org.hortonmachine.gears.utils.time.ETimeUtilities;

/**
 * Scenario: a 3 m homogeneous-by-layer (3 x 1 m Van Genuchten layers) soil
 * column starting from a hydrostatic initial condition (water table at the
 * bottom boundary), with the top interface held at saturation (Dirichlet,
 * psi=0) and free drainage at the bottom. Inputs come from
 * {@code /input/gpkg/VanGenuchtenDirichlet.gpkg}.
 * 
 * <p>Note that the boundary conditions are as follows:
 * <ul>
 * <li>top: all 0 values and per Dirichlet condition that means psi always = 0 (saturated). Once the
 * simulation starts, psi at the top will remain 0 constant. That is like a wetting front moving downward,
 * due to the drier soil below. That creates the gradient, driving the flow downward. </li>
 * <li>bottom: free drainage (Neumann condition with zero gradient).</li>
 * </ul>
 *
 * <p>
 * With the top boundary held saturated and the bottom free to drain, the
 * column can only wet up over time relative to its initial dry-at-the-top
 * hydrostatic profile: suction must not increase (get drier) anywhere,
 * and the mass balance must close at every step.
 *
 * @author Andrea Antonello
 * @author Niccolo' Tubini
 */
public class TestVanGenuchtenDirichletGpkg extends WGTestCase {

	public void testVanGenuchtenDirichlet() throws Exception {

		String startDate = "2015-01-15 00:00";
		String endDate = "2015-12-15 00:00";

		String inputsPath = getRes("/input/gpkg/VanGenuchtenDirichlet.gpkg");

		var inputsHandler = new Whetgeo1DInputsHandler(inputsPath);
		inputsHandler.read();
		int KMAX = inputsHandler.KMAX;

		var topBC = RichardsBoundaryConditionType.TOP_DIRICHLET;
		var bottomBC = RichardsBoundaryConditionType.BOTTOM_FREE_DRAINAGE;

		RichardsSolver1D solver = new RichardsSolver1D();
		solver.z = inputsHandler.z;
		solver.spaceDeltaZ = inputsHandler.spaceDelta;
		solver.psiIC = inputsHandler.psi;
		solver.temperature = inputsHandler.temperatureIC;
		solver.controlVolume = inputsHandler.controlVolume;
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
		solver.maxPonding = 0.0;
		solver.typeClosureEquation = new String[] { "Van Genuchten" };
		solver.typeEquationState = new String[] { "Van Genuchten" };
		solver.typeUHCModel = new String[] { "Mualem Van Genuchten" };
		solver.typeUHCTemperatureModel = "notemperature";
		solver.interfaceHydraulicConductivityModel = "max";
		solver.topBCType = topBC;
		solver.bottomBCType = bottomBC;
		solver.delta = 0;
		solver.tTimeStep = 3600;
		solver.timeDelta = 1800;
		solver.newtonTolerance = 0.00000000001;
		solver.nestedNewton = 1;
		solver.picardIteration = 1;
		solver.stationID = 0;

		double[] initialPsi = inputsHandler.psi.clone();
		double[] finalPsi = null;
		double maxAbsWaterVolumeError = 0.0;

		String outputPath = getTmpPath("VanGenuchtenDirichletOutput", ".gpkg");
		var progressMonitor = new LogProgressMonitor("TestVanGenuchtenDirichletGpkg");
		try (var topBCIterator = inputsHandler.iterateTimeseries("timeseries_topBC", startDate, endDate, 1000);
				var bottomBCIterator = inputsHandler.iterateTimeseries("timeseries_bottomBC", startDate, endDate,
						1000);
				var writer = new Whetgeo1DOutputsHandler(outputPath, 500)) {
			progressMonitor.beginTask(" -> Running test", -1);

			writer.writeIntervalMinutes = 60 * 24; // write once a day, keep the output manageable
			writer.eta = inputsHandler.eta;
			writer.etaDual = inputsHandler.etaDual;
			writer.controlVolume = inputsHandler.controlVolume;
			writer.psi = inputsHandler.psi;
			writer.temperatureIC = inputsHandler.temperatureIC;
			// snapshot of the per-layer SWRC parameters, so the output gpkg is
			// self-contained (e.g. for chart annotations) without needing the input gpkg
			writer.parameterID = inputsHandler.parameterID;
			writer.swrcThetaS = inputsHandler.thetaS;
			writer.swrcThetaR = inputsHandler.thetaR;
			writer.swrcKs = inputsHandler.Ks;
			writer.swrcN = inputsHandler.par1SWRC;
			writer.swrcAlpha = inputsHandler.par2SWRC;
			// so a chart of the output can tell a fixed head apart from a flux or a
			// coupled forcing, not just see the raw numeric value
			writer.topBCType = topBC.name();
			writer.bottomBCType = bottomBC.name();

			while (topBCIterator.next() && bottomBCIterator.next()) {
				long timestamp = topBCIterator.timestamp();
				solver.inTopBC = new HashMap<>(Map.of(solver.stationID, topBCIterator.values()));
				solver.inBottomBC = new HashMap<>(Map.of(solver.stationID, bottomBCIterator.values()));
				solver.inCurrentDate = ETimeUtilities.INSTANCE.TIME_FORMATTER_UTC.format(new Date(timestamp));

				solver.solve();

				maxAbsWaterVolumeError = Math.max(maxAbsWaterVolumeError, Math.abs(solver.outErrorVolume));
				finalPsi = solver.outWaterSuctions.clone();

				// RichardsSolver1D doesn't solve temperature (notemperature model): the
				// input array is written back unchanged at every step.
				writer.timestamp = timestamp;
				writer.temperature = solver.temperature;
				writer.theta = solver.outWaterContent;
				writer.waterSuction = solver.outWaterSuctions;
				writer.darcyVelocity = solver.outDarcyVelocity;
				writer.errorVolume = solver.outErrorVolume;
				writer.topBC = solver.outTopBCValue;
				writer.bottomBC = solver.outBottomBCValue;
				writer.write();
			}
			progressMonitor.done();
		}

		// the water balance residual must stay negligible at every step
		assertTrue("water volume balance residual too large: " + maxAbsWaterVolumeError,
				maxAbsWaterVolumeError < 1e-6);

		// with the top interface held saturated (psi=0) and free drainage at the
		// bottom, the column can only wet up relative to the initial dry-at-the-top
		// hydrostatic profile: suction must not have increased anywhere
		for (int i = 0; i < KMAX; i++) {
			assertTrue("cell " + i + " (eta=" + inputsHandler.eta[i] + ") got drier instead of wetting up: initial psi="
					+ initialPsi[i] + " final psi=" + finalPsi[i], finalPsi[i] >= initialPsi[i] - 1e-9);
		}

		// the top-most cell, closest to the saturated Dirichlet boundary, must have
		// wet up the most in absolute terms compared to the bottom-most cell
		double topWetting = finalPsi[KMAX - 1] - initialPsi[KMAX - 1];
		assertTrue("top cell should have wet up over the course of the run: delta psi=" + topWetting,
				topWetting > 0.0);
	}
}
