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
 * Scenario: same 3 m, 3-layer Van Genuchten soil column as
 * {@link TestVanGenuchtenDirichletGpkg}, but driven by an actual rainfall
 * timeseries through a <b>coupled</b> top boundary condition instead of a
 * fixed Dirichlet head, with free drainage at the bottom. Inputs come from
 * {@code /input/gpkg/RichardsCoupled_VG.gpkg}.
 *
 * <p>
 * {@code TOP_COUPLED} reuses the same flux (Neumann-style) equation as
 * {@code TOP_NEUMANN}, but the grid carries one extra pseudo control-volume
 * stacked on top of the 150 real soil cells: a "Water Depth" surface-ponding
 * reservoir (equation state {@code theta(psi) = max(psi, 0)}, i.e. its state
 * variable *is* the ponding depth in metres). This lets rainfall in excess of
 * the soil's infiltration capacity accumulate as a real solved state instead
 * of being discarded as bookkeeping-only runoff, and couples back into the
 * Richards column through the pond cell's own (soil-derived) conductivity at
 * the interface.
 *
 * <p>
 * Unlike the Dirichlet scenario, rainfall here is intermittent, so the
 * column isn't guaranteed to monotonically wet up. Instead this test checks
 * physical invariants that must hold regardless of the exact rainfall
 * sequence: the water balance must close at every step, and the water
 * content of every real soil cell must stay within the physical bounds of
 * its own Van Genuchten curve ({@code thetaR..thetaS}) throughout the run.
 *
 * @author Niccolo' Tubini
 * @author Andrea Antonello
 */
public class TestVanGenuchtenGpkg extends WGTestCase {

	public void testVanGenuchten() throws Exception {

		String startDate = "2015-01-15 00:00";
		String endDate = "2015-12-15 00:00";

		String inputsPath = getRes("/input/gpkg/RichardsCoupled_VG.gpkg");

		var inputsHandler = new Whetgeo1DInputsHandler(inputsPath);
		inputsHandler.read();
		int KMAX = inputsHandler.KMAX;
		// the last cell is the "Water Depth" pseudo control-volume (see class
		// javadoc); its water content isn't a saturation fraction, so it's
		// excluded from the thetaR..thetaS physical-bounds check below.
		int topSoilIdx = KMAX - 2;

		var topBC = RichardsBoundaryConditionType.TOP_COUPLED;
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
		solver.typeClosureEquation = new String[] { "Water Depth", "Van Genuchten" };
		solver.typeEquationState = new String[] { "Water Depth", "Van Genuchten" };
		solver.typeUHCModel = new String[] { "", "Mualem Van Genuchten" };
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
		double maxThetaBoundsViolation = 0.0;

		String outputPath = getTmpPath("VanGenuchtenCoupledOutput", ".gpkg");
		var progressMonitor = new LogProgressMonitor("TestVanGenuchtenGpkg");
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

				for (int i = 0; i <= topSoilIdx; i++) {
					int paramId = inputsHandler.parameterID[i];
					double theta = solver.outWaterContent[i];
					double violation = Math.max(inputsHandler.thetaR[paramId] - theta,
							theta - inputsHandler.thetaS[paramId]);
					// negative violation means theta is within its physical bounds, so we only care about positive
					maxThetaBoundsViolation = Math.max(maxThetaBoundsViolation, violation);
				}

				// RichardsSolver1D doesn't solve temperature (notemperature model): leave
				// writer.temperature unset rather than write the IC array back unchanged
				// at every step, which would just be a flat, non-informative column.
				writer.timestamp = timestamp;
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

		// every real soil cell's water content must stay within its own Van
		// Genuchten curve's physical bounds (thetaR..thetaS) at every step
		assertTrue("water content left its physical thetaR..thetaS bounds by " + maxThetaBoundsViolation,
				maxThetaBoundsViolation < 1e-6);

		// sanity check that the coupled rainfall boundary actually drove the
		// simulation (not a vacuously-passing no-op run)
		double topSoilChange = Math.abs(finalPsi[topSoilIdx] - initialPsi[topSoilIdx]);
		assertTrue("top soil cell's suction did not change at all over the run: delta psi=" + topSoilChange,
				topSoilChange > 1e-6);
	}
}
