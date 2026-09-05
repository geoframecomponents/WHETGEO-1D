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
import org.geoframe.whetgeo1d.io.Whetgeo1DInputsHandler;
import org.geoframe.whetgeo1d.io.Whetgeo1DOutputsHandler;
import org.geoframe.whetgeo1d.solvers.RichardsSolverWithRootWaterUptake1D;
import org.hortonmachine.gears.libs.monitor.LogProgressMonitor;
import org.hortonmachine.gears.utils.time.ETimeUtilities;

/**
 * {@code RichardsSolverWithRootWaterUptake1D}'s own native test - not routed
 * through any downstream project's coupling (find that in GEOSPACE-1D's
 * {@code TestGeospacePriestleyTaylorGpkg} covers). Same real 3-layer Van
 * Genuchten soil column and real rainfall top BC {@link TestVanGenuchtenGpkg}
 * already uses, extended with {@code thetaWP}/{@code thetaFC} (see
 * {@code tools.BuildVanGenuchtenRootWaterUptakeFixture} - those two values
 * are synthetic placeholders, not real site data, since this is a new
 * solver-mechanics test with no historical scenario to match).
 *
 * <p>
 * {@code stressedETs} here is a simple, self-limiting per-step demand - not a
 * live GEOET coupling, but the same linear stress-factor shape
 * {@code JarvisStressFactorSolverWithNetRadiation}/{@code LinearWaterStressFactor}
 * use ({@code (theta-thetaWp)/(thetaFc-thetaWp)}, clamped to {@code [0,1]}),
 * inlined here so this test stays WHETGEO-1D-native with no GEOET dependency.
 *
 * <p>
 * Applied to the middle layer (parameterID 2, cells 50-99), not the top one:
 * this fixture's real initial condition (borrowed as-is from
 * {@code RichardsCoupled_VG.gpkg}, {@code psi0=-3m} at the surface) starts the
 * top layer already ~0.08 below this fixture's own synthetic {@code thetaWP}
 * - a mismatch between an invented threshold and a borrowed IC, not a solver
 * issue, but it means the top layer can't meaningfully test "does uptake stay
 * above thetaWP" since it starts below it regardless of uptake. The middle
 * layer starts with a healthy ~0.2 margin above its own {@code thetaWP},
 * giving the clamp real room to be exercised.
 *
 * @author Niccolo' Tubini
 * @author Andrea Antonello
 */
public class TestVanGenuchtenRootWaterUptakeGpkg extends WGTestCase {

	public void testVanGenuchtenRootWaterUptake() throws Exception {

		String startDate = "2015-01-15 00:00";
		String endDate = "2015-02-15 00:00";

		String inputsPath = getRes("/input/gpkg/RichardsCoupled_VG_RootUptake.gpkg");

		var inputsHandler = new Whetgeo1DInputsHandler(inputsPath);
		inputsHandler.read();
		int KMAX = inputsHandler.KMAX;
		// the last cell is the "Water Depth" pond pseudo control-volume (see
		// TestVanGenuchtenGpkg's javadoc); excluded from root uptake and from the
		// thetaR..thetaS/thetaWP physical-bounds checks below
		int topSoilIdx = KMAX - 2;

		var topBC = RichardsBoundaryConditionType.TOP_COUPLED;
		var bottomBC = RichardsBoundaryConditionType.BOTTOM_FREE_DRAINAGE;

		RichardsSolverWithRootWaterUptake1D solver = new RichardsSolverWithRootWaterUptake1D();
		solver.z = inputsHandler.z;
		solver.spaceDeltaZ = inputsHandler.spaceDelta;
		solver.psiIC = inputsHandler.psi;
		solver.temperature = inputsHandler.temperatureIC;
		solver.controlVolume = inputsHandler.controlVolume;
		solver.ks = inputsHandler.Ks;
		solver.thetaS = inputsHandler.thetaS;
		solver.thetaR = inputsHandler.thetaR;
		solver.thetaWP = inputsHandler.thetaWP;
		solver.thetaFC = inputsHandler.thetaFC;
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

		// base (unstressed) demand: 0.02 mm/step, requested from cells 80-99 (part of
		// the middle layer, which starts with a healthy margin above its own
		// thetaWP - see the class javadoc), scaled down by each cell's own linear
		// stress factor - tapers to 0 right at thetaWP, same as the real
		// Jarvis-driven feedback would
		int demandRangeStart = 80;
		int demandRangeEnd = 99;
		double baseDemandMm = 0.02;
		double[] stressedETs = new double[KMAX]; // all zero for the first solve() call (see below)

		double maxAbsWaterVolumeError = 0.0;
		double maxThetaBoundsViolation = 0.0;
		double maxThetaWpViolation = 0.0;

		String outputPath = getTmpPath("VanGenuchtenRootWaterUptakeOutput", ".gpkg");
		var progressMonitor = new LogProgressMonitor("TestVanGenuchtenRootWaterUptakeGpkg");
		try (inputsHandler;
				var topBCIterator = inputsHandler.iterateTimeseries("timeseries_topBC", startDate, endDate, 1000);
				var bottomBCIterator = inputsHandler.iterateTimeseries("timeseries_bottomBC", startDate, endDate,
						1000);
				var writer = new Whetgeo1DOutputsHandler(outputPath, 500)) {
			progressMonitor.beginTask(" -> Running test", -1);

			writer.eta = inputsHandler.eta;
			writer.etaDual = inputsHandler.etaDual;
			writer.controlVolume = inputsHandler.controlVolume;
			writer.psi = inputsHandler.psi;
			writer.temperatureIC = inputsHandler.temperatureIC;
			writer.parameterID = inputsHandler.parameterID;
			writer.swrcThetaS = inputsHandler.thetaS;
			writer.swrcThetaR = inputsHandler.thetaR;
			writer.swrcKs = inputsHandler.Ks;
			writer.swrcN = inputsHandler.par1SWRC;
			writer.swrcAlpha = inputsHandler.par2SWRC;
			writer.topBCType = topBC.name();
			writer.bottomBCType = bottomBC.name();

			while (topBCIterator.next() && bottomBCIterator.next()) {
				long timestamp = topBCIterator.timestamp();
				solver.inTopBC = new HashMap<>(Map.of(solver.stationID, topBCIterator.values()));
				solver.inBottomBC = new HashMap<>(Map.of(solver.stationID, bottomBCIterator.values()));
				solver.inCurrentDate = ETimeUtilities.INSTANCE.TIME_FORMATTER_UTC.format(new Date(timestamp));
				solver.stressedETs = stressedETs;

				solver.solve();

				maxAbsWaterVolumeError = Math.max(maxAbsWaterVolumeError, Math.abs(solver.outErrorVolume));

				// thetaR..thetaS is a general physical-validity invariant of the equation of
				// state, checked over every real cell regardless of demand
				for (int i = 0; i <= topSoilIdx; i++) {
					int paramId = inputsHandler.parameterID[i];
					double theta = solver.outWaterContent[i];
					double boundsViolation = Math.max(inputsHandler.thetaR[paramId] - theta,
							theta - inputsHandler.thetaS[paramId]);
					maxThetaBoundsViolation = Math.max(maxThetaBoundsViolation, boundsViolation);
				}

				// thetaWP is only checked where demand is actually applied - it's a
				// synthetic threshold for this fixture, not guaranteed to sit below every
				// cell's real (borrowed) initial condition elsewhere in the profile (see
				// the class javadoc)
				for (int i = demandRangeStart; i <= demandRangeEnd; i++) {
					int paramId = inputsHandler.parameterID[i];
					double theta = solver.outWaterContent[i];
					maxThetaWpViolation = Math.max(maxThetaWpViolation, inputsHandler.thetaWP[paramId] - theta);
				}

				// self-limiting demand for the NEXT solve() call, from THIS step's theta
				// (lagged, same structure as GEOSPACE-1D's coupled loop)
				for (int i = demandRangeStart; i <= demandRangeEnd; i++) {
					int paramId = inputsHandler.parameterID[i];
					double theta = solver.outWaterContent[i];
					double stressFactor = (theta - inputsHandler.thetaWP[paramId])
							/ (inputsHandler.thetaFC[paramId] - inputsHandler.thetaWP[paramId]);
					stressFactor = Math.max(0, Math.min(1, stressFactor));
					stressedETs[i] = baseDemandMm * stressFactor;
				}

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

		assertTrue("water volume balance residual too large: " + maxAbsWaterVolumeError,
				maxAbsWaterVolumeError < 1e-6);
		assertTrue("water content left its physical thetaR..thetaS bounds by " + maxThetaBoundsViolation,
				maxThetaBoundsViolation < 1e-6);
		assertTrue("water content dropped below thetaWP (root-uptake clamp failed) by " + maxThetaWpViolation,
				maxThetaWpViolation < 1e-6);
	}
}
