/*
 * GNU GPL v3 License
 *
 * Copyright 2026 Andrea Antonello
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
import java.util.Map;

import org.geoframe.whetgeo.WGTestCase;
import org.geoframe.whetgeo1d.core.boundaryconditions.IBoundaryCondition;
import org.geoframe.whetgeo1d.solvers.HeatDiffusionSolver1D;
import org.hortonmachine.gears.io.geoframe.whetgeo.Whetgeo1DInputsHandler;
import org.hortonmachine.gears.io.geoframe.whetgeo.Whetgeo1DOutputsHandler;
import org.hortonmachine.gears.libs.monitor.LogProgressMonitor;
import org.hortonmachine.gears.utils.time.ETimeUtilities;

/**
 * Scenario: a 1 m homogeneous, saturated soil column starting from
 * 10 degC, with a fixed temperature at the top (Dirichlet, 20 degC) and an
 * insulated bottom (Neumann, prescribed flux = 0, i.e. no heat exchange with
 * the deep subsurface). Inputs come from
 * {@code /input/gpkg/HeatDiffusionInsulatedBottomEquilibration.gpkg}.
 *
 * <p>
 * With the bottom sealed, the column cannot lose energy anywhere: physically
 * it must eventually warm to a uniform 20 degC everywhere, and the bottom
 * boundary flux must stay at exactly the prescribed 0 the whole time.
 * Simulated for 2 months (hourly steps), long enough to fully equilibrate. 
 * Therefore:
 * <ol>
 * <li>the bottom-interface diffusive flux must be exactly 0 at every step
 * (this is a direct check of the Neumann wiring itself);</li>
 * <li>the magnitude of the top-interface flux must never increase over time
 * </li>
 * <li>by the end of the run every cell's temperature must have reached the
 * boundary value of 20 degC.</li>
 * </ol>
 *
 * @author Andrea Antonello
 */
public class TestHeatDiffusionInsulatedBottomEquilibration extends WGTestCase {

	public void testInsulatedBottomEquilibratesToTopTemperature() throws Exception {

		String startDate = "2000-01-01 00:00";
		String endDate = "2000-04-30 00:00"; // 2 months

		double topTemperatureCelsius = 20.0;

		String inputsPath = getRes("/input/gpkg/HeatDiffusionInsulatedBottomEquilibration.gpkg");

		var inputsHandler = new Whetgeo1DInputsHandler(inputsPath);
		inputsHandler.read();
		int KMAX = inputsHandler.KMAX;

		var topBC = IBoundaryCondition.DiffusionBoundaryConditionType.TOP_DIRICHLET;
		var bottomBC = IBoundaryCondition.DiffusionBoundaryConditionType.BOTTOM_NEUMANN;

		HeatDiffusionSolver1D solver = new HeatDiffusionSolver1D();
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

		double maxAbsEnergyError = 0.0;
		double maxAbsBottomFlux = 0.0;
		double previousAbsTopFlux = Double.POSITIVE_INFINITY;
		boolean topFluxMagnitudeNeverIncreases = true;

		String outputPath = getTmpPath("HeatDiffusionInsulatedBottomEquilibrationOutput", ".gpkg");
		var progressMonitor = new LogProgressMonitor("TestHeatDiffusionInsulatedBottomEquilibration");
		try (var topBCIterator = inputsHandler.iterateTimeseries("temperature_top_interface", startDate, endDate,
				1000);
				var bottomBCIterator = inputsHandler.iterateTimeseries("temperature_bottom_interface", startDate,
						endDate, 1000);
				var writer = new Whetgeo1DOutputsHandler(outputPath, 500)) {
			progressMonitor.beginTask(" -> Running test", -1);

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
			writer.topBCType = topBC.name();
			writer.bottomBCType = bottomBC.name();

			while (topBCIterator.next() && bottomBCIterator.next()) {
				long timestamp = topBCIterator.timestamp();
				solver.inTopBC = Map.of(solver.stationID, topBCIterator.values());
				solver.inBottomBC = Map.of(solver.stationID, bottomBCIterator.values());
				solver.inCurrentDate = ETimeUtilities.INSTANCE.TIME_FORMATTER_UTC.format(new Date(timestamp));

				solver.solve();

				maxAbsEnergyError = Math.max(maxAbsEnergyError, Math.abs(solver.outErrorInternalEnergy));
				maxAbsBottomFlux = Math.max(maxAbsBottomFlux, Math.abs(solver.outDiffusionHeatFlux[0]));

				double absTopFlux = Math.abs(solver.outDiffusionHeatFlux[KMAX]);
				if (absTopFlux > previousAbsTopFlux + 1e-9) {
					topFluxMagnitudeNeverIncreases = false;
				}
				previousAbsTopFlux = absTopFlux;

				writer.timestamp = timestamp;
				writer.temperature = solver.outTemperature;
				// outTheta is left unset by this solver/equation-state combo (plain
				// SoilInternalEnergy has no liquid-water-content computation), so
				// writing it would just be a column of meaningless constant zeros
				writer.internalEnergy = solver.outInternalEnergy;
				writer.heatFlux = solver.outDiffusionHeatFlux;
				writer.error = solver.outErrorInternalEnergy;
				writer.topBC = solver.outHeatFluxTop;
				writer.bottomBC = solver.outHeatFluxBottom;
				writer.write();
			}
			progressMonitor.done();
		}

		// check that the energy error is small enough to be considered negligible
		assertTrue("energy balance residual too large: " + maxAbsEnergyError, maxAbsEnergyError < 1e-3);

		// 1) the insulated bottom must carry exactly zero flux throughout: this is a
		//    direct check that the BOTTOM_NEUMANN wiring applies the prescribed 0 flux
		assertEquals("bottom (insulated) boundary is not carrying zero flux", 0.0, maxAbsBottomFlux, 1e-6);

		// 2) heat can only flow into a column that starts colder than its one active
		//    boundary and cannot escape anywhere else: the top-interface flux magnitude
		//    must never increase over time
		assertTrue("top boundary flux magnitude increased at some point during the run",
				topFluxMagnitudeNeverIncreases);

		// 3) with no way to lose energy, the whole column must have equilibrated to
		//    the top boundary temperature by the end of the 120-day run
		double targetKelvin = topTemperatureCelsius + 273.15;
		for (int i = 0; i < KMAX; i++) {
			assertEquals("cell " + i + " (eta=" + inputsHandler.eta[i] + ") did not equilibrate to top temperature",
					targetKelvin, solver.outTemperature[i], 1e-3);
		}
	}
}
