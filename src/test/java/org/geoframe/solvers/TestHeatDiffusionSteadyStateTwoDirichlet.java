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
import org.hortonmachine.gears.utils.time.ETimeUtilities;

/**
 * Scenario: a 1 m homogeneous, saturated soil column with a fixed temperature
 * at both ends (Dirichlet top = 20 degC, Dirichlet bottom = 5 degC), starting
 * from a 12 degC. Inputs come from
 * {@code /input/gpkg/HeatDiffusionSteadyStateTwoDirichlet.gpkg}.
 *
 * <p>
 * Simulated for 2 months so the column is expected to have reached
 * conductive steady state well before the run ends. Therefore:
 * <ol>
 * <li>the temperature profile must match the analytical linear steady-state
 * gradient between the two boundary temperatures;</li>
 * <li>the heat flux must be uniform across every interface, since it is steady state.</li>
 * </ol>
 *
 * <p>
 * Note: {@link HeatDiffusionSolver1D} does not move water. The {@code psiIC}
 * (initial pressure head) and SWRC inputs are only used once, at setup, to
 * derive each cell's water content and hence its thermal properties (heat
 * capacity, thermal conductivity); that water content is then held fixed at
 * its initial value for the rest of the run &mdash; not to be confused with
 * ice/freezing physics, which is a separate feature ({@code WithFreezingThawing}
 * solver variants) not used here. Only temperature evolves via
 * {@code solve()}. A run where
 * water actually moves would use {@code HeatAdvectionDiffusionSolver1D},
 * which runs its own {@code Richards1DKernel} every step.
 *
 * @author Andrea Antonello
 */
public class TestHeatDiffusionSteadyStateTwoDirichlet extends WGTestCase {

	public void testSteadyStateLinearProfile() throws Exception {

		String startDate = "2000-01-01 00:00";
		String endDate = "2000-03-01 00:00"; // 2 months, enough to reach steady state

		double topTemperatureCelsius = 20.0;
		double bottomTemperatureCelsius = 5.0;
		double etaTop = 0.0;
		double etaBottom = -1.0;

		String inputsPath = getRes("/input/gpkg/HeatDiffusionSteadyStateTwoDirichlet.gpkg");

		var inputsHandler = new Whetgeo1DInputsHandler(inputsPath);
		inputsHandler.read();
		int KMAX = inputsHandler.KMAX;

		var topBC = IBoundaryCondition.DiffusionBoundaryConditionType.TOP_DIRICHLET;
		var bottomBC = IBoundaryCondition.DiffusionBoundaryConditionType.BOTTOM_DIRICHLET;

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

		try (var topBCIterator = inputsHandler.iterateTimeseries("temperature_top_interface", startDate, endDate,
				1000);
				var bottomBCIterator = inputsHandler.iterateTimeseries("temperature_bottom_interface", startDate,
						endDate, 1000)) {

			while (topBCIterator.next() && bottomBCIterator.next()) {
				long timestamp = topBCIterator.timestamp();
				solver.inTopBC = Map.of(solver.stationID, topBCIterator.values());
				solver.inBottomBC = Map.of(solver.stationID, bottomBCIterator.values());
				solver.inCurrentDate = ETimeUtilities.INSTANCE.TIME_FORMATTER_UTC.format(new Date(timestamp));

				solver.solve();

				maxAbsEnergyError = Math.max(maxAbsEnergyError, Math.abs(solver.outErrorInternalEnergy));
			}
		}
		// check that the energy error is small enough to be considered negligible
		assertTrue("energy error too large: " + maxAbsEnergyError, maxAbsEnergyError < 1e-4);

		// 1) after 2 months the column must match the analytical linear steady-state
		//    profile between the two fixed boundary temperatures
		double topTemperatureKelvin = topTemperatureCelsius + 273.15;
		double bottomTemperatureKelvin = bottomTemperatureCelsius + 273.15;
		double domainLength = etaTop - etaBottom;
		for (int i = 0; i < KMAX; i++) {
			double eta = inputsHandler.eta[i];
			double expectedKelvin = bottomTemperatureKelvin
					+ (topTemperatureKelvin - bottomTemperatureKelvin) * (eta - etaBottom) / domainLength;
			assertEquals("steady-state temperature mismatch at eta=" + eta, expectedKelvin, solver.outTemperature[i],
					1e-4);
		}

		// 2) since we are in a steady state, the heat flux must be uniform across the interfaces.
		double fluxBottom = solver.outDiffusionHeatFlux[0];
		double fluxTop = solver.outDiffusionHeatFlux[KMAX];
		assertEquals("flux not uniform between bottom and top interfaces at steady state", fluxBottom, fluxTop,
				1e-4);
		for (int k = 1; k < KMAX; k++) {
			assertEquals("flux not uniform at internal interface " + k, fluxBottom, solver.outDiffusionHeatFlux[k],
					1e-4);
		}
	}
}
