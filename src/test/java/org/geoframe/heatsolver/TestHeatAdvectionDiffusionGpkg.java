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
import org.geoframe.whetgeo1d.boundaryconditions.IBoundaryCondition.RichardsBoundaryConditionType;
import org.geoframe.whetgeo1d.heatsolver.HeatAdvectionDiffusionSolver1DMain;
import org.hortonmachine.gears.io.geoframe.whetgeo.Whetgeo1DInputsHandler;
import org.hortonmachine.gears.io.geoframe.whetgeo.Whetgeo1DOutputsHandler;
import org.hortonmachine.gears.utils.time.ETimeUtilities;

/**
 * Test the {@link TestHeatAdvectionDiffusionGpkg} module.
 * 
 * 
 * @author Niccolo' Tubini
 */
public class TestHeatAdvectionDiffusionGpkg extends WGTestCase {

	public void testHeatAdvectionDiffusion() throws Exception {

		String startDate = "2003-01-01 00:00";
		String endDate = "2007-01-01 00:00";
		String folder = "/home/hydrologis/TMP/UNITN/whetgeo1d/TestHeatAdvectionDiffusion/";
		String inputsPath = folder + "HeatAdvectionDiffusion.gpkg";
		String outputsPath = folder + "HeatAdvectionDiffusion_output.gpkg";
		Files.deleteIfExists(Path.of(outputsPath));

		var topInternalEnergyBC = DiffusionBoundaryConditionType.TOP_DIRICHLET;
		var bottomInternalEnergyBC = DiffusionBoundaryConditionType.BOTTOM_DIRICHLET;
		var topRichardsBC = RichardsBoundaryConditionType.TOP_COUPLED;
		var bottomRichardsBC = RichardsBoundaryConditionType.BOTTOM_FREE_DRAINAGE;

		var precipTable = "timeseries_Precip";
		var airTTable = "timeseries_airT";
		var bottomTTable = "timeseries_bottomT_T0135";

		var inputsHandler = new Whetgeo1DInputsHandler(inputsPath);
		inputsHandler.read();

		HeatAdvectionDiffusionSolver1DMain solver = new HeatAdvectionDiffusionSolver1DMain();

		solver.z = inputsHandler.z;
		solver.spaceDeltaZ = inputsHandler.spaceDelta;
		solver.psiIC = inputsHandler.psi;
		solver.temperature = inputsHandler.temperatureIC;
		solver.controlVolume = inputsHandler.controlVolume;
		solver.soilParticlesDensity = inputsHandler.soilParticlesDensity;
		solver.thermalConductivitySoilParticles = inputsHandler.soilParticlesThermalConductivity;
		solver.specificThermalCapacitySoilParticles = inputsHandler.soilParticlesSpecificHeatCapacity;
		solver.ks = inputsHandler.Ks;
		solver.thetaS = inputsHandler.thetaS;
		solver.thetaR = inputsHandler.thetaR;
		solver.par1SWRC = inputsHandler.par1SWRC;
		solver.par2SWRC = inputsHandler.par2SWRC;
		solver.par3SWRC = inputsHandler.par3SWRC;
		solver.par4SWRC = inputsHandler.par4SWRC;
		solver.par5SWRC = inputsHandler.par5SWRC;
		solver.meltingTemperature = inputsHandler.meltingTemperature;
		solver.alphaSpecificStorage = inputsHandler.alphaSS;
		solver.betaSpecificStorage = inputsHandler.betaSS;
		solver.inEquationStateID = inputsHandler.equationStateID;
		solver.inParameterID = inputsHandler.parameterID;
		solver.typeClosureEquation = new String[] { "Water depth", "Van Genuchten" };

		solver.typeRichardsEquationState = new String[] { "Water depth", "Van Genuchten" };
		solver.typeUHCModel = new String[] { "", "Mualem Van Genuchten" };
		solver.interfaceHydraulicConductivityModel = "max";
		solver.typeUHCTemperatureModel = "notemperature";

		solver.typeInternalEnergyEquationState = new String[] { "Water heat capacity", "SoilHeatCapacity" };
		solver.typeThermalConductivity = new String[] { "Water", "Cosenza" };
		solver.interfaceThermalConductivityModel = "max";
		solver.topRichardsBCType = topRichardsBC;
		solver.bottomRichardsBCType = bottomRichardsBC;
		solver.topInternalEnergyBCType = topInternalEnergyBC;
		solver.bottomInternalEnergyBCType = bottomInternalEnergyBC;
		solver.tTimeStep = 3600;
		solver.timeDelta = 3600;
		solver.newtonTolerance = Math.pow(10, -12);
		solver.nestedNewton = 1;
		solver.picardIteration = 1;
		solver.stationID = 135;

		try (var topInternalEnergyBCIter = inputsHandler.iterateTimeseries(airTTable, startDate, endDate, 1000);
				var bottomInternalEnergyBCIter = inputsHandler.iterateTimeseries(bottomTTable, startDate, endDate,
						1000);
				var top_bot_RichardsBCIter = inputsHandler.iterateTimeseries(precipTable, startDate, endDate, 1000);
				var writer = new Whetgeo1DOutputsHandler(outputsPath, 500)) {

			writer.writeIntervalMinutes = 60 * 6; // write every 6 hours

			writer.eta = inputsHandler.eta;
			writer.etaDual = inputsHandler.etaDual;
			writer.controlVolume = inputsHandler.controlVolume;
			writer.psi = inputsHandler.psi;
			writer.temperatureIC = inputsHandler.temperatureIC;

			int iterCount = 0;
			while (topInternalEnergyBCIter.next()) {
				bottomInternalEnergyBCIter.next();
				top_bot_RichardsBCIter.next();

				long timestamp = topInternalEnergyBCIter.timestamp();
				solver.inInternalEnergyTopBC = Map.of(solver.stationID, topInternalEnergyBCIter.values());
				solver.inInternalEnergyBottomBC = Map.of(solver.stationID, bottomInternalEnergyBCIter.values());
				solver.inRichardsTopBC = Map.of(solver.stationID, top_bot_RichardsBCIter.values());
				solver.inRichardsBottomBC = Map.of(solver.stationID, top_bot_RichardsBCIter.values());
				solver.inCurrentDate = ETimeUtilities.INSTANCE.TIME_FORMATTER_UTC.format(new Date(timestamp));
				iterCount++;
				if (iterCount % 100 == 0) {
					System.out.println(iterCount + ") Solving for timestamp: " + solver.inCurrentDate);
				}
				solver.solve();

				writer.timestamp = timestamp;
				writer.temperature = solver.outTemperatures;
				writer.theta = solver.outWaterContent;
				writer.heatFlux = solver.outHeatFlux;
				writer.error = solver.outErrorInternalEnergy;
				writer.waterSuction = solver.outWaterSuctions;
				writer.darcyVelocity = solver.outDarcyVelocity;
				writer.errorVolume = solver.outErrorVolume;
				writer.write();
			}

		}
	}

}
