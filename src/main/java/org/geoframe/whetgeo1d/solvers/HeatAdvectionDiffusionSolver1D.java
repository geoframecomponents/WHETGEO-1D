/*
 * GNU GPL v3 License
 *
 * Copyright 2019 Niccolo` Tubini
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

package org.geoframe.whetgeo1d.solvers;

import java.util.ArrayList;
import java.util.Map;
import java.util.List;

import org.geoframe.closureequation.closureequation.Parameters;
import org.geoframe.closureequation.equationstate.EquationState;
import org.geoframe.whetgeo1d.core.boundaryconditions.IBoundaryCondition;
import org.geoframe.whetgeo1d.core.boundaryconditions.IBoundaryCondition.DiffusionBoundaryConditionType;
import org.geoframe.whetgeo1d.core.boundaryconditions.IBoundaryCondition.RichardsBoundaryConditionType;
import org.geoframe.whetgeo1d.core.derivedquantities.ComputeQuantitiesHeatAdvectionDiffusion;
import org.geoframe.whetgeo1d.core.derivedquantities.ComputeQuantitiesRichards;
import org.geoframe.whetgeo1d.core.finitevolume.AdvectionDiffusion1DKernel;
import org.geoframe.whetgeo1d.core.finitevolume.Richards1DKernel;
import org.geoframe.whetgeo1d.core.state.WGGeometry;
import org.geoframe.whetgeo1d.core.state.WGProblemQuantities;
import org.hortonmachine.gears.libs.modules.HMModel;

import oms3.annotations.Author;
import oms3.annotations.Bibliography;
import oms3.annotations.Description;
import oms3.annotations.Documentation;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Keywords;
import oms3.annotations.License;
import oms3.annotations.Out;
import oms3.annotations.Unit;

@Description("Solve the heat advection diffusion equation in the conservative form for the 1D domain.")
@Documentation("")
@Author(name = "Niccolo' Tubini, and Riccardo Rigon", contact = "tubini.niccolo@gmail.com")
@Keywords("Hydrology, Richards, Infiltration")
@Bibliography("Casulli (2010)")
//@Label()
//@Name()
//@Status()
@License("General Public License Version 3 (GPLv3)")
public class HeatAdvectionDiffusionSolver1D extends HMModel {

	/*
	 * WATER THERMAL PROPERTIES
	 */

	@Description("Water density. Default value 1000.0 [kg m-3].")
	@In
	@Unit("kg m-3")
	public double waterDensity = 1000.0;

	@Description("Ice density. Default value 920.0 [kg m-3].")
	@In
	@Unit("kg m-3")
	public double iceDensity = 920.0;

	@Description("Specific thermal capacity of water. Default value 4188.0 [J kg-1 K-1].")
	@In
	@Unit("J kg-1 K-1")
	public double specificThermalCapacityWater = 4188.0;

	@Description("Specific thermal capacity of ice. Default value 2117.0 [J kg-1 K-1].")
	@In
	@Unit("J kg-1 K-1")
	public double specificThermalCapacityIce = 2117.0;

	@Description("Thermal conductivity of water. Default value 0.6 [W m-1 K-1].")
	@In
	@Unit("W m-1 K-1")
	public double thermalConductivityWater = 0.6;

	@Description("Thermal conductivity of ice. Default value 2.29 [W m-1 K-1].")
	@In
	@Unit("W m-1 K-1")
	public double thermalConductivityIce = 2.29;

	@Description("Latent heat of fusion. Default value 333700 [J kg-1].")
	@In
	@Unit("J kg-1")
	public double latentHeatFusion = 333700;

	/*
	 * SOIL PARAMETERS
	 */
	@Description("The hydraulic conductivity at saturation")
	@In
	@Unit("m/s")
	public double[] ks;

	@Description("Saturated water content")
	@In
	@Unit("-")
	public double[] thetaS;

	@Description("Residual water content")
	@In
	@Unit("-")
	public double[] thetaR;

	@Description("First parameter of SWRC")
	@In
	@Unit("-")
	public double[] par1SWRC;

	@Description("Second parameter of SWRC")
	@In
	@Unit("-")
	public double[] par2SWRC;

	@Description("Third parameter of SWRC")
	@In
	@Unit("-")
	public double[] par3SWRC;

	@Description("Fourth parameter of SWRC")
	@In
	@Unit("-")
	public double[] par4SWRC;

	@Description("Fifth parameter of SWRC")
	@In
	@Unit("-")
	public double[] par5SWRC;

	@Description("Aquitard compressibility")
	@In
	@Unit("1/Pa")
	public double[] alphaSpecificStorage;

	@Description("Water compressibility")
	@In
	@Unit("1/Pa")
	public double[] betaSpecificStorage;

	@Description("Coefficient for water suction dependence on temperature")
	@In
	@Unit("K")
	public double beta0 = -776.45;

	@Description("Reference temperature for soil water content")
	@In
	@Unit("K")
	public double referenceTemperatureSWRC = 278.15;

	@Description("Reference temperature to compute internal energy")
	@In
	@Unit("K")
	public double referenceTemperatureInternalEnergy = 273.15;

	@Description("Soil particles density")
	@In
	@Unit("kg m-3")
	public double[] soilParticlesDensity;

	@Description("Specific thermal capacity of soil particles")
	@In
	@Unit("J kg-1 K-1")
	public double[] specificThermalCapacitySoilParticles;

	@Description("Thermal conductivity of soil particles")
	@In
	@Unit("W m-1 K-1")
	public double[] thermalConductivitySoilParticles;

	@Description("Melting temperature")
	@In
	@Unit("K")
	public double[] meltingTemperature;

	@Description("Control volume label defining the equation state")
	@In
	@Unit("-")
	public int[] inEquationStateID;

	@Description("Control volume label defining the set of the paramters")
	@In
	@Unit("-")
	public int[] inParameterID;

	/*
	 * MODELS - closure equation - conductivity model - interface conductivity model
	 */

	// Richards equation

	@Description("It is possibile to chose between 3 different models to compute "
			+ "the soil hydraulic properties: Van Genuchten; Brooks and Corey; Kosugi unimodal")
	@In
	public String[] typeClosureEquation;

	@Description("It is possibile to chose between 3 different models to compute "
			+ "the soil hydraulic properties: Van Genuchten; Brooks and Corey; Kosugi unimodal")
	@In
	public String[] typeRichardsEquationState;

	@Description("It is possible to choose among these models:" + "Mualem Van Genuchten, Mualem Brooks Corey, ....")
	@In
	public String[] typeUHCModel;

	@Description("It is possible to choose among these models:" + "notemperature, ....")
	@In
	public String typeUHCTemperatureModel;

	@Description("Hydraulic conductivity at control volume interface can be evaluated as"
			+ " the average of kappas[i] and kappas[i+1]" + " the maximum between kappas[i] and kappas[i+1]"
			+ " the minimum between kappas[i] and kappas[i+1]"
			+ " a weighted average of kappas[i] and kappas[i+1] where weights are dx[i] and dx[i+1]")
	@In
	public String interfaceHydraulicConductivityModel;

	// Heat equation
	@Description("Equation state")
	@In
	public String[] typeInternalEnergyEquationState;

	@Description("Thermal conductivity models")
	@In
	public String[] typeThermalConductivity;

	@Description("Thermal conductivity at control volume interface can be evaluated as"
			+ " the average of kappas[i] and kappas[i+1]" + " the maximum between kappas[i] and kappas[i+1]"
			+ " the minimum between kappas[i] and kappas[i+1]"
			+ " a weighted average of kappas[i] and kappas[i+1] where weights are dx[i] and dx[i+1]")
	@In
	public String interfaceThermalConductivityModel;

	/*
	 * INITIAL CONDITION
	 */
	@Description("Initial condition for water suction read from grid NetCDF file")
	@In
	@Unit("m")
	public double[] psiIC;

	@Description("Initial condition for temperature read from grid NetCDF file")
	@In
	@Unit("K")
	public double[] temperature;

	/*
	 * GEOMETRY
	 */
	@Description("z coordinate read from grid NetCDF file")
	@In
	@Unit("m")
	public double[] z;

	@Description("Space delta to compute gradients read from grid NetCDF file")
	@In
	@Unit("m")
	public double[] spaceDeltaZ;

	@Description("Length of control volumes read from grid NetCDF file")
	@In
	@Unit("m")
	public double[] controlVolume;

	@Description("Maximum ponding depth")
	@In
	@Unit("m")
	public double maxPonding;

	/*
	 * TIME STEP
	 */
	@Description("Time amount at every time-loop")
	@In
	@Unit("s")
	public double tTimeStep;

	@Description("Time step of integration")
	@In
	@Unit("s")
	public double timeDelta;

	/*
	 * ITERATION PARAMETERS
	 */
	@Description("Tolerance for Newton iteration")
	@In
	public double newtonTolerance;

	@Description("Control parameter for nested Newton algorithm:" + "0 --> simple Newton method"
			+ "1 --> nested Newton method")
	@In
	public int nestedNewton;

	@Description("Damped factor for Newton algorithm")
	@In
	public double delta = 0.0;

	@Description("Number of Picard iteration to update the diffusive flux matrix")
	@In
	public int picardIteration = 1;

	/*
	 * BOUNDARY CONDITIONS
	 */
	@Description("The station ID in the timeseries file")
	@In
	@Unit("-")
	public int stationID;

	@Description("The HashMap with the time series of the boundary condition at the top of soil column")
	@In
	@Unit("m")
	public Map<Integer, double[]> inRichardsTopBC;

	@Description("It is possibile to chose between 2 different kind "
			+ "of boundary condition at the top of the domain: " + "- Dirichlet boundary condition --> Top Dirichlet"
			+ "- Neumann boundary condition --> Top Neumann")
	@In
	public RichardsBoundaryConditionType topRichardsBCType;

	@Description("The Map with the time series of the boundary condition at the bottom of soil column")
	@In
	@Unit("m")
	public Map<Integer, double[]> inRichardsBottomBC;

	@Description("It is possibile to chose among 2 different kind "
			+ "of boundary condition at the bottom of the domain: "
			+ "- Dirichlet boundary condition --> Bottom Dirichlet" + "- Neumann boundary condition --> Bottom Neumann")
	@In
	public RichardsBoundaryConditionType bottomRichardsBCType;

	@Description("The Map with the time series of the boundary condition at the top of soil column")
	@In
	@Unit("m")
	public Map<Integer, double[]> inInternalEnergyTopBC;

	@Description("It is possibile to chose between 2 different kind "
			+ "of boundary condition at the top of the domain: " + "- Dirichlet boundary condition --> Top Dirichlet"
			+ "- Neumann boundary condition --> Top Neumann")
	@In
	public DiffusionBoundaryConditionType topInternalEnergyBCType;

	@Description("The Map with the time series of the boundary condition at the bottom of soil column")
	@In
	@Unit("")
	public Map<Integer, double[]> inInternalEnergyBottomBC;

	@Description("It is possibile to chose among 2 different kind "
			+ "of boundary condition at the bottom of the domain: "
			+ "- Dirichlet boundary condition --> Bottom Dirichlet" + "- Neumann boundary condition --> Bottom Neumann")
	@In
	public DiffusionBoundaryConditionType bottomInternalEnergyBCType;

	@Description("")
	@In
	@Unit("")
	public Map<Integer, double[]> inSaveDate;

	@Description("The current date of the simulation.")
	@In
	@Out
	public String inCurrentDate;

	/*
	 * OUTPUT
	 */
	
	@Description("Water suction")
	@Out
	public double[] outWaterSuctions;
	
	@Description("Temperature")
	@Out
	public double[] outTemperatures;
	
	@Description("Water content")
	@Out
	public double[] outWaterContent;
	
	@Description("Heat flux")
	@Out
	public double[] outHeatFlux;
	
	@Description("Darcy velocity")
	@Out
	public double[] outDarcyVelocity;
	
	@Description("Error in internal energy")
	@Out
	public double outErrorInternalEnergy;
	
	@Description("Error in water volume")
	@Out
	public double outErrorVolume;
	
	@Description("ArrayList of variable to be stored in the buffer writer")
	@Out
	public ArrayList<double[]> outputToBuffer;

	@Description("Control variable")
	@Out
	public boolean doProcessBuffer;

	//////////////////////////////////////////
	//////////////////////////////////////////

	@Description("Maximun number of Newton iterations")
	private final int MAXITER_NEWT = 50;

	@Description("Number of control volume for domain discetrization")
	@Unit(" ")
	private int KMAX;

	@Description("It is needed to iterate on the date")
	private int step;

	@Description("Control value to save output:" + "- 1 save the current time step output" + "- 0 do not save")
	private double saveDate;

	private Richards1DKernel richardsSolver;
	private AdvectionDiffusion1DKernel advectionDiffusionSolver;
	private WGProblemQuantities variables;
	private WGGeometry geometry;
	private Parameters parameters;
	private ComputeQuantitiesRichards computeQuantitiesRichards;
	private ComputeQuantitiesHeatAdvectionDiffusion computeQuantitiesHeatAdvectionDiffusion;
	private IBoundaryCondition topRichardsBoundaryCondition;
	private IBoundaryCondition bottomRichardsBoundaryCondition;
	private IBoundaryCondition topInternalEnergyBoundaryCondition;
	private IBoundaryCondition bottomInternalEnergyBoundaryCondition;

	@Execute
	public void solve() {

		if (step == 0) {
			KMAX = psiIC.length;

			variables = new WGProblemQuantities(psiIC, temperature, inEquationStateID, inParameterID);
			geometry = new WGGeometry(z, spaceDeltaZ, controlVolume);
			parameters = new Parameters(waterDensity, iceDensity, specificThermalCapacityWater,
					specificThermalCapacityIce, thermalConductivityWater, thermalConductivityIce, latentHeatFusion,
					referenceTemperatureInternalEnergy, referenceTemperatureSWRC, beta0, thetaS, thetaR,
					soilParticlesDensity, specificThermalCapacitySoilParticles, thermalConductivitySoilParticles,
					meltingTemperature, par1SWRC, par2SWRC, par3SWRC, par4SWRC, par5SWRC, ks, alphaSpecificStorage,
					betaSpecificStorage);

			computeQuantitiesRichards = new ComputeQuantitiesRichards(typeClosureEquation, typeRichardsEquationState,
					typeUHCModel, typeUHCTemperatureModel, interfaceHydraulicConductivityModel, topRichardsBCType,
					bottomRichardsBCType, variables, geometry, parameters);

			computeQuantitiesHeatAdvectionDiffusion = new ComputeQuantitiesHeatAdvectionDiffusion(typeClosureEquation,
					typeInternalEnergyEquationState, typeThermalConductivity, interfaceThermalConductivityModel,
					topInternalEnergyBCType, bottomInternalEnergyBCType, variables, geometry, parameters);

			outputToBuffer = new ArrayList<double[]>();

			List<EquationState> richardsEquationState = computeQuantitiesRichards.getRichardsStateEquation();

			topRichardsBoundaryCondition = IBoundaryCondition.createRichardsSimpleBoundaryCondition(topRichardsBCType);
			bottomRichardsBoundaryCondition = IBoundaryCondition
					.createRichardsSimpleBoundaryCondition(bottomRichardsBCType);

			richardsSolver = new Richards1DKernel(topRichardsBoundaryCondition,
					bottomRichardsBoundaryCondition, KMAX, nestedNewton, newtonTolerance, delta, MAXITER_NEWT,
					richardsEquationState);

			topInternalEnergyBoundaryCondition = IBoundaryCondition
					.createDiffusionBoundaryCondition(topInternalEnergyBCType);
			bottomInternalEnergyBoundaryCondition = IBoundaryCondition
					.createDiffusionBoundaryCondition(bottomInternalEnergyBCType);

			advectionDiffusionSolver = new AdvectionDiffusion1DKernel(topInternalEnergyBoundaryCondition,
					bottomInternalEnergyBoundaryCondition, KMAX);

		} // close step==0

		doProcessBuffer = false;
		variables.richardsTopBCValue = 0.0;
		if (topRichardsBCType == RichardsBoundaryConditionType.TOP_NEUMANN
				|| topRichardsBCType == RichardsBoundaryConditionType.TOP_COUPLED) {
			variables.richardsTopBCValue = (inRichardsTopBC.get(stationID)[0] / 1000) / tTimeStep;
		} else {
			variables.richardsTopBCValue = inRichardsTopBC.get(stationID)[0] / 1000;
		}

		variables.richardsBottomBCValue = 0.0;
		variables.richardsBottomBCValue = inRichardsBottomBC.get(stationID)[0];

		variables.internalEnergyTopBCValue = 0.0;
		if (topInternalEnergyBCType == DiffusionBoundaryConditionType.TOP_NEUMANN) {
			variables.internalEnergyTopBCValue = inInternalEnergyTopBC.get(stationID)[0] / tTimeStep;
		} else {
			variables.internalEnergyTopBCValue = inInternalEnergyTopBC.get(stationID)[0] + 273.15;
		}

		variables.internalEnergyBottomBCValue = 0.0;
		if (bottomInternalEnergyBCType == DiffusionBoundaryConditionType.BOTTOM_NEUMANN) {
			variables.internalEnergyBottomBCValue = inInternalEnergyBottomBC.get(stationID)[0] / tTimeStep;
		} else {
			variables.internalEnergyBottomBCValue = inInternalEnergyBottomBC.get(stationID)[0] + 273.15;
		}

		if(inSaveDate != null) {
			saveDate = -1.0;
			saveDate = inSaveDate.get(stationID)[0];
		}
		outputToBuffer.clear();

		double sumTimeDelta = 0;

		while (sumTimeDelta < tTimeStep) {

			if (sumTimeDelta + timeDelta > tTimeStep) {
				timeDelta = tTimeStep - sumTimeDelta;
			}
			sumTimeDelta = sumTimeDelta + timeDelta;

			/*
			 * Compute water volumes
			 */
			computeQuantitiesRichards.computeWaterVolume(KMAX);

			/*
			 * Compute heat capacity
			 */
			computeQuantitiesHeatAdvectionDiffusion.computeHeatCapacity(KMAX);
			computeQuantitiesHeatAdvectionDiffusion.computeInternalEnergy(KMAX);

			/*
			 * Compute thermal conductivity
			 */
			computeQuantitiesHeatAdvectionDiffusion.computeThermalConductivity(KMAX);
			computeQuantitiesHeatAdvectionDiffusion.computeInterfaceThermalConductivity(KMAX);
			variables.lambdasInterface[KMAX] = 0.6;

			computeQuantitiesHeatAdvectionDiffusion.computeTransportedQuantity(KMAX);

			/*
			 * Compute xStar
			 */
			computeQuantitiesRichards.computeXStar(KMAX);

			/*
			 * Solve Richards equation
			 */
			for (int picard = 0; picard < picardIteration; picard++) {

				/*
				 * Compute hydraulic conductivity
				 * 
				 */
				computeQuantitiesRichards.computeHydraulicConductivity(KMAX);

				computeQuantitiesRichards.computeInterfaceHydraulicConductivity(KMAX);

				/*
				 * Solve PDE
				 */
				variables.waterSuctions = richardsSolver.solve(timeDelta, variables.richardsBottomBCValue,
						variables.richardsTopBCValue, KMAX, variables.kappasInterface, variables.volumes,
						geometry.spaceDeltaZ, variables.ETs, variables.waterSuctions, variables.temperatures,
						variables.parameterID, variables.equationStateID);

			} // close Picard iteration

			/*
			 * compute run-off
			 */
			computeQuantitiesRichards.computeRunOff(KMAX, maxPonding);

			/*
			 * Compute - water volume and total water volume - water content
			 */
			computeQuantitiesRichards.computeWaterVolumeNew(KMAX);
			computeQuantitiesRichards.computeThetasNew(KMAX);

			/*
			 * Fluxes
			 */
			computeQuantitiesRichards.computeDarcyVelocities(KMAX);

			/*
			 * compute error
			 */
			computeQuantitiesRichards.computeError(KMAX, timeDelta);

			/*
			 * New heat capacity
			 */
			computeQuantitiesHeatAdvectionDiffusion.computeHeatCapacityNew(KMAX);

			/*
			 * Solve heat advection-diffusion equation
			 */

			for (int k = 0; k < KMAX; k++) {
				variables.temperatures[k] = variables.temperatures[k] - 273.15;
			}
			variables.internalEnergyTopBCValue = variables.internalEnergyTopBCValue - 273.15;
			variables.internalEnergyBottomBCValue = variables.internalEnergyBottomBCValue - 273.15;

			if (variables.thetasNew[variables.thetasNew.length - 1] <= 0) {
				KMAX = KMAX - 1;
				variables.internalEnergy -= variables.internalEnergys[variables.thetasNew.length - 1];
			}

			variables.temperatures = advectionDiffusionSolver.solve(timeDelta, variables.internalEnergyBottomBCValue,
					variables.internalEnergyTopBCValue, KMAX, variables.lambdasInterface, variables.heatCapacitysNew,
					variables.heatCapacitys, geometry.spaceDeltaZ, variables.heatSourcesSinksTerm,
					variables.temperatures, variables.waterSuctions, variables.darcyVelocities,
					variables.waterCapacityTransported, variables.parameterID, variables.equationStateID);

			for (int k = 0; k < KMAX; k++) {
				variables.temperatures[k] = variables.temperatures[k] + 273.15;
			}
			variables.internalEnergyTopBCValue = variables.internalEnergyTopBCValue + 273.15;
			variables.internalEnergyBottomBCValue = variables.internalEnergyBottomBCValue + 273.15;

			computeQuantitiesHeatAdvectionDiffusion.computeInternalEnergyNew(KMAX);

			computeQuantitiesHeatAdvectionDiffusion.computeConductionHeatFlux(KMAX);
			computeQuantitiesHeatAdvectionDiffusion.computeAdvectionHeatFlux(KMAX);
			computeQuantitiesHeatAdvectionDiffusion.computeHeatFlux(KMAX);

			/*
			 * Compute error heat equation
			 */
			computeQuantitiesHeatAdvectionDiffusion.computeError(KMAX, timeDelta);

			if (variables.thetasNew[variables.thetasNew.length - 1] <= 0) {
				KMAX = KMAX + 1;
				variables.temperatures[KMAX - 1] = variables.temperatures[KMAX - 2];
			}

		}

		
		outWaterSuctions = variables.waterSuctions;
		outTemperatures = variables.temperatures;
		outWaterContent = variables.thetasNew;
		outHeatFlux = variables.heatFluxs;
		outDarcyVelocity = variables.darcyVelocities;
		outErrorInternalEnergy = variables.errorInternalEnergy;
		outErrorVolume = variables.errorVolume;
		/*
		 * aggiugere output da salvare e modificare il buffer e il writer
		 */
		if (inSaveDate != null && saveDate == 1) {
			outputToBuffer.add(variables.waterSuctions);
			outputToBuffer.add(variables.temperatures);
			outputToBuffer.add(variables.thetasNew);
			outputToBuffer.add(variables.heatFluxs);
			outputToBuffer.add(variables.darcyVelocities);
			outputToBuffer.add(new double[] { variables.errorInternalEnergy });
			outputToBuffer.add(new double[] { variables.errorVolume });
			doProcessBuffer = true;
		}
		step++;

	} //// MAIN CYCLE END ////

} /// CLOSE ///
