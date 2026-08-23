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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.geoframe.closureequation.closureequation.Parameters;
import org.geoframe.closureequation.equationstate.EquationState;
import org.geoframe.whetgeo1d.core.boundaryconditions.IBoundaryCondition;
import org.geoframe.whetgeo1d.core.boundaryconditions.IBoundaryCondition.DiffusionBoundaryConditionType;
import org.geoframe.whetgeo1d.core.derivedquantities.ComputeQuantitiesInternalEnergy;
import org.geoframe.whetgeo1d.core.finitevolume.Diffusion1DKernel;
import org.geoframe.whetgeo1d.utils.GFGeometry;
import org.geoframe.whetgeo1d.utils.ProblemQuantities;
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

@Description("Solve the pure heat diffusion equation in the conservative form for the 1D domain.")
@Documentation("")
@Author(name = "Niccolo' Tubini, and Riccardo Rigon", contact = "tubini.niccolo@gmail.com")
@Keywords("Hydrology, heat diffusion")
@Bibliography("Casulli (2010)")
//@Label()
//@Name()
//@Status()
@License("General Public License Version 3 (GPLv3)")
public class HeatDiffusionSolver1D extends HMModel {

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
	@Description("Closure equation models")
	@In
	public String[] typeClosureEquation;

	@Description("Equation state")
	@In
	public String[] typeEquationState;

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
	public Map<Integer, double[]> inTopBC;

	@Description("It is possibile to chose between 2 different kind "
			+ "of boundary condition at the top of the domain: "
			+ "- Dirichlet boundary condition --> IBoundaryCondition.DiffusionBoundaryConditionType.TOP_DIRICHLET"
			+ "- Neumann boundary condition --> IBoundaryCondition.DiffusionBoundaryConditionType.TOP_NEUMANN")
	@In
	public DiffusionBoundaryConditionType topBCType;

	@Description("The HashMap with the time series of the boundary condition at the bottom of soil column")
	@In
	@Unit("m")
	public Map<Integer, double[]> inBottomBC;

	@Description("It is possibile to chose among 2 different kind "
			+ "of boundary condition at the bottom of the domain: "
			+ "- Dirichlet boundary condition --> Bottom Dirichlet" + "- Neumann boundary condition --> Bottom Neumann")
	@In
	public DiffusionBoundaryConditionType bottomBCType;

	@Description("The current date of the simulation.")
	@In
	@Out
	public String inCurrentDate;

	/*
	 * OUTPUT
	 */
	
	@Description("Temperature at the current time step")
	@Unit("K")
	@Out
	public double[] outTemperature;
	
	@Description("Water content at the current time step")
	@Unit("-")
	@Out
	public double[] outTheta;
	
	@Description("Internal energy at the current time step")
	@Unit("J kg-1")
	@Out
	public double[] outInternalEnergy;
	
	@Description("Diffusion heat flux at the current time step")
	@Unit("W m-2")
	@Out
	public double[] outDiffusionHeatFlux;
	
	@Description("Error on internal energy at the current time step")
	@Unit("J kg-1")
	@Out
	public double outErrorInternalEnergy;
	
	@Description("Heat flux at the top of the domain at the current time step")
	@Unit("W m-2")	
	@Out
	public double outHeatFluxTop;
	
	@Description("Heat flux at the bottom of the domain at the current time step")
	@Unit("W m-2")
	@Out
	public double outHeatFluxBottom;
	

	@Description("ArrayList of variable to be stored in the buffer writer")
	@Out
	public List<double[]> outputToBuffer;


	//////////////////////////////////////////
	//////////////////////////////////////////

	@Description("Maximun number of Newton iterations")
	private final int MAXITER_NEWT = 50;

	@Description("Number of control volume for domain discetrization")
	@Unit(" ")
	private int KMAX;

	@Description("It is needed to iterate on the date")
	private int step;

	private Diffusion1DKernel diffusionSolver;
	private ProblemQuantities variables;
	private GFGeometry geometry;
	private Parameters parameters;
	private ComputeQuantitiesInternalEnergy computeQuantitiesInternalEnergy;
	private IBoundaryCondition topBoundaryCondition;
	private IBoundaryCondition bottomBoundaryCondition;

	@Execute
	public void solve() {

		if (step == 0) {
			KMAX = psiIC.length;

			variables = new ProblemQuantities(psiIC, temperature, inEquationStateID, inParameterID);
			geometry = new GFGeometry(z, spaceDeltaZ, controlVolume);
			parameters = new Parameters(waterDensity, iceDensity, specificThermalCapacityWater,
					specificThermalCapacityIce, thermalConductivityWater, thermalConductivityIce, latentHeatFusion,
					referenceTemperatureInternalEnergy, referenceTemperatureSWRC, beta0, thetaS, thetaR,
					soilParticlesDensity, specificThermalCapacitySoilParticles, thermalConductivitySoilParticles,
					meltingTemperature, par1SWRC, par2SWRC, par3SWRC, par4SWRC, par5SWRC, ks, alphaSpecificStorage,
					betaSpecificStorage);

			computeQuantitiesInternalEnergy = new ComputeQuantitiesInternalEnergy(typeClosureEquation,
					typeEquationState, typeThermalConductivity, interfaceThermalConductivityModel, topBCType,
					bottomBCType, variables, geometry, parameters);

			outputToBuffer = new ArrayList<double[]>();

			List<EquationState> equationState = computeQuantitiesInternalEnergy.getInternalEnergyStateEquation();

			topBoundaryCondition = IBoundaryCondition.createDiffusionBoundaryCondition(topBCType);
			bottomBoundaryCondition = IBoundaryCondition.createDiffusionBoundaryCondition(bottomBCType);

			diffusionSolver = new Diffusion1DKernel(topBoundaryCondition, bottomBoundaryCondition, KMAX,
					nestedNewton, newtonTolerance, delta, MAXITER_NEWT, equationState);

		} // close step==0

		variables.internalEnergyTopBCValue = 0.0;
		if (topBCType == IBoundaryCondition.DiffusionBoundaryConditionType.TOP_NEUMANN) {
			variables.internalEnergyTopBCValue = inTopBC.get(stationID)[0] / tTimeStep;
		} else {
			variables.internalEnergyTopBCValue = inTopBC.get(stationID)[0] + 273.15;
		}

		variables.internalEnergyBottomBCValue = 0.0;
		if (bottomBCType == IBoundaryCondition.DiffusionBoundaryConditionType.BOTTOM_NEUMANN) {
			variables.internalEnergyBottomBCValue = inBottomBC.get(stationID)[0] / tTimeStep;
		} else {
			variables.internalEnergyBottomBCValue = inBottomBC.get(stationID)[0] + 273.15;
		}

		outputToBuffer.clear();

		double sumTimeDelta = 0;

		while (sumTimeDelta < tTimeStep) {

			if (sumTimeDelta + timeDelta > tTimeStep) {
				timeDelta = tTimeStep - sumTimeDelta;
			}
			sumTimeDelta = sumTimeDelta + timeDelta;

			/*
			 * Compute internal energy
			 */
			computeQuantitiesInternalEnergy.computeInternalEnergy(KMAX);

			/*
			 * Compute xStar
			 */
			computeQuantitiesInternalEnergy.computeXStar(KMAX);

			/*
			 * Solve PDE
			 */
			for (int picard = 0; picard < picardIteration; picard++) {

				/*
				 * Compute thermal conductivity
				 * 
				 */
				computeQuantitiesInternalEnergy.computeThermalConductivity(KMAX);

				computeQuantitiesInternalEnergy.computeInterfaceThermalConductivity(KMAX);

				/*
				 * Solve PDE
				 */
				variables.temperatures = diffusionSolver.solve(timeDelta, variables.internalEnergyBottomBCValue,
						variables.internalEnergyTopBCValue, KMAX, variables.lambdasInterface, variables.internalEnergys,
						geometry.spaceDeltaZ, variables.heatSourcesSinksTerm, variables.temperatures,
						variables.waterSuctions, variables.parameterID, variables.equationStateID);

			} // close Picard iteration

			/*
			 * Compute - internal energy and total internal energy
			 */
			computeQuantitiesInternalEnergy.computeInternalEnergyNew(KMAX);

			/*
			 * Compute fluxes
			 */
			computeQuantitiesInternalEnergy.computeConductionHeatFlux(KMAX);

			// TODO: check if this addition makes indeed sense. Added the following
			// line because top/bottom values are not set anywere (hence resulting
			// in always being 0 per java default)
			variables.heatFluxTop = variables.conductionHeatFluxs[KMAX];
			variables.heatFluxBottom = variables.conductionHeatFluxs[0];

			/*
			 * Compute error
			 */
			computeQuantitiesInternalEnergy.computeErrorDiffusion(KMAX, timeDelta);

		}
		
		outTemperature = variables.temperatures.clone();
		outTheta = variables.thetas.clone();
		outInternalEnergy = variables.internalEnergys.clone();
		outDiffusionHeatFlux = variables.conductionHeatFluxs.clone();
		outErrorInternalEnergy = variables.errorInternalEnergy;
		outHeatFluxTop = variables.heatFluxTop;
		outHeatFluxBottom = variables.heatFluxBottom;

		// TODO remove this at some point
		outputToBuffer.add(variables.temperatures);
		outputToBuffer.add(variables.thetas);
		outputToBuffer.add(variables.internalEnergys);
		outputToBuffer.add(variables.conductionHeatFluxs);
		outputToBuffer.add(new double[] { variables.errorInternalEnergy });
		outputToBuffer.add(new double[] { variables.heatFluxTop });
		outputToBuffer.add(new double[] { variables.heatFluxBottom });
		step++;

	} //// MAIN CYCLE END ////

} /// CLOSE ///
