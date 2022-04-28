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

package it.geoframe.blogspot.whetgeo1d.solutetransport;


import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import it.geoframe.blogspot.closureequation.closureequation.Parameters;
import it.geoframe.blogspot.closureequation.equationstate.EquationState;

import it.geoframe.blogspot.whetgeo1d.boundaryconditions.BoundaryCondition;
import it.geoframe.blogspot.whetgeo1d.boundaryconditions.DiffusionSimpleBoundaryConditionFactory;
import it.geoframe.blogspot.whetgeo1d.boundaryconditions.RichardsSimpleBoundaryConditionFactory;
import it.geoframe.blogspot.whetgeo1d.data.*;

import it.geoframe.blogspot.whetgeo1d.data.ComputeQuantitiesRichards;
import it.geoframe.blogspot.whetgeo1d.data.Geometry;
import it.geoframe.blogspot.whetgeo1d.data.ProblemQuantities;
import it.geoframe.blogspot.whetgeo1d.pdefinitevolume.AdvectionDiffusion1DFiniteVolumeSolver;

import it.geoframe.blogspot.whetgeo1d.pdefinitevolume.Richards1DFiniteVolumeSolver;
import oms3.annotations.*;


@Description("Solve the solute advection dispersion equation in the conservative form for the 1D domain")
@Documentation("")
@Author(name = "Concetta D'Amato, Niccolo' Tubini, and Riccardo Rigon", contact = "concetta.damato@unitn.it")
@Keywords("Hydrology, Richards, Solute transport, Infiltration")
@Bibliography("Casulli (2010), Stumpp et al., 2012")
//@Label()
//@Name()
//@Status()
@License("General Public License Version 3 (GPLv3)")
public class RichardsRootConservativeSoluteADESolver1DMain2 {

	
	/* 
	 * SOLUTE TRANSPORT PARAMETERS
	 */
	@Description("Molecular Diffusion in free water. Default value 10-9[m2 s-1].")
	@In 
	@Unit ("m2 s-1")
	public double[] molecularDiffusion;
	
	@Description("")
	@In 
	@Unit ("m")
	public double[] longitudinalDispersivity; 
	
	//@Description("")
	//@In 
	//@Unit ("-")
	//public double tortuosityFactor; 
	

	/*
	 * SOIL PARAMETERS
	 */
	@Description("The hydraulic conductivity at saturation")
	@In 
	@Unit ("m/s")
	public double[] ks;

	@Description("Saturated water content")
	@In 
	@Unit ("-")
	public double[] thetaS;

	@Description("Residual water content")
	@In 
	@Unit ("-")
	public double[] thetaR;

	@Description("First parameter of SWRC")
	@In 
	@Unit ("-")
	public double[] par1SWRC;

	@Description("Second parameter of SWRC")
	@In 
	@Unit ("-")
	public double[] par2SWRC;

	@Description("Third parameter of SWRC")
	@In 
	@Unit ("-")
	public double[] par3SWRC;

	@Description("Fourth parameter of SWRC")
	@In 
	@Unit ("-")
	public double[] par4SWRC;

	@Description("Fifth parameter of SWRC")
	@In 
	@Unit ("-")
	public double[] par5SWRC;

	@Description("Aquitard compressibility")
	@In 
	@Unit ("1/Pa")
	public double[] alphaSpecificStorage;

	@Description("Water compressibility")
	@In 
	@Unit ("1/Pa")
	public double[] betaSpecificStorage;

	@Description("Coefficient for water suction dependence on temperature")
	@In 
	@Unit ("K")
	public double beta0 = -776.45;

	@Description("Reference temperature for soil water content")
	@In 
	@Unit ("K")
	public double referenceTemperatureSWRC = 278.15;
	
	/*@Description("Reference temperature to compute internal energy")
	@In 
	@Unit ("K")
	public double referenceTemperatureInternalEnergy = 273.15;
	
	@Description("Soil particles density")
	@In 
	@Unit ("kg m-3")
	public double[] soilParticlesDensity;

	@Description("Specific thermal capacity of soil particles")
	@In 
	@Unit ("J kg-1 K-1")
	public double[] specificThermalCapacitySoilParticles;

	@Description("Thermal conductivity of soil particles")
	@In 
	@Unit ("W m-1 K-1")
	public double[] thermalConductivitySoilParticles;

	@Description("Melting temperature")
	@In 
	@Unit ("K")
	public double[] meltingTemperature;*/
	
	@Description("Control volume label defining the equation state")
	@In 
	@Unit("-")
	public int[] inEquationStateID;

	@Description("Control volume label defining the set of the paramters")
	@In 
	@Unit("-")
	public int[] inParameterID;
	
	/*
	 * MODELS
	 * - closure equation
	 * - conductivity model
	 * - interface conductivity model
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
	
	@Description("It is possible to choose among these models:"
			+ "Mualem Van Genuchten, Mualem Brooks Corey, ....")
	@In 
	public String[] typeUHCModel;

	@Description("It is possible to choose among these models:"
			+ "notemperature, ....")
	@In 
	public String typeUHCTemperatureModel;


	@Description("Hydraulic conductivity at control volume interface can be evaluated as"
			+ " the average of kappas[i] and kappas[i+1]"
			+ " the maximum between kappas[i] and kappas[i+1]"
			+ " the minimum between kappas[i] and kappas[i+1]"
			+ " a weighted average of kappas[i] and kappas[i+1] where weights are dx[i] and dx[i+1]")
	@In
	public String interfaceHydraulicConductivityModel;
	
	// Heat equation	
	/*@Description("Equation state")
	@In 
	public String[] typeInternalEnergyEquationState;

	@Description("Thermal conductivity models")
	@In 
	public String[] typeThermalConductivity;

	@Description("Thermal conductivity at control volume interface can be evaluated as"
			+ " the average of kappas[i] and kappas[i+1]"
			+ " the maximum between kappas[i] and kappas[i+1]"
			+ " the minimum between kappas[i] and kappas[i+1]"
			+ " a weighted average of kappas[i] and kappas[i+1] where weights are dx[i] and dx[i+1]")
	@In
	public String interfaceThermalConductivityModel;*/
	
	@Description("Dispersion Coefficient at control volume interface can be evaluated as"
			+ " the average of kappas[i] and kappas[i+1]"
			+ " the maximum between kappas[i] and kappas[i+1]"
			+ " the minimum between kappas[i] and kappas[i+1]"
			+ " a weighted average of kappas[i] and kappas[i+1] where weights are dx[i] and dx[i+1]")
	@In
	public String interfaceDispersionModel;


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
	public double[] temperatureIC;
	
	@Description("Initial condition for concentration read from grid NetCDF file")
	@In
	@Unit("-")
	public double[] concentrationIC;

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
	@Unit ("s")
	public double tTimeStep;

	@Description("Time step of integration")
	@In
	@Unit ("s")
	public double timeDelta;

	/*
	 * ITERATION PARAMETERS
	 */
	@Description("Tolerance for Newton iteration")
	@In
	public double newtonTolerance;

	@Description("Control parameter for nested Newton algorithm:"
			+"0 --> simple Newton method"
			+"1 --> nested Newton method")
	@In
	public int nestedNewton; 

	@Description("Damped factor for Newton algorithm")
	@In
	public double delta = 0.0; 
	
	@Description("Number of Picard iteration to update the diffusive flux matrix")
	@In
	public int picardIteration=1;

	/*
	 *  BOUNDARY CONDITIONS
	 */
	@Description("The station ID in the timeseries file")
	@In
	@Unit ("-")
	public int stationID;
	
	@Description("The HashMap with the time series of the boundary condition at the top of soil column")
	@In
	@Unit ("m")
	public HashMap<Integer, double[]> inRichardsTopBC;

	@Description("It is possibile to chose between 2 different kind "
			+ "of boundary condition at the top of the domain: "
			+ "- Dirichlet boundary condition --> Top Dirichlet"
			+ "- Neumann boundary condition --> Top Neumann")
	@In 
	public String topRichardsBCType;

	@Description("The HashMap with the time series of the boundary condition at the bottom of soil column")
	@In
	@Unit ("m")
	public HashMap<Integer, double[]> inRichardsBottomBC;
	
	@Description("It is possibile to chose among 2 different kind "
			+ "of boundary condition at the bottom of the domain: "
			+ "- Dirichlet boundary condition --> Bottom Dirichlet"
			+ "- Neumann boundary condition --> Bottom Neumann")
	@In 
	public String bottomRichardsBCType;
	
	/*@Description("The HashMap with the time series of the boundary condition at the top of soil column")
	@In
	@Unit ("m")
	public HashMap<Integer, double[]> inInternalEnergyTopBC;

	@Description("It is possibile to chose between 2 different kind "
			+ "of boundary condition at the top of the domain: "
			+ "- Dirichlet boundary condition --> Top Dirichlet"
			+ "- Neumann boundary condition --> Top Neumann")
	@In 
	public String topInternalEnergyBCType;
	

	@Description("The HashMap with the time series of the boundary condition at the bottom of soil column")
	@In
	@Unit ("")
	public HashMap<Integer, double[]> inInternalEnergyBottomBC;
	
	@Description("It is possibile to chose among 2 different kind "
			+ "of boundary condition at the bottom of the domain: "
			+ "- Dirichlet boundary condition --> Bottom Dirichlet"
			+ "- Neumann boundary condition --> Bottom Neumann")
	@In 
	public String bottomInternalEnergyBCType;*/
	
	@Description("The HashMap with the time series of the boundary condition at the top of soil column")
	@In
	@Unit ("m")
	public HashMap<Integer, double[]> inSoluteTopBC;
	
	@Description("It is possibile to chose between 2 different kind "
			+ "of boundary condition at the top of the domain: "
			+ "- Dirichlet boundary condition --> Top Dirichlet"
			+ "- Neumann boundary condition --> Top Neumann")
	@In 
	public String topSoluteBCType;
	
	@Description("The HashMap with the time series of the boundary condition at the bottom of soil column")
	@In
	@Unit ("")
	public HashMap<Integer, double[]> inSoluteBottomBC;
	
	@Description("It is possibile to chose among 2 different kind "
			+ "of boundary condition at the bottom of the domain: "
			+ "- Dirichlet boundary condition --> Bottom Dirichlet"
			+ "- Neumann boundary condition --> Bottom Neumann")
	@In 
	public String bottomSoluteBCType;

	@Description("")
	@In
	@Unit ("")
	public HashMap<Integer, double[]> inSaveDate;

	@Description("The current date of the simulation.")
	@In
	@Out
	public String inCurrentDate;


	/*
	 * OUTPUT
	 */

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
	@Unit (" ")
	private int KMAX; 

	@Description("It is needed to iterate on the date")
	private int step;

	@Description("Control value to save output:"
			+ "- 1 save the current time step output"
			+ "- 0 do not save")
	private double saveDate;

	private Richards1DFiniteVolumeSolver richardsSolver;
	private AdvectionDiffusion1DFiniteVolumeSolver advectionDispersionSolver;
	private ProblemQuantities variables;
	private Geometry geometry;
	private Parameters parameters;
	private ComputeQuantitiesRichards computeQuantitiesRichards;
	private ComputeQuantitiesSoluteAdvectionDispersion computeQuantitiesSoluteAdvectionDispersion;
	private BoundaryCondition topRichardsBoundaryCondition;
	private BoundaryCondition bottomRichardsBoundaryCondition;
	private BoundaryCondition topSoluteBoundaryCondition; //cambiato
	private BoundaryCondition bottomSoluteBoundaryCondition;  //cambiato
	private RichardsSimpleBoundaryConditionFactory boundaryRichardsConditionFactory;
	private DiffusionSimpleBoundaryConditionFactory boundarySoluteConditionFactory;  //cambiato

	@Execute
	public void solve() {



		if(step==0){
			KMAX = psiIC.length;

			variables = ProblemQuantities.getInstance(psiIC, temperatureIC, concentrationIC, inEquationStateID, inParameterID);
			geometry = Geometry.getInstance(z, spaceDeltaZ, controlVolume);
			parameters = Parameters.getInstance(molecularDiffusion,longitudinalDispersivity,referenceTemperatureSWRC, beta0,
					thetaS, thetaR, par1SWRC, par2SWRC, par3SWRC, par4SWRC, par5SWRC, ks, alphaSpecificStorage, betaSpecificStorage); // HO FATTO UN NUOVO getInstance su closure equation 

			computeQuantitiesRichards = new ComputeQuantitiesRichards(typeClosureEquation, typeRichardsEquationState, typeUHCModel, typeUHCTemperatureModel, interfaceHydraulicConductivityModel, topRichardsBCType, bottomRichardsBCType);

			computeQuantitiesSoluteAdvectionDispersion = new ComputeQuantitiesSoluteAdvectionDispersion(typeClosureEquation, interfaceDispersionModel, topSoluteBCType, bottomSoluteBCType); //	CAPIRE SE SI DEVE LASCIARE typeInternalEnergyEquationState
			
			outputToBuffer = new ArrayList<double[]>();

			List<EquationState> richardsEquationState = computeQuantitiesRichards.getRichardsStateEquation();

			
			boundaryRichardsConditionFactory = new RichardsSimpleBoundaryConditionFactory();
			topRichardsBoundaryCondition = boundaryRichardsConditionFactory.createBoundaryCondition(topRichardsBCType);
			bottomRichardsBoundaryCondition = boundaryRichardsConditionFactory.createBoundaryCondition(bottomRichardsBCType);
			
			richardsSolver = new Richards1DFiniteVolumeSolver(topRichardsBoundaryCondition, bottomRichardsBoundaryCondition, KMAX, nestedNewton, newtonTolerance, delta, MAXITER_NEWT, richardsEquationState);

			
			boundarySoluteConditionFactory = new DiffusionSimpleBoundaryConditionFactory();
			topSoluteBoundaryCondition = boundarySoluteConditionFactory.createBoundaryCondition(topSoluteBCType);
			bottomSoluteBoundaryCondition = boundarySoluteConditionFactory.createBoundaryCondition(bottomSoluteBCType);	
			
			advectionDispersionSolver = new AdvectionDiffusion1DFiniteVolumeSolver(topSoluteBoundaryCondition, bottomSoluteBoundaryCondition, KMAX);

			for(int element = 0; element < KMAX; element++) {
			variables.soluteQuantitiesTransported [element] = 1 ;} //questo lo abbiamo aggiunto perchè il metodo per la ADE vuole in input una waterCapacityTransported

		} // close step==0
		

		doProcessBuffer = false;
		System.out.println(inCurrentDate);
		
		variables.richardsTopBCValue = 0.0;
		if(topRichardsBCType.equalsIgnoreCase("Top Neumann") || topRichardsBCType.equalsIgnoreCase("TopNeumann") || topRichardsBCType.equalsIgnoreCase("Top Coupled") || topRichardsBCType.equalsIgnoreCase("TopCoupled")) {
			variables.richardsTopBCValue = (inRichardsTopBC.get(stationID)[0]/1000)/tTimeStep;
		} else {
			variables.richardsTopBCValue = inRichardsTopBC.get(stationID)[0]/1000;
		}
		

		variables.richardsBottomBCValue = 0.0;
		variables.richardsBottomBCValue = inRichardsBottomBC.get(stationID)[0];


		/*variables.internalEnergyTopBCValue = 0.0;
		if(topInternalEnergyBCType.equalsIgnoreCase("Top Neumann") || topInternalEnergyBCType.equalsIgnoreCase("TopNeumann")) {
			variables.internalEnergyTopBCValue = inInternalEnergyTopBC.get(stationID)[0]/tTimeStep;
		} else {
			variables.internalEnergyTopBCValue = inInternalEnergyTopBC.get(stationID)[0]+273.15;
		}
		
		variables.internalEnergyBottomBCValue = 0.0;
		if(bottomInternalEnergyBCType.equalsIgnoreCase("Bottom Neumann") || bottomInternalEnergyBCType.equalsIgnoreCase("BottomNeumann")) {
			variables.internalEnergyBottomBCValue = inInternalEnergyBottomBC.get(stationID)[0]/tTimeStep;
		} else {
			variables.internalEnergyBottomBCValue = inInternalEnergyBottomBC.get(stationID)[0]+273.15;
		}
		*/
		
		variables.soluteTopBCValue = 0.0;
		if(topSoluteBCType.equalsIgnoreCase("Top Neumann") || topSoluteBCType.equalsIgnoreCase("TopNeumann")) {
			variables.soluteTopBCValue = inSoluteTopBC.get(stationID)[0]/tTimeStep;
		} else {
			variables.soluteTopBCValue = inSoluteTopBC.get(stationID)[0];
		}
		

		variables.soluteBottomBCValue = 0.0;
		if(bottomSoluteBCType.equalsIgnoreCase("Bottom Neumann") || bottomSoluteBCType.equalsIgnoreCase("BottomNeumann")) {
			variables.soluteBottomBCValue = inSoluteBottomBC.get(stationID)[0]/tTimeStep;
		} else {
			variables.soluteBottomBCValue = inSoluteBottomBC.get(stationID)[0];
		}

		saveDate = -1.0;
		saveDate = inSaveDate.get(stationID)[0];
		outputToBuffer.clear();

		double sumTimeDelta = 0;

		
		while(sumTimeDelta < tTimeStep) {

			
			if(sumTimeDelta + timeDelta>tTimeStep) {
				timeDelta = tTimeStep - sumTimeDelta;
			}
			sumTimeDelta = sumTimeDelta + timeDelta;

	
			
			/*
			 * Compute water volumes
			 */
			computeQuantitiesRichards.computeWaterVolume(KMAX);
			computeQuantitiesRichards.computeThetas(KMAX);
			
			/*
			 * Compute heat capacity
			 */
			//computeQuantitiesSoluteAdvectionDispersion.computeHeatCapacity(KMAX); //Non mi serve 
			
			computeQuantitiesSoluteAdvectionDispersion.computeWaterVolumeConcentrations(KMAX);
			/*
			 * Compute dispersion coefficient
			 */
			computeQuantitiesSoluteAdvectionDispersion.computeThetasInterface(KMAX);
			
			computeQuantitiesSoluteAdvectionDispersion.computeTortuosityFactorsInterface(KMAX);
			
			computeQuantitiesSoluteAdvectionDispersion.computeDispersionCoefficients(KMAX); 
			computeQuantitiesSoluteAdvectionDispersion.computeDispersionFactors(KMAX);
			
			
			
			//variables.lambdasInterface[KMAX] = 0.6;
			//variables.dispersionFactorsInterface[KMAX] = variables.dispersionFactorsInterface[KMAX-1]; già lo faccio dentro il metodo

			//computeQuantitiesSoluteAdvectionDispersion.computeTransportedQuantity(KMAX); //Non mi serve
			
			/*
			 * Compute xStar
			 */
			computeQuantitiesRichards.computeXStar(KMAX);
			
			

			/*
			 * Solve Richards equation
			 */
			for(int picard=0; picard<picardIteration; picard++) {

				/*
				 * Compute hydraulic conductivity
				 * 
				 */	
				computeQuantitiesRichards.computeHydraulicConductivity(KMAX);

				computeQuantitiesRichards.computeInterfaceHydraulicConductivity(KMAX);
				
				/*
				 * Solve PDE
				 */
				variables.waterSuctions = richardsSolver.solve(timeDelta, variables.richardsBottomBCValue, variables.richardsTopBCValue, KMAX, variables.kappasInterface,
						variables.volumes, geometry.spaceDeltaZ, variables.ETs, variables.waterSuctions, variables.temperatures, variables.parameterID, variables.equationStateID);

			} // close Picard iteration
			
			
			
			/*
			 * compute run-off
			 */
			computeQuantitiesRichards.computeRunOff(KMAX, maxPonding);
			

			/*
			 * Compute 
			 * - water volume and total water volume
			 * - water content
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
			
		
			/*
			 * Solve solute advection-dispersion equation
			 */

			
			/*for(int k=0; k<KMAX; k++) {
				variables.temperatures[k] = variables.temperatures[k]-273.15;
			}
			variables.internalEnergyTopBCValue = variables.internalEnergyTopBCValue-273.15; 
			variables.internalEnergyBottomBCValue = variables.internalEnergyBottomBCValue-273.15; 
			*/
			if(variables.thetasNew[variables.thetasNew.length-1]<=0) {
				KMAX = KMAX-1;
				variables.waterVolumeConcentration-=variables.waterVolumeConcentrations[variables.thetasNew.length-1];
			}
			

			/*variables.temperatures = advectionDiffusionSolver.solve(timeDelta, variables.internalEnergyBottomBCValue, variables.internalEnergyTopBCValue, KMAX, variables.lambdasInterface,
						variables.heatCapacitysNew, variables.heatCapacitys, geometry.spaceDeltaZ, variables.heatSourcesSinksTerm, variables.temperatures, variables.waterSuctions, variables.darcyVelocities, 
						variables.waterCapacityTransported, variables.parameterID, variables.equationStateID);
			*/
			
			
			variables.concentrations = advectionDispersionSolver.solve(timeDelta, variables.soluteBottomBCValue, variables.soluteTopBCValue, KMAX, variables.dispersionFactors,
					variables.volumesNew, variables.volumes, geometry.spaceDeltaZ, variables.soluteSourcesSinksTerm, variables.concentrations, variables.waterSuctions, variables.darcyVelocities, 
					variables.soluteQuantitiesTransported, variables.parameterID, variables.equationStateID);
		
			
			
			/*for(int k=0; k<KMAX; k++) {
				variables.temperatures[k] = variables.temperatures[k]+273.15;
			}
			variables.internalEnergyTopBCValue = variables.internalEnergyTopBCValue+273.15; 
			variables.internalEnergyBottomBCValue = variables.internalEnergyBottomBCValue+273.15; 
*/
			
			computeQuantitiesSoluteAdvectionDispersion.computeWaterVolumeConcentrationsNew(KMAX);

	
			//computeQuantitiesSoluteAdvectionDispersion.computeConductionHeatFlux(KMAX); // 
			//computeQuantitiesSoluteAdvectionDispersion.computeAdvectionHeatFlux(KMAX);
			//computeQuantitiesSoluteAdvectionDispersion.computeHeatFlux(KMAX);
			
			computeQuantitiesSoluteAdvectionDispersion.computeDispersionSoluteFluxes(KMAX);
			computeQuantitiesSoluteAdvectionDispersion.computeAdvectionSoluteFluxes(KMAX);
			computeQuantitiesSoluteAdvectionDispersion.computeSoluteFluxes(KMAX);
			
			computeQuantitiesSoluteAdvectionDispersion.computeAverageSoluteConcentration(KMAX);
			computeQuantitiesSoluteAdvectionDispersion.computeAverageWaterVolumeSoluteConcentration(KMAX);
			
			/*
			 * Compute error advection dispersion equation 
			 */
			computeQuantitiesSoluteAdvectionDispersion.computeError(KMAX, timeDelta); 

			
			/*if(variables.thetasNew[variables.thetasNew.length-1]<=0) {
				KMAX = KMAX+1;
				variables.temperatures[KMAX-1] = variables.temperatures[KMAX-2]; // 
			}*/
			
			if(variables.thetasNew[variables.thetasNew.length-1]<=0) {
				KMAX = KMAX+1;
				variables.concentrations[KMAX-1] = variables.concentrations[KMAX-2];
			}

		
		}

		if(saveDate == 1) {
			outputToBuffer.add(variables.waterSuctions);
			outputToBuffer.add(variables.thetasNew);
			outputToBuffer.add(variables.volumesNew);
			outputToBuffer.add(variables.darcyVelocities);
			outputToBuffer.add(variables.ETs);
			
			outputToBuffer.add(variables.concentrations);
			outputToBuffer.add(variables.waterVolumeConcentrationsNew);
			
			outputToBuffer.add(variables.soluteFluxes);
			outputToBuffer.add(variables.dispersionSoluteFluxes);
			outputToBuffer.add(variables.advectionSoluteFluxes);
		
			outputToBuffer.add(new double[] {variables.errorWaterVolumeConcentration});
			outputToBuffer.add(new double[] {variables.errorVolume});
			outputToBuffer.add(new double[] {variables.averageSoluteConcentration});
			outputToBuffer.add(new double[] {variables.averageWaterVolumeSoluteConcentration});
			
			doProcessBuffer = true;
		} else {
			//			System.out.println("SaveDate = " + saveDate);
		}
		step++;

	} //// MAIN CYCLE END ////

}  /// CLOSE ///








































