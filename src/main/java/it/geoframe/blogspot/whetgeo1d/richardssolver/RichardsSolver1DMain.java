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

package it.geoframe.blogspot.whetgeo1d.richardssolver;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import it.geoframe.blogspot.conductivitymodel.ConductivityEquation;
import it.geoframe.blogspot.conductivitymodel.ConductivityEquationFactory;
import it.geoframe.blogspot.conductivitymodel.UnsaturatedHydraulicConductivityTemperatureFactory;
import it.geoframe.blogspot.interfaceconductivity.InterfaceConductivity;
import it.geoframe.blogspot.interfaceconductivity.SimpleInterfaceConductivityFactory;
import it.geoframe.blogspot.whetgeo1d.data.Geometry;
import it.geoframe.blogspot.whetgeo1d.data.ProblemQuantities;
import it.geoframe.blogspot.whetgeo1d.equationstate.EquationStateFactory;
import it.geoframe.blogspot.whetgeo1d.equationstate.WaterDepth;
import oms3.annotations.*;
import it.geoframe.blogspot.closureequation.ClosureEquation;
import it.geoframe.blogspot.closureequation.Parameters;
import it.geoframe.blogspot.closureequation.SoilWaterRetentionCurveFactory;
import it.geoframe.blogspot.equationstate.EquationState;




@Description("Solve the Richards equation for the 1D domain.")
@Documentation("")
@Author(name = "Niccolo' Tubini, and Riccardo Rigon", contact = "tubini.niccolo@gmail.com")
@Keywords("Hydrology, Richards, Infiltration")
@Bibliography("Casulli (2010)")
//@Label()
//@Name()
//@Status()
@License("General Public License Version 3 (GPLv3)")
public class RichardsSolver1DMain {

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

	@Description("Water content at the whilting point")
	@In 
	@Unit ("-")
	public double[] thetaWp;

	@Description("Water content at field capacity")
	@In 
	@Unit ("-")
	public double[] thetaFc;

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
	public double temperatureR = 278.15;
	
	@Description("Control volume label defining the rheology")
	@In 
	@Unit("-")
	public int[] inRheologyID;

	@Description("Control volume label defining the set of the paramters")
	@In 
	@Unit("-")
	public int[] inParameterID;
	
	/*
	 * MODELS
	 */

	@Description("It is possibile to chose between 3 different models to compute "
			+ "the soil hydraulic properties: Van Genuchten; Brooks and Corey; Kosugi unimodal")
	@In 
	public String soilHydraulicModel;

	@Description("It is possible to choose among these models:"
			+ "notemperature, ....")
	@In 
	public String typeUHCTemperatureModel;

	@Description("It is possible to choose among these models:"
			+ "Mualem Van Genuchten, Mualem Brooks Corey, ....")
	@In 
	public String typeUHCModel;

	@Description("Hydraulic conductivity at control volume interface can be evaluated as"
			+ " the average of kappas[i] and kappas[i+1]"
			+ " the maximum between kappas[i] and kappas[i+1]"
			+ " the minimum between kappas[i] and kappas[i+1]"
			+ " a weighted average of kappas[i] and kappas[i+1] where weights are dx[i] and dx[i+1]")
	@In
	public String interfaceHydraulicConductivityModel;


	/////////////////////////////////////////////

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
	 * Time step
	 */
	@Description("Time amount at every time-loop")
	@In
	@Unit ("s")
	public double tTimestep;

	@Description("Time step of integration")
	@In
	@Unit ("s")
	public double timeDelta;

	/*
	 * Iteration parameter
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
	 * Lysimeter
	 */
	@Description("Stressed Evapotranspiration for each layer")
	@In 
	@Unit ("mm/s")
	public double [] StressedETs;

	@Description("Boolean value to use lysimeter")
	@In
	public boolean lysimeter = false;

	@Description("Stressed Evapotranspiration for each layer")
	@Out
	@Unit ("m")
	public double[] ETs;

	@Description("Sum of Stressed Evapotranspiration")
	@Out 
	@Unit ("m")
	public double sumETs;	

	@Description("Sum of Stressed Evapotranspiration")
	@Out 
	@Unit ("-")
	public double[] thetasNew;	

	/*
	 *  BOUNDARY CONDITIONS
	 */

	@Description("The HashMap with the time series of the boundary condition at the top of soil column")
	@In
	@Unit ("m")
	public HashMap<Integer, double[]> inTopBC;

	@Description("It is possibile to chose between 2 different kind "
			+ "of boundary condition at the top of the domain: "
			+ "- Dirichlet boundary condition --> Top Dirichlet"
			+ "- Neumann boundary condition --> Top Neumann")
	@In 
	public String topBCType;

	@Description("The HashMap with the time series of the boundary condition at the bottom of soil column")
	@In
	@Unit ("m")
	public HashMap<Integer, double[]> inBottomBC;

	@Description("")
	@In
	@Unit ("")
	public HashMap<Integer, double[]> inSaveDate;

	@Description("It is possibile to chose among 3 different kind "
			+ "of boundary condition at the bottom of the domain: "
			+ "- Dirichlet boundary condition --> Bottom Dirichlet"
			+ "- Neumann boundary condition --> Bottom Neumann"
			+ "- Impervious boundary condition --> Bottom Impervious")
	@In 
	public String bottomBCType;

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
	final int MAXITER_NEWT = 50;

	@Description("Top boundary condition according with topBCType")
	@Unit ("")
	double topBC;

	@Description("Bottom boundary condition according with bottomBCType")
	@Unit ("")
	double bottomBC;

	@Description("Number of control volume for domain discetrization")
	@Unit (" ")
	int KMAX; 

	@Description("It is needed to iterate on the date")
	int step;

	double saveDate;
	///////////////////////////////
	double volume1;
	double volume2;
	double volumeLost;

	int[] rheologyID;
	int[] parameterID;

	PDE1DSolver richardsSolver;
	ProblemQuantities variables;
	Geometry geometry;
	Parameters rehologyParameters;
	SoilWaterRetentionCurveFactory soilWaterRetentionCurveFactory;

	@Description("This list contains the objects that describes the state equations of the problem")
	List<EquationState> equationState;

	@Description("Object for the SWRC model")
	ClosureEquation soilWaterRetentionCurve;

	@Description("Object dealing with the state equation")
//	EquationState internalEnergy;
	EquationStateFactory equationStateFactory;


	@Description("Object dealing with the hydraulic conductivity model")
	ConductivityEquation hydraulicConductivity;
	ConductivityEquationFactory conductivityEquationFactory;
	UnsaturatedHydraulicConductivityTemperatureFactory unsaturatedHydraulicConductivityTemperatureFactory;


	@Description("This object compute the interface hydraulic conductivity accordingly with the prescribed method.")
	InterfaceConductivity interfaceConductivity;
	SimpleInterfaceConductivityFactory interfaceConductivityFactory;

	@Execute
	public void solve() {

		/// da cancellare
		//		if(inCurrentDate.equalsIgnoreCase("2015-11-01 02:00")) {
		//			System.out.println("RICHARDS 1D "+inCurrentDate);
		//		}


		if(step==0){
			KMAX = psiIC.length;

			variables = ProblemQuantities.getInstance(psiIC, temperature);
			geometry = Geometry.getInstance(z, spaceDeltaZ, controlVolume);
			rehologyParameters = Parameters.getInstance(1000.0, 970, 4188,
					2117, 0.6, 2.29, 333700, 273.15,
					thetaS, thetaR, new double[] {-9999.0}, new double[] {-9999.0}, new double[] {-9999.0},
					new double[] {-9999.0}, par1SWRC, par2SWRC, par3SWRC, par4SWRC, par5SWRC, ks, alphaSpecificStorage, betaSpecificStorage);

			rheologyID = inRheologyID.clone();
			parameterID = inParameterID.clone();

			outputToBuffer = new ArrayList<double[]>();


			soilWaterRetentionCurveFactory = new SoilWaterRetentionCurveFactory();
			soilWaterRetentionCurve = soilWaterRetentionCurveFactory.create(soilHydraulicModel);

			equationStateFactory = new EquationStateFactory();

			equationState = new ArrayList<EquationState>();
			equationState.add(new WaterDepth(null));
			equationState.add(equationStateFactory.create(soilHydraulicModel, soilWaterRetentionCurve));


			conductivityEquationFactory = new ConductivityEquationFactory();
			hydraulicConductivity = conductivityEquationFactory.create(typeUHCModel, soilWaterRetentionCurve);

			unsaturatedHydraulicConductivityTemperatureFactory = new UnsaturatedHydraulicConductivityTemperatureFactory();
			hydraulicConductivity = unsaturatedHydraulicConductivityTemperatureFactory.create(typeUHCTemperatureModel, soilWaterRetentionCurve, hydraulicConductivity);


			interfaceConductivityFactory = new SimpleInterfaceConductivityFactory();
			interfaceConductivity = interfaceConductivityFactory.createInterfaceConductivity(interfaceHydraulicConductivityModel);

			richardsSolver = new PDE1DSolver(topBCType, bottomBCType, KMAX, nestedNewton, newtonTolerance, delta, MAXITER_NEWT,
					equationState, rheologyID, parameterID);

			if(topBCType.equalsIgnoreCase("Top Dirichlet")||topBCType.equalsIgnoreCase("TopDirichlet")) {
				KMAX = KMAX-1;
			}

		} // close step==0

		doProcessBuffer = false;


		topBC = 0.0;
		if(topBCType.equalsIgnoreCase("Top Neumann") || topBCType.equalsIgnoreCase("TopNeumann")) {
			topBC = (inTopBC.get(0)[0]/1000)/tTimestep;
		} else {
			topBC = inTopBC.get(0)[0]/1000;
		}
		

		bottomBC = 0.0;
		if(inBottomBC != null)
			bottomBC = inBottomBC.get(0)[0];
		if(bottomBCType.equalsIgnoreCase("Bottom Neumann") || bottomBCType.equalsIgnoreCase("BottomNeumann")) {
			bottomBC = bottomBC/tTimestep;
		}

		saveDate = -1.0;
		saveDate = inSaveDate.get(0)[0];
		outputToBuffer.clear();

		double sumTimeDelta = 0;


		while(sumTimeDelta < tTimestep) {

			variables.waterVolume = 0.0;
			variables.waterVolumeNew = 0.0;
			variables.runOff = 0.0;
			volumeLost = 0.0;

			if(sumTimeDelta + timeDelta>tTimestep) {
				timeDelta = tTimestep - sumTimeDelta;
			}
			sumTimeDelta = sumTimeDelta + timeDelta;



			/*
			 * Compute water volumes
			 */
			for(int element = 0; element < KMAX; element++) {
				variables.volumes[element] = equationState.get(rheologyID[element]).stateEquation(variables.waterSuctions[element], variables.temperatures[element], parameterID[element], element);
				variables.waterVolume += variables.volumes[element];
			}

			
			/*
			 * Compute xStar
			 */
			for(int element=0; element<KMAX; element++) {
				equationState.get(rheologyID[element]).computeXStar(variables.temperatures[element], parameterID[element], element);

			}
			
			
			/*
			 * Compute source term for evapotranspiration
			 * 
			 */
			if(lysimeter==true) {
				for(int i = 0; i < KMAX-2; i++) {
					variables.ETs[i] = (StressedETs[i]/1000)/tTimestep*timeDelta;
				} 
			} else {
				for(int i = 0; i < KMAX-2; i++) {
					variables.ETs[i] = 0.0;
				} 
			}
			

			if(lysimeter==true) {
				variables.sumETs = 0;
				for(int element = 0; element < KMAX-1; element++) {
					if (variables.ETs[element] > (variables.volumes[element] - thetaWp[parameterID[element]]*geometry.controlVolume[element])){
						variables.ETs[element] = variables.volumes[element] - thetaWp[parameterID[element]]*geometry.controlVolume[element];
						System.out.println("Errore nel calcolo di ETs. E' maggiore di volumes[i] - thetaR[i]*geometry.controlVolumex[i] ");} 
					else if (variables.ETs[element] <= (variables.volumes[element] - thetaWp[parameterID[element]]*geometry.controlVolume[element])){
						variables.ETs[element] = variables.ETs[element];}
					variables.sumETs = variables.sumETs + variables.ETs[element];
				}
			}
			/*
			 * Solve PDE
			 */
			for(int picard=0; picard<picardIteration; picard++) {

				/*
				 * Compute hydraulic conductivity
				 * 
				 */
				for(int element = 0; element < KMAX; element++) {
					if(element==KMAX-1) {
						if(this.topBCType.equalsIgnoreCase("Top Dirichlet") || this.topBCType.equalsIgnoreCase("TopDirichlet")){
							double kappaTop = hydraulicConductivity.k(topBC, variables.temperatures[element], parameterID[element], element); 
							variables.kappasInterface[element+1] =  interfaceConductivity.compute(variables.kappas[element], kappaTop, geometry.controlVolume[element], geometry.controlVolume[element]);
							variables.kappasInterface[element] =  interfaceConductivity.compute(variables.kappas[element-1],variables.kappas[element],geometry.controlVolume[element-1], geometry.controlVolume[element]);
						} else {
							variables.kappas[element] = hydraulicConductivity.k(variables.waterSuctions[element], variables.temperatures[element], parameterID[element-1], element-1);
						}
					} else {
						variables.kappas[element] = hydraulicConductivity.k(variables.waterSuctions[element], variables.temperatures[element], parameterID[element], element);
					}
					if(element==0) {
						if(this.bottomBCType.equalsIgnoreCase("Bottom Free Drainage") || this.bottomBCType.equalsIgnoreCase("BottomFreeDrainage")){
							variables.kappasInterface[element] =  variables.kappas[element];
						} else if (this.bottomBCType.equalsIgnoreCase("Bottom Impervious") || this.bottomBCType.equalsIgnoreCase("BottomImpervious")) {
							variables.kappasInterface[element] = + 0.0;
						} else {
							double kappaBottom = hydraulicConductivity.k(bottomBC, variables.temperatures[element], parameterID[element], element); 
							variables.kappasInterface[element] =  interfaceConductivity.compute(variables.kappas[element], kappaBottom, geometry.controlVolume[element], geometry.controlVolume[element]);
						}
					}
					else {
						variables.kappasInterface[element] = interfaceConductivity.compute(variables.kappas[element-1],variables.kappas[element],geometry.controlVolume[element-1], geometry.controlVolume[element]);
					}

				}			



				richardsSolver.solve(topBC, bottomBC, timeDelta, KMAX);

			} // close Picard iteration

			/*
			 * compute run-off
			 */
			if(this.topBCType.equalsIgnoreCase("Top Neumann") || this.topBCType.equalsIgnoreCase("TopNeumann")) {
				if(maxPonding>0 && variables.waterSuctions[KMAX -1]>maxPonding) {
					volumeLost = (variables.waterSuctions[KMAX -1] - maxPonding);
					variables.waterSuctions[KMAX -1] = maxPonding;
					variables.runOff += volumeLost;
				}
			}

			/*
			 * Compute 
			 * - water content
			 * - water volume
			 * - total water volume
			 */
			for(int element = 0; element < KMAX; element++) {		
				variables.volumes[element] = equationState.get(rheologyID[element]).stateEquation(variables.waterSuctions[element], variables.temperatures[element], parameterID[element], element);
				variables.waterVolumeNew += variables.volumes[element];
				if(element<KMAX-1) {
					variables.thetas[element] = soilWaterRetentionCurve.f(variables.waterSuctions[element], parameterID[element]);
				}
			}


			/*
			 * Fluxes
			 */
			for(int k = 0; k < KMAX; k++) {
				if( k == 0 ) {
					if(this.bottomBCType.equalsIgnoreCase("Bottom Free Drainage") || this.bottomBCType.equalsIgnoreCase("BottomFreeDrainage")){
						variables.darcyVelocities[k] = -variables.kappasInterface[k];
						variables.darcyVelocitiesCapillary[k] = 0.0; 
						variables.darcyVelocitiesGravity[k] = -variables.kappasInterface[k];
						variables.poreVelocities[k] = variables.darcyVelocities[k]/(variables.thetas[k]-rehologyParameters.thetaR[parameterID[k]]);
						variables.celerities[k] = -9999.0;
						variables.kinematicRatio[k] = variables.celerities[k]/variables.poreVelocities[k];
					} else if (this.bottomBCType.equalsIgnoreCase("Bottom Impervious") || this.bottomBCType.equalsIgnoreCase("BottomImpervious")) {
						variables.darcyVelocities[k] = + 0.0;
						variables.darcyVelocitiesCapillary[k] = 0.0; 
						variables.darcyVelocitiesGravity[k] = 0.0;
						variables.poreVelocities[k] = 0.0;
						variables.celerities[k] = 0.0;
						variables.kinematicRatio[k] = Double.NaN;

					} else {
						variables.darcyVelocities[k] = -variables.kappasInterface[k] * ( (variables.waterSuctions[k]-bottomBC)/geometry.spaceDeltaZ[k] + 1 );
						variables.darcyVelocitiesCapillary[k] = -variables.kappasInterface[k] * (variables.waterSuctions[k]-bottomBC)/geometry.spaceDeltaZ[k]; 
						variables.darcyVelocitiesGravity[k] = -variables.kappasInterface[k];
						variables.poreVelocities[k] = variables.darcyVelocities[k]/(variables.thetas[k]-rehologyParameters.thetaR[parameterID[k]]);
						variables.celerities[k] = -9999.0;
						variables.kinematicRatio[k] = variables.celerities[k]/variables.poreVelocities[k];

					}

				} else if( k == KMAX-1){
					if(this.topBCType.equalsIgnoreCase("Top Dirichlet") || this.topBCType.equalsIgnoreCase("TopDirichlet")){
						variables.darcyVelocities[k+1] = -variables.kappasInterface[k+1] * ( (topBC-variables.waterSuctions[k])/geometry.spaceDeltaZ[k+1] +1 );
						variables.darcyVelocitiesCapillary[k+1] = -variables.kappasInterface[k+1] * (topBC-variables.waterSuctions[k])/geometry.spaceDeltaZ[k+1];
						variables.darcyVelocitiesGravity[k+1] = -variables.kappasInterface[k+1];
						variables.poreVelocities[k+1] = variables.darcyVelocities[k+1]/(variables.thetas[k]-rehologyParameters.thetaR[parameterID[k]]);
						variables.celerities[k+1] = -9999.0;
						variables.kinematicRatio[k+1] = variables.celerities[k+1]/variables.poreVelocities[k+1];
						
						variables.darcyVelocities[k] = -variables.kappasInterface[k] * ( (variables.waterSuctions[k]-variables.waterSuctions[k-1])/geometry.spaceDeltaZ[k] +1 );
						variables.darcyVelocitiesCapillary[k] = -variables.kappasInterface[k] * (variables.waterSuctions[k]-bottomBC)/geometry.spaceDeltaZ[k];
						variables.darcyVelocitiesGravity[k] = -variables.kappasInterface[k];
						variables.poreVelocities[k] = variables.darcyVelocities[k]/interfaceConductivity.compute(variables.thetas[k-1]-rehologyParameters.thetaR[parameterID[k-1]],variables.thetas[k]-rehologyParameters.thetaR[parameterID[k]]
								,geometry.controlVolume[k-1], geometry.controlVolume[k]);
						variables.celerities[k] = -9999.0;
						variables.kinematicRatio[k] = variables.celerities[k]/variables.poreVelocities[k];
					} else {
						variables.darcyVelocities[k] = -variables.kappasInterface[k] * ( (variables.waterSuctions[k]-variables.waterSuctions[k-1])/geometry.spaceDeltaZ[k] +1 );
						variables.darcyVelocitiesCapillary[k] = -variables.kappasInterface[k] * (variables.waterSuctions[k]-bottomBC)/geometry.spaceDeltaZ[k];
						variables.darcyVelocitiesGravity[k] = -variables.kappasInterface[k];
						variables.poreVelocities[k] = variables.darcyVelocities[k]/interfaceConductivity.compute(variables.thetas[k-1]-rehologyParameters.thetaR[parameterID[k-1]],variables.thetas[k]-rehologyParameters.thetaR[parameterID[k]]
								,geometry.controlVolume[k-1], geometry.controlVolume[k]);
						variables.celerities[k] = -9999.0;
						variables.kinematicRatio[k] = variables.celerities[k]/variables.poreVelocities[k];
					}
				}else {
					variables.darcyVelocities[k] = -variables.kappasInterface[k] * ( (variables.waterSuctions[k]-variables.waterSuctions[k-1])/geometry.spaceDeltaZ[k] +1 );
					variables.darcyVelocitiesCapillary[k] = -variables.kappasInterface[k] * (variables.waterSuctions[k]-bottomBC)/geometry.spaceDeltaZ[k];
					variables.darcyVelocitiesGravity[k] = -variables.kappasInterface[k];
					variables.poreVelocities[k] = variables.darcyVelocities[k]/interfaceConductivity.compute(variables.thetas[k-1]-rehologyParameters.thetaR[parameterID[k-1]],variables.thetas[k]-rehologyParameters.thetaR[parameterID[k]]
							,geometry.controlVolume[k-1], geometry.controlVolume[k]);
					variables.celerities[k] = -9999.0;
					variables.kinematicRatio[k] = variables.celerities[k]/variables.poreVelocities[k];
				}
			}

			/*
			 * compute error
			 */
			if(this.topBCType.equalsIgnoreCase("Top Dirichlet") || this.topBCType.equalsIgnoreCase("TopDirichlet")) {
				variables.errorVolume = variables.waterVolumeNew - variables.waterVolume - timeDelta*(-variables.darcyVelocities[KMAX-1] + variables.darcyVelocities[0]) + variables.sumETs + volumeLost;
			} else {
				variables.errorVolume = variables.waterVolumeNew - variables.waterVolume - timeDelta*(topBC + variables.darcyVelocities[0]) + variables.sumETs + volumeLost;
			}
		}


		if(saveDate == 1) {
			outputToBuffer.add(variables.waterSuctions);
			outputToBuffer.add(variables.thetas);
			outputToBuffer.add(variables.volumes);
			outputToBuffer.add(variables.darcyVelocities);
			outputToBuffer.add(variables.darcyVelocitiesCapillary);
			outputToBuffer.add(variables.darcyVelocitiesGravity);
			outputToBuffer.add(variables.poreVelocities);
			outputToBuffer.add(variables.celerities);
			outputToBuffer.add(variables.kinematicRatio);
			outputToBuffer.add(variables.ETs);
			outputToBuffer.add(new double[] {variables.errorVolume});
			outputToBuffer.add(new double[] {topBC*tTimestep*1000}); // I want to have rainfall height instead of water flux
			if(bottomBCType.equalsIgnoreCase("Bottom Neumann") || bottomBCType.equalsIgnoreCase("BottomNeumann")) {
				bottomBC = bottomBC*tTimestep;
			}
			outputToBuffer.add(new double[] {bottomBC});
			outputToBuffer.add(new double[] {variables.runOff/tTimestep}); // surface runoff
			doProcessBuffer = true;
		} else {
			//			System.out.println("SaveDate = " + saveDate);
		}
		step++;

	} //// MAIN CYCLE END ////

}  /// CLOSE Richards1d ///



