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

package it.geoframe.blogspot.whetgeo1d.data;

import java.util.ArrayList;
import java.util.List;

import it.geoframe.blogspot.closureequation.closureequation.ClosureEquation;
import it.geoframe.blogspot.closureequation.closureequation.Parameters;
import it.geoframe.blogspot.closureequation.closureequation.SoilWaterRetentionCurveFactory;
import it.geoframe.blogspot.closureequation.conductivitymodel.ConductivityEquation;
import it.geoframe.blogspot.closureequation.conductivitymodel.ConductivityEquationFactory;
import it.geoframe.blogspot.closureequation.conductivitymodel.UnsaturatedHydraulicConductivityTemperatureFactory;
import it.geoframe.blogspot.closureequation.equationstate.EquationState;
import it.geoframe.blogspot.closureequation.interfaceconductivity.InterfaceConductivity;
import it.geoframe.blogspot.closureequation.interfaceconductivity.SimpleInterfaceConductivityFactory;
import it.geoframe.blogspot.whetgeo1d.equationstate.EquationStateFactory;
import it.geoframe.blogspot.whetgeo1d.equationstate.SoilInternalEnergy;
import it.geoframe.blogspot.whetgeo1d.surfaceproperties.*;
import oms3.annotations.Author;
import oms3.annotations.Bibliography;
import oms3.annotations.Description;
import oms3.annotations.Documentation;
import oms3.annotations.Keywords;
import oms3.annotations.License;

@Description("This class compute all the quantities to solve the energy equation.")
@Documentation("")
@Author(name = "Niccolo' Tubini and Riccardo Rigon", contact = "tubini.niccolo@gmail.com")
@Keywords("Richards equation, numerical solver, finite volume ")
@Bibliography("")
//@Label()
//@Name()
//@Status()
@License("General Public License Version 3 (GPLv3)")


public class ComputeQuantitiesSoluteAdvectionDispersion {
	
	
	private ProblemQuantities variables;
	private Geometry geometry;
	private Parameters parameters;
	

	@Description("List containing the closure equations")
	private List<ClosureEquation> soilWaterRetentionCurve;
	
	@Description("Factory for the closure equations")
	private SoilWaterRetentionCurveFactory soilWaterRetentionCurveFactory;
	
	@Description("List containig the objects that describes the state equations of the problem")
	private List<EquationState> equationState;

	@Description("Factory for the equation state")
	private EquationStateFactory equationStateFactory;

	@Description("Object dealing with the thermal conductivity model")
	private List<ConductivityEquation> thermalConductivity;
	private ConductivityEquationFactory conductivityEquationFactory;

	/*@Description("This object compute the interface thermal conductivity accordingly with the prescribed method.")
	private InterfaceConductivity interfaceConductivity;
	private SimpleInterfaceConductivityFactory interfaceConductivityFactory;*/
	
	@Description("This object compute the interface thermal conductivity accordingly with the prescribed method.")
	private InterfaceConductivity interfaceDispersion;
	private SimpleInterfaceConductivityFactory interfaceDispersionFactory;


	private String topBCType;
	private String bottomBCType;
	
//	private double tmp;
	
	public ComputeQuantitiesSoluteAdvectionDispersion(String[] typeClosureEquation, String interfaceDispersionModel, String topBCType, String bottomBCType) {
		
		variables = ProblemQuantities.getInstance();
		geometry = Geometry.getInstance();
		parameters = Parameters.getInstance();


		soilWaterRetentionCurveFactory = new SoilWaterRetentionCurveFactory();
		
		soilWaterRetentionCurve = new ArrayList<ClosureEquation>();
		for(int i=0; i<typeClosureEquation.length; i++) {
			soilWaterRetentionCurve.add(soilWaterRetentionCurveFactory.create(typeClosureEquation[i]));
		}

		equationStateFactory = new EquationStateFactory();

		/*equationState = new ArrayList<EquationState>();
		for(int i=0; i<typeEquationState.length; i++) {
			equationState.add(equationStateFactory.create(typeEquationState[i], soilWaterRetentionCurve.get(i)));
		}*/

		/*conductivityEquationFactory = new ConductivityEquationFactory();
		thermalConductivity = new ArrayList<ConductivityEquation>();
		for(int i=0; i<typeThermalConductivity.length; i++) {
			thermalConductivity.add(conductivityEquationFactory.create(typeThermalConductivity[i], soilWaterRetentionCurve.get(i)));
		}*/
		
		interfaceDispersionFactory = new SimpleInterfaceConductivityFactory();
		interfaceDispersion = interfaceDispersionFactory.createInterfaceConductivity(interfaceDispersionModel);
		
		this.topBCType = topBCType;
		this.bottomBCType = bottomBCType;
		
	}
	
	
	public List<EquationState> getInternalEnergyStateEquation(){
		return equationState;
	}
	
	
	
	/*public void computeHeatCapacity(int KMAX) {
		
		for(int element = 0; element < KMAX; element++) {
			variables.heatCapacitys[element] = equationState.get(variables.equationStateID[element]).equationState(variables.temperatures[element], variables.waterSuctions[element], variables.parameterID[element], element);
		}
	}*/
	
	public void computeWaterVolumeConcentrations(int KMAX) {
		
		variables.waterVolumeConcentration = 0.0;
		for(int element = 0; element < KMAX; element++) {

			variables.waterVolumeConcentrations[element] = variables.volumes[element] * variables.concentrations[element];
			variables.waterVolumeConcentration += variables.waterVolumeConcentrations[element];
		}
	}
	
	/*public void computeHeatCapacityNew(int KMAX) {
		
		for(int element = 0; element < KMAX; element++) {
			variables.heatCapacitysNew[element] = equationState.get(variables.equationStateID[element]).equationState(variables.temperatures[element], variables.waterSuctions[element], variables.parameterID[element], element);
		}
	}*/
	
	public void computeSoluteSourcesSinksTerm(int KMAX) {

	variables.sumSoluteSourceSinkTerm = 0;
		
		for(int element = 0; element < KMAX; element++) {

			variables.soluteSourcesSinksTerm[element] = variables.ETs[element] * variables.concentrations[element];;
			variables.sumSoluteSourceSinkTerm = variables.sumSoluteSourceSinkTerm + variables.soluteSourcesSinksTerm[element];
		}
	}
	
	/*public void computeTransportedQuantity(int KMAX) {
		
		for(int element = 0; element < KMAX; element++) {
			variables.waterCapacityTransported[element] =  parameters.waterDensity*parameters.specificThermalCapacityWater;
		}
	}*/
	
	
public void computeWaterVolumeConcentrationsNew(int KMAX) {
		
		variables.waterVolumeConcentrationNew = 0.0;
		for(int element = 0; element < KMAX; element++) {

			variables.waterVolumeConcentrationsNew[element] = variables.volumesNew[element] * variables.concentrations[element];;
			variables.waterVolumeConcentrationNew += variables.waterVolumeConcentrationsNew[element];
		}
	}
	
	
	public void computeThetasInterface(int KMAX) {
		
		for(int k = 1; k <= KMAX-1; k++) {
			variables.thetasInterface[k] = interfaceDispersion.compute(variables.thetas[k-1],variables.thetas[k], geometry.controlVolume[k-1], geometry.controlVolume[k]);
		}			
		
		variables.thetasInterface[0] = variables.thetas[0];
		variables.thetasInterface[KMAX] = variables.thetas[KMAX-1];
		
		/*// bottom interface 
		if(this.bottomBCType.equalsIgnoreCase("Bottom Dirichlet") || this.bottomBCType.equalsIgnoreCase("BottomDirichlet")){
			variables.lambdasInterface[0] = thermalConductivity.get(variables.equationStateID[0]).k(variables.soluteBottomBCValue, variables.waterSuctions[0], variables.parameterID[0], 0);
		} else {
			variables.lambdasInterface[0] = - 9999.0;
		}
		
		if(this.topBCType.equalsIgnoreCase("Top Dirichlet") || this.bottomBCType.equalsIgnoreCase("TopDirichlet")){
			variables.lambdasInterface[KMAX] = thermalConductivity.get(variables.equationStateID[0]).k(variables.soluteBottomBCValue, variables.waterSuctions[0], variables.parameterID[0], 0);
		} else {
			variables.lambdasInterface[KMAX] = - 9999.0;
		}*/

		
	}
	public void computeTortuosityFactorsInterface(int KMAX) {
		
		for(int element = 0; element <= KMAX-1; element++) {
			
			
			//Computing tortuosity factor according to Millington and Quirk 1961.
			
			variables.tortuosityFactors[element] = Math.pow(variables.thetas[element], 2.33333333333333)/ Math.pow(parameters.thetaS[variables.parameterID[element]],2);}
		
		
			//Computing tortuosity factor at the interface.
		for(int k = 1; k <= KMAX-1; k++) {
			variables.tortuosityFactorsInterface[k] = interfaceDispersion.compute(variables.tortuosityFactors[k-1],variables.tortuosityFactors[k], geometry.controlVolume[k-1], geometry.controlVolume[k]);
		}			
		
		variables.tortuosityFactorsInterface[0] = variables.tortuosityFactors[0];
		variables.tortuosityFactorsInterface[KMAX] = variables.tortuosityFactors[KMAX-1];
		
		//for(int element = 0; element < KMAX; element++) {
			//variables.tortuosityFactorsInterface [element] = 1 ;} 

	}
	
	public void computeDispersionCoefficients(int KMAX) {
		
		
		for(int element = 1; element <= KMAX-1; element++) {
			
			
			//CALCOLO DEL COEFFICIENTE DI DISPERSIONE secondo Bear(1972) da Stumpp et all., (2012)
			//The computation of the Dispersion Coefficient is in the interface of the control volume.
			
			variables.dispersionCoefficients[element] = (parameters.longitudinalDispersivity[variables.parameterID[element]] * Math.abs(variables.darcyVelocities[element]))/variables.thetasInterface[element] + parameters.molecularDiffusion[variables.parameterID[element]] * variables.tortuosityFactorsInterface[element];
			if (variables.thetasInterface[element]==0) {variables.dispersionCoefficients[element] = 0;}
			if (Double.isNaN(variables.dispersionCoefficients[element])) {variables.dispersionCoefficients[element] = 0;}	
			//System.out.println("variables.dispersionCoefficients is  = "+ variables.dispersionCoefficients[element]);
		}	
		
		variables.dispersionCoefficients[0] = variables.dispersionCoefficients[1];
		variables.dispersionCoefficients[KMAX] = variables.dispersionCoefficients[KMAX-1];
	}

public void computeDispersionFactors(int KMAX) {
		
		for(int element = 0; element <= KMAX; element++) {
			
			variables.dispersionFactors[element] = variables.dispersionCoefficients[element] * variables.thetasInterface[element];
		}			

	}
	
	/*public void computeInterfaceThermalConductivity(int KMAX) {
		
		for(int k = 1; k <= KMAX-1; k++) {
			variables.lambdasInterface[k] = interfaceConductivity.compute(variables.lambdas[k-1],variables.lambdas[k], geometry.controlVolume[k-1], geometry.controlVolume[k]);
		}			
		
		// bottom interface 
		if(this.bottomBCType.equalsIgnoreCase("Bottom Dirichlet") || this.bottomBCType.equalsIgnoreCase("BottomDirichlet")){
			variables.lambdasInterface[0] = thermalConductivity.get(variables.equationStateID[0]).k(variables.internalEnergyBottomBCValue, variables.waterSuctions[0], variables.parameterID[0], 0);
		} else {
			variables.lambdasInterface[0] = - 9999.0;
		}
		
		if(this.topBCType.equalsIgnoreCase("Top Dirichlet") || this.bottomBCType.equalsIgnoreCase("TopDirichlet")){
			variables.lambdasInterface[KMAX] = thermalConductivity.get(variables.equationStateID[0]).k(variables.internalEnergyBottomBCValue, variables.waterSuctions[0], variables.parameterID[0], 0);
		} else {
			variables.lambdasInterface[KMAX] = - 9999.0;
		}

		
	}*/
	

	
	/*public void computeConductionHeatFlux(int KMAX) {
		
		for(int k = 1; k <= KMAX-1; k++) {
			variables.conductionHeatFluxs[k] = -variables.lambdasInterface[k] * (variables.temperatures[k]-variables.temperatures[k-1])/geometry.spaceDeltaZ[k];
		}
		
		// bottom interface
		if (this.bottomBCType.equalsIgnoreCase("Bottom Dirichlet") || this.bottomBCType.equalsIgnoreCase("BottomDirichlet")) {
			variables.conductionHeatFluxs[0] = -variables.lambdasInterface[0] * (variables.temperatures[0]-variables.internalEnergyBottomBCValue)/geometry.spaceDeltaZ[0];
		} else {
			variables.conductionHeatFluxs[0] = -variables.internalEnergyBottomBCValue;
		}
		
		// top interface
		if (this.topBCType.equalsIgnoreCase("Top Dirichlet") || this.topBCType.equalsIgnoreCase("TopDirichlet")) {
			variables.conductionHeatFluxs[KMAX] = -variables.lambdasInterface[KMAX] * (variables.internalEnergyTopBCValue-variables.temperatures[KMAX-1])/geometry.spaceDeltaZ[KMAX];
		} else {
			variables.conductionHeatFluxs[KMAX] = -variables.internalEnergyTopBCValue;
		}
		
	}*/
	
public void computeDispersionSoluteFluxes(int KMAX) {
		
		for(int k = 1; k <= KMAX-1; k++) {
			variables.dispersionSoluteFluxes[k] = -variables.dispersionFactors[k] * (variables.concentrations[k]-variables.concentrations[k-1])/geometry.spaceDeltaZ[k+1];
		}
		
		// bottom interface
		if (this.bottomBCType.equalsIgnoreCase("Bottom Dirichlet") || this.bottomBCType.equalsIgnoreCase("BottomDirichlet")) {
			variables.dispersionSoluteFluxes[0] = -variables.dispersionFactors[0] * (variables.concentrations[0]-variables.soluteBottomBCValue)/geometry.spaceDeltaZ[0];
		} else {
			variables.dispersionSoluteFluxes[0] = -variables.soluteBottomBCValue;
		}
		
		// top interface
		if (this.topBCType.equalsIgnoreCase("Top Dirichlet") || this.topBCType.equalsIgnoreCase("TopDirichlet")) {
			variables.dispersionSoluteFluxes[KMAX] = -variables.dispersionFactors[KMAX] * (variables.soluteTopBCValue-variables.concentrations[KMAX-1])/geometry.spaceDeltaZ[KMAX];
 		} else {
			variables.dispersionSoluteFluxes[KMAX] = -variables.soluteTopBCValue;
		}
		
	}
	
	/*public void computeAdvectionHeatFlux(int KMAX) {
		
		for(int k = 1; k <= KMAX-1; k++) {
			variables.advectionHeatFluxs[k] = variables.waterCapacityTransported[k]*( 0.5*variables.darcyVelocities[k]*(variables.temperatures[k]-parameters.referenceTemperatureInternalEnergy+variables.temperatures[k-1]-parameters.referenceTemperatureInternalEnergy)
					- 0.5*Math.abs(variables.darcyVelocities[k])*(variables.temperatures[k]-parameters.referenceTemperatureInternalEnergy-variables.temperatures[k-1]+parameters.referenceTemperatureInternalEnergy));
		}
		
		// bottom interface
		variables.advectionHeatFluxs[0] = variables.waterCapacityTransported[0]*( 0.5*variables.darcyVelocities[0]*(variables.temperatures[0]-parameters.referenceTemperatureInternalEnergy + variables.internalEnergyBottomBCValue-parameters.referenceTemperatureInternalEnergy)
				- 0.5*Math.abs(variables.darcyVelocities[0])*(variables.temperatures[0]-parameters.referenceTemperatureInternalEnergy - variables.internalEnergyBottomBCValue+parameters.referenceTemperatureInternalEnergy));

//		variables.advectionHeatFluxs[0] = variables.waterCapacityTransported[0]*( 0.5*variables.darcyVelocities[0]*(variables.temperatures[0] + variables.internalEnergyBottomBCValue)
//				- 0.5*Math.abs(variables.darcyVelocities[0])*(variables.temperatures[0] - variables.internalEnergyBottomBCValue));
		/*
		 * FIXME: check the case for Neumann boundary condition: T_BC = Flux*DeltaZ/lambda_0 - T_0
		 */
		// top interface
/*		variables.advectionHeatFluxs[KMAX] = variables.waterCapacityTransported[KMAX-1]*( 0.5*variables.darcyVelocities[KMAX]*(variables.internalEnergyTopBCValue-parameters.referenceTemperatureInternalEnergy + variables.temperatures[KMAX-1]-parameters.referenceTemperatureInternalEnergy)
				- 0.5*Math.abs(variables.darcyVelocities[KMAX])*(variables.internalEnergyTopBCValue-parameters.referenceTemperatureInternalEnergy - variables.temperatures[KMAX-1]+parameters.referenceTemperatureInternalEnergy));

//		variables.advectionHeatFluxs[KMAX] = variables.waterCapacityTransported[KMAX-1]*( 0.5*variables.darcyVelocities[KMAX]*(variables.internalEnergyTopBCValue + variables.temperatures[KMAX-1])
//				- 0.5*Math.abs(variables.darcyVelocities[KMAX])*(variables.internalEnergyTopBCValue - variables.temperatures[KMAX-1]));
		
	}*/
	
public void computeAdvectionSoluteFluxes(int KMAX) {
		
		for(int k = 1; k <= KMAX-1; k++) {
			variables.advectionSoluteFluxes[k] = 0.5*variables.darcyVelocities[k]*(variables.concentrations[k]+variables.concentrations[k-1])-0.5*Math.abs(variables.darcyVelocities[k])*(variables.concentrations[k]-variables.concentrations[k-1]);
		}
		
		// bottom interface
		variables.advectionSoluteFluxes[0] = 0.5*variables.darcyVelocities[0]*(variables.concentrations[0]+variables.soluteBottomBCValue)-0.5*Math.abs(variables.darcyVelocities[0])*(variables.concentrations[0]- variables.soluteBottomBCValue); 

		/*
		 * FIXME: check the case for Neumann boundary condition: T_BC = Flux*DeltaZ/lambda_0 - T_0
		 */
		
		// top interface
		variables.advectionSoluteFluxes[KMAX] = 0.5*variables.darcyVelocities[KMAX]*(variables.soluteTopBCValue+ variables.concentrations[KMAX-1])-0.5*Math.abs(variables.darcyVelocities[KMAX])*(variables.soluteTopBCValue - variables.concentrations[KMAX-1]);
	
	}

	public void computeSoluteFluxes(int KMAX) {
		
		for(int k = 0; k <= KMAX; k++) {
			variables.soluteFluxes[k] = variables.dispersionSoluteFluxes[k]+variables.advectionSoluteFluxes[k];
		}
		
	}
	
	public void computeAverageSoluteConcentration(int KMAX) {
		
		variables.averageSoluteConcentration=0;
		
		for(int k = 0; k < KMAX-1; k++) {
			variables.averageSoluteConcentration += variables.concentrations[k];}
		variables.averageSoluteConcentration = variables.averageSoluteConcentration/(KMAX-1);
		}
	
	

	
	public void computeAverageWaterVolumeSoluteConcentration(int KMAX) {
		
		variables.averageWaterVolumeSoluteConcentration=0;
		
		for(int k = 0; k < KMAX-1; k++) {
			
			variables.averageWaterVolumeSoluteConcentration += variables.concentrations[k]*variables.volumes[k];}
		variables.averageWaterVolumeSoluteConcentration = variables.averageWaterVolumeSoluteConcentration/(KMAX-1);
		}
	

	/*public void computeError(int KMAX, double timeDelta) {
		variables.errorInternalEnergy = variables.internalEnergyNew - variables.internalEnergy - timeDelta*(-variables.conductionHeatFluxs[KMAX]-variables.advectionHeatFluxs[KMAX] + variables.conductionHeatFluxs[0] + variables.advectionHeatFluxs[0]);

	}*/
	
	public void computeError(int KMAX, double timeDelta) {
		variables.errorWaterVolumeConcentration = variables.waterVolumeConcentrationNew - variables.waterVolumeConcentration - timeDelta*(-variables.dispersionSoluteFluxes[KMAX]-variables.advectionSoluteFluxes[KMAX] + variables.dispersionSoluteFluxes[0] + variables.advectionSoluteFluxes[0]) + variables.sumSoluteSourceSinkTerm;
		//variables.errorWaterVolumeConcentration = variables.waterVolumeConcentrationNew - variables.waterVolumeConcentration - timeDelta*(-variables.dispersionSoluteFluxes[KMAX] + variables.dispersionSoluteFluxes[0]);
	}
	
	public void computeTimeVariationWaterVolumesConcentration(int KMAX, double timeDelta) {
		
		for(int k = 0; k <= KMAX-1; k++) {
			variables.timeVariationWaterVolumesConcentration[k] = -timeDelta*(variables.dispersionSoluteFluxes[k+1]+variables.advectionSoluteFluxes[k+1] - variables.dispersionSoluteFluxes[k] - variables.advectionSoluteFluxes[k]) + variables.soluteSourcesSinksTerm[k];
			//variables.waterVolumesConcentrationNew - variables.waterVolumesConcentration = -timeDelta*(variables.dispersionSoluteFluxes[k+1]+variables.advectionSoluteFluxes[k+1] + variables.dispersionSoluteFluxes[k] + variables.advectionSoluteFluxes[k])+variables.soluteSourcesSinksTerm[k];	
		}
	
	}
}




