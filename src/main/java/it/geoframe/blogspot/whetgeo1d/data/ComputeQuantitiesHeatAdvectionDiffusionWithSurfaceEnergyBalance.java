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


public class ComputeQuantitiesHeatAdvectionDiffusionWithSurfaceEnergyBalance {
	
	
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

	@Description("This object compute the interface thermal conductivity accordingly with the prescribed method.")
	private InterfaceConductivity interfaceConductivity;
	private SimpleInterfaceConductivityFactory interfaceConductivityFactory;
	
	@Description("This object evaluate the surface albedo.")
	private SurfaceAlbedo surfaceAlbedoModel;
	private SurfaceAlbedoFactory surfaceAlbedoFactory;
	
	@Description("This object evaluate the surface emissivity.")
	private SurfaceEmissivity surfaceEmissivityModel;
	private SurfaceEmissivityFactory surfaceEmissivityFactory;
	
	@Description("This object evaluate the surface aereodynamic resistance.")
	private SurfaceAereodynamicResistance surfaceAereodynamicResistanceModel;
	private SurfaceAereodynamicResistanceFactory surfaceAereodynamicResistanceFactory;


	private String topBCType;
	private String bottomBCType;
	
//	private double tmp;
	
	public ComputeQuantitiesHeatAdvectionDiffusionWithSurfaceEnergyBalance(String[] typeClosureEquation, String[] typeEquationState, String[] typeThermalConductivity,
			String interfaceHydraulicConductivityModel, String bottomBCType, String surfaceAlbedoType, String surfaceEmissivityType,
			String surfaceAereodynamicResistanceType) {
		
		variables = ProblemQuantities.getInstance();
		geometry = Geometry.getInstance();
		parameters = Parameters.getInstance();


		soilWaterRetentionCurveFactory = new SoilWaterRetentionCurveFactory();
		
		soilWaterRetentionCurve = new ArrayList<ClosureEquation>();
		for(int i=0; i<typeClosureEquation.length; i++) {
			soilWaterRetentionCurve.add(soilWaterRetentionCurveFactory.create(typeClosureEquation[i]));
		}

		equationStateFactory = new EquationStateFactory();

		equationState = new ArrayList<EquationState>();
		for(int i=0; i<typeEquationState.length; i++) {
			equationState.add(equationStateFactory.create(typeEquationState[i], soilWaterRetentionCurve.get(i)));
		}

		conductivityEquationFactory = new ConductivityEquationFactory();
		thermalConductivity = new ArrayList<ConductivityEquation>();
		for(int i=0; i<typeThermalConductivity.length; i++) {
			thermalConductivity.add(conductivityEquationFactory.create(typeThermalConductivity[i], soilWaterRetentionCurve.get(i)));
		}
		
		interfaceConductivityFactory = new SimpleInterfaceConductivityFactory();
		interfaceConductivity = interfaceConductivityFactory.createInterfaceConductivity(interfaceHydraulicConductivityModel);
		
		surfaceAlbedoFactory = new SurfaceAlbedoFactory();
		surfaceAlbedoModel = surfaceAlbedoFactory.create(surfaceAlbedoType);
		
		surfaceEmissivityFactory = new SurfaceEmissivityFactory();
		surfaceEmissivityModel = surfaceEmissivityFactory.create(surfaceEmissivityType);
		
		surfaceAereodynamicResistanceFactory = new SurfaceAereodynamicResistanceFactory();
		surfaceAereodynamicResistanceModel = surfaceAereodynamicResistanceFactory.create(surfaceAereodynamicResistanceType);
		
//		this.topBCType = topBCType;
		this.bottomBCType = bottomBCType;
		
	}
	
	
	public List<EquationState> getInternalEnergyStateEquation(){
		return equationState;
	}
	
	
	
	public void computeHeatCapacity(int KMAX) {
		
		for(int element = 0; element < KMAX; element++) {
			variables.heatCapacitys[element] = equationState.get(variables.equationStateID[element]).equationState(variables.temperatures[element], variables.waterSuctions[element], variables.parameterID[element], element);
//			System.out.println("\t\t\tHeatCapacitys "+element+" "+variables.heatCapacitys[element]);
		}
	}
	
	public void computeInternalEnergy(int KMAX) {
		
		variables.internalEnergy = 0.0;
		for(int element = 0; element < KMAX; element++) {
//			variables.internalEnergys[element] = equationState.get(variables.equationStateID[element]).equationState(variables.temperatures[element], variables.waterSuctions[element], variables.parameterID[element], element);
			variables.internalEnergys[element] = variables.heatCapacitys[element]*(variables.temperatures[element]-parameters.referenceTemperatureInternalEnergy);
			variables.internalEnergy += variables.internalEnergys[element];
		}
	}
	
	public void computeHeatCapacityNew(int KMAX) {
		
		for(int element = 0; element < KMAX; element++) {
			variables.heatCapacitysNew[element] = equationState.get(variables.equationStateID[element]).equationState(variables.temperatures[element], variables.waterSuctions[element], variables.parameterID[element], element);
		}
	}
	
	public void computeInternalEnergyNew(int KMAX) {
		
		variables.internalEnergyNew = 0.0;
		for(int element = 0; element < KMAX; element++) {
//			variables.internalEnergys[element] = equationState.get(variables.equationStateID[element]).equationState(variables.temperatures[element], variables.waterSuctions[element], variables.parameterID[element], element);
			variables.internalEnergysNew[element] = variables.heatCapacitysNew[element]*(variables.temperatures[element]-parameters.referenceTemperatureInternalEnergy);
			variables.internalEnergyNew += variables.internalEnergysNew[element];
		}
	}
	
	public void computeTransportedQuantity(int KMAX) {
		
		for(int element = 0; element < KMAX; element++) {
			variables.waterCapacityTransported[element] =  parameters.waterDensity*parameters.specificThermalCapacityWater;
		}
	}
	
	public void computeHeatSourcesSinksTerm(int KMAX) {

	variables.sumSoluteSourceSinkTerm = 0;
		
		for(int element = 0; element < KMAX; element++) {

//			variables.heatSourcesSinksTerm[element] = variables.ETs[element] * parameters.waterDensity*parameters.specificThermalCapacityWater*(variables.temperatures[element]-parameters.referenceTemperatureInternalEnergy);
			variables.heatSourcesSinksTerm[element] = variables.ETs[element] * 2.45*Math.pow(10,6);
			variables.sumHeatSourceSinkTerm = variables.sumHeatSourceSinkTerm + variables.heatSourcesSinksTerm[element];
		}
	}
	
	public void computeThermalConductivity(int KMAX) {
		
		for(int element = 0; element < KMAX; element++) {
			variables.lambdas[element] = thermalConductivity.get(variables.equationStateID[element]).k(variables.temperatures[element], variables.waterSuctions[element], variables.parameterID[element], element);
		}			

	}
	
	public void computeInterfaceThermalConductivity(int KMAX) {
		
		for(int k = 1; k <= KMAX-1; k++) {
			variables.lambdasInterface[k] = interfaceConductivity.compute(variables.lambdas[k-1],variables.lambdas[k], geometry.controlVolume[k-1], geometry.controlVolume[k]);
		}			
		
		// bottom interface 
		if(this.bottomBCType.equalsIgnoreCase("Bottom Dirichlet") || this.bottomBCType.equalsIgnoreCase("BottomDirichlet")){
			variables.lambdasInterface[0] = thermalConductivity.get(variables.equationStateID[0]).k(variables.internalEnergyBottomBCValue, variables.waterSuctions[0], variables.parameterID[0], 0);
		} else {
			variables.lambdasInterface[0] = - 9999.0;
		}
		
//		if(this.topBCType.equalsIgnoreCase("Top Dirichlet") || this.bottomBCType.equalsIgnoreCase("TopDirichlet")){
//			variables.lambdasInterface[KMAX] = thermalConductivity.get(variables.equationStateID[0]).k(variables.internalEnergyBottomBCValue, variables.waterSuctions[0], variables.parameterID[0], 0);
//		} else {
		variables.lambdasInterface[KMAX] = - 9999.0;
//		}

		
	}
	
	public void computeConductionHeatFlux(int KMAX) {
		
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
//		if (this.topBCType.equalsIgnoreCase("Top Dirichlet") || this.topBCType.equalsIgnoreCase("TopDirichlet")) {
//			variables.conductionHeatFluxs[KMAX] = -variables.lambdasInterface[KMAX] * (variables.internalEnergyTopBCValue-variables.temperatures[KMAX-1])/geometry.spaceDeltaZ[KMAX];
//		} else {
//			variables.conductionHeatFluxs[KMAX] = -variables.internalEnergyTopBCValue;
//		}
		variables.conductionHeatFluxs[KMAX] = (1-variables.surfaceAlbedo)*variables.shortWaveIn + variables.surfaceEmissivity*variables.longWaveIn - variables.linearizedLongWaveOut*variables.temperatures[KMAX-1] 
				+ variables.sensibleHeatFlux;
		
	}
	
	public void computeAdvectionHeatFlux(int KMAX) {
		
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
		variables.advectionHeatFluxs[KMAX] = variables.waterCapacityTransported[KMAX-1]*( 0.5*variables.darcyVelocities[KMAX]*(variables.rainT-parameters.referenceTemperatureInternalEnergy + variables.temperatures[KMAX-1]-parameters.referenceTemperatureInternalEnergy)
				- 0.5*Math.abs(variables.darcyVelocities[KMAX])*(variables.rainT-parameters.referenceTemperatureInternalEnergy - variables.temperatures[KMAX-1]+parameters.referenceTemperatureInternalEnergy));

//		variables.advectionHeatFluxs[KMAX] = variables.waterCapacityTransported[KMAX-1]*( 0.5*variables.darcyVelocities[KMAX]*(variables.internalEnergyTopBCValue + variables.temperatures[KMAX-1])
//				- 0.5*Math.abs(variables.darcyVelocities[KMAX])*(variables.internalEnergyTopBCValue - variables.temperatures[KMAX-1]));
		
	}

	public void computeHeatFlux(int KMAX) {
		
		for(int k = 0; k <= KMAX; k++) {
			variables.heatFluxs[k] = variables.conductionHeatFluxs[k]+variables.advectionHeatFluxs[k];
		}
		
	}
	
	public void computeLinearizedLongWaveOut(int KMAX) {
		
		variables.linearizedLongWaveOut = variables.constantStefanBoltzmann*variables.surfaceEmissivity*Math.pow(variables.temperatures[KMAX-1],3);
	}
	
	public void computeSensibleHeatCoefficient() {
		
		variables.sensibleHeatCoefficient = variables.airDensity*variables.airHeatCapacity/surfaceAereodynamicResistanceModel.evaluate(variables.referenceHeight, variables.surfaceRoughness, variables.windVelocity, variables.surfaceZeroHeightDisplacement);
	}
	
	public void computeSurfaceAlbedo(int KMAX) {
		
		variables.surfaceAlbedo = surfaceAlbedoModel.evaluate(variables.thetas[KMAX-1], variables.surfaceAlbedo);
	}
	
	public void computeSurfaceEmissivity(int KMAX) {
		
		variables.surfaceEmissivity = surfaceEmissivityModel.evaluate(variables.thetas[KMAX-1], variables.surfaceEmissivity);
	}
	
	public void computeLongWaveOut(int KMAX) {
		
		variables.longWaveOut = -(variables.surfaceEmissivity-1)*variables.longWaveIn - variables.linearizedLongWaveOut*variables.temperatures[KMAX-1];
	}
	
	public void computeShortWaveOut() {
		
		variables.shortWaveOut = -variables.surfaceAlbedo*variables.shortWaveIn;
	}
	
	public void computeSensibleHeatFlux(int KMAX) {
		
		variables.sensibleHeatFlux = variables.sensibleHeatCoefficient*(variables.airT-variables.temperatures[KMAX-1]);
	}
	
	
	/*
	 * To be Fixed
	 */

	public void computeError(int KMAX, double timeDelta) {
		variables.errorInternalEnergy = variables.internalEnergyNew - variables.internalEnergy - timeDelta*(-variables.conductionHeatFluxs[KMAX]-variables.advectionHeatFluxs[KMAX] + variables.conductionHeatFluxs[0] + variables.advectionHeatFluxs[0]);

	}
	
}
