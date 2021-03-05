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
import oms3.annotations.Author;
import oms3.annotations.Bibliography;
import oms3.annotations.Description;
import oms3.annotations.Documentation;
import oms3.annotations.Keywords;
import oms3.annotations.License;

@Description("This class compute all the quantities to solve the Richards' equation.")
@Documentation("")
@Author(name = "Niccolo' Tubini and Riccardo Rigon", contact = "tubini.niccolo@gmail.com")
@Keywords("Richards equation, numerical solver, finite volume ")
@Bibliography("")
//@Label()
//@Name()
//@Status()
@License("General Public License Version 3 (GPLv3)")
public class ComputeQuantitiesRichards {
	
	
	private ProblemQuantities variables;
	private Geometry geometry;
	private Parameters parameters;
	

	@Description("List containing the closure equations")
	private List<ClosureEquation> closureEquation;
	
	@Description("Factory for the closure equations")
	private SoilWaterRetentionCurveFactory soilWaterRetentionCurveFactory;
	
	@Description("List containig the objects that describes the state equations of the problem")
	private List<EquationState> equationState;

	@Description("Factory for the equation state")
	private EquationStateFactory equationStateFactory;

	@Description("Object dealing with the hydraulic conductivity model")
	private List<ConductivityEquation> hydraulicConductivity;
	private ConductivityEquationFactory conductivityEquationFactory;
	private UnsaturatedHydraulicConductivityTemperatureFactory unsaturatedHydraulicConductivityTemperatureFactory;

	@Description("This object compute the interface hydraulic conductivity accordingly with the prescribed method.")
	private InterfaceConductivity interfaceConductivity;
	private SimpleInterfaceConductivityFactory interfaceConductivityFactory;

	private String topBCType;
	private String bottomBCType;
	
//	private double tmp;
	
	public ComputeQuantitiesRichards(String[] typeClosureEquation, String[] typeEquationState, String[] typeUHCModel, String typeUHCTemperatureModel,
			String interfaceHydraulicConductivityModel, String topBCType, String bottomBCType) {
		
		variables = ProblemQuantities.getInstance();
		geometry = Geometry.getInstance();
		parameters = Parameters.getInstance();


		soilWaterRetentionCurveFactory = new SoilWaterRetentionCurveFactory();
		
		closureEquation = new ArrayList<ClosureEquation>();
//		soilWaterRetentionCurve.add(soilWaterRetentionCurveFactory.create("Water Depth"));
//		soilWaterRetentionCurve.add(soilWaterRetentionCurveFactory.create(soilHydraulicModel));
		for(int i=0; i<typeClosureEquation.length; i++) {
			closureEquation.add(soilWaterRetentionCurveFactory.create(typeClosureEquation[i]));
		}

		equationStateFactory = new EquationStateFactory();

		equationState = new ArrayList<EquationState>();
//		equationState.add(equationStateFactory.create("Water Depth", closureEquation.get(0)));
//		equationState.add(equationStateFactory.create(soilHydraulicModel, soilWaterRetentionCurve.get(1)));

		for(int i=0; i<typeEquationState.length; i++) {
			equationState.add(equationStateFactory.create(typeEquationState[i], closureEquation.get(i)));
		}


		hydraulicConductivity = new ArrayList<ConductivityEquation>();

		conductivityEquationFactory = new ConductivityEquationFactory();
		if(typeClosureEquation[0].equalsIgnoreCase("Water Depth") || typeClosureEquation[0].equalsIgnoreCase("Water Depth")) {
//			typeUHCModel[0] = typeUHCModel[variables.parameterID[variables.parameterID.length-2]];
//			hydraulicConductivity.add(conductivityEquationFactory.create(typeUHCModel[0], closureEquation.get(variables.equationStateID[variables.equationStateID.length-2])));
			typeUHCModel[0] = typeUHCModel[variables.equationStateID[variables.equationStateID.length-2]];
			hydraulicConductivity.add(conductivityEquationFactory.create(typeUHCModel[0], closureEquation.get(variables.equationStateID[variables.equationStateID.length-2])));
			for(int i=1; i<typeUHCModel.length; i++) {
				hydraulicConductivity.add(conductivityEquationFactory.create(typeUHCModel[i], closureEquation.get(i)));
			}
		} else {
			for(int i=0; i<typeUHCModel.length; i++) {
				hydraulicConductivity.add(conductivityEquationFactory.create(typeUHCModel[i], closureEquation.get(i)));
			}
		}

		unsaturatedHydraulicConductivityTemperatureFactory = new UnsaturatedHydraulicConductivityTemperatureFactory();
		if(typeClosureEquation[0].equalsIgnoreCase("Water Depth") || typeClosureEquation[0].equalsIgnoreCase("Water Depth")) {
			hydraulicConductivity.add(0, unsaturatedHydraulicConductivityTemperatureFactory.create(typeUHCTemperatureModel, closureEquation.get(variables.equationStateID[variables.equationStateID.length-2]), hydraulicConductivity.get(0)) );
			for(int i=1; i<typeUHCModel.length; i++) {
				hydraulicConductivity.add(i, unsaturatedHydraulicConductivityTemperatureFactory.create(typeUHCTemperatureModel, closureEquation.get(i), hydraulicConductivity.get(i)) );
			}
		} else {
			for(int i=0; i<typeUHCModel.length; i++) {
				hydraulicConductivity.add(i, unsaturatedHydraulicConductivityTemperatureFactory.create(typeUHCTemperatureModel, closureEquation.get(i), hydraulicConductivity.get(i)) );
			}
		}

		interfaceConductivityFactory = new SimpleInterfaceConductivityFactory();
		interfaceConductivity = interfaceConductivityFactory.createInterfaceConductivity(interfaceHydraulicConductivityModel);
		
		this.topBCType = topBCType;
		this.bottomBCType = bottomBCType;
		
	}
	
	
	public List<EquationState> getRichardsStateEquation(){
		return equationState;
	}
	
	
	public void resetRunOff() {
		variables.runOff = 0.0;
	}
	
	
	public void computeWaterVolume(int KMAX) {
		
		variables.waterVolume = 0.0;
		for(int element = 0; element < KMAX; element++) {
			variables.volumes[element] = equationState.get(variables.equationStateID[element]).equationState(variables.waterSuctions[element], variables.temperatures[element], variables.parameterID[element], element);
			variables.waterVolume += variables.volumes[element];
		}
	}
	
	public void computeWaterVolumeNew(int KMAX) {
		
		variables.waterVolumeNew = 0.0;
		for(int element = 0; element < KMAX; element++) {
			variables.volumesNew[element] = equationState.get(variables.equationStateID[element]).equationState(variables.waterSuctions[element], variables.temperatures[element], variables.parameterID[element], element);
			variables.waterVolumeNew += variables.volumesNew[element];
		}
	}
	
	public void computeThetas(int KMAX) {
		
		for(int element = 0; element < KMAX; element++) {
			variables.thetas[element] = closureEquation.get(variables.equationStateID[element]).f(variables.waterSuctions[element], variables.temperatures[element], variables.parameterID[element]);
		}
	}
	
	public void computeThetasNew(int KMAX) {
		
		for(int element = 0; element < KMAX; element++) {
			variables.thetasNew[element] = closureEquation.get(variables.equationStateID[element]).f(variables.waterSuctions[element], variables.temperatures[element], variables.parameterID[element]);
		}
	}
	
	
	public void computeXStar(int KMAX) {
		
		for(int element=0; element<KMAX; element++) {
			equationState.get(variables.equationStateID[element]).computeXStar(variables.temperatures[element], variables.parameterID[element], element);
		}
		
	}
	
	
	public void computeHydraulicConductivity(int KMAX) {
		
		for(int element = 0; element < KMAX; element++) {
			variables.kappas[element] = hydraulicConductivity.get(variables.equationStateID[element]).k(variables.waterSuctions[element], variables.temperatures[element], variables.parameterID[element], element);
		}			

//		// element == KMAX-1
//		if(topBCType.equalsIgnoreCase("Top Coupled") || this.topBCType.equalsIgnoreCase("TopCoupled")){
//			variables.kappas[KMAX-1] = hydraulicConductivity.k(variables.waterSuctions[KMAX-1], variables.temperatures[KMAX-1], variables.parameterID[KMAX-1], KMAX-1);
//		} else {
//			variables.kappas[KMAX-1] = hydraulicConductivity.k(variables.waterSuctions[KMAX-1], variables.temperatures[KMAX-1], variables.parameterID[KMAX-2], KMAX-2);
//		}
	}
	
	public void computeInterfaceHydraulicConductivity(int KMAX) {
		
		for(int k = 1; k < KMAX-1; k++) {
			variables.kappasInterface[k] = interfaceConductivity.compute(variables.kappas[k-1],variables.kappas[k], geometry.controlVolume[k-1], geometry.controlVolume[k]);
		}			
		
		// element == 0
		if(this.bottomBCType.equalsIgnoreCase("Bottom Free Drainage") || this.bottomBCType.equalsIgnoreCase("BottomFreeDrainage")){
			variables.kappasInterface[0] =  variables.kappas[0];
		} else if (this.bottomBCType.equalsIgnoreCase("Bottom Impervious") || this.bottomBCType.equalsIgnoreCase("BottomImpervious")) {
			variables.kappasInterface[0] = + 0.0;
		} else {
			variables.kappasInterface[0] = hydraulicConductivity.get(variables.equationStateID[0]).k(variables.richardsBottomBCValue, variables.temperatures[0], variables.parameterID[0], 0);
		} 
		
		// element == KMAX-1
		if(this.topBCType.equalsIgnoreCase("Top Dirichlet") || this.topBCType.equalsIgnoreCase("TopDirichlet")){
			variables.kappasInterface[KMAX-1] = interfaceConductivity.compute(variables.kappas[KMAX-2],variables.kappas[KMAX-1], geometry.controlVolume[KMAX-2], geometry.controlVolume[KMAX-1]); 
			variables.kappasInterface[KMAX] = hydraulicConductivity.get(variables.equationStateID[KMAX-1]).k(variables.richardsTopBCValue, variables.temperatures[KMAX-1], variables.parameterID[KMAX-1], KMAX-1); 
		} else if(this.topBCType.equalsIgnoreCase("Top Neumann") || this.topBCType.equalsIgnoreCase("TopNeumann")){
			variables.kappasInterface[KMAX-1] = interfaceConductivity.compute(variables.kappas[KMAX-2],variables.kappas[KMAX-1], geometry.controlVolume[KMAX-2], geometry.controlVolume[KMAX-1]);
			variables.kappasInterface[KMAX] = -9999.0;
		} else if(this.topBCType.equalsIgnoreCase("Top coupled") || this.topBCType.equalsIgnoreCase("TopCoupled")){
			variables.kappasInterface[KMAX-1] = variables.kappas[KMAX-1]; 
			variables.kappasInterface[KMAX] = -9999.0; 
		}
	}
	
	public void computeDarcyVelocities(int KMAX) {
		
		for(int k = 1; k < KMAX; k++) {
			variables.darcyVelocities[k] = -variables.kappasInterface[k] * ( (variables.waterSuctions[k]-variables.waterSuctions[k-1])/geometry.spaceDeltaZ[k] + 1 );
		}
		
		// element == 0
		if(this.bottomBCType.equalsIgnoreCase("Bottom Free Drainage") || this.bottomBCType.equalsIgnoreCase("BottomFreeDrainage")){
			variables.darcyVelocities[0] = -variables.kappasInterface[0];
		} else if (this.bottomBCType.equalsIgnoreCase("Bottom Impervious") || this.bottomBCType.equalsIgnoreCase("BottomImpervious")) {
			variables.darcyVelocities[0] = + 0.0;
		} else if (this.bottomBCType.equalsIgnoreCase("Bottom Dirichlet") || this.bottomBCType.equalsIgnoreCase("BottomDirichlet")) {
			variables.darcyVelocities[0] = -variables.kappasInterface[0] * ( (variables.waterSuctions[0]-variables.richardsBottomBCValue)/geometry.spaceDeltaZ[0] + 1 );
		} else if (this.bottomBCType.equalsIgnoreCase("Bottom Neumann") || this.bottomBCType.equalsIgnoreCase("BottomNeumann")) {
			variables.darcyVelocities[0] = variables.richardsBottomBCValue;
		}
		
		// element == KMAX-1
		if(this.topBCType.equalsIgnoreCase("Top Dirichlet") || this.topBCType.equalsIgnoreCase("TopDirichlet")){
//			variables.darcyVelocities[KMAX-1] = -variables.kappasInterface[KMAX-1] * ( (variables.waterSuctions[KMAX-1]-variables.waterSuctions[KMAX-2])/geometry.spaceDeltaZ[KMAX-1] +1 );
			variables.darcyVelocities[KMAX] = -variables.kappasInterface[KMAX] * ( (variables.richardsTopBCValue-variables.waterSuctions[KMAX-1])/geometry.spaceDeltaZ[KMAX] +1 );
		} else if (this.topBCType.equalsIgnoreCase("Top Neumann") || this.topBCType.equalsIgnoreCase("TopNeumann")){
//			variables.darcyVelocities[KMAX-1] = -variables.kappasInterface[KMAX-1] * ( (variables.waterSuctions[KMAX-1]-variables.waterSuctions[KMAX-2])/geometry.spaceDeltaZ[KMAX-1] +1 );
			variables.darcyVelocities[KMAX] = -variables.richardsTopBCValue;
		} else if (this.topBCType.equalsIgnoreCase("Top coupled") || this.topBCType.equalsIgnoreCase("TopCoupled")){
//			variables.darcyVelocities[KMAX-1] = -variables.kappasInterface[KMAX-1] * ( (variables.waterSuctions[KMAX-1]-variables.waterSuctions[KMAX-2])/geometry.spaceDeltaZ[KMAX-1] +1 );
			variables.darcyVelocities[KMAX] = -variables.richardsTopBCValue;
		}
		
	}
	
	public void computeDarcyVelocitiesCapillary(int KMAX) {
		
		for(int k = 1; k < KMAX; k++) {
			variables.darcyVelocitiesCapillary[k] = -variables.kappasInterface[k] * (variables.waterSuctions[k]-variables.waterSuctions[k-1])/geometry.spaceDeltaZ[k];
		}
		
		// element == 0
		if(this.bottomBCType.equalsIgnoreCase("Bottom Free Drainage") || this.bottomBCType.equalsIgnoreCase("BottomFreeDrainage")){
			variables.darcyVelocitiesCapillary[0] = + 0.0;
		} else if (this.bottomBCType.equalsIgnoreCase("Bottom Impervious") || this.bottomBCType.equalsIgnoreCase("BottomImpervious")) {
			variables.darcyVelocitiesCapillary[0] = + 0.0;
		} else if (this.bottomBCType.equalsIgnoreCase("Bottom Dirichlet") || this.bottomBCType.equalsIgnoreCase("BottomDirichlet")) {
			variables.darcyVelocitiesCapillary[0] = -variables.kappasInterface[0] * (variables.waterSuctions[0]-variables.richardsBottomBCValue)/geometry.spaceDeltaZ[0];
		}
		
		// element == KMAX-1
		if(this.topBCType.equalsIgnoreCase("Top Dirichlet") || this.topBCType.equalsIgnoreCase("TopDirichlet")){
//			variables.darcyVelocitiesCapillary[KMAX-1] = -variables.kappasInterface[KMAX-1] * (variables.waterSuctions[KMAX-1]-variables.waterSuctions[KMAX-2])/geometry.spaceDeltaZ[KMAX-1];
			variables.darcyVelocitiesCapillary[KMAX] = -variables.kappasInterface[KMAX] * (variables.richardsTopBCValue-variables.waterSuctions[KMAX-1])/geometry.spaceDeltaZ[KMAX];
		} else if (this.topBCType.equalsIgnoreCase("Top Neumann") || this.topBCType.equalsIgnoreCase("TopNeumann")){
//			variables.darcyVelocitiesCapillary[KMAX-1] = -variables.kappasInterface[KMAX-1] * (variables.waterSuctions[KMAX-1]-variables.waterSuctions[KMAX-2])/geometry.spaceDeltaZ[KMAX-1];
			variables.darcyVelocitiesCapillary[KMAX] = -9999.0;
		} else if (this.topBCType.equalsIgnoreCase("Top coupled") || this.topBCType.equalsIgnoreCase("TopCoupled")){
//			variables.darcyVelocitiesCapillary[KMAX-1] = -variables.kappasInterface[KMAX-1] * (variables.waterSuctions[KMAX-1]-variables.waterSuctions[KMAX-2])/geometry.spaceDeltaZ[KMAX-1];
			variables.darcyVelocitiesCapillary[KMAX] = -9999.0;
		}
		
	}
	
	public void computeDarcyVelocitiesGravity(int KMAX) {
		
		for(int k = 1; k < KMAX; k++) {
			variables.darcyVelocitiesGravity[k] = -variables.kappasInterface[k];
		}
		
		// element == 0
		if(this.bottomBCType.equalsIgnoreCase("Bottom Free Drainage") || this.bottomBCType.equalsIgnoreCase("BottomFreeDrainage")){
			variables.darcyVelocitiesGravity[0] = -variables.kappasInterface[0];
		} else if (this.bottomBCType.equalsIgnoreCase("Bottom Impervious") || this.bottomBCType.equalsIgnoreCase("BottomImpervious")) {
			variables.darcyVelocitiesGravity[0] = + 0.0;
		} else if (this.bottomBCType.equalsIgnoreCase("Bottom Dirichlet") || this.bottomBCType.equalsIgnoreCase("BottomDirichlet")) {
			variables.darcyVelocitiesGravity[0] = -variables.kappasInterface[0];
		}
		
		// element == KMAX-1
		if(this.topBCType.equalsIgnoreCase("Top Dirichlet") || this.topBCType.equalsIgnoreCase("TopDirichlet")){
//			variables.darcyVelocitiesGravity[KMAX-1] = -variables.kappasInterface[KMAX-1];
			variables.darcyVelocitiesGravity[KMAX] = -variables.kappasInterface[KMAX];
		} else if (this.topBCType.equalsIgnoreCase("Top Neumann") || this.topBCType.equalsIgnoreCase("TopNeumann")){
//			variables.darcyVelocitiesGravity[KMAX-1] = -variables.kappasInterface[KMAX-1];
			variables.darcyVelocitiesGravity[KMAX] = -9999.0;
		} else if (this.topBCType.equalsIgnoreCase("Top coupled") || this.topBCType.equalsIgnoreCase("TopCoupled")){
//			variables.darcyVelocitiesGravity[KMAX-1] = -variables.kappasInterface[KMAX-1];
			variables.darcyVelocitiesGravity[KMAX] = -9999.0;
		}	
		
	}
	
	public void computePoreVelocities(int KMAX) {
		
		for(int k = 1; k < KMAX-1; k++) {
			variables.poreVelocities[k] = variables.darcyVelocities[k]/interfaceConductivity.compute(variables.thetas[k-1]-parameters.thetaR[variables.parameterID[k-1]],variables.thetas[k]-parameters.thetaR[variables.parameterID[k]]
					,geometry.controlVolume[k-1], geometry.controlVolume[k]);
		}
		
		// element == 0
		if(this.bottomBCType.equalsIgnoreCase("Bottom Free Drainage") || this.bottomBCType.equalsIgnoreCase("BottomFreeDrainage")){
			variables.poreVelocities[0] = variables.darcyVelocities[0]/(variables.thetas[0]-parameters.thetaR[variables.parameterID[0]]);
		} else if (this.bottomBCType.equalsIgnoreCase("Bottom Impervious") || this.bottomBCType.equalsIgnoreCase("BottomImpervious")) {
			variables.poreVelocities[0] = + 0.0;
		} else if (this.bottomBCType.equalsIgnoreCase("Bottom Dirichlet") || this.bottomBCType.equalsIgnoreCase("BottomDirichlet"))  {
			variables.poreVelocities[0] = variables.darcyVelocities[0]/interfaceConductivity.compute(closureEquation.get(variables.equationStateID[0]).f(variables.richardsBottomBCValue, variables.temperatures[0], variables.parameterID[0])-parameters.thetaR[variables.parameterID[0]], variables.thetas[0]-parameters.thetaR[variables.parameterID[0]]
					,geometry.controlVolume[0], geometry.controlVolume[0]);
		} else if (this.bottomBCType.equalsIgnoreCase("Bottom Dirichlet") || this.bottomBCType.equalsIgnoreCase("BottomDirichlet"))  {
			variables.poreVelocities[0] = variables.darcyVelocities[0]/(variables.thetas[0]-parameters.thetaR[variables.parameterID[0]]);
		}
		
		// element == KMAX-1
		if(this.topBCType.equalsIgnoreCase("Top Dirichlet") || this.topBCType.equalsIgnoreCase("TopDirichlet")){
			variables.poreVelocities[KMAX-1] = variables.darcyVelocities[KMAX-1]/interfaceConductivity.compute(variables.thetas[KMAX-2]-parameters.thetaR[variables.parameterID[KMAX-2]],variables.thetas[KMAX-1]-parameters.thetaR[variables.parameterID[KMAX-1]]
					,geometry.controlVolume[KMAX-2], geometry.controlVolume[KMAX-1]);
			variables.poreVelocities[KMAX] = variables.darcyVelocities[KMAX]/interfaceConductivity.compute(closureEquation.get(variables.equationStateID[KMAX-1]).f(variables.richardsTopBCValue, variables.temperatures[KMAX-1], variables.parameterID[KMAX-1])-parameters.thetaR[variables.parameterID[KMAX-1]], variables.thetas[KMAX-1]-parameters.thetaR[variables.parameterID[KMAX-1]]
					,geometry.controlVolume[KMAX-1], geometry.controlVolume[KMAX-1]);
		} else if (this.topBCType.equalsIgnoreCase("Top Neumann") || this.topBCType.equalsIgnoreCase("TopNeumann")){
			variables.poreVelocities[KMAX-1] = variables.darcyVelocities[KMAX-1]/interfaceConductivity.compute(variables.thetas[KMAX-2]-parameters.thetaR[variables.parameterID[KMAX-2]],variables.thetas[KMAX-1]-parameters.thetaR[variables.parameterID[KMAX-1]]
					,geometry.controlVolume[KMAX-2], geometry.controlVolume[KMAX-1]);
			variables.poreVelocities[KMAX] = variables.darcyVelocities[KMAX-1]/(variables.thetas[KMAX-1]-parameters.thetaR[variables.parameterID[KMAX-1]]);
		} else if (this.topBCType.equalsIgnoreCase("Top coupled") || this.topBCType.equalsIgnoreCase("TopCoupled")){
			variables.poreVelocities[KMAX-1] = variables.darcyVelocities[KMAX-1]/interfaceConductivity.compute(closureEquation.get(variables.equationStateID[KMAX-1]).f(variables.waterSuctions[KMAX-1], variables.temperatures[KMAX-1], variables.parameterID[KMAX-1])-parameters.thetaR[variables.parameterID[KMAX-1]], variables.thetas[KMAX-2]-parameters.thetaR[variables.parameterID[KMAX-2]]
					,geometry.controlVolume[KMAX-2], geometry.controlVolume[KMAX-2]);
			variables.poreVelocities[KMAX] = -9999.0;
		}	
		
	}
	
	public void computeCelerities(int KMAX) {
		
		for(int k = 1; k < KMAX-1; k++) {
			variables.celerities[k] = -9999.0;
		}
		
		// element == 0
		if(this.bottomBCType.equalsIgnoreCase("Bottom Free Drainage") || this.bottomBCType.equalsIgnoreCase("BottomFreeDrainage")){
			variables.celerities[0] = -9999.0;
		} else if (this.bottomBCType.equalsIgnoreCase("Bottom Impervious") || this.bottomBCType.equalsIgnoreCase("BottomImpervious")) {
			variables.celerities[0] = + 0.0;
		} else if (this.bottomBCType.equalsIgnoreCase("Bottom Dirichlet") || this.bottomBCType.equalsIgnoreCase("BottomDirichlet"))  {
			variables.celerities[0] = -9999.0;
		} else if (this.bottomBCType.equalsIgnoreCase("Bottom Dirichlet") || this.bottomBCType.equalsIgnoreCase("BottomDirichlet"))  {
			variables.celerities[0] = -9999.0;
		}
		
		// element == KMAX-1
		if(this.topBCType.equalsIgnoreCase("Top Dirichlet") || this.topBCType.equalsIgnoreCase("TopDirichlet")){
			variables.celerities[KMAX-1] = -9999.0;
			variables.celerities[KMAX] = -9999.0;
		} else {
			variables.celerities[KMAX-1] = -9999.0;
			variables.celerities[KMAX] = -9999.0;
		}	
		
	}

	public void computeKinematicRatio(int KMAX) {
		
		for(int k = 0; k < KMAX+1; k++) {
			variables.kinematicRatio[k] = variables.celerities[k]/variables.poreVelocities[k];
		}

	}
	
	public void computeRunOff(int KMAX, double maxPonding) {
		
		if(this.topBCType.equalsIgnoreCase("Top Neumann") || this.topBCType.equalsIgnoreCase("TopNeumann")) {
			
			if(maxPonding>0 && variables.waterSuctions[KMAX -1]>maxPonding) {
				variables.volumeLost = (variables.waterSuctions[KMAX -1] - maxPonding);
				variables.waterSuctions[KMAX -1] = maxPonding;
				variables.runOff += variables.volumeLost;
			}
			
		}
	}
	
	public void computeError(int KMAX, double timeDelta) {

		variables.errorVolume = variables.waterVolumeNew - variables.waterVolume - timeDelta*(-variables.darcyVelocities[KMAX] + variables.darcyVelocities[0]) + variables.sumETs + variables.volumeLost;

	}
	
}
