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

package richardssolver;
import java.util.List;

import data.*;
import stateequation.*;
//import interfaceConductivity.*;
import oms3.annotations.*;
import richardsboundaryconditions.*;
import newtonalgorithm.NestedNewtonThomas;

@Description("This code solve the mixed form of Richards equation."
		+ "A semi-implicit finite volume method is used to discretize the equation, and the non-linear system is solved using the nested Newton algorithm.")
@Documentation("")
@Author(name = "Niccolo' Tubini and Riccardo Rigon", contact = "tubini.niccolo@gmail.com")
@Keywords("Richards equation, numerical solver, finite volume ")
@Bibliography("")
//@Label()
//@Name()
//@Status()
@License("General Public License Version 3 (GPLv3)")
public class PDE1DSolver {

	@Description("Time step of integration")
	@Unit ("s")
	public double timeDelta;


	//	// BOUNDARY CONDITIONS

	@Description("It is possibile to chose between 3 different kind "
			+ "of boundary condition at the top of the domain: "
			+ "- Dirichlet boundary condition --> Top Dirichlet"
			+ "- Neumann boundary condition --> Top Neumann"
			+ "- Newton's law for heat transfer boundary condition --> Top Newton")
	public String topBCType;

	@Description("It is possibile to chose among 2 different kind "
			+ "of boundary condition at the bottom of the domain: "
			+ "- Dirichlet boundary condition --> Bottom Dirichlet"
			+ "- Neumann boundary condition --> Bottom Neumann")
	public String bottomBCType;

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

	
	@Description("Vector collecting the lower diagonal entries of the coefficient matrix")
	@Unit ("?")
	double[] lowerDiagonal;

	@Description("Vector collecting the main diagonal entries of the coefficient matrix")
	@Unit ("?")
	double[] mainDiagonal;

	@Description("Vector collecting the upper diagonal entries of the coefficient matrix")
	@Unit ("?")
	double[] upperDiagonal;

	@Description("Right hand side vector of the scalar equation to solve")
	@Unit ("-")
	double[] rhss;
	
	@Description("")
	@Unit ("-")
	int[] rheologyID;

	@Description("")
	@Unit ("-")
	int[] parameterID;

	@Description("Hydraulic conductivity at the cell interface i+1/2")
	@Unit ("m/s)")
	double kP;

	@Description("Hydraulic conductivity at the cell interface i-1/2")
	@Unit ("m/s")
	double kM;

	@Description("Number of control volume for domain discetrization")
	@Unit (" ")
	int KMAX;

	@Description("Distance between two control volumes. Used to compute gradients.")
	@Unit ("m")
	double[] spaceDelta;

	boolean checkData = false;

	@Description("Object to perform the nested Newton algortithm")
	NestedNewtonThomas nestedNewtonAlg;

	@Description("This list contains the objects that describes the state equations of the problem")
	List<StateEquation> stateEquation;

	@Description("This object compute the diagonal and right hand side entries for the uppermost cell accordingly with the prescribed top boundary condition.")
	BoundaryCondition topBoundaryCondition;

	@Description("This object compute the diagonal and right hand side entries for the lowermost cell accordingly with the prescribed bottom boundary condition.")
	BoundaryCondition bottomBoundaryCondition;
	
	ProblemQuantities variables;
	Geometry geometry;

    //////////////////////////////



	public PDE1DSolver( String topBCType, String bottomBCType, int KMAX, int nestedNewton, double newtonTolerance, double delta,
			int MAXITER_NEWT, List<StateEquation>  stateEquation, int[] rheologyID, int[] parameterID) {


		SimpleBoundaryConditionFactory boundCondFactory = new SimpleBoundaryConditionFactory();
		this.topBCType = topBCType;
		this.bottomBCType = bottomBCType;
		topBoundaryCondition = boundCondFactory.createBoundaryCondition(topBCType);		
		bottomBoundaryCondition = boundCondFactory.createBoundaryCondition(bottomBCType);	

		nestedNewtonAlg = new NestedNewtonThomas(nestedNewton, newtonTolerance, MAXITER_NEWT, KMAX, stateEquation, delta, parameterID, rheologyID);

		this.KMAX = KMAX;

		lowerDiagonal = new double[KMAX];
		mainDiagonal  = new double[KMAX];
		upperDiagonal = new double[KMAX];
		rhss 		  = new double[KMAX];
		kP 		      = 0.0;
		kM	          = 0.0;

		variables = ProblemQuantities.getInstance();
		geometry = Geometry.getInstance();
		
		this.parameterID = parameterID;
		this.rheologyID = rheologyID;
		
	}


	/**
	 * 
	 * @param topBC
	 * @param bottomBC
	 * @param inCurrentDate
	 * @param timeDelta
	 */
	public void solve(double topBC, double bottomBC, double timeDelta, int KMAX) {


		this.timeDelta = timeDelta;
		this.KMAX = KMAX;


		/* COEFFICIENT MATRIX IS BUILD BY THREE VECTORS COLLECTING ELEMENTS OF THE THREE DIAGONAL:
				   	 - lower diagonal x_(i+1)
				   	 - main diagonal  x_i
				   	 - upper diagonal x_(i-1)
				   	 - right-hand-side 
		*/
		for(int i = 0; i < KMAX; i++) {
			if( i == 0 ) {
				
				kP = variables.kappasInterface[i+1];
				kM = variables.kappasInterface[i];
				lowerDiagonal[i] =  bottomBoundaryCondition.lowerDiagonal(-999.0, kP, kM, geometry.spaceDeltaZ[i+1], geometry.spaceDeltaZ[i], timeDelta);
				mainDiagonal[i] = bottomBoundaryCondition.mainDiagonal(-999.0, kP, kM, geometry.spaceDeltaZ[i+1], geometry.spaceDeltaZ[i], timeDelta);
				upperDiagonal[i] = bottomBoundaryCondition.upperDiagonal(-999.0, kP, kM, geometry.spaceDeltaZ[i+1], geometry.spaceDeltaZ[i], timeDelta);
				rhss[i] = variables.volumes[i] + bottomBoundaryCondition.rightHandSide(bottomBC, kP, kM, geometry.spaceDeltaZ[i+1], geometry.spaceDeltaZ[i], timeDelta) - variables.ETs[i];

			} else if(i == KMAX -1) {

				kP = variables.kappasInterface[i+1];//0.0;
				kM = variables.kappasInterface[i];
				lowerDiagonal[i] = topBoundaryCondition.lowerDiagonal(-999.0, kP, kM, geometry.spaceDeltaZ[i+1], geometry.spaceDeltaZ[i], timeDelta); 
				mainDiagonal[i] = topBoundaryCondition.mainDiagonal(-999.0, kP, kM, geometry.spaceDeltaZ[i+1], geometry.spaceDeltaZ[i], timeDelta);
				upperDiagonal[i] = topBoundaryCondition.upperDiagonal(-999.0, kP, kM,  geometry.spaceDeltaZ[i+1], geometry.spaceDeltaZ[i], timeDelta);
				rhss[i] = variables.volumes[i] + topBoundaryCondition.rightHandSide(topBC, kP, kM, geometry.spaceDeltaZ[i+1], geometry.spaceDeltaZ[i-1], timeDelta);

			} else {

				kP = variables.kappasInterface[i+1];
				kM = variables.kappasInterface[i];
				lowerDiagonal[i] = -kM*timeDelta/geometry.spaceDeltaZ[i];
				mainDiagonal[i] = kM*timeDelta/geometry.spaceDeltaZ[i]  + kP*timeDelta/geometry.spaceDeltaZ[i+1];
				upperDiagonal[i] = -kP*timeDelta/geometry.spaceDeltaZ[i+1];
				rhss[i] = variables.volumes[i] + timeDelta*(kP - kM) - variables.ETs[i]; 

			}

		}
		if(checkData == true) {

			System.out.println("Lower:");
			for(int k=0;k<KMAX;k++){
				System.out.println(lowerDiagonal[k]);
			}

			System.out.println("\n\nMain:");
			for(int k=0;k<KMAX;k++){
				System.out.println(mainDiagonal[k]);
			}

			System.out.println("\n\nUpper:");
			for(int k=0;k<KMAX;k++){
				System.out.println(upperDiagonal[k]);
			}

			System.out.println("\n\nRhs:");
			for(int k=0;k<KMAX;k++){
				System.out.println(rhss[k]);
			}
		}
		
		/* 
		 * NESTED NEWTON ALGORITHM /
		 */
	
		nestedNewtonAlg.set(variables.waterSuctions, variables.temperatures, variables.xStar1, mainDiagonal, upperDiagonal, lowerDiagonal, rhss, KMAX);
		variables.waterSuctions = nestedNewtonAlg.solver();



	} //// MAIN CYCLE END ////


} 



