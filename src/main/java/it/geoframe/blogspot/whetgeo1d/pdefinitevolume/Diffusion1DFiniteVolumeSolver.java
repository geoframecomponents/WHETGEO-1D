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

package it.geoframe.blogspot.whetgeo1d.pdefinitevolume;
import java.util.List;

import it.geoframe.blogspot.whetgeo1d.boundaryconditions.*;
import it.geoframe.blogspot.whetgeo1d.data.*;
import it.geoframe.blogspot.closureequation.equationstate.EquationState;
import it.geoframe.blogspot.numerical.newtonalgorithm.NestedNewtonThomas;
import oms3.annotations.*;

@Description("This code solve the diffusione equation."
		+ "A semi-implicit finite volume method is used to discretize the equation, and the non-linear system is solved using the nested Newton algorithm.")
@Documentation("")
@Author(name = "Niccolo' Tubini and Riccardo Rigon", contact = "tubini.niccolo@gmail.com")
@Keywords("Richards equation, numerical solver, finite volume ")
@Bibliography("")
//@Label()
//@Name()
//@Status()
@License("General Public License Version 3 (GPLv3)")

public class Diffusion1DFiniteVolumeSolver {

	@Description("Time step of integration")
	@Unit ("s")
	private double timeDelta;

	@Description("Maximun number of Newton iterations")
	final int MAXITER_NEWT = 50;

	@Description("Vector collecting the lower diagonal entries of the coefficient matrix")
	@Unit ("?")
	private double[] lowerDiagonal;

	@Description("Vector collecting the main diagonal entries of the coefficient matrix")
	@Unit ("?")
	private double[] mainDiagonal;

	@Description("Vector collecting the upper diagonal entries of the coefficient matrix")
	@Unit ("?")
	private double[] upperDiagonal;

	@Description("Right hand side vector of the scalar equation to solve")
	@Unit ("-")
	private double[] rhss;
	
	@Description("Conductivity at the cell interface i+1/2")
	@Unit ("m/s)")
	private double kP;

	@Description("Conductivity at the cell interface i-1/2")
	@Unit ("m/s")
	private double kM;

	@Description("Number of control volume for domain discetrization")
	@Unit (" ")
	private int KMAX;

	private boolean checkData = true;

	@Description("Object to perform the nested Newton algortithm")
	private NestedNewtonThomas nestedNewtonAlg;

	@Description("This list contains the objects that describes the state equations of the problem")
	private List<EquationState> equationState;

	@Description("This object compute the diagonal and right hand side entries for the uppermost cell accordingly with the prescribed top boundary condition.")
	private BoundaryCondition topBoundaryCondition;

	@Description("This object compute the diagonal and right hand side entries for the lowermost cell accordingly with the prescribed bottom boundary condition.")
	private BoundaryCondition bottomBoundaryCondition;
	
    //////////////////////////////



	public Diffusion1DFiniteVolumeSolver( BoundaryCondition topBoundaryCondition, BoundaryCondition bottomBoundaryCondition, int KMAX, int nestedNewton, double newtonTolerance, double delta,
			int MAXITER_NEWT, List<EquationState> equationState) {

		this.topBoundaryCondition = topBoundaryCondition;		
		this.bottomBoundaryCondition = bottomBoundaryCondition;	

		nestedNewtonAlg = new NestedNewtonThomas(nestedNewton, newtonTolerance, MAXITER_NEWT, KMAX, equationState, delta);

		this.KMAX = KMAX;

		lowerDiagonal = new double[KMAX];
		mainDiagonal  = new double[KMAX];
		upperDiagonal = new double[KMAX];
		rhss 		  = new double[KMAX];
		kP 		      = 0.0;
		kM	          = 0.0;

	}


	/**
	 * 
	 * @param topBC
	 * @param bottomBC
	 * @param inCurrentDate
	 * @param timeDelta
	 */
	public double[] solve(double timeDelta, double bottomBCValue, double topBCValue, int KMAX, double[] kappasInterface, double[] conservedQuantities, double[] spaceDeltaZ,
			double[] sourceSinkTerm, double[] x, double[] y, int[] parameterID, int[] rheologyID) {


		this.timeDelta = timeDelta;
		this.KMAX = KMAX;


		/* COEFFICIENT MATRIX IS BUILD BY THREE VECTORS COLLECTING ELEMENTS OF THE THREE DIAGONAL:
				   	 - lower diagonal x_(i+1)
				   	 - main diagonal  x_i
				   	 - upper diagonal x_(i-1)
				   	 - right-hand-side 
		*/
		for(int i = 1; i < KMAX-1; i++) {
			
				kP = kappasInterface[i+1];
				kM = kappasInterface[i];
				lowerDiagonal[i] = -kM*timeDelta/spaceDeltaZ[i];
				mainDiagonal[i] = kM*timeDelta/spaceDeltaZ[i] + kP*timeDelta/spaceDeltaZ[i+1];
				upperDiagonal[i] = -kP*timeDelta/spaceDeltaZ[i+1];
				rhss[i] = conservedQuantities[i] - sourceSinkTerm[i]; 
				
		}
		
		
		// i == 0			
		kP = kappasInterface[1];
		kM = kappasInterface[0];
		lowerDiagonal[0] =  bottomBoundaryCondition.lowerDiagonal(-999.0, kP, kM, spaceDeltaZ[1], spaceDeltaZ[0], timeDelta);
		mainDiagonal[0] = bottomBoundaryCondition.mainDiagonal(-999.0, kP, kM, spaceDeltaZ[1], spaceDeltaZ[0], timeDelta);
		upperDiagonal[0] = bottomBoundaryCondition.upperDiagonal(-999.0, kP, kM, spaceDeltaZ[1], spaceDeltaZ[0], timeDelta);
		rhss[0] = conservedQuantities[0] + bottomBoundaryCondition.rightHandSide(bottomBCValue, kP, kM, spaceDeltaZ[1], spaceDeltaZ[0], timeDelta) - sourceSinkTerm[0];

		// i == KMAX -1
		kP = kappasInterface[KMAX];
		kM = kappasInterface[KMAX-1];
		lowerDiagonal[KMAX-1] = topBoundaryCondition.lowerDiagonal(-999.0, kP, kM, spaceDeltaZ[KMAX], spaceDeltaZ[KMAX-1], timeDelta); 
		mainDiagonal[KMAX-1] = topBoundaryCondition.mainDiagonal(-999.0, kP, kM, spaceDeltaZ[KMAX], spaceDeltaZ[KMAX-1], timeDelta);
		upperDiagonal[KMAX-1] = topBoundaryCondition.upperDiagonal(-999.0, kP, kM, spaceDeltaZ[KMAX], spaceDeltaZ[KMAX-1], timeDelta);
		rhss[KMAX -1] = conservedQuantities[KMAX-1] + topBoundaryCondition.rightHandSide(topBCValue, kP, kM, spaceDeltaZ[KMAX], spaceDeltaZ[KMAX-1], timeDelta) - sourceSinkTerm[KMAX-1];

		 
//		if(checkData == true) {
//
//			System.out.println("Lower:");
//			for(int k=0;k<KMAX;k++){
//				System.out.println(lowerDiagonal[k]);
//			}
//
//			System.out.println("\n\nMain:");
//			for(int k=0;k<KMAX;k++){
//				System.out.println(mainDiagonal[k]);
//			}
//
//			System.out.println("\n\nUpper:");
//			for(int k=0;k<KMAX;k++){
//				System.out.println(upperDiagonal[k]);
//			}
//
//			System.out.println("\n\nRhs:");
//			for(int k=0;k<KMAX;k++){
//				System.out.println(rhss[k]);
//			}
//		}
		
		/* 
		 * NESTED NEWTON ALGORITHM /
		 */
	
		nestedNewtonAlg.set(x, y, mainDiagonal, upperDiagonal, lowerDiagonal, rhss, KMAX, parameterID, rheologyID);
		return x = nestedNewtonAlg.solver();



	} //// MAIN CYCLE END ////


} 



