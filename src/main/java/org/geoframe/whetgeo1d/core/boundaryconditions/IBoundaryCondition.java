/*
 * GNU GPL v3 License
 *
 * Copyright 2019  Niccolo` Tubini
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
package org.geoframe.whetgeo1d.core.boundaryconditions;

import org.geoframe.whetgeo1d.core.boundaryconditions.diffusion.DiffusionBottomBoundaryConditionDirichlet;
import org.geoframe.whetgeo1d.core.boundaryconditions.diffusion.DiffusionBottomBoundaryConditionNeumann;
import org.geoframe.whetgeo1d.core.boundaryconditions.diffusion.DiffusionBottomBoundaryConditionNoGradient;
import org.geoframe.whetgeo1d.core.boundaryconditions.diffusion.DiffusionTopBoundaryConditionDirichlet;
import org.geoframe.whetgeo1d.core.boundaryconditions.diffusion.DiffusionTopBoundaryConditionNeumann;
import org.geoframe.whetgeo1d.core.boundaryconditions.richards.RichardsBottomBoundaryConditionDirichlet;
import org.geoframe.whetgeo1d.core.boundaryconditions.richards.RichardsBottomBoundaryConditionFreeDrainage;
import org.geoframe.whetgeo1d.core.boundaryconditions.richards.RichardsBottomBoundaryConditionImpervious;
import org.geoframe.whetgeo1d.core.boundaryconditions.richards.RichardsBottomBoundaryConditionNeumann;
import org.geoframe.whetgeo1d.core.boundaryconditions.richards.RichardsTopBoundaryConditionDirichlet;
import org.geoframe.whetgeo1d.core.boundaryconditions.richards.RichardsTopBoundaryConditionNeumann;
import org.geoframe.whetgeo1d.utils.EnumUtils;

/**
 * The boundary condition interface.
 * 
 * For the 1D problem the matrix is tridiagonal and thus it is necessary to
 * compute and store only the main, the lower, and the upper diagonal of the
 * matrix.
 * 
 * @author Niccolo' Tubini
 *
 */
public interface IBoundaryCondition {
	enum DiffusionBoundaryConditionType {
		TOP_DIRICHLET, TOP_NEUMANN, BOTTOM_DIRICHLET, BOTTOM_NEUMANN, BOTTOM_NO_GRADIENT;
		
		public static DiffusionBoundaryConditionType fromString(String value) {
	        return EnumUtils.fromString(DiffusionBoundaryConditionType.class, value);
	    }
	}

	enum RichardsBoundaryConditionType {
		TOP_DIRICHLET, TOP_NEUMANN, TOP_COUPLED, BOTTOM_DIRICHLET, BOTTOM_IMPERVIOUS, BOTTOM_NEUMANN,
		BOTTOM_FREE_DRAINAGE, BOTTOM_SEEPAGE;
		
		public static RichardsBoundaryConditionType fromString(String value) {
	        return EnumUtils.fromString(RichardsBoundaryConditionType.class, value);
	    }
	}

	/**
	 * 
	 * This method computes the upper diagonal term of the coefficient matrix T of
	 * the system
	 * 
	 * @param bC          value of the boundary condition
	 * @param kP          value of the conductivity at the interface of volumes (i)
	 *                    and (i+1)
	 * @param kM          value of the conductivity at the interface of volumes (i)
	 *                    and (i-1)
	 * @param spaceDeltaP distance between nodes of volumes (i) and (i+1)
	 * @param spaceDeltaM distance between nodes of volumes (i) and (i-1)
	 * @param tTimestep   time step of the simulation
	 * @return
	 */
	double upperDiagonal(double bC, double kP, double kM, double spaceDeltaP, double spaceDeltaM, double timeDelta);

	/**
	 * 
	 * This method computes the main diagonal term of the coefficient matrix T of
	 * the system
	 * 
	 * @param bC          value of the boundary condition
	 * @param kP          value of the conductivity at the interface of volumes (i)
	 *                    and (i+1)
	 * @param kM          value of the conductivity at the interface of volumes (i)
	 *                    and (i-1)
	 * @param spaceDeltaP distance between nodes of volumes (i) and (i+1)
	 * @param spaceDeltaM distance between nodes of volumes (i) and (i-1)
	 * @param tTimestep   time step of the simulation
	 * @return
	 */
	double mainDiagonal(double bC, double kP, double kM, double spaceDeltaP, double spaceDeltaM, double timeDelta);

	/**
	 *
	 * This method computes the lower diagonal term of the coefficient matrix T of
	 * the system
	 * 
	 * @param bC          value of the boundary condition
	 * @param kP          value of the conductivity at the interface of volumes (i)
	 *                    and (i+1)
	 * @param kM          value of the conductivity at the interface of volumes (i)
	 *                    and (i-1)
	 * @param spaceDeltaP distance between nodes of volumes (i) and (i+1)
	 * @param spaceDeltaM distance between nodes of volumes (i) and (i-1)
	 * @param tTimestep   time step of the simulation
	 * @return
	 */
	double lowerDiagonal(double bC, double kP, double kM, double spaceDeltaP, double spaceDeltaM, double timeDelta);

	/**
	 *
	 * This method computes the right-hand side term of the system
	 * 
	 * @param bC          value of the boundary condition
	 * @param kP          value of the conductivity at the interface of volumes (i)
	 *                    and (i+1)
	 * @param kM          value of the conductivity at the interface of volumes (i)
	 *                    and (i-1)
	 * @param spaceDeltaP distance between nodes of volumes (i) and (i+1)
	 * @param spaceDeltaM distance between nodes of volumes (i) and (i-1)
	 * @param tTimestep   time step of the simulation
	 * @return
	 */
	double rightHandSide(double bC, double kP, double kM, double spaceDeltaP, double spaceDeltaM, double timeDelta);

	/**
	 * Factory method to create the diffusion boundary condition object based on the
	 * type of boundary condition.
	 * 
	 * @param type type of boundary condition to create.
	 * @return the created boundary condition object.
	 */
	static IBoundaryCondition createDiffusionBoundaryCondition(DiffusionBoundaryConditionType type) {
		IBoundaryCondition boundaryCondition = null;
		if (type == DiffusionBoundaryConditionType.TOP_DIRICHLET) {
			boundaryCondition = new DiffusionTopBoundaryConditionDirichlet();
		} else if (type == DiffusionBoundaryConditionType.TOP_NEUMANN) {
			boundaryCondition = new DiffusionTopBoundaryConditionNeumann();
		} else if (type == DiffusionBoundaryConditionType.BOTTOM_DIRICHLET) {
			boundaryCondition = new DiffusionBottomBoundaryConditionDirichlet();
		} else if (type == DiffusionBoundaryConditionType.BOTTOM_NEUMANN) {
			boundaryCondition = new DiffusionBottomBoundaryConditionNeumann();
		} else if (type == DiffusionBoundaryConditionType.BOTTOM_NO_GRADIENT) {
			boundaryCondition = new DiffusionBottomBoundaryConditionNoGradient();
		} else {
			throw new IllegalArgumentException("Boundary condition type not recognized");
		}
		return boundaryCondition;
	}

	/**
	 * Factory method to create the Richards boundary condition object based on the type
	 * of boundary condition.
	 * 
	 * @param type type of boundary condition to create.
	 * @return the created boundary condition object.
	 */
	static IBoundaryCondition createRichardsSimpleBoundaryCondition(RichardsBoundaryConditionType type) {
		IBoundaryCondition boundaryCondition = null;
		if (type == RichardsBoundaryConditionType.TOP_DIRICHLET) {
			boundaryCondition = new RichardsTopBoundaryConditionDirichlet();
		} else if (type == RichardsBoundaryConditionType.TOP_NEUMANN
				|| type == RichardsBoundaryConditionType.TOP_COUPLED) {
			boundaryCondition = new RichardsTopBoundaryConditionNeumann();
		} else if (type == RichardsBoundaryConditionType.BOTTOM_DIRICHLET) {
			boundaryCondition = new RichardsBottomBoundaryConditionDirichlet();
		} else if (type == RichardsBoundaryConditionType.BOTTOM_IMPERVIOUS) {
			boundaryCondition = new RichardsBottomBoundaryConditionImpervious();
		} else if (type == RichardsBoundaryConditionType.BOTTOM_NEUMANN) {
			boundaryCondition = new RichardsBottomBoundaryConditionNeumann();
		} else if (type == RichardsBoundaryConditionType.BOTTOM_FREE_DRAINAGE
				|| type == RichardsBoundaryConditionType.BOTTOM_SEEPAGE) {
			boundaryCondition = new RichardsBottomBoundaryConditionFreeDrainage();
		} else {
			throw new IllegalArgumentException("Boundary condition type not recognized");
		}
		return boundaryCondition;
	}
}
