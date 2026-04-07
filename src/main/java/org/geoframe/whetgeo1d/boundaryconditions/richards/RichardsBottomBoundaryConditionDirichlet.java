/*
 * GNU GPL v3 License
 *
 * Copyright 2017 Niccolo` Tubini
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

package org.geoframe.whetgeo1d.boundaryconditions.richards;

import org.geoframe.whetgeo1d.boundaryconditions.IBoundaryCondition;

/**
 * This class compute the element of the coefficient matrix and the right-hand side
 * when a Dirichlet boundary condition is applied at the bottom of the domain.
 * @author Niccolo' Tubini
 *
 */
public class RichardsBottomBoundaryConditionDirichlet implements IBoundaryCondition {

	@Override
	public double upperDiagonal(double bC, double kP, double kM, double spaceDeltaP, double spaceDeltaM, double timeDelta) {
		return -kP * timeDelta / spaceDeltaP;
	}
	
	
	@Override
	public double mainDiagonal(double bC, double kP, double kM, double spaceDeltaP, double spaceDeltaM, double timeDelta) {
		return kM * timeDelta / spaceDeltaM + kP * timeDelta / spaceDeltaP;
	}
	
	
	@Override
	public double lowerDiagonal(double bC, double kP, double kM, double spaceDeltaP, double spaceDeltaM, double timeDelta) {
		return 0.0;
	}

	
	@Override
	public double rightHandSide(double bC, double kP, double kM, double spaceDeltaP, double spaceDeltaM, double timeDelta) {
		return timeDelta * (kP - kM) + kM * timeDelta / spaceDeltaM * bC;
	}
}
