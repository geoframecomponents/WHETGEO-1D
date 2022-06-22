/*
 * GNU GPL v3 License
 *
 * Copyright 2017  Niccolo` Tubini
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
package it.geoframe.blogspot.whetgeo1d.boundaryconditions;

/**
 * A simple design factory to create a BoundaryCondition objects.
 */

public class DiffusionSimpleBoundaryConditionFactory {
	
	/**
	 * Creates a new BoundaryCondition object.
	 * 
	 * @param type of boundary condition
	 * @return boundCond
	 */
	
	public BoundaryCondition createBoundaryCondition (String type) {

		BoundaryCondition boundaryCondition = null;
		if(type.equalsIgnoreCase("Top Dirichlet") || type.equalsIgnoreCase("TopDirichlet")){
			boundaryCondition = new DiffusionTopBoundaryConditionDirichlet();
		} else if(type.equalsIgnoreCase("Top Neumann") || type.equalsIgnoreCase("TopNeumann")){
			boundaryCondition = new DiffusionTopBoundaryConditionNeumann();
		} else if(type.equalsIgnoreCase("Bottom Dirichlet") || type.equalsIgnoreCase("BottomDirichlet")){
			boundaryCondition = new DiffusionBottomBoundaryConditionDirichlet();
		} else if(type.equalsIgnoreCase("Bottom Neumann") || type.equalsIgnoreCase("BottomNeumann")){
			boundaryCondition = new DiffusionBottomBoundaryConditionNeumann();
		} else if(type.equalsIgnoreCase("Bottom No Gradient") || type.equalsIgnoreCase("Bottom No Gradient")){
			boundaryCondition = new DiffusionBottomBoundaryConditionNoGradient();
		} else {
			System.out.println("\n\n\tERROR: please check the boundary condition name.");
		}


		return boundaryCondition;
		}
}
