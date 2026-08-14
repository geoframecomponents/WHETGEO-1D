/*
 * GNU GPL v3 License
 *
 * Copyright 2026 Andrea Antonello
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

/**
 * GEOframe components directly usable in simulation environments.
 * <p>
 * Each class is an equation solver: it configures boundary conditions and
 * an equation state and runs the time-stepping loop. The
 * actual numerics are handled by a kernel in {@link org.geoframe.whetgeo1d.core.finitevolume}.
 * Classes in {@code org.geoframe.whetgeo1d.core} (kernels, equation states, boundary conditions,
 * derived quantities, surface properties) are internal building blocks used by the solvers
 * below and are not meant to be picked directly.
 * <p>
 * Class names follow one fixed grammar, so the name alone tells you what a solver does:
 * <pre>
 *     &lt;Equation&gt;Solver(With&lt;Feature&gt;)*1D
 * </pre>
 * <ul>
 *   <li><b>{@code <Equation>}</b> &mdash; the base equation being solved, e.g. {@code HeatDiffusion},
 *       {@code HeatAdvectionDiffusion}, {@code Richards}, {@code RichardsConservativeSoluteADE}.</li>
 *   <li><b>{@code Solver}</b> &mdash; always present immediately after the equation name.</li>
 *   <li><b>{@code With<Feature>}</b> &mdash; zero or more optional physics toggles, each appended
 *       independently, e.g. {@code WithFreezingThawing}, {@code WithSurfaceEnergyBalance},
 *       {@code WithRootWaterUptake}. A bare name with no {@code With} clause is the plain base
 *       case; every {@code With} clause is an add-on to that base equation, never part of it.
 *       When more than one feature applies, they are chained in a fixed order: physics/source-term
 *       extensions (e.g. {@code WithFreezingThawing}, {@code WithRootWaterUptake}) before
 *       boundary-treatment extensions (e.g. {@code WithSurfaceEnergyBalance}).</li>
 *   <li><b>{@code 1D}</b> &mdash; always the final suffix, denoting the domain dimensionality.</li>
 * </ul>
 * For example, {@code HeatDiffusionSolverWithFreezingThawingWithSurfaceEnergyBalance1D} solves
 * the heat diffusion equation, with freezing/thawing of water and a surface-energy-balance top
 * boundary condition, over a 1D domain.
 */
package org.geoframe.whetgeo1d.solvers;
