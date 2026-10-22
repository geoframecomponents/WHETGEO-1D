import org.geoframe.whetgeo1d.io.Whetgeo1DInputsHandler
import static org.geoframe.whetgeo1d.core.boundaryconditions.IBoundaryCondition.DiffusionBoundaryConditionType.*
import org.geoframe.whetgeo1d.solvers.HeatDiffusionSolver1D
import org.geoframe.whetgeo1d.io.Whetgeo1DOutputsHandler
import org.hortonmachine.HM
import org.hortonmachine.gears.utils.chart.Scatter

// run HeatDiffusionSolver1D in hourly mode

def workspace = "/home/hydrologis/storage/lavori_tmp/UNITN/WHETGEO1D/"
def dataFolder = workspace + "heatdiffusion/" as File
def inGpkg = dataFolder.absolutePath + "/heatdiffusion.gpkg" 
def outGpkg = dataFolder.absolutePath + "/heatdiffusion_output.gpkg" 
def startDate = "2000-01-01 00:00"
def endDate = "2000-03-01 00:00"
def stationId = 0

def tableTopBC = "timeseries_topBC";
def tableBotBC = "timeseries_bottomBC";

// cleanup first
deleteGpkg(outGpkg)

def inputHandler = new Whetgeo1DInputsHandler(inGpkg)
inputHandler.read()

def kmax = inputHandler.KMAX
println "Simulation on column with " + kmax + 
        " cells [" + inputHandler.etaDual[kmax] + " to " +
        inputHandler.etaDual[0] + "]"

def topBC = TOP_DIRICHLET
def bottomBC = BOTTOM_DIRICHLET

////////////////////////////////////
// DEFINE THE SOLVER
////////////////////////////////////
def solver = new HeatDiffusionSolver1D()
// grid, initial condition and soil parameters, as read from the input gpkg
inputHandler.applyTo(solver)
solver.beta0 = -766.45
solver.referenceTemperatureSWRC = 278.15
// models
solver.typeClosureEquation = ["Van Genuchten"] as String[]
solver.typeEquationState = ["SoilInternalEnergy"] as String[]
solver.typeThermalConductivity = ["Cosenza"] as String[]
solver.interfaceThermalConductivityModel = "max"
solver.topBCType = topBC
solver.bottomBCType = bottomBC
// numerics: 1 hour steps
solver.tTimeStep = 3600
solver.timeDelta = 3600
solver.delta = 0
solver.newtonTolerance = 1e-6
solver.nestedNewton = 1
solver.picardIteration = 1
solver.stationID = stationId

////////////////////////////////////
// get data iterators and output writer
////////////////////////////////////

// keep a few profiles along the way, to see the approach to steady state
def snapshotHours = [1, 6, 24, 72]
def snapshots = [:]
snapshots["start"] = inputHandler.temperatureIC

def stepCount = 0
try (inputHandler
    def topIterator = inputHandler.iterateTimeseries(tableTopBC, startDate, endDate, 1000)
    def bottomIterator = inputHandler.iterateTimeseries(tableBotBC, startDate, endDate, 1000)
    def outWriter = new Whetgeo1DOutputsHandler(outGpkg, 500)
    ){
    // grid, initial condition, soil parameters and boundary condition types, written once
    outWriter.describe(inputHandler, topBC, bottomBC)

    while (topIterator.next() && bottomIterator.next()) {
        def timestamp = topIterator.timestamp()
        def tsStr = HM.ts2str(timestamp)
        println "Processing timestamp: " + tsStr
        solver.inTopBC = [(stationId): topIterator.values()]
        solver.inBottomBC = [(stationId): bottomIterator.values()]
        solver.inCurrentDate = tsStr

        solver.solve()

        stepCount++
        if (stepCount in snapshotHours) {
            snapshots["after ${stepCount} h"] = solver.outTemperature.clone()
        }

        outWriter.timestamp = timestamp
        outWriter.temperature = solver.outTemperature
        outWriter.internalEnergy = solver.outInternalEnergy
        outWriter.heatFlux = solver.outDiffusionHeatFlux
        outWriter.error = solver.outErrorInternalEnergy
        outWriter.topBC = solver.outHeatFluxTop
        outWriter.bottomBC = solver.outHeatFluxBottom
        outWriter.write()
    } 
}

snapshots["after ${stepCount} hours"] = solver.outTemperature.clone()
println "Simulation geopackage written: " + outGpkg


////////////////////////////////////
// PLOT RESULTS
////////////////////////////////////
def chart = new Scatter("Temperature profile approaching steady state")
chart.setXLabel("temperature [degC]")
chart.setYLabel("eta [m]")
// a profile: keep the points in depth order instead of sorting them by temperature
chart.setSortByX(false)
snapshots.each { name, temperatures ->
    chart.addSeries(name, temperatures.collect { it - 273.15 } as double[], inputHandler.eta)
}
chart.setShowLines(snapshots.collect { true })
chart.setShowShapes(snapshots.collect { false })

chart.getImage(900, 600)






////////////////////////////////////////////
// HELPER FUNCTIONS
////////////////////////////////////////////
def deleteGpkg(path){
    if ((path as File).exists()) {
        (path as File).delete()
    }    
}





