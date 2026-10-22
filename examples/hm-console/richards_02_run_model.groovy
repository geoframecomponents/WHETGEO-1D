import org.geoframe.whetgeo1d.io.Whetgeo1DInputsHandler
import org.geoframe.whetgeo1d.core.boundaryconditions.IBoundaryCondition
import org.geoframe.whetgeo1d.solvers.RichardsSolver1D
import org.geoframe.whetgeo1d.io.Whetgeo1DOutputsHandler
import org.hortonmachine.HM
import org.hortonmachine.gears.utils.chart.Scatter

// run HeatDiffusionSolver1D in hourly mode

def workspace = "/home/hydrologis/storage/lavori_tmp/UNITN/WHETGEO1D/"
def dataFolder = workspace + "richards/" as File
def inGpkg = dataFolder.absolutePath + "/richards.gpkg" 
def outGpkg = dataFolder.absolutePath + "/richards_output.gpkg" 
def startDate = "2015-01-15 00:00"
def endDate = "2015-12-15 00:00"
def stationId = 0

// cleanup first
deleteGpkg(outGpkg)

def inputHandler = new Whetgeo1DInputsHandler(inGpkg)
inputHandler.read()

def kmax = inputHandler.KMAX
def topBC = IBoundaryCondition.RichardsBoundaryConditionType.TOP_DIRICHLET
def bottomBC = IBoundaryCondition.RichardsBoundaryConditionType.BOTTOM_FREE_DRAINAGE

////////////////////////////////////
// DEFINE THE SOLVER
////////////////////////////////////
def solver = new RichardsSolver1D()
// grid, initial condition and soil parameters, as read from the input gpkg
inputHandler.applyTo(solver)
solver.beta0 = -766.45
solver.referenceTemperatureSWRC = 278.15
solver.maxPonding = 0.0
// models
solver.typeClosureEquation = ["Van Genuchten"] as String[]
solver.typeEquationState = ["Van Genuchten"] as String[]
solver.typeUHCModel = ["Mualem Van Genuchten"] as String[]
solver.typeUHCTemperatureModel = "notemperature"
solver.interfaceHydraulicConductivityModel = "max"
solver.topBCType = topBC
solver.bottomBCType = bottomBC
// numerics: 1 hour steps
solver.tTimeStep = 3600
solver.timeDelta = 3600
solver.delta = 0
solver.newtonTolerance = 1e-11
solver.nestedNewton = 1
solver.picardIteration = 1
solver.stationID = stationId

////////////////////////////////////
// get data iterators and output writer
////////////////////////////////////

// keep a few profiles along the way, to see the approach to steady state
def snapshotHours = [3, 6, 12, 18, 24, 30]
def snapshots = [:]
//snapshots["start"] = inputHandler.temperatureIC

def stepCount = 0
try (inputHandler
    def topIterator = inputHandler.iterateTimeseries("timeseries_topBC", startDate, endDate, 1000)
    def bottomIterator = inputHandler.iterateTimeseries("timeseries_bottomBC", startDate, endDate, 1000)
    def outWriter = new Whetgeo1DOutputsHandler(outGpkg, 500)
    ){
    outWriter.writeIntervalMinutes = 60 * 24 // one output per day is enough here
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
            snapshots["after ${stepCount} h"] = solver.outWaterContent.clone()
        }

        outWriter.timestamp = timestamp
        outWriter.theta = solver.outWaterContent
        outWriter.waterSuction = solver.outWaterSuctions
        outWriter.darcyVelocity = solver.outDarcyVelocity
        outWriter.errorVolume = solver.outErrorVolume
        outWriter.topBC = solver.outTopBCValue
        outWriter.bottomBC = solver.outBottomBCValue
        outWriter.write()
    } 
}

snapshots["after ${stepCount} hours"] = solver.outWaterContent.clone()
println "Simulation geopackage written: " + outGpkg


////////////////////////////////////
// PLOT RESULTS
////////////////////////////////////
def chart = new Scatter("Water content: the wetting front moving down")
chart.setXLabel("theta [-]")
chart.setYLabel("eta [m]")
// a profile: keep the points in depth order instead of sorting them by theta
chart.setSortByX(false)
snapshots.each { name, thetas -> chart.addSeries(name, thetas, inputHandler.eta) }
chart.setShowShapes(snapshots.collect { false })
chart.addHorizontalAnnotation("layer 1 | layer 2", -1.0)
chart.addHorizontalAnnotation("layer 2 | layer 3", -2.0)
chart.getImage(900, 700)



////////////////////////////////////////////
// HELPER FUNCTIONS
////////////////////////////////////////////
def deleteGpkg(path){
    if ((path as File).exists()) {
        (path as File).delete()
    }    
}





