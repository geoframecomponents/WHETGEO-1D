import org.geoframe.whetgeo1d.io.Whetgeo1DInputsHandler
import org.geoframe.whetgeo1d.core.boundaryconditions.IBoundaryCondition
import org.geoframe.whetgeo1d.solvers.RichardsSolver1D
import org.geoframe.whetgeo1d.io.Whetgeo1DOutputsHandler
import org.hortonmachine.HM
import org.hortonmachine.gears.utils.chart.Scatter
import org.hortonmachine.gears.utils.chart.StackedTimeSeries

// run HeatDiffusionSolver1D in hourly mode

def workspace = "/home/hydrologis/storage/lavori_tmp/UNITN/WHETGEO1D/"
def dataFolder = workspace + "richards_rain/" as File
def inGpkg = dataFolder.absolutePath + "/richards_rain.gpkg" 
def outGpkg = dataFolder.absolutePath + "/richards_rain_output.gpkg" 
def startDate = "2015-01-15 00:00"
def endDate = "2015-12-15 00:00"
def stationId = 0

// cleanup first
deleteGpkg(outGpkg)

def inputHandler = new Whetgeo1DInputsHandler(inGpkg)
inputHandler.read()
def kmax = inputHandler.KMAX
def pondIndex = kmax - 1 // the pseudo cell is the last one, cells go bottom to top
println "Column with ${kmax - 1} soil cells + 1 pond cell"

def topBC = IBoundaryCondition.RichardsBoundaryConditionType.TOP_COUPLED
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
// models, by equationStateID: 0 = pond, 1 = soil
solver.typeClosureEquation = ["Water Depth", "Van Genuchten"] as String[]
solver.typeEquationState = ["Water Depth", "Van Genuchten"] as String[]
solver.typeUHCModel = ["", "Mualem Van Genuchten"] as String[]
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

// what to chart: rain, pond depth and water content near the top of each layer
def watchDepths = [-0.1, -1.1, -2.1]
// find cells for that depths
def watchCells = watchDepths.collect { d -> (0..<pondIndex).min { Math.abs(inputHandler.eta[it] - d) } }
def times = []
def rainValues = []
def pondValues = []
def thetaValues = watchDepths.collect { [] }

def cumulativeRain = 0.0
def cumulativeDrainage = 0.0

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
        
        if (stepCount % 1000 == 0) {
            println "Processing timestamp: " + tsStr
        }
        
        def rain = topIterator.values()[0]
        solver.inTopBC = [(stationId): topIterator.values()]
        solver.inBottomBC = [(stationId): bottomIterator.values()]
        solver.inCurrentDate = tsStr
        solver.solve()

        stepCount++
        
        
        def thetas = solver.outWaterContent
        
        cumulativeRain += rain / 1000.0
        // darcy velocity at the bottom interface [m/s], negative downwards
        cumulativeDrainage += -solver.outDarcyVelocity[0] * solver.tTimeStep
        
        // collect data for info/chart
        times << timestamp
        rainValues << rain
        pondValues << thetas[pondIndex] * 1000.0
        watchCells.eachWithIndex { cell, k -> thetaValues[k] << thetas[cell] }

        outWriter.timestamp = timestamp
        outWriter.theta = thetas
        outWriter.waterSuction = solver.outWaterSuctions
        outWriter.darcyVelocity = solver.outDarcyVelocity
        outWriter.errorVolume = solver.outErrorVolume
        outWriter.topBC = solver.outTopBCValue
        outWriter.bottomBC = solver.outBottomBCValue
        outWriter.write()
    } 
}

println "Simulation geopackage written: " + outGpkg
println String.format("Rain:              %8.1f mm", cumulativeRain * 1000)
println String.format("Drained at base:   %8.1f mm", cumulativeDrainage * 1000)

////////////////////////////////////
// PLOT RESULTS
////////////////////////////////////
def timesArray = times as long[]
def chart = new StackedTimeSeries("Rain, ponding and soil water content")
def rainPanel = chart.addPanel("rain [mm/h]").addSeries("rain", timesArray, rainValues as double[])
rainPanel.asBars = true
rainPanel.inverted = true
chart.addPanel("pond [mm]").addSeries("pond depth", timesArray, pondValues as double[])
thetaPanel = chart.addPanel("theta [-]")
thetaPanel.weight = 3
watchDepths.eachWithIndex { depth, k ->
    thetaPanel.addSeries("theta at ${depth} m", timesArray, thetaValues[k] as double[])
}
chart.getImage(900, 600)



////////////////////////////////////////////
// HELPER FUNCTIONS
////////////////////////////////////////////
def deleteGpkg(path){
    if ((path as File).exists()) {
        (path as File).delete()
    }    
}





