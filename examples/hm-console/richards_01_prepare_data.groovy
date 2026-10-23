import org.hortonmachine.HM
import org.geoframe.whetgeo1d.io.Whetgeo1DInputsBuilder
import org.geoframe.whetgeo1d.io.Whetgeo1DInputsHandler

//   A 3 m column made of three 1 m soil layers (coarse on top, finer below)
//   starts in hydrostatic equilibrium with the water table at its base, so it is
//   dry near the surface. Then the surface is kept saturated (psi = 0) while
//   the base drains freely. How fast does the water move down, and how do the
//   layers change it?

def workspace = "/home/hydrologis/storage/lavori_tmp/UNITN/WHETGEO1D/"
def dataFolder = workspace + "richards/" as File
def gpkg = dataFolder.absolutePath + "/richards.gpkg" 
def startDate = "2015-01-15 00:00"
def endDate = "2015-12-15 00:00"

def tableTopBC = "timeseries_topBC";

checkFolder(dataFolder)
deleteGpkg(gpkg)

// create the db and use the input builder to create the db tables
HM.connectGeopackage(gpkg).withCloseable{ db ->
    def builder = new Whetgeo1DInputsBuilder(db)
    builder
    // three 1 m layers of 50 cells each, every layer with its own parameter set
        .layer(0.0, 50, 0, 1)
        .layer(-1.0, 50, 0, 2)
        .layer(-2.0, 50, 0, 3)
        .domainBottom(-3.0)
        // psi = 0 at the water table (the base), -3 m at the surface
        .initialCondition(0.0, -3.0, 273.15)
        .initialCondition(-3.0, 0.0, 273.15)
        // Van Genuchten parameters, one set per layer: water content at saturation and
        // residual [-], alpha [1/m], n [-], saturated hydraulic conductivity Ks [m/s].
        // Layer 1 drains fastest, layer 3 is the finest.
        .parameterSet(thetaS: 0.36689, thetaR: 0.05385, alpha: 2.991, n: 2.54723,
                alphaSpecificStorage: 0.0, betaSpecificStorage: 0.0, Ks: 5.40583e-05)
        .parameterSet(thetaS: 0.46291, thetaR: 0.06707, alpha: 0.445, n: 1.65275,
                alphaSpecificStorage: 0.0, betaSpecificStorage: 0.0, Ks: 5.18521e-06)
        .parameterSet(thetaS: 0.48664, thetaR: 0.13335, alpha: 0.846, n: 1.34174,
                alphaSpecificStorage: 0.0, betaSpecificStorage: 0.0, Ks: 9.73284e-07)
        .vanGenuchtenSlots()
        // top: psi = 0 m constant
        .constantTimeseries(tableTopBC, startDate, endDate, 60, 0.0)
        // bottom: free drainage, which needs no data
        // write record to gpkg
        .write()
}

println "Input data written to: " + gpkg










////////////////////////////////////////////
// HELPER FUNCTIONS
////////////////////////////////////////////
def checkFolder(path){
    if (!(path as File).exists()) {
        path.mkdirs()
    }
}

def deleteGpkg(path){
    if ((path as File).exists()) {
        (path as File).delete()
    }    
}

