import org.hortonmachine.HM
import org.geoframe.whetgeo1d.io.Whetgeo1DInputsBuilder
import org.geoframe.whetgeo1d.io.Whetgeo1DInputsHandler

//   A 1 m saturated soil column sits at 12 degC. From one day to the next its
//   surface is held at 20 degC and its base at 5 degC. How does the
//   temperature profile evolve, and where does it settle?

def workspace = "/home/hydrologis/storage/lavori_tmp/UNITN/WHETGEO1D/"
def dataFolder = workspace + "heatdiffusion/" as File
def gpkg = dataFolder.absolutePath + "/heatdiffusion.gpkg" 

def startDate = "2000-01-01 00:00"
def endDate = "2000-03-01 00:00"
def topTempC = 20.0 
def bottomTempC = 5.0 
def initTempC = 285.15

def tableTopBC = "timeseries_topBC";
def tableBotBC = "timeseries_bottomBC";


checkFolder(dataFolder)
deleteGpkg(gpkg)

// create the db and use the input builder to create the db tables
HM.connectGeopackage(gpkg).withCloseable{ db ->
    def builder = new Whetgeo1DInputsBuilder(db)
    // create the grid. 
    // * one single layer 1 meter from 0 to bottom split into 50 cells
    // * closure model 0 and parameter set 1
    builder.layer(0.0, 50, 0, 1)
    // close the grid
    .domainBottom(-1.0)
    // initial conditions: eta, psi and t
    .initialCondition(0.0, 0.0, initTempC)
    .initialCondition(-1.0, 0.0, initTempC)
    // parameter set 1
    .parameterSet(
            "spDensity": 2500, 
            "spConductivity": 2.5, 
            "spSpecificHeatCapacity": 1000,
            "thetaS": 0.38, 
            "thetaR": 0.068, 
            "meltingT": 273.15, 
            "alpha": 0.8, 
            "n": 1.09,
            "alphaSpecificStorage": 0.0,
            "betaSpecificStorage": 0.0
        )
    .vanGenuchtenSlots() // n and alpha 
    // boundary temperatures, hourly
    .constantTimeseries(tableTopBC, startDate, endDate, 60, topTempC)
    .constantTimeseries(tableBotBC, startDate, endDate, 60, bottomTempC)
    // write the gpkg
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

