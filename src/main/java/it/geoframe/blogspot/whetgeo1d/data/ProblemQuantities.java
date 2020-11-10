package it.geoframe.blogspot.whetgeo1d.data;


public class ProblemQuantities {
	
	private static ProblemQuantities uniqueInstance;
	
	public static ProblemQuantities getInstance() {
		return uniqueInstance;
	}
	
	public static ProblemQuantities getInstance(double[] icWaterSuction, double[] icTemperature, int[] rheologyID, int[] parameterID) {
		if (uniqueInstance == null) {
			uniqueInstance = new ProblemQuantities(icWaterSuction, icTemperature, rheologyID, parameterID);
		}
		return uniqueInstance;
	}
	
	
	public double[] waterSuctions;
	public double[] temperatures;
	public double[] thetas;
    public double[] thetasNew;    
	public double[] kappas;
	public double[] kappasInterface;
	public double[] volumes;
	public double[] volumesNew;
	public double[] darcyVelocities;
	public double[] darcyVelocitiesCapillary;
	public double[] darcyVelocitiesGravity;
	public double[] poreVelocities;
	public double[] celerities;        // Rasmussen et al. 2000
	public double[] kinematicRatio;  // Rasmussen et al. 2000
	public double[] waterSuctionStar1;
	public double[] waterSuctionStar2;
	public double[] waterSuctionStar3;
	
	public int[] rheologyID;
	public int[] parameterID;
	
	public double runOff;
	public double kappaBottom;
	public double waterVolume;
	public double waterVolumeNew;
	public double errorVolume;
	public double volumeLost;
	public double richardsTopBCValue;
	public double richardsBottomBCValue;

	public double sumETs;
	public double[] ETs;
	
	private ProblemQuantities(double[] icWaterSuction, double[] icTemperature, int[] rheologyID, int[] parameterID) {
		
		waterSuctions = icWaterSuction.clone();
		temperatures = icTemperature.clone();
		thetas = new double[icWaterSuction.length];
		thetasNew = new double[icWaterSuction.length];
		kappas = new double[icWaterSuction.length+1];
		kappasInterface = new double[icWaterSuction.length+1];
		volumes = new double[icWaterSuction.length];
		volumesNew = new double[icWaterSuction.length];
		darcyVelocities = new double[icWaterSuction.length+1];
		darcyVelocitiesCapillary = new double[icWaterSuction.length+1];
		darcyVelocitiesGravity = new double[icWaterSuction.length+1];
		poreVelocities = new double[icWaterSuction.length+1];
		celerities = new double[icWaterSuction.length+1];
		kinematicRatio = new double[icWaterSuction.length+1];
		ETs = new double[icWaterSuction.length];
		waterSuctionStar1 = new double[icWaterSuction.length];
		waterSuctionStar2 = new double[icWaterSuction.length];
		waterSuctionStar3 = new double[icWaterSuction.length];
		
		this.rheologyID = rheologyID.clone();
		this.parameterID = parameterID.clone();
		
		
	}


}
