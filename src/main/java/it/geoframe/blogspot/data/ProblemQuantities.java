package it.geoframe.blogspot.data;

import java.util.ArrayList;
import java.util.List;

public class ProblemQuantities {
	
	private static ProblemQuantities uniqueInstance;
	
	public static ProblemQuantities getInstance() {
		/*if (uniqueInstance == null) {
			uniqueInstance = new Variables(waterSuction, temperature);
		}*/
		return uniqueInstance;
	}
	
	public static ProblemQuantities getInstance(double[] icWaterSuction, double[] icTemperature) {
		if (uniqueInstance == null) {
			uniqueInstance = new ProblemQuantities(icWaterSuction, icTemperature);
		}
		return uniqueInstance;
	}
	
	
	public double[] waterSuctions;
	public double[] temperatures;
	public double[] thetasOld;
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
	public double[] ETs;
	public double[] xStar1;
	public double[] xStar2;
	public double[] xStar3;
	
	public double sumETs;
	public double runOff;
	public double kappaBottom;
	public double waterVolume;
	public double waterVolumeNew;
	public double errorVolume;
	

	
	private ProblemQuantities(double[] icWaterSuction, double[] icTemperature) {
		
		waterSuctions = icWaterSuction.clone();
		temperatures = icTemperature.clone();
		thetasOld = new double[icWaterSuction.length];
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
		xStar1 = new double[icWaterSuction.length];
		xStar2 = new double[icWaterSuction.length];
		xStar3 = new double[icWaterSuction.length];
		
	}


}
