package it.geoframe.blogspot.whetgeo1d.data;

import java.util.ArrayList;
import java.util.List;

public class Geometry {
	
	private static Geometry uniqueInstance;
	
	public static Geometry getInstance() {
		/*if (uniqueInstance == null) {
			uniqueInstance = new Variables(waterSuction, temperature);
		}*/
		return uniqueInstance;
	}
	
	public static Geometry getInstance(double[] z, double[] spaceDeltaZ, double[] controlVolume) {
		if (uniqueInstance == null) {
			uniqueInstance = new Geometry(z, spaceDeltaZ, controlVolume);
		}
		return uniqueInstance;
	}
	
	
	public double[] z;
	public double[] spaceDeltaZ;
	public double[] controlVolume;
	
	private Geometry(double[] z, double[] spaceDeltaZ, double[] controlVolume) {
		
		this.z = z.clone();
		this.spaceDeltaZ = spaceDeltaZ.clone();
		this.controlVolume = controlVolume.clone();
		
	}


}
