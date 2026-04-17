package org.geoframe.whetgeo1d.utils;

public class GFGeometry {
	
	
	public double[] z;
	public double[] spaceDeltaZ;
	public double[] controlVolume;
	
	public GFGeometry(double[] z, double[] spaceDeltaZ, double[] controlVolume) {
		
		this.z = z.clone();
		this.spaceDeltaZ = spaceDeltaZ.clone();
		this.controlVolume = controlVolume.clone();
		
	}


}
