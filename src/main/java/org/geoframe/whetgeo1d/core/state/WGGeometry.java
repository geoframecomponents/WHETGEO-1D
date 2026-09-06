package org.geoframe.whetgeo1d.core.state;

public class WGGeometry {
	
	
	public double[] z;
	public double[] spaceDeltaZ;
	public double[] controlVolume;
	
	public WGGeometry(double[] z, double[] spaceDeltaZ, double[] controlVolume) {
		
		this.z = z.clone();
		this.spaceDeltaZ = spaceDeltaZ.clone();
		this.controlVolume = controlVolume.clone();
		
	}


}
