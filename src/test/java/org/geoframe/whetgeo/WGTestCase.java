package org.geoframe.whetgeo;

import java.io.File;
import java.net.URISyntaxException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Paths;

import org.hortonmachine.gears.io.timedependent.OmsTimeSeriesIteratorReader;

import junit.framework.TestCase;

public abstract class WGTestCase extends TestCase {
	protected String getRes(String name) throws Exception {
		URL url = this.getClass().getResource(name);
		if (url == null) {
			throw new Exception("Resource not found: " + name);
		}
		return Paths.get(url.toURI()).toString();
	}
	
	protected String getTmpPath(String prefix, String ext) throws Exception {
		File tempFile = Files.createTempFile(prefix, ext).toFile();
		String pathOutput = tempFile.getAbsolutePath();
		return pathOutput;
	}
	
	
	protected OmsTimeSeriesIteratorReader getTimeseriesReader( String inPath, String id, String startDate, String endDate,
			int timeStepMinutes ) throws URISyntaxException {
		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file = inPath;
		reader.idfield = "ID";
		reader.tStart = startDate;
		reader.tTimestep = timeStepMinutes;
		reader.tEnd = endDate;
		reader.fileNovalue = "-9999";
		reader.initProcess();
		return reader;
	}
}
