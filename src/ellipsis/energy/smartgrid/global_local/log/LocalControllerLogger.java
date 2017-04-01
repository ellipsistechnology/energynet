package ellipsis.energy.smartgrid.global_local.log;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.complex.Complex;

import ellipsis.energy.grid.Unit;
import ellipsis.energy.smartgrid.global_local.LocalGlobalOptimiser.CostToGo;
import ellipsis.util.CSVLogger;

public class LocalControllerLogger extends CSVLogger
{
	//// Instance ////
	
	private static ThreadLocal<Map<String,LocalControllerLogger>> instances = new ThreadLocal<Map<String,LocalControllerLogger>>();
	private static int defaulLogLevel;
	
	public static LocalControllerLogger instance(String name)
	{
		Map<String, LocalControllerLogger> loggerMap = instances.get();
		if(loggerMap == null)
		{
			loggerMap = new HashMap<String, LocalControllerLogger>();
			instances.set(loggerMap);
		}
		
		LocalControllerLogger logger = loggerMap.get(name);
		if(logger == null)
		{
			logger = new LocalControllerLogger(name);
			logger.resetData();
			loggerMap.put(name, logger);
		}
		return logger;
	}
	
	private LocalControllerLogger(String name)
	{
		setLogLevel(defaulLogLevel);
		
		try
		{
			File file = new File("/opt/energynet/local_controller_"+name+".csv");
			file.createNewFile();
			setOut(new PrintStream(file));
		}
		catch (IOException e)
		{
			throw new RuntimeException(e);
		}
		resetData();
	}
	
	public static void setDefaultLogLevel(int logLevel)
	{
		defaulLogLevel = logLevel;
	}
	
	
	//// Logging ////
	
	private String currentSubnet;

	public void logState(String subnet, int centralIteration, int localIteration, int time, Complex[] state, int logLevel)
	{
		if(!subnet.equals(currentSubnet))
		{
			printDataAndReset(logLevel);
			printHeading(subnet, logLevel);
			currentSubnet = subnet;
		}
		
		addData("Central", centralIteration);
		addData("Local", localIteration);
		addData("Time", time);
		
		int i = 0;
		for (Complex c : state)
		{
			addData("State "+i+" (Re)", c.getReal());
			addData("State "+i+" (Im)", c.getImaginary());
			++i;
		}
	}

	public void logControls(CostToGo estimatedControls, CostToGo appliedControls, int logLevelDebug)
	{
		for (Unit unit : estimatedControls.controls.keySet())
		{
			Complex control = estimatedControls.controls.get(unit);
			addData("Control1 "+unit.getName()+" (Re)", control.getReal());
			addData("Control1 "+unit.getName()+" (Im)", control.getImaginary());
		}
		
		for (Unit unit : appliedControls.controls.keySet())
		{
			Complex control = appliedControls.controls.get(unit);
			addData("Control2 "+unit.getName()+" (Re)", control.getReal());
			addData("Control2 "+unit.getName()+" (Im)", control.getImaginary());
		}
	}
}
