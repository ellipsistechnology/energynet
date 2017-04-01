package ellipsis.energy.smartgrid.global_local.log;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.complex.Complex;

import ellipsis.energy.grid.Unit;
import ellipsis.util.CSVLogger;

public class ControlLogger extends CSVLogger
{
	//// Instance ////
	
	private static ThreadLocal<Map<String, ControlLogger>> instances = new ThreadLocal<Map<String, ControlLogger>>();
	
	public static ControlLogger instance()
	{
		return instance("controls");
	}
	
	public static ControlLogger instance(String key)
	{
		if(instances.get() == null)
			instances.set(new HashMap<String, ControlLogger>());
		
		ControlLogger instance = instances.get().get(key);
		if(instance == null)
		{
			instance = new ControlLogger(key);
			instances.get().put(key, instance);
		}
		return instance;
	}
	
	private ControlLogger(String key)
	{
		setOut("/opt/energynet/"+key+".csv");
		resetData();
	}
	
	
	//// Logging ////
	
	private boolean headingsPrinted = false;
	String subnet;
	int centralIteration; 
	int adpIteration;
	int timeIndex;
	public double highestVoltage[] = new double[12];

	public void logControls(
			String subnet,
			int centralIteration, 
			int adpIteration, 
			int timeIndex, 
			Map<Unit, Complex> originalControls, 
			Map<Unit, Complex> appliedControls,
			int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		this.subnet = subnet;
		this.centralIteration = centralIteration; 
		this.adpIteration = adpIteration;
		this.timeIndex = timeIndex;
		
		if(!shouldLog())
			return;
		
		addData("Central Iteration", centralIteration);
		addData("Local Iteration", adpIteration);
		addData("Time", timeIndex);
		
		for (Unit unit : originalControls.keySet())
		{
			String name = unit.getName();
			addData(name+"(original)", originalControls.get(unit).getReal());
			addData(name+"(applied)", appliedControls.get(unit).getReal());
		}
		
		if(!headingsPrinted)
		{
			printDataAndContinue(logLevel, true);
			headingsPrinted = true;
		}
	}

	private boolean shouldLog()
	{
		return (this.centralIteration == 2 && this.adpIteration > 26 && this.subnet.equals("Local (Bus 5)"));
	}
	
	@Override
	public void printDataAndContinue(int logLevel, boolean includeHeadings)
	{
		if(!shouldLog())
			return;
		
		addData("Highest Voltage", highestVoltage[timeIndex]);
		
		super.printDataAndContinue(logLevel, includeHeadings);
	}

	public void setHighestVoltage(double highestVoltage, int timeIndex, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		if(timeIndex < 12)
			this.highestVoltage[timeIndex] = highestVoltage;
	}
}
