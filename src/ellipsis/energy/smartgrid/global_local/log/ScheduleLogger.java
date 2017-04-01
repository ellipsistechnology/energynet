package ellipsis.energy.smartgrid.global_local.log;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import org.apache.commons.math3.complex.Complex;

import ellipsis.energy.grid.Unit;
import ellipsis.energy.smartgrid.global_local.LocalGlobalOptimiser;
import ellipsis.energy.smartgrid.global_local.LocalGlobalOptimiser.CostToGo;
import ellipsis.util.CSVLogger;

public class ScheduleLogger extends CSVLogger
{
	//// Instance ////
	
	private static ThreadLocal<ScheduleLogger> instances = new ThreadLocal<ScheduleLogger>();
	
	public static ScheduleLogger instance()
	{
		ScheduleLogger instance = instances.get();
		if(instance == null)
		{
			instance = new ScheduleLogger();
			instances.set(instance);
		}
		return instance;
	}
	
	private ScheduleLogger()
	{
		try
		{
			File file = new File("/opt/energynet/schedule.csv");
			file.createNewFile();
			setOut(new PrintStream(file));
			resetData();
		}
		catch (IOException e)
		{
			throw new RuntimeException(e);
		}
	}
	
	
	//// Logging ////
	
	public void logCentralIteration(int centralIteration, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		addData("Central Iteration", centralIteration);
	}
	
	public void logLocalIteration(int localIteration, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		addData("Local Iteration", localIteration);
	}
	
	public void logTime(int time, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		addData("Time", time);
	}
	
	public void logControl(
			LocalGlobalOptimiser localGlobalOptimiser, 
			CostToGo ctg,
			int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		if(ctg.controls == null)
		{
			fillData("-");
		}
		else
		{
			addData("Controller", localGlobalOptimiser.toString());
			for (Unit unit : ctg.controls.keySet())
			{
				Complex control = ctg.controls.get(unit);
				addData(unit.getName(), control.getReal());
			}
		}
	}
}
