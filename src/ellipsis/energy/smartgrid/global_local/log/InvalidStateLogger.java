package ellipsis.energy.smartgrid.global_local.log;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import org.apache.commons.math3.complex.Complex;

import ellipsis.energy.grid.Bus;
import ellipsis.energy.grid.Unit;
import ellipsis.energy.smartgrid.global_local.LocalGlobalOptimiser.CostToGo;
import ellipsis.util.CSVLogger;

public class InvalidStateLogger extends CSVLogger
{
	//// Instance ////
	
	private static ThreadLocal<InvalidStateLogger> instances = new ThreadLocal<InvalidStateLogger>();
	
	public static InvalidStateLogger instance()
	{
		InvalidStateLogger instance = instances.get();
		if(instance == null)
		{
			instance = new InvalidStateLogger();
			instances.set(instance);
		}
		return instance;
	}
	
	private InvalidStateLogger()
	{
		try
		{
			File file = new File("/opt/energynet/invalid_state.csv");
			file.createNewFile();
			setOut(new PrintStream(file));
		}
		catch (IOException e)
		{
			throw new RuntimeException(e);
		}
		resetData();
	}
	
	
	//// Logging ////
	
	public void logInvlaidState(int centralIteration, int localIteration, int timeIndex, CostToGo ctgEstimate, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		addData("Central Iteration", centralIteration);
		addData("Local Iteration", localIteration);
		addData("Time", timeIndex);
		for (Unit unit : ctgEstimate.controls.keySet())
		{
			String name = unit.getName();
			Complex control = ctgEstimate.controls.get(unit);
			addData("R{"+name+"}", control.getReal());
			addData("I{"+name+"}", control.getImaginary());
		}
		fillData("");
	}

	public void logInvalidCentralVoltage(int centralIteration, int timeIndex, Bus bus, double v, int logLevel)
	{
		printHeading("Error - Invalid voltage found calculating global cost: " +
				"Central iteration = " + centralIteration + 
				", Time = " + timeIndex +
				", Bus = '" + bus.getName() +
				"', Voltage = " + v,
				logLevel);
	}
}
