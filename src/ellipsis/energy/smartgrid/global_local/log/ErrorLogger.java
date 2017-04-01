package ellipsis.energy.smartgrid.global_local.log;

import org.apache.commons.math3.linear.RealVector;

import ellipsis.util.CSVLogger;

public class ErrorLogger extends CSVLogger
{
	//// Instance ////
	
	private static ThreadLocal<ErrorLogger> instances = new ThreadLocal<ErrorLogger>();
	
	public static ErrorLogger instance()
	{
		ErrorLogger instance = instances.get();
		if(instance == null)
		{
			instance = new ErrorLogger();
			instances.set(instance);
		}
		return instance;
	}
	
	private ErrorLogger()
	{
		setOut("/opt/energynet/error-deltas.csv");
		resetData();
	}
	
	
	//// Logging ////

	public void logVoltageError(int centralIteration, String subnet, int t, RealVector voltageError, RealVector voltageErrorDelta, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		if(!(subnet.equals("Local (Bus 5)") && t == 2))
			return;
		
		addData("Central Iteration", centralIteration);
		addData("Sub-network", subnet);
		addData("Time", t);
		for (int i = 0; i < voltageErrorDelta.getDimension(); i++)
		{
			addData("E"+i, voltageError.getEntry(i));
		}

		for (int i = 0; i < voltageErrorDelta.getDimension(); i++)
		{
			addData("D"+i, voltageErrorDelta.getEntry(i));
		}
	}
}
