package ellipsis.energy.smartgrid.global_local.log;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import ellipsis.util.CSVLogger;

public class ValidSolutionLogger extends CSVLogger
{
	//// Instance ////
	
	private static ThreadLocal<ValidSolutionLogger> instances = new ThreadLocal<ValidSolutionLogger>();
	
	public static ValidSolutionLogger instance()
	{
		ValidSolutionLogger instance = instances.get();
		if(instance == null)
		{
			instance = new ValidSolutionLogger();
			instances.set(instance);
		}
		return instance;
	}
	
	private ValidSolutionLogger()
	{
		try
		{
			File file = new File("/opt/energynet/valid_solutions.csv");
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
	
	private int lastCentral = -1;
	private int lastLocal = -1;
	
	public void logLocalValid(int centralIteration, int localIteration, String subnet, boolean valid, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		if(centralIteration > lastCentral)
		{
			lastCentral = centralIteration;
			lastLocal = -1;
		}
		
		if(localIteration > lastLocal)
		{
			addData("Central", centralIteration);
			addData("Local", localIteration);
			lastLocal = localIteration;
		}
		
		addData(subnet, valid);
	}
	
	public void logCentralValid(boolean valid, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		fillData("");
		addData("Central Valid", valid);
		fillData("");
	}
}
