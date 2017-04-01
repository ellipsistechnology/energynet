package ellipsis.energy.smartgrid.global_local.log;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Map;

import ellipsis.energy.grid.Bus;
import ellipsis.util.CSVLogger;

public class PartitioningLogger extends CSVLogger
{
	//// Instance ////
	
	private static ThreadLocal<PartitioningLogger> instances = new ThreadLocal<PartitioningLogger>();
	
	public static PartitioningLogger instance()
	{
		PartitioningLogger instance = instances.get();
		if(instance == null)
		{
			instance = new PartitioningLogger();
			instances.set(instance);
		}
		return instance;
	}
	
	private PartitioningLogger()
	{
		try
		{
			File file = new File("/opt/energynet/partitions.csv");
			file.createNewFile();
			setOut(new PrintStream(file));
		}
		catch (IOException e)
		{
			throw new RuntimeException(e);
		}
		resetDataWithColumns("Bus", "Sensitivity(P)","Delta P","Sensitivity(Q)","Delta Q","Value");
	}

	
	//// Logging ////

	public void logTerms(
			Map<String, Integer> busNumbers, 
			int slackIndex, 
			String subnet,
			double[] vArgRow, 
			double[] vAbsRow, 
			double[] deltaPQMaxArray, 
			double[] vArgTerms, 
			double[] vAbsTerms,
			int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		int size = vAbsTerms.length/2;
		String[] busNames = getBusNames(busNumbers, slackIndex);
		
		printHeading(subnet, logLevel);
		printHeading("Arguments", logLevel);
		for (int i = 0; i < size; i++)
		{
			String name = busNames[i];
			addData(name, vArgRow[i], deltaPQMaxArray[i], vArgRow[i+size], deltaPQMaxArray[i+size], vArgTerms[i]);
		}
		printDataAndContinue(logLevel, true);
		
		printHeading("Absolutes", logLevel);
		for (int i = 0; i < size; i++)
		{
			String name = busNames[i];
			addData(name, vAbsRow[i], deltaPQMaxArray[i], vAbsRow[i+size], deltaPQMaxArray[i+size], vAbsTerms[i]);
		}
		printDataAndContinue(logLevel, true);
	}

	private String[] getBusNames(Map<String, Integer> busNumbers, int slackIndex)
	{
		String names[] = new String[busNumbers.size()];
		for (String name : busNumbers.keySet())
		{
			int index = busNumbers.get(name);
			if(index == slackIndex)
				continue;
			if(index > slackIndex)
				--index;
			
			names[index] = name;
		}
		
		return names;
	}

	public void logSubnet(HashSet<Bus> subnet, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		printHeading("Initial Composition:", logLevel);
		for (Bus bus : subnet)
		{
			printCell(bus.getName(), logLevel);
		}
		endRow(logLevel);
	}
}
