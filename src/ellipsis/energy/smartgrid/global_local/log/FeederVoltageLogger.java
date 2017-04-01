package ellipsis.energy.smartgrid.global_local.log;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexUtils;
import org.apache.commons.math3.linear.RealVector;

import ellipsis.energy.calculation.AnalysisResults;
import ellipsis.energy.grid.Bus;
import ellipsis.util.CSVLogger;

public class FeederVoltageLogger extends CSVLogger
{
	private String targetControlledBus;
	private int targetTime;
	private Complex lastApproximateLocalVoltage;

	
	//// Instance ////
	
	private static ThreadLocal<FeederVoltageLogger> instances = new ThreadLocal<FeederVoltageLogger>();
	
	public static FeederVoltageLogger instance()
	{
		FeederVoltageLogger instance = instances.get();
		if(instance == null)
		{
			instance = new FeederVoltageLogger();
			instances.set(instance);
		}
		return instance;
	}
	
	private FeederVoltageLogger()
	{
		setOut("/opt/energynet/feeder_voltage.csv");
		reset();
	}
	
	
	//// Logging ////

	public void reset()
	{
		resetDataWithColumns("Central Iteration", "True Voltage", "Global Approximate Voltage", "Local Approximate Voltage", "Error");
	}

	public void logHeading(int logLevel)
	{
		printHeading("Controlled bus: '"+targetControlledBus+"; Time: "+targetTime, logLevel);
	}

	public void logRow(int centralIteration, int timeIndex, Bus controlledBus, AnalysisResults results, RealVector errors, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		if(targetTime != timeIndex || !targetControlledBus.equals(controlledBus.getName()))
			return;
		
		addData("Central Iteration", centralIteration);
//		addData("True Voltage", lastTrueVoltage.abs());
		addData("Local Approximate Voltage", lastApproximateLocalVoltage.abs());
		
		int size = errors.getDimension()/2;
		int index = indexOf(controlledBus.getName(), results);
		double abs = errors.getEntry(size + index);
//		double arg = errors.getEntry(index);
		if(abs < 0)
		{
			abs = -abs;
//			arg += Math.PI;
		}
//		Complex error = ComplexUtils.polar2Complex(abs, arg);
		addData("Error", abs);//error.abs());
	}

	private int indexOf(String controlledBus, AnalysisResults results)
	{
		int index = results.getBusNumbers().get(controlledBus);
		if(index > results.getSlackIndex())
			index -= 1;
		return index;
	}
	
	public void logTrueVoltage(int timeIndex, AnalysisResults results, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		if(timeIndex != targetTime)
			return;
		
		Complex voltage = results.getBusVoltage(targetControlledBus);
		addData("True Voltage", voltage.abs());
	}
	
	public void logLocalVoltage(int timeIndex, Bus controlledBus, Complex voltage, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		if(timeIndex != targetTime || !targetControlledBus.equals(controlledBus.getName()))
			return;
		
		this.lastApproximateLocalVoltage = voltage;
	}

	public void logApproximateCentralVoltages(int t, AnalysisResults results, RealVector approximateVoltages, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		if(t != targetTime)
			return;
		
		int index = indexOf(targetControlledBus, results);
		int size = approximateVoltages.getDimension()/2;
		Complex voltage = ComplexUtils.polar2Complex(approximateVoltages.getEntry(index+size), approximateVoltages.getEntry(index));
		addData("Global Approximate Voltage", voltage.abs());
	}
	
	
	//// Accessors ////

	public String getTargetControlledBus()
	{
		return targetControlledBus;
	}

	public void setTargetControlledBus(String targetControlledBus)
	{
		this.targetControlledBus = targetControlledBus;
	}

	public int getTargetTime()
	{
		return targetTime;
	}

	public void setTargetTime(int targetTime)
	{
		this.targetTime = targetTime;
	}
}
