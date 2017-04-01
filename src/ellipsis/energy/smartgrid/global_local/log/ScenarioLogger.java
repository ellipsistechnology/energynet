package ellipsis.energy.smartgrid.global_local.log;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Collection;
import java.util.Map;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.RealVector;

import ellipsis.energy.calculation.AnalysisResults;
import ellipsis.energy.grid.Bus;
import ellipsis.energy.grid.Unit;
import ellipsis.energy.smartgrid.global_local.LocalController;
import ellipsis.energy.smartgrid.global_local.LocalGlobalOptimiser;
import ellipsis.energy.smartgrid.global_local.LocalGlobalOptimiser.CostToGo;
import ellipsis.util.CSVLogger;

public class ScenarioLogger extends CSVLogger
{
//	private static final String TIME = "Time";
	private static final String LOCAL_ITERATION = "Local iteration";
	private static final String BUS = "Bus";
	private static final String POWER_R = "R{S}";
	private static final String POWER_I = "I{S}";
	private static final String VOLTAGE_ABS = "V_abs";
	private static final String VOLTAGE_ARG = "V_arg";
	
	//// Instance ////
	
	private static ThreadLocal<ScenarioLogger> instances = new ThreadLocal<ScenarioLogger>();
	
	public static ScenarioLogger instance()
	{
		ScenarioLogger instance = instances.get();
		if(instance == null)
		{
			instance = new ScenarioLogger();
			instances.set(instance);
		}
		return instance;
	}
	
	private ScenarioLogger()
	{
		try
		{
			File file = new File("/opt/energynet/scenario.csv");
			file.createNewFile();
			setOut(new PrintStream(file));
			reset();
		}
		catch (IOException e)
		{
			throw new RuntimeException(e);
		}
	}
	
	
	//// Logging ////
	
	private int time;

	public void initialState(AnalysisResults results, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		reset();
		
		Map<String, Integer> busNumbers = results.getBusNumbers();
		for (String busName : busNumbers.keySet())
		{
			Complex power = results.getBusPower(busName);
			Complex voltage = results.getBusVoltage(busName);
			
			addData(busName, power.getReal(), power.getImaginary(), voltage.getArgument(), voltage.abs());
		}
		
		printDataAndReset(logLevel);
	}

	private void reset()
	{
		resetDataWithColumns(BUS, POWER_R, POWER_I, VOLTAGE_ARG, VOLTAGE_ABS);
	}

	public void startLocalOptimisation(int centralIteration, LocalController localController, int logLevel, int horizon)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		// Central iteration/local optimisation header:
		printHeading("Central "+centralIteration, logLevel);
		printCell(localController.toString(), logLevel);
		endRow(logLevel);
		
		// Initialise columns:
		addColumn(LOCAL_ITERATION);
		
		for (Unit unit : localController.getControllableUnits())
		{
			for(int t = 0; t < horizon+1; ++t)
				addColumn("~"+t+"_CTG");
			for(int t = 0; t < horizon+1; ++t)
				addColumn(t+"_CTG");
			for(int t = 0; t < horizon+1; ++t)
				addColumn(t+"_control("+unit.getName()+")");
		}
		
		Collection<Bus> busses = localController.getGrid().getBusses();
		for (Bus bus : busses)
		{
			String busName = bus.getName();
			for(int t = 0; t < horizon+1; ++t)
				addColumn(busName, t, "_Delta R{S} (");
			for(int t = 0; t < horizon+1; ++t)
				addColumn(busName, t, "_Delta I{S} (");
			for(int t = 0; t < horizon+1; ++t)
				addColumn(busName, t, "_Delta V_arg (");
			for(int t = 0; t < horizon+1; ++t)
				addColumn(busName, t, "_Delta V_abs (");
			for(int t = 0; t < horizon+1; ++t)
				addColumn(busName, t, "_R{S} (");
			for(int t = 0; t < horizon+1; ++t)
				addColumn(busName, t, "_I{S} (");
			for(int t = 0; t < horizon+1; ++t)
				addColumn(busName, t, "_V_arg (");
			for(int t = 0; t < horizon+1; ++t)
				addColumn(busName, t, "_V_abs (");
		}
		for(int t = 0; t < horizon+1; ++t)
			addColumn(t+"_Cost");
	}

	private void addColumn(String busName, int t, String middleText)
	{
		addColumn(t+middleText+busName+")");
	}


	public void startLocalIteration(int localIteration, LocalGlobalOptimiser localGlobalOptimiser, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		addData(LOCAL_ITERATION, localIteration);
	}


	public void endLocalIteration(int localIteration, LocalGlobalOptimiser localGlobalOptimiser, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		fillData("-");
	}

	public void logDeltaS(RealVector deltaS, Collection<Bus> busses, AnalysisResults results, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		int size = deltaS.getDimension()/2;
		for (Bus bus : busses)
		{
			String busName = bus.getName();
			Integer index = LocalGlobalOptimiser.indexOf(bus, results);
			double busDeltaS_r = deltaS.getEntry(index);
			double busDeltaS_i = deltaS.getEntry(index+size);
			addData(time+"_Delta R{S} ("+busName+")", busDeltaS_r);
			addData(time+"_Delta I{S} ("+busName+")", busDeltaS_i);
		}
	}


	public void logDeltaV(RealVector deltaV, Collection<Bus> busses, AnalysisResults results, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		int size = deltaV.getDimension()/2;
		for (Bus bus : busses)
		{
			String busName = bus.getName();
			Integer index = LocalGlobalOptimiser.indexOf(bus, results);
			double busDeltaV_arg = deltaV.getEntry(index);
			double busDeltaV_abs = deltaV.getEntry(index+size);
			addData(time+"_Delta V_arg ("+busName+")", busDeltaV_arg);
			addData(time+"_Delta V_abs ("+busName+")", busDeltaV_abs);
		}
	}


	public void logNewState(Complex[] voltages, Collection<Bus> busses, Map<String, Integer> busNumbers, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
//		Map<String, Integer> busNumbers = results.getBusNumbers();
//		for (String busName : busNumbers.keySet())
		for (Bus bus : busses)
		{
			String busName = bus.getName();
			int index = busNumbers.get(busName);
			
			Complex power = bus.netPower();//results.getBusPower(busName);
			Complex voltage = voltages[index];//results.getBusVoltage(busName);
			
			addData(time+"_R{S} ("+busName+")", power.getReal());
			addData(time+"_I{S} ("+busName+")", power.getImaginary());
			
			addData(time+"_V_arg ("+busName+")", voltage.getArgument());
			addData(time+"_V_abs ("+busName+")", voltage.abs());
		}
	}


	public void logCost(double cost, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		addData(time+"_Cost", cost);
	}

	public void logCTG(CostToGo ctg, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		addData("~"+time+"_CTG", ctg.cost);
		if(ctg.controls != null)
			for (Unit unit : ctg.controls.keySet())
			{
				double control = ctg.controls.get(unit).getReal();
				addData(time+"_control("+unit.getName()+")", control);
			}
	}

	public void logActualCTG(CostToGo ctg, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		addData(time+"_CTG", ctg.cost);
	}

	public void endLocalOptimisation(int logLevel)
	{
		printDataAndReset(logLevel);
	}

	
	//// Accessors ////

	public int getTime()
	{
		return time;
	}

	public void setTime(int time)
	{
		this.time = time;
	}
}
