package ellipsis.energy.smartgrid.global_local.log;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;
import java.util.Queue;
import java.util.Stack;

import org.apache.commons.math3.complex.Complex;

import ellipsis.energy.grid.Unit;
import ellipsis.energy.grid.res.SolarPanel;
import ellipsis.energy.smartgrid.ControllableDemand;
import ellipsis.energy.smartgrid.global_local.LocalGlobalOptimiser;
import ellipsis.energy.smartgrid.global_local.LocalGlobalOptimiser.CostToGo;
import ellipsis.util.CSVLogger;

public class LocalGlobalOptimiserLogger extends CSVLogger
{
	public static final String LOGGER_ITERATIONS = "Training Costs";
	public static final String LOGGER_LOCAL_ANALYSIS = "Local Analysis";
	
	
	//// Instance ////
	
	private static ThreadLocal<Map<String, LocalGlobalOptimiserLogger>> instance = new ThreadLocal<Map<String, LocalGlobalOptimiserLogger>>();
	
	private LocalGlobalOptimiserLogger(){}
	
	public static LocalGlobalOptimiserLogger instance(String name)
	{
		Map<String, LocalGlobalOptimiserLogger> loggerMap = instance.get();
		if(loggerMap == null)
		{
			loggerMap = new HashMap<String, LocalGlobalOptimiserLogger>();
			instance.set(loggerMap);
		}
		
		LocalGlobalOptimiserLogger logger = loggerMap.get(name);
		if(logger == null)
		{
			logger = new LocalGlobalOptimiserLogger();
			logger.resetData();
			loggerMap.put(name, logger);
		}
		return logger;
	}
	
	
	//// Logging - Iteratins ////

	public void logCost(int timeIndex, double cost, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		addData(""+timeIndex, cost);
	}
	
	
	////Logging - Local Analysis ////
	
	private Queue<CostToGo> ctgEstimates = new LinkedList<CostToGo>();
	private Stack<CostToGo> ctgActuals = new Stack<CostToGo>();
	private LocalGlobalOptimiser controller;
	private String targetController;

	public void startIteration(int k, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		if(controller.getGrid().getName().equals(targetController))
		{
			printHeading("Iteration,"+k, logLevel);
			
			// Setup columns:
			addColumn("Time");
			addColumn("Cost");
			addData("Rand", "");
			addData("CTG Estimate", "");
			addData("CTG Actual", "");
			addData("CTG Error", "");
			
			// Initial state:
			for (SolarPanel sp : controller.getGrid().get(SolarPanel.class))
				addData(sp.getName(), sp.getPowerOutput().getReal());
			
			for (ControllableDemand cd : controller.getGrid().get(ControllableDemand.class))
				addData(cd.getName(), cd.getChargeRate().getReal());
		}
	}

	public void logRandom(boolean b, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		addData("Rand", b);
	}

	public void logEstimate(CostToGo ctgEstimate, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		ctgEstimates.add(ctgEstimate);
	}

	public void logActual(CostToGo ctg, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		ctgActuals.push(ctg);
	}

	public void setController(LocalGlobalOptimiser controller)
	{
		this.controller = controller;
	}

	public void setTargetController(String controller)
	{
		this.targetController = controller;
	}

	public void printIterations(int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		if(!controller.getGrid().getName().equals(targetController))
		{
			resetData();
			return;
		}
		
		// Setup data:
		while(!ctgEstimates.isEmpty())
		{
			CostToGo estimate = ctgEstimates.poll();
			CostToGo actual = ctgActuals.pop();
			addData("CTG Estimate", estimate.cost);
			addData("CTG Actual", actual.cost);
			addData("CTG Error", actual.cost-estimate.cost);
			
			Map<Unit, Complex> controls = estimate.controls;
			if(controls == null)
			{
				fillData("N/A");
			}
			else
			{
				for (Unit unit : controls.keySet())
					addData(unit.getName(), controls.get(unit).getReal());
			}
		}
		
		// Print:
		printDataAndReset(logLevel);
	}

	@Override
	public void resetData()
	{
		super.resetData();
		ctgEstimates = new LinkedList<CostToGo>();
		ctgActuals = new Stack<CostToGo>();
	}

	public void popEstimate(int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		ctgEstimates.poll();
	}
}
