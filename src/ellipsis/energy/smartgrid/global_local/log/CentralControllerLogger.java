package ellipsis.energy.smartgrid.global_local.log;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.List;
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

public class CentralControllerLogger extends CSVLogger
{
	public static final String LOGGER_COST_AND_ERROR_ANALYSIS = "Cost and Error Analysis";
	public static final String LOGGER_COST_TO_GO_ANALYSIS = "Cost-to-go Analysis";
	
	
	//// Instance ////
	
	private static ThreadLocal<Map<String, CentralControllerLogger>> instance = new ThreadLocal<Map<String, CentralControllerLogger>>();
	
	private CentralControllerLogger(){}
	
	public static CentralControllerLogger instance(String name)
	{
		Map<String, CentralControllerLogger> loggerMap = instance.get();
		if(loggerMap == null)
		{
			loggerMap = new HashMap<String, CentralControllerLogger>();
			instance.set(loggerMap);
		}
		
		CentralControllerLogger logger = loggerMap.get(name);
		if(logger == null)
		{
			logger = new CentralControllerLogger();
			logger.resetData();
			loggerMap.put(name, logger);
			if(name.equals(LOGGER_COST_AND_ERROR_ANALYSIS))
			{
				try
				{
					File file = new File("/opt/energynet/controls_and_errors.csv");
					file.createNewFile();
					logger.setOut(new PrintStream(file));
				}
				catch (IOException e)
				{
					throw new RuntimeException(e);
				}
			}
		}
		return logger;
	}
	

	//// Logging - Cost and Error Analysis ////
	
	public void logError(Bus controlledBus, AnalysisResults results, RealVector error, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		String controlledBusName = controlledBus.getName();
		Integer controlledBusIndex = results.getBusNumbers().get(controlledBusName);
		int slackIndex = results.getSlackIndex();
		if(controlledBusIndex == slackIndex)
			throw new RuntimeException("Controlling the slack bus!");
		if(controlledBusIndex > slackIndex)
			--controlledBusIndex;
		double busError = error.getEntry(controlledBusIndex);
		addData(controlledBusName+" Error", busError);
	}

	/**
	 * Used for Cost and Error Analysis and Voltage Analysis.
	 * @param k
	 * @param t
	 */
	public void logIterationAndTime(int k, int t, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		addData("Iteration", k);
		addData("Time", t);
	}

	public void logControls(Bus controlledBus, List<CostToGo> schedule, int t, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		Map<Unit, Complex> controls;
		if(schedule.size() <= t)
			controls = null;
		else
			controls = schedule.get(t).controls;
		for (Unit unit : controlledBus.getChildren())
		{
			String unitName = unit.getName();
			if(controls == null)
			{
				if(columnData.containsKey(unitName))
					addData(unitName, "Invalid");
			}
			else
			{
				Complex control = controls.get(unit);
				if(control != null)
					addData(unitName, control.getReal());
			}
		}
	}
	
	
	//// Logging - Cost-to-go Analysis ////
	
	private int centralIteration;
	private int minLocalIteration;

	public void startCentralOptimisation(int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		printHeading("Iteration,,Cost-to-go", logLevel);
		resetDataWithColumns("Central", "Local", "Global");
	}

	public void startCentralIteration(int k, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		centralIteration = k;
		minLocalIteration = -1;
	}

	public void endCentralIteration(double costToGo, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		fillData(""); // Finish local iterations
		addData("Central", centralIteration);
		addData("Local", "");
		addData("Global", costToGo);
		fillData("");
	}

	public void startLocalOptimisation(LocalController localController, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
	}

	public void logLocalIteration(LocalGlobalOptimiser localGlobalOptimiser, int k, double costToGo, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		if(k > minLocalIteration) // log iteration numbers only once per central iteration
		{
			minLocalIteration = k;
			addData("Central", centralIteration);
			addData("Local", k);
		}
		
		addData(localGlobalOptimiser.getGrid().getName(), costToGo);
	}

	public void logOptions(
			int adpIterations, 
			double gamma,
			double bandwidth, 
			double forgettingFactor, 
			int history, 
			boolean centralIterationsWithMemory,
			boolean individualExpectationApproximations, 
			boolean random,
			double alpha, 
			boolean resetExpectationApproximations, 
			boolean localIterationsWithMemory, 
			boolean errorsWithMemory, 
			double beta, 
			int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		printHeading("Options:", logLevel);
		
		addData("Local Iterations", adpIterations);
		addData("Alpha", alpha);
		addData("Beta", beta);
		addData("Gamma", gamma);
		addData("Bandwidth", bandwidth);
		addData("Forgetting Factor", forgettingFactor);
		addData("History", history);
		addData("Central Memory", centralIterationsWithMemory);
		addData("Local Memory", localIterationsWithMemory);
		addData("Error Memory", errorsWithMemory);
		addData("Individual Approximations", individualExpectationApproximations); 
		addData("Reset Approximations", resetExpectationApproximations);
		addData("Random", random);
		printDataAndReset(logLevel);
	}
	
	
	//// Logging - Voltage Analysis ////
	
	
}
