package ellipsis.energy.smartgrid.global_local;

import static ellipsis.energy.smartgrid.global_local.log.LocalGlobalOptimiserLogger.LOGGER_ITERATIONS;
import static ellipsis.energy.smartgrid.global_local.log.LocalGlobalOptimiserLogger.LOGGER_LOCAL_ANALYSIS;
import static ellipsis.util.CSVLogger.LOG_LEVEL_DEBUG;
import static ellipsis.util.CSVLogger.LOG_LEVEL_ERROR;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Random;
import java.util.Set;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexUtils;
import org.apache.commons.math3.linear.ArrayRealVector;

import com.mls.util.Util;

import ellipsis.energy.calculation.AnalysisResults;
import ellipsis.energy.calculation.LoadFlowAnalyser;
import ellipsis.energy.grid.Bus;
import ellipsis.energy.grid.Capacitor;
import ellipsis.energy.grid.Grid;
import ellipsis.energy.grid.Load;
import ellipsis.energy.grid.Unit;
import ellipsis.energy.grid.res.RenewableGenerator;
import ellipsis.energy.grid.res.SolarPanel;
import ellipsis.energy.smartgrid.ControllableDemand;
import ellipsis.energy.smartgrid.ControllableDemand.CDLoad;
import ellipsis.energy.smartgrid.global_local.log.CentralControllerLogger;
import ellipsis.energy.smartgrid.global_local.log.ControlLogger;
import ellipsis.energy.smartgrid.global_local.log.FeederVoltageLogger;
import ellipsis.energy.smartgrid.global_local.log.InvalidStateLogger;
import ellipsis.energy.smartgrid.global_local.log.LocalControllerLogger;
import ellipsis.energy.smartgrid.global_local.log.LocalGlobalOptimiserLogger;
import ellipsis.energy.smartgrid.global_local.log.ScenarioLogger;
import ellipsis.energy.smartgrid.global_local.log.ScheduleLogger;
import ellipsis.energy.smartgrid.global_local.log.ValidSolutionLogger;
import ellipsis.util.CSVLogger;
import ellipsis.util.NonParametricApproximator;
import ellipsis.util.TreeCache;
import ellipsis.util.TreeCache.LazyCallback;

public class LocalGlobalOptimiser
{
	public static final int ADP_ITERATIONS = 30;
	protected static final boolean CACHE = false;
	protected static final boolean INDIVIDUAL_EXPECTATION_APPROXIMATIONS = true;
	protected static final boolean RESET_EXPECTATION_APPROXIMATIONS = true; // this needs resetting if we are applying controls each central iteration
	protected static final boolean CENTRAL_ITERATIONS_WITH_MEMORY = true;
	protected static final boolean LOCAL_ITERATIONS_WITH_MEMORY = false;
	protected static final boolean ERRORS_WITH_MEMORY = false;
	protected static final boolean ERROR_DELTAS_WITH_MEMORY = true;
	protected static final double BETA = 0.2; // Voltage error memory rate
	protected static final double ALPHA = 0.3; // Central iteration memory rate
	protected static final boolean RANDOM = true;
	protected static final boolean RESET_RAND2 = false;
	protected static final boolean STATS = false;
	protected static final boolean STATS_VERBOSE = false;
	protected static final boolean LOG_TIMING = false;

	private static final int CONTROL_RESOLUTION = 6;
	
	// Non-parametric approximation parameters:
	private static final double MAX_CTG = 1;//0.05;
	protected static final double GAMMA = 0.8;
	protected static final double FORGETTING_FACTOR = 0;//0.5;//.2; FIXME FORGETTING_FACTOR
	protected static final double BANDWIDTH = 1.3e6/4;
	protected static final int HISTORY = ADP_ITERATIONS;
	
	protected Grid grid;
	protected double basePower;
	protected double baseVoltage;
	protected double minVoltage = 0.95;
	protected double maxVoltage = 1.05;
	
	protected List<Unit> controllableUnits;
    protected Collection<Bus> busses;
    protected Collection<ControllableDemand> cds;
	ArrayList<Load> loads;
	
	protected TreeCache<Double> costCache = new TreeCache<Double>();
	
    private double lastCostToGo;
    
    protected int centralIteration;

	protected boolean viableSolution;
    
    // Approximate DP:
    private NonParametricApproximator[] expectationApproximations;
    private boolean expectationApproximationInitialised = false;
    protected int adpIteration;
	protected boolean finalADPCostToGo = false;
	protected List<CostToGo> previousCostToGoSchedule;
	private List<CostToGo> previousADPControlSchedule;
	private List<CostToGo> costToGoSchedule;

    // STATS:
    private int ctgCount = 0;
    private int minCount = 0;
    private int intCount = 0;
    private int costCount = 0;
    private int ctgDepth = 0;
    private int maxCtgDepth = 0;
    protected Timer timer = new Timer();

    protected final String TIME_COST = "cost";
    protected final String TIME_CTG = "ctg";
    protected final String TIME_MIN = "min";
    protected final String TIME_INT = "int";
    
    String timerID_costoToGo_cost() {return "costToGo-"+grid.getName()+"-cost";}
    String timerID_costoToGo_controls() {return "costToGo-"+grid.getName()+"-controls";}
    String timerID_costoToGo_ctgApproximation() {return "costToGo-"+grid.getName()+"-ctgApproximation";}
    String timerID_costoToGo_approximation() {return "costToGo-"+grid.getName()+"-ctgApproximation";}
    String timerID_costoToGo_approximation_findMinimum() {return "costToGo-"+grid.getName()+"-ctgApproximation-findMinimum";}
    
    String timerID_costoToGo_approximation_findMinimum_value() {return "costToGo-"+grid.getName()+"-ctgApproximation-findMinimum-value";}
    String timerID_costoToGo_approximation_random() {return "costToGo-"+grid.getName()+"-ctgApproximation-random";}
    String timerID_iterations_costToGo() {return "iterations-"+grid.getName()+"-costToGo";}
    String timerID_adp_initApproximations() {return "adp-"+grid.getName()+"-initApproximations";}
    String timerID_storeState() {return "storeState";}
    String timerID_restoreState() {return "restoreState";}
    
    protected void startTimer(String id)
    {
    	com.mls.util.Timer.getTimer(id).start();
    }
    
    protected void pauseTimer(String id)
    {
    	com.mls.util.Timer.getTimer(id).stop();
    }
    
    protected void logTime(String id, boolean seconds)
    {
    	if(!LOG_TIMING)
    		return;
    	
    	com.mls.util.Timer t = com.mls.util.Timer.getTimer(id);
		long time = t.getTimeElapsed();
    	if(seconds)
    		time /= 1000;
		System.out.println("Time for "+id+": "+time+(seconds?"s":"ms"));
		t.reset();
    }
    
    protected static class Timer
    {
        protected Map<String, Long> times = new HashMap<String, Long>();
        protected Map<String, Long> starts = new HashMap<String, Long>();
        
    	public void start(String key)
    	{
    		starts.put(key, System.currentTimeMillis());
    	}
    	
    	public void stop(String key)
    	{
    		long start = starts.put(key, Long.valueOf(0));
    		Long previousTime = times.get(key);
    		if(previousTime == null)
    			previousTime = Long.valueOf(0);
    		times.put(key, System.currentTimeMillis()-start + previousTime);
    	}
    	
    	public long getTime(String key)
    	{
    		return times.get(key);
    	}
    }
    
	public LocalGlobalOptimiser(Grid grid, double basePower, double baseVoltage)
	{
		this.grid = grid;
		this.basePower = basePower;
		this.baseVoltage = baseVoltage;
		
		resetControllableUnits();
		
		busses = grid.getBusses();
		cds = grid.get(ControllableDemand.class);
		
		loads = new ArrayList<Load>();
		for (Load load : grid.getLoads()) 
		{
			if(!(load instanceof CDLoad || load instanceof Capacitor))
				loads.add(load);
		}
	}

	protected void resetControllableUnits() 
	{
		controllableUnits = new ArrayList<Unit>();
		controllableUnits.addAll(grid.get(RenewableGenerator.class));
		controllableUnits.addAll(grid.get(ControllableDemand.class));
	}
	
	public void optimiseControls(List<Forecast> forecasts)
    {
		String timerID = "optimiseControls-"+grid.getName();
		startTimer(timerID);
		
		if(STATS)
		{
			ctgCount = 0;
		    minCount = 0;
		    intCount = 0;
		    costCount = 0;
		    ctgDepth = 0;
		    maxCtgDepth = 0;
		    timer = new Timer();
		}
		
		CostToGo costToGo;
		try
		{
			costToGo = adp(forecasts);
		}
		catch(NoSolutionException e)
		{
			e.printStackTrace();
			this.viableSolution = false;
			costToGoSchedule = previousCostToGoSchedule;
			return;
		}
		catch(InvalidStateException e)
		{
			e.printStackTrace();
			this.viableSolution = false;
			costToGoSchedule = previousCostToGoSchedule;
			return;
		}
    	
    	if(costToGo == null)
    		return;

    	lastCostToGo = costToGo.cost;
    	
    	if(STATS)
    		logStats();
    	
    	logTime(timerID, false);

		previousCostToGoSchedule = costToGoSchedule;
    }

	protected CostToGo adp(List<Forecast> forecasts)
	{
		LocalGlobalOptimiserLogger iterationsLogger = LocalGlobalOptimiserLogger.instance(LOGGER_ITERATIONS);
		LocalGlobalOptimiserLogger bus5Logger = LocalGlobalOptimiserLogger.instance(LOGGER_LOCAL_ANALYSIS);
		bus5Logger.setController(this);
		bus5Logger.setTargetController("Local (Bus 5)");
		ScenarioLogger scenarioLogger = ScenarioLogger.instance();
		
		if(RESET_RAND2)
			rand2 = new Random(13);
		
		// Initialise value function approximations:
		if(!expectationApproximationInitialised || RESET_EXPECTATION_APPROXIMATIONS)
		{
			startTimer(timerID_adp_initApproximations());
			expectationApproximations = new NonParametricApproximator[forecasts.size()];
			
			for (int t = 0; t < expectationApproximations.length; t++)
			{
				expectationApproximations[t] = new NonParametricApproximator();
				expectationApproximations[t].setDefaultValue(MAX_CTG);
				expectationApproximations[t].setForgettingFactor(FORGETTING_FACTOR);
//				expectationApproximations[t].setBandwidth(BANDWIDTH);
				expectationApproximations[t].setHistorySize(HISTORY);
			}
			expectationApproximationInitialised = true;
			pauseTimer(timerID_adp_initApproximations());
		}

		Map<ControllableDemand, Double> capacities = new HashMap<ControllableDemand, Double>();
		Map<Unit, Complex> powers = new HashMap<Unit, Complex>();
		Map<SolarPanel, Double> irradiance = new HashMap<SolarPanel, Double>();
		Map<Load, Complex> loadMap = new HashMap<Load, Complex>();
		storeState(capacities, powers, irradiance, loadMap);
		
		bus5Logger.resetData();
		
		// Train value function approximation:
		CostToGo costToGo;
		System.out.print("\tLocal Iterations for "+grid.getName()+": ");
		for (this.adpIteration = 0; this.adpIteration < ADP_ITERATIONS; this.adpIteration++)
		{
			System.out.print(this.adpIteration+",");
			bus5Logger.startIteration(this.adpIteration, CSVLogger.LOG_LEVEL_DEBUG);
			iterationsLogger.addData("Iteration\\Time", this.adpIteration);
			scenarioLogger.startLocalIteration(this.adpIteration, this, CSVLogger.LOG_LEVEL_DEBUG);

			this.viableSolution = true;
			previousADPControlSchedule = costToGoSchedule;
			costToGoSchedule = new ArrayList<CostToGo>();
			
	    	try
	    	{
	    		startTimer(timerID_iterations_costToGo());
	    		costToGo = costToGo(forecasts, 0);
	    		pauseTimer(timerID_iterations_costToGo());
	    	}
	    	catch(NoSolutionException e)
	    	{
	    		this.viableSolution = false;
	    		startTimer(timerID_iterations_costToGo());
	    		costToGo = costToGo(forecasts, 0);
	    		pauseTimer(timerID_iterations_costToGo());
	    	}

		    restoreState(capacities, powers, irradiance, loadMap);

			bus5Logger.printIterations(LOG_LEVEL_DEBUG);
			scenarioLogger.endLocalIteration(this.adpIteration, this, CSVLogger.LOG_LEVEL_DEBUG);
			ValidSolutionLogger.instance().logLocalValid(centralIteration, adpIteration, getGrid().getName(), this.viableSolution, CSVLogger.LOG_LEVEL_DEBUG);
		}
		iterationsLogger.printHeading(grid.getName(), LOG_LEVEL_DEBUG);
		iterationsLogger.printData(LOG_LEVEL_DEBUG);
		
		// Final run:
		scenarioLogger.startLocalIteration(ADP_ITERATIONS, this, CSVLogger.LOG_LEVEL_DEBUG);

		this.viableSolution = true;		
		previousADPControlSchedule = costToGoSchedule;
		costToGoSchedule = new ArrayList<CostToGo>();

		this.finalADPCostToGo = true;
		System.out.println(ADP_ITERATIONS);
		startTimer(timerID_iterations_costToGo());
		costToGo = costToGo(forecasts, 0);
		pauseTimer(timerID_iterations_costToGo());

		restoreState(capacities, powers, irradiance, loadMap);
		this.finalADPCostToGo  = false;

		// Final Logging:
		scenarioLogger.endLocalIteration(ADP_ITERATIONS, this, CSVLogger.LOG_LEVEL_DEBUG);
		iterationsLogger.resetData();
		ScheduleLogger.instance().printDataAndContinue(LOG_LEVEL_DEBUG, true);
		ValidSolutionLogger.instance().logLocalValid(centralIteration, ADP_ITERATIONS, getGrid().getName(), this.viableSolution, CSVLogger.LOG_LEVEL_DEBUG);
		LocalControllerLogger.instance(grid.getName()).printDataAndReset(LOG_LEVEL_DEBUG);
		
		logTime(timerID_iterations_costToGo(), false);
		logTime(timerID_costoToGo_cost(), false);
		logTime(timerID_costoToGo_controls(), false);
		logTime(timerID_costoToGo_ctgApproximation(), false);
		logTime(timerID_costoToGo_approximation(), false);
		logTime(timerID_costoToGo_approximation_findMinimum(), false);

		if(LOG_TIMING)
			System.out.println(
					"nextCount: "+nextCount+"\n" +
					"findMinimumCount: "+findMinimumCount+"\n" +
					"whileTrueCount: "+whileTrueCount+"\n" +
					"invalidSolutionCount: "+invalidSolutionCount);
		nextCount = 0;
		findMinimumCount = 0;
		whileTrueCount = 0;
		invalidSolutionCount = 0;
		
		logTime(timerID_costoToGo_approximation_findMinimum_value(), false);
		logTime(timerID_costoToGo_approximation_random(), false);
		logTime(timerID_adp_initApproximations(), false);
		logTime(timerID_storeState(), false);
		logTime(timerID_restoreState(), false);
		
		return costToGo;
	}

	protected void restoreState(
			Map<ControllableDemand, Double> capacities,
			Map<Unit, Complex> powers, 
			Map<SolarPanel, Double> irradiance, 
			Map<Load, Complex> loadMap)
	{
		startTimer(timerID_restoreState());
		
		// Reset CD capacities:
		for (ControllableDemand cd : capacities.keySet()) 
		{
			cd.setCapacity(capacities.get(cd));
		}
		
		// Reset powers:
		for (Unit unit : powers.keySet())
		{
			Complex power = powers.get(unit);
			if(unit instanceof SolarPanel)
			{
				SolarPanel solarPanel = (SolarPanel)unit;
				solarPanel.setPowerOutput(power);
				Double E = irradiance.get(solarPanel);
				solarPanel.setIrradiance(E);
			}
			else if(unit instanceof ControllableDemand)
				((ControllableDemand)unit).setChargeRate(power);
		}
		
		// Reset loads:
		for (Load load : loadMap.keySet())
		{
			load.setLoad(loadMap.get(load));
		}
		
		pauseTimer(timerID_restoreState());
	}

	protected void storeState(
			Map<ControllableDemand, Double> capacities,
			Map<Unit, Complex> powers, 
			Map<SolarPanel, Double> irradiance, 
			Map<Load, Complex> loadMap)
	{
		startTimer(timerID_storeState());
		
		// Store current CD capacities:
		for (Unit unit : controllableUnits) 
		{
			if(unit instanceof ControllableDemand)
			{
				ControllableDemand cd = (ControllableDemand)unit;
				capacities.put(cd, cd.getCapacity());
			}
		}
		
		// Store current powers and irradiance:
		for (Unit unit : controllableUnits)
		{
			Complex power;
			if(unit instanceof SolarPanel)
			{
				SolarPanel solarPanel = (SolarPanel)unit;
				power = solarPanel.getPowerOutput();
				irradiance.put(solarPanel, solarPanel.getIrradiance());
			}
			else if(unit instanceof ControllableDemand)
				power = ((ControllableDemand)unit).getChargeRate();
			else
				throw new RuntimeException("Invalid unit found: "+unit);
			
			powers.put(unit, power);
		}
		
		// Store loads:
		for (Load load : loads)
		{
			loadMap.put(load, load.getLoad());
		}
		
		pauseTimer(timerID_storeState());
	}

	public static class CostToGo
    {
    	public double cost;
    	public Map<Unit, Complex> controls;
    	
		public void setControls(List<Unit> controllableUnits, List<Complex> currentControl) 
		{
			controls = new LinkedHashMap<Unit, Complex>();
			for(int i = 0; i < currentControl.size(); ++i)
			{
				controls.put(controllableUnits.get(i), currentControl.get(i));
			}
		}
		
		@Override
		public boolean equals(Object obj)
		{
		    if(!(obj instanceof CostToGo) || obj == null)
		        return false;

		    CostToGo ctg2 = (CostToGo)obj;
		    
		    for (Unit unit : controls.keySet())
            {
                if(!controls.get(unit).equals(ctg2.controls.get(unit)))
                    return false;
            }
		    
		    return true;
		}
		
		@Override
		public int hashCode()
		{
		    int hash = 0;
            
            for (Complex c : controls.values())
            {
                hash ^= c.hashCode();
            }
            
            return hash;
		}
		
		@Override
		public String toString()
		{
			StringBuffer sb = new StringBuffer();
			for (Unit unit : controls.keySet())
			{
				sb.append(unit.getName());
				sb.append("=>");
				sb.append(controls.get(unit));
				sb.append(", ");
			}
			return sb.toString();
		}

		public CostToGo copy()
		{
			CostToGo newCTG = new CostToGo();
			newCTG.cost = this.cost;
			newCTG.controls = new HashMap<Unit, Complex>(this.controls);
			
			return newCTG;
		}
    }

//	FIXME 1 costToGo(): Whatever is coming back from this is giving me the increasing local CTGs.
	private CostToGo costToGo(final List<Forecast> forecasts, final int timeIndex) throws InvalidStateException
	{
		ScenarioLogger scenarioLogger = ScenarioLogger.instance();
		scenarioLogger.setTime(timeIndex);
		
		if(STATS)
		{
			++ctgCount;
			++ctgDepth;
			if(ctgDepth > maxCtgDepth)
				maxCtgDepth = ctgDepth;
			timer.start(TIME_CTG);
		}
		
// FIXME 1 Debug
ControllableDemand storage = (ControllableDemand) grid.getUnit("Bus 21 CD");
if(timeIndex == 12 && storage != null)
{
	double rate = storage.getChargeRate().getReal();
	if(rate == -650e3)
		Util.nullop();
}

		startTimer(timerID_costoToGo_cost());
		double cost = cost(timeIndex);
		pauseTimer(timerID_costoToGo_cost());

		scenarioLogger.logCost(cost, CSVLogger.LOG_LEVEL_DEBUG);
		
		LocalGlobalOptimiserLogger bus5Logger = LocalGlobalOptimiserLogger.instance(LOGGER_LOCAL_ANALYSIS);
		bus5Logger.addData("Time", timeIndex);
		bus5Logger.addData("Cost", cost);

		LocalGlobalOptimiserLogger.instance(LOGGER_ITERATIONS).logCost(timeIndex, cost, CSVLogger.LOG_LEVEL_DEBUG);
		
		// If at the end of the forecasts then return final cost:
		if(timeIndex == forecasts.size())
		{
			CostToGo ctg = new CostToGo();
			ctg.cost = cost;

			if(STATS)
			{
				--ctgDepth;
				timer.stop(TIME_CTG);
			}
			return ctg;
		}
		
		// Get a list of controls (changes in power output) for each controllable unit:
    	final Map<Unit, List<Complex>> availableControls = new HashMap<Unit, List<Complex>>();
    	Forecast forecast = forecasts.get(timeIndex);
    	startTimer(timerID_costoToGo_controls());
    	for (Unit unit : controllableUnits) 
    	{
    		// Find control range:
    		Complex maxPower;
			Complex minPower;
    		if(unit instanceof RenewableGenerator)
    		{
    			RenewableGenerator dg = (RenewableGenerator)unit;
				maxPower = getExpectedMaximumPower(dg, forecast);
				minPower = getExpectedMinimumPower(dg, forecast);
    		}
    		else if(unit instanceof ControllableDemand)
    		{
    			ControllableDemand cd = (ControllableDemand)unit;
				maxPower = getExpectedMaximumPower(cd);
// FIXME 1 Debug
if(timeIndex == 11 && adpIteration == 9 && maxPower.getReal() <= 0)
	Util.nullop();
				minPower = getExpectedMinimumPower(cd);
    		}
    		else
    		{
    			throw new RuntimeException(unit+" is not a RenewableGenerator or ControllableDemand");
    		}
    			
			// Calculate control options:
			List<Complex> options = new ArrayList<Complex>(CONTROL_RESOLUTION);
			Complex stepSize = maxPower.subtract(minPower).divide(CONTROL_RESOLUTION);
			Complex power = minPower;
			for(int i = 0; i <= CONTROL_RESOLUTION; ++i) 
			{
				if(power.getReal() > maxPower.getReal() || power.getImaginary() > maxPower.getImaginary())
					power = maxPower; // this is required to overcome rounding issues
				options.add(power);
				power = power.add(stepSize);
			}
			
			availableControls.put(unit, options);
		}
    	pauseTimer(timerID_costoToGo_controls());
		
		if(STATS)
		{
			++minCount;
			timer.stop(TIME_CTG);
			timer.start(TIME_MIN);
		}

		CostToGo minCTG;
		startTimer(timerID_costoToGo_ctgApproximation());
		minCTG = costToGoExpectationApproximation(availableControls, timeIndex, forecasts);
		pauseTimer(timerID_costoToGo_ctgApproximation());
		
		if(STATS)
		{
			timer.stop(TIME_MIN);
			timer.start(TIME_CTG);
		}
		
		minCTG.cost = cost + GAMMA*minCTG.cost;

		if(STATS)
		{
			if(STATS_VERBOSE)
				logStats();
			
			--ctgDepth;
			timer.stop(TIME_CTG);
		}
		
		if(timeIndex == 11)
		{
			CentralControllerLogger costToGoAnalysisLogger = CentralControllerLogger.instance(CentralControllerLogger.LOGGER_COST_TO_GO_ANALYSIS);
			costToGoAnalysisLogger.logLocalIteration(this, this.adpIteration, minCTG.cost, CSVLogger.LOG_LEVEL_DEBUG);
		}

		return minCTG;
	}

	private void logStats() {
		System.out.print("Depth="+ctgDepth+",\tMax Depth="+maxCtgDepth+"," +
			"\tctgCount="+ctgCount+",\tminCount="+minCount+",\tintCount="+intCount+",");
		for (String key : timer.times.keySet()) 
		{
			System.out.print("\t"+key+"="+timer.getTime(key)+",");
		}
		System.out.println();
	}

	protected Complex getExpectedMinimumPower(ControllableDemand cd) 
	{
		double capacity = cd.getCapacity();
		
		// No charge left to discharge:
		if(capacity <= 0)
			return Complex.ZERO;
		
		// Almost empty, so limit allowed discharge rate:
		double maxDischargeRate = cd.getMaxDischargeRate();
		if(capacity < maxDischargeRate)
			return new Complex(-capacity);
		
		// Enough capacity left to discharge:
		return new Complex(-maxDischargeRate);
	}

	protected Complex getExpectedMaximumPower(ControllableDemand cd) 
	{
		double capacity = cd.getCapacity();
		double maxCapacity = cd.getMaxCapacity();
		
		// Charged fully, can't charge any more:
		if(capacity >= maxCapacity)
			return Complex.ZERO;
		
		// Almost fully charged, limit allowed charge rate:
		double roomLeft = maxCapacity-capacity;
		double maxChargeRate = cd.getMaxChargeRate();
		if(maxChargeRate > roomLeft)
			return new Complex(roomLeft);
		
		// Enough room left to charge at max charge rate:
		return new Complex(maxChargeRate);
	}

	protected Complex getExpectedMinimumPower(RenewableGenerator dg, Forecast forecast) 
	{
		return Complex.ZERO;
	}

	protected Complex getExpectedMaximumPower(RenewableGenerator dg, Forecast forecast) 
	{
		return dg.getExpectedMaximumPower(forecast);
	}

	/**
	 * 
	 * @return The cost of the grid in its current state (i.e. slack bus abs power) 
	 * or 0 if the state is invalid.
	 */
    public double cost(final int timeIndex)
    {
    	if(STATS) timer.start(TIME_COST);
		if(STATS) ++costCount;
		LazyCallback<Double> callback = new LazyCallback<Double>() {
			@Override
			public Double value() {
				LoadFlowAnalyser lfa = new LoadFlowAnalyser(grid);
				lfa.setBasePower(basePower);
				lfa.setBaseVoltage(baseVoltage);
				lfa.setIterations(100);
				lfa.setTargetError(0.000001);
				
				AnalysisResults results = lfa.analyse();
				
				FeederVoltageLogger.instance().logTrueVoltage(timeIndex, results, CSVLogger.LOG_LEVEL_DEBUG);
				
				// Check voltages: 
				for (Bus bus : busses) 
				{
					double v = results.getBusVoltage(bus.getName()).abs();
					if(v < minVoltage || v > maxVoltage)
					{
						InvalidStateLogger.instance().logInvalidCentralVoltage(centralIteration, timeIndex, bus, v, LOG_LEVEL_ERROR);
					}
				}

				String slackBusName = grid.getSlackBus().getName();
				Complex busPower = results.getBusPower(slackBusName);
				double slackPower = busPower.abs()*Math.signum(busPower.getReal());

				return slackPower;
			}
		};
		
		if(CACHE)
		{
			Object[] state = state();
			Double slackPower = costCache.lazyGet(callback, state);
			if(STATS) timer.stop(TIME_COST);
			
			return slackPower;
		}
		else
		{
			if(STATS) timer.stop(TIME_COST);
			
			return callback.value();
		}
	}
	
	
	//////////////////////////////////////
	//// Value Function Approximation ////
	//////////////////////////////////////

	private static Random rand = new Random(0);
	double randomGaussian()
	{
		if(RANDOM)
			return rand.nextGaussian();
		else
			return 0;
	}
	
static int whileTrueCount;
static int invalidSolutionCount;
	private CostToGo costToGoExpectationApproximation(
    		Map<Unit, List<Complex>> availableControls,  
    		int timeIndex, 
    		List<Forecast> forecasts)
	{
		LocalGlobalOptimiserLogger bus5Logger = LocalGlobalOptimiserLogger.instance(LOGGER_LOCAL_ANALYSIS);
		ScenarioLogger scenarioLogger = ScenarioLogger.instance();
		LocalControllerLogger.instance(grid.getName()).logState(grid.getName(), centralIteration, adpIteration, timeIndex, state(), LOG_LEVEL_DEBUG);
		
		Forecast forecast = forecasts.get(timeIndex);

		// Loop until a valid control is found:
		Set<CostToGo> invalidControls = new HashSet<CostToGo>();
		List<CostToGo> validControls = new ArrayList<CostToGo>();

		while(true)
		{
++whileTrueCount;
			// Variables used in catch statement:
			Map<Load, Complex> oldLoads = null;
			double oldIrradiance = 0;
			Map<ControllableDemand, Double> oldCapacities = null;
			CostToGo ctgEstimate = null;
			
			// Get approximate cost-to-go and corresponding controls for the current state:
			CostToGo originalEstimate;
			try
			{
				startTimer(timerID_costoToGo_approximation());
				originalEstimate = approximateExpectation(availableControls, invalidControls, validControls, forecasts, timeIndex);
				pauseTimer(timerID_costoToGo_approximation());
			}
			catch(NoSolutionException e)
			{
				if(finalADPCostToGo)
					throw e;
				
				viableSolution = false;
				CostToGo ctg = new CostToGo();
				ctg.cost = MAX_CTG;
				
				// Fill in cost schedule for remaining times from previous schedule:
				ctg.controls = previousADPControlSchedule.get(timeIndex).controls;
				this.costToGoSchedule.add(ctg);
				for (int t = timeIndex+1; t < costToGoSchedule.size(); t++)
				{
					costToGoSchedule.add(previousADPControlSchedule.get(t));
				}
				
				return ctg;
			}
			ctgEstimate = originalEstimate.copy();
			combineControlsWithPrevious(timeIndex, ctgEstimate);
			bus5Logger.logEstimate(ctgEstimate, CSVLogger.LOG_LEVEL_DEBUG);
		
			// Apply the chosen controls:
			applyControls(ctgEstimate, timeIndex);

			this.costToGoSchedule.add(ctgEstimate);
				
			// Apply a random sample of the stochastic variables:
			double irradianceMean = forecast.getIrradianceMean();
//			double irradianceSD = forecast.getIrradianceStandardDeviation();
			double loadMean = forecast.getLoadMean();
//			double loadSD = forecast.getLoadStandardDeviation();
			oldLoads = new HashMap<Load, Complex>();
			Map<Load, Complex> newLoads = new HashMap<Load, Complex>();
			for (Load load : this.loads) 
			{
				oldLoads.put(load, load.getLoad());
				
				double abs = loadMean;//FIXME 1 + randomGaussian()*loadSD;
				double arg = load.getLoad().getArgument();
				Complex power = ComplexUtils.polar2Complex(abs, arg);
				newLoads.put(load, power);
			}
			
			double irradiance = irradianceMean;//FIXME 1 + randomGaussian()*irradianceSD;
			oldIrradiance = setStochasticState(newLoads, irradiance);
	
			// Tick over CD capacities (and store current values for restore below):
			oldCapacities = passTime();

			try
			{
				// Estimate the cost-to-go for the next state:
				CostToGo ctg = costToGo(forecasts, timeIndex+1);
				double ctgNext;
				if(ctg == null || Double.isNaN(ctg.cost)) // invalid state so give it a high cost
					ctgNext = MAX_CTG;
				else
					ctgNext = ctg.cost;
				
				ControlLogger.instance().logControls(grid.getName(), centralIteration, adpIteration, timeIndex, originalEstimate.controls, ctgEstimate.controls, LOG_LEVEL_DEBUG);
				ControlLogger.instance().printDataAndContinue(LOG_LEVEL_DEBUG, false);
				
				// Logging:
				bus5Logger.logActual(ctg, CSVLogger.LOG_LEVEL_DEBUG);
				ScheduleLogger.instance().logCentralIteration(centralIteration, CSVLogger.LOG_LEVEL_DEBUG);
				ScheduleLogger.instance().logLocalIteration(adpIteration, CSVLogger.LOG_LEVEL_DEBUG);
				ScheduleLogger.instance().logTime(timeIndex, CSVLogger.LOG_LEVEL_DEBUG);
				ScheduleLogger.instance().logControl(this, ctgEstimate, CSVLogger.LOG_LEVEL_DEBUG);
				scenarioLogger.setTime(timeIndex);
				scenarioLogger.logActualCTG(ctg, CSVLogger.LOG_LEVEL_DEBUG);
				LocalControllerLogger.instance(grid.getName()).logControls(originalEstimate, ctgEstimate, LOG_LEVEL_DEBUG);
				
				// Reset the state to after controls have been applied, and before 
				// stochastic variables have been applied and update the estimate:
				setStochasticState(oldLoads, oldIrradiance);
				applyCapacities(oldCapacities);
				applyControls(ctgEstimate, timeIndex);
// FIXME 1 logging
//if(timeIndex == 11 && (adpIteration == 0 || adpIteration == 30))
//	Util.nullop();
//if(timeIndex == 11/* && adpIteration == 30*/)
//{
//	Complex control = ctgEstimate.controls.get(cds.iterator().next());
//	double real = control.getReal();
//	double im = control.getImaginary();
//	System.out.println(real+","+im+","+ctgNext);
//}
				updateApproximateCostToGo(ctgNext, timeIndex);
				
				CostToGo ctgSelected = new CostToGo();
				ctgSelected.controls = ctgEstimate.controls;
				ctgSelected.cost = ctgNext;

				return ctgSelected;
			}
			catch(InvalidStateException e)
			{
++invalidSolutionCount;
				// Remember the invalid control:
				invalidControls.add(originalEstimate);
				if(!validControls.isEmpty())
					validControls.remove(0);

				InvalidStateLogger.instance().logInvlaidState(centralIteration, adpIteration, timeIndex, ctgEstimate, CSVLogger.LOG_LEVEL_DEBUG);
				
				// Unlog:
				bus5Logger.popEstimate(CSVLogger.LOG_LEVEL_DEBUG);
				
				// Remove control from schedule:
				this.costToGoSchedule.remove(this.costToGoSchedule.size()-1);
				
				// Reset state and try again:
				setStochasticState(oldLoads, oldIrradiance);
				applyCapacities(oldCapacities);
				continue;
			}
		}
	}

	/**
	 * 
	 * @return The capacities prior to passing time.
	 */
	protected Map<ControllableDemand, Double> passTime()
	{
		Map<ControllableDemand, Double> capacities = new HashMap<ControllableDemand, Double>();
		for (ControllableDemand cd : this.cds)
		{
			capacities.put(cd, cd.getCapacity());
			cd.passTime(1);
		}
		return capacities;
	}

	protected void applyCapacities(Map<ControllableDemand, Double> capacities)
	{
		for (ControllableDemand cd : this.cds)
		{
			cd.setCapacity(capacities.get(cd));
		}
	}

	private double setStochasticState(Map<Load, Complex> loads, double irradiance)
	{
		double oldIrradiance = -1;
		
		for (Load load : loads.keySet())
		{
			Complex power = loads.get(load);
			load.setLoad(power);
		}
		int size = controllableUnits.size();
		for (int i = 0; i < size; ++i) 
		{
			Unit unit = controllableUnits.get(i);
			if(unit instanceof SolarPanel)
			{
				SolarPanel solarPanel = (SolarPanel)unit;
				if(oldIrradiance == -1)
					oldIrradiance = solarPanel.getIrradiance();
				solarPanel.setIrradiance(irradiance);
			}
		}
		
		return oldIrradiance;
	}

	static Random rand2 = new Random(13);
	private CostToGo approximateExpectation(final Map<Unit, List<Complex>> availableControls, Set<CostToGo> invalidControls, List<CostToGo> validControls, final List<Forecast> forecasts, final int timeIndex)
	{
		@SuppressWarnings("unused")
		double alpha = Math.exp(-5*(this.adpIteration)/ADP_ITERATIONS); // From 100% probability to almost 0%.
				//1.0/(0.5*this.adpIteration+1.0);
		LocalGlobalOptimiserLogger bus5Logger = LocalGlobalOptimiserLogger.instance(LOGGER_LOCAL_ANALYSIS);
		ScenarioLogger scenarioLogger = ScenarioLogger.instance();
		CostToGo minCTG;

		startTimer(timerID_costoToGo_approximation_findMinimum());
		minCTG = approximateExpectation_findMinimum(availableControls, invalidControls, validControls, timeIndex, bus5Logger);
		pauseTimer(timerID_costoToGo_approximation_findMinimum());

		// FIXME 1 Removed random control addition
//		startTimer(timerID_costoToGo_approximation_random());
//		minCTG = randomiseMinCTG(minCTG, 1-alpha, availableControls, timeIndex);
//		pauseTimer(timerID_costoToGo_approximation_random());

		scenarioLogger.logCTG(minCTG, CSVLogger.LOG_LEVEL_DEBUG);

		return minCTG;
	}

	CostToGo randomiseMinCTG(CostToGo minCTG, double alpha, Map<Unit, List<Complex>> availableControls, int timeIndex)
	{
		// Random CTG:
		CostToGo randCTG = new CostToGo();
		randCTG.controls = new HashMap<Unit, Complex>();
		for (Unit unit : availableControls.keySet())
		{
			List<Complex> controls = availableControls.get(unit);
			int controlCount = controls.size();
			int r = (int)(rand2.nextDouble()*controlCount);
			if(r == controlCount)
				r = controlCount - 1;
			Complex control = controls.get(r);
			randCTG.controls.put(unit, control);
		}
		
		// Combine CTGs:
		CostToGo newCTG = new CostToGo();
		newCTG.controls = new HashMap<Unit, Complex>();
		for (Unit unit : minCTG.controls.keySet())
		{
			Complex minControls = minCTG.controls.get(unit);
			Complex randControls = randCTG.controls.get(unit);
			Complex newControls = minControls.multiply(alpha).add(randControls.multiply(1-alpha));
			newCTG.controls.put(unit, newControls);
		}
		
		// Calculate cost:
		applyControls(newCTG, timeIndex);
		
		// Tick over CD capacities:
		Map<ControllableDemand, Double> capacities = passTime();
		
		// Get cost to go:
		Complex[] state = state();
		ArrayRealVector vX = trainingInput(state);
		NonParametricApproximator expectationApproximation = expectationApproximation(timeIndex);
		newCTG.cost = expectationApproximation.value(vX);
		
		// Restore capacities:
		applyCapacities(capacities);
		
		return newCTG;
	}

static int nextCount = 0;
static int findMinimumCount = 0;
	private CostToGo approximateExpectation_findMinimum(
			final Map<Unit, List<Complex>> availableControls,
			Set<CostToGo> invalidControls, 
			List<CostToGo> validControls,
			final int timeIndex,
			LocalGlobalOptimiserLogger bus5Logger)
	{
++findMinimumCount;
		CostToGo minCTG;
		bus5Logger.logRandom(false, CSVLogger.LOG_LEVEL_DEBUG);
		
		if(!validControls.isEmpty())
		{
			return validControls.get(0);
		}
		
		// Choose the best by iterating over all control combinations:
		Iterator<CostToGo> iterator = new Iterator<CostToGo>()
		{
			List<Unit> units = new ArrayList<Unit>(availableControls.keySet());
			int valueIndeces[] = new int[availableControls.size()];
			boolean hasNext = true;
			
			@Override
			public void remove() {}
			
			@Override
			public CostToGo next()
			{
				if(!hasNext)
					throw new NoSuchElementException();
++nextCount;
				// Setup and apply controls:
				CostToGo ctg = new CostToGo();
				ctg.controls = new HashMap<Unit, Complex>();
				for (int i = 0; i < units.size(); i++)
				{
					Unit unit = units.get(i);
					List<Complex> controls = availableControls.get(unit);
					Complex control = controls.get(valueIndeces[i]);
					ctg.controls.put(unit, control);
				}
if(timeIndex == 11 && adpIteration == 30)
	Util.nullop();
		        applyControls(ctg, timeIndex);

		        // Tick over capacities:
				Map<ControllableDemand, Double> capacities = passTime();
				
				// Ensure there is a bandwidth:
				Complex[] state = state();
				NonParametricApproximator expectationApproximation = expectationApproximation(timeIndex);
				double[] bw = expectationApproximation.getBandwidth();
				if(bw == null)
				{
					bw = Util.doubleArray(state.length*2, BANDWIDTH);
					expectationApproximation.setBandwidth(bw);
				}
				
				// Get cost to go:
				ArrayRealVector vX = trainingInput(state);
				startTimer(timerID_costoToGo_approximation_findMinimum_value());
				ctg.cost = expectationApproximation.value(vX);
				pauseTimer(timerID_costoToGo_approximation_findMinimum_value());
				
				// Restore capacities:
				applyCapacities(capacities);
				
				if(Double.isNaN(ctg.cost))
					throw new RuntimeException("Invalid cost: NaN");
				
				// Increment:
				increment(0);
				
				return ctg;
			}

			private void increment(int i)
			{
				if(i == valueIndeces.length)
				{
					hasNext = false;
					return;
				}
				
				Unit unit = units.get(i);
				List<Complex> controls = availableControls.get(unit);
				int controlsSize = controls.size();
				
				if(valueIndeces[i] == controlsSize-1) // at the end
				{
					valueIndeces[i] = 0;
					increment(i+1);
				}
				else
				{
					valueIndeces[i] += 1;
				}
			}
			
			@Override
			public boolean hasNext()
			{
				return hasNext;
			}
		};
		
		minCTG = new CostToGo();
		minCTG.cost = Double.MAX_VALUE;
// FIXME 1 Logging
if(timeIndex == 11 && adpIteration%5 == 4)
	Util.nullop();
		while (iterator.hasNext())
		{
			CostToGo ctg = (CostToGo) iterator.next();
		
			if(contains(invalidControls, ctg))
				continue;
			
			validControls.add(ctg);

			if(ctg.cost < minCTG.cost)
			{
				minCTG.cost = ctg.cost;
				minCTG.controls = new HashMap<Unit, Complex>(ctg.controls);
			}
		}
// FIXME 1 Logging
if(timeIndex == 11 && adpIteration%5 == 4)
	Util.nullop();
		
		Collections.sort(validControls, new Comparator<CostToGo>()
		{
			@Override
			public int compare(CostToGo o1, CostToGo o2)
			{
				return Double.compare(o1.cost, o2.cost);
			}
		});
		
		if(minCTG.controls == null)
			throw new NoSolutionException();
		
		return minCTG;
	}

	private boolean contains(Set<CostToGo> invalidControls, CostToGo ctg)
	{
//if(adpIteration == 30)
//{
//		for (CostToGo costToGo : invalidControls)
//		{
//			if(costToGo.equals(ctg))
//				return true;
//		}
//		return false;
//}
//else
//{
	return invalidControls.contains(ctg);
//}
	}

    public void combineControlsWithPrevious(final int timeIndex, CostToGo ctg)
    {
    	if(ctg.controls == null)
    		return;
    	
    	if(LOCAL_ITERATIONS_WITH_MEMORY)
    	{
    		if(adpIteration > 0)
    		{
    			Map<Unit, Complex> previousControls = this.previousADPControlSchedule.get(timeIndex).controls;
    			combineControlsWithPrevious(ctg, previousControls, 1/(1+adpIteration/4));//1-adpIteration/((double)ADP_ITERATIONS+1)); 
    		}
    	}
    	
        if(CENTRAL_ITERATIONS_WITH_MEMORY)
		{
			if (centralIteration > 0)
			{
				Map<Unit, Complex> previousControls = previousCostToGoSchedule.get(timeIndex).controls;
				combineControlsWithPrevious(ctg, previousControls, ALPHA);
			}
		}
    }

	private void combineControlsWithPrevious(CostToGo ctg, Map<Unit, Complex> previousControls, double gamma)
	{
		if(previousControls != null)
		{
			for (Unit unit : ctg.controls.keySet())
			{
				Complex previousControl = previousControls.get(unit);
				Complex currentControl = ctg.controls.get(unit);
				ctg.controls.put(unit, 
						previousControl.multiply((1-gamma)).
						add(
						currentControl.multiply(gamma)));
			}
		}
	}

    private NonParametricApproximator expectationApproximation(int timeIndex)
	{
		int index;
		if(INDIVIDUAL_EXPECTATION_APPROXIMATIONS)
			index = timeIndex;
		else
			index = 0;
		NonParametricApproximator expectationApproximation = expectationApproximations[index];
		return expectationApproximation;
	}

	protected void applyControls(CostToGo ctg, int time)
	{
		if(time == -1)
			throw new RuntimeException("Invalid time -1");
		
		if(ctg.controls == null)
			return;
		
		for (int i = 0; i < controllableUnits.size(); ++i) 
        {
            Unit unit = controllableUnits.get(i);
            Complex control = ctg.controls.get(unit);
			if(unit instanceof RenewableGenerator)
            	((RenewableGenerator)unit).setPowerOutput(control);
            else if(unit instanceof ControllableDemand)
            	((ControllableDemand)unit).setChargeRate(control);
        }
	}

	private void updateApproximateCostToGo(double updatedCostToGo, int timeIndex)
	{
        NonParametricApproximator expectationApproximation = expectationApproximation(timeIndex);
        
		// Add sample:
        ArrayRealVector vX = trainingInput(state());
		expectationApproximation.addSample(vX, updatedCostToGo);
        
        // Tune bandwidth:
//		double heuristicBW = 1.0/expectationApproximation.sampleCount();
//		expectationApproximation.setBandwidth(Util.doubleArray(vX.getDimension(), heuristicBW));
//		BandwidthSelector bs = new BandwidthSelector(expectationApproximation);
//		bs.crossValidate();
	}

	private ArrayRealVector trainingInput(Complex[] state)
	{
		int stateCount = state.length*2;
		int stateSize = stateCount;
		double trainingState[] = new double[stateSize];

		for (int i = 0; i < state.length; ++i)
		{
			double real = state[i].getReal();
			double imaginary = state[i].getImaginary();
			trainingState[2*i] = real;///basePower;
			trainingState[2*i+1] = imaginary;///basePower;
		}
		
		return new ArrayRealVector(trainingState);
	}

	/**
	 * Builds a representation of state from the power at each bus and the capacity of all {@link ControllableDemand}s.
	 * {@link ControllableDemand}s.
	 */
	private Complex[] state() 
	{
		Complex[] state = new Complex[busses.size()+cds.size()];
		int i = 0;
		for (Bus bus : busses)
		{
			state[i++] = bus.netPower();
		}
		for (ControllableDemand cd : cds)
		{
			state[i++] = new Complex(cd.getCapacity());
		}
		return state;
	}
	
	
	
	//// Helpers ////

	public static int indexOf(Bus bus, AnalysisResults r)
	{
		int n = r.getBusNumbers().get(bus.getName());
        int i;
        if(n > r.getSlackIndex())
            i = n - 1;
        else
            i = n;
        
        return i;
	}
	

	//// Accessors ////

	public Grid getGrid() {
		return grid;
	}

	public double getBasePower() {
		return basePower;
	}

	public double getBaseVoltage() {
		return baseVoltage;
	}
	
	public double getLastCostToGo()
    {
        return lastCostToGo;
    }
	
	public List<CostToGo> getCostToGoSchedule()
	{
		return costToGoSchedule;
	}

	public int getCentralIteration()
	{
		return centralIteration;
	}

	public void setCentralIteration(int centralIteration)
	{
		this.centralIteration = centralIteration;
	}

	@Override
	public String toString()
	{
		return "Controller["+grid.getName()+"]";
	}

	public List<Unit> getControllableUnits()
	{
		return controllableUnits;
	}

	public boolean isViableSolution()
	{
		return viableSolution;
	}
}
