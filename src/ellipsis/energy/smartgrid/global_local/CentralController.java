package ellipsis.energy.smartgrid.global_local;

import static ellipsis.util.CSVLogger.LOG_LEVEL_DEBUG;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexUtils;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.FieldVector;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import com.mls.util.ThreadPool;

import ellipsis.energy.calculation.AnalysisResults;
import ellipsis.energy.calculation.LoadFlowAnalyser;
import ellipsis.energy.grid.Bus;
import ellipsis.energy.grid.Capacitor;
import ellipsis.energy.grid.DistributedSource;
import ellipsis.energy.grid.Grid;
import ellipsis.energy.grid.Line;
import ellipsis.energy.grid.Load;
import ellipsis.energy.grid.Source;
import ellipsis.energy.grid.Unit;
import ellipsis.energy.grid.res.RenewableGenerator;
import ellipsis.energy.grid.res.SolarPanel;
import ellipsis.energy.grid.res.WindTurbine;
import ellipsis.energy.smartgrid.ControllableDemand;
import ellipsis.energy.smartgrid.global_local.log.CentralControllerLogger;
import ellipsis.energy.smartgrid.global_local.log.ErrorLogger;
import ellipsis.energy.smartgrid.global_local.log.FeederVoltageLogger;
import ellipsis.energy.smartgrid.global_local.log.PartitioningLogger;
import ellipsis.energy.smartgrid.global_local.log.ScenarioLogger;
import ellipsis.energy.smartgrid.global_local.log.ScheduleLogger;
import ellipsis.energy.smartgrid.global_local.log.ValidSolutionLogger;
import ellipsis.util.CSVLogger;

public class CentralController extends LocalGlobalOptimiser
{
    @SuppressWarnings("all")
	private static final boolean TEST = false;

//    private static final double MAX_V_ARG_ERROR = Math.PI/16;
//    private static final double MAX_V_ABS_ERROR = 0.009;//0.012;
    private static final int MAX_PARTITION_UNITS = 4;
    private static final int MAX_PARTITION_SIZE = 7;

	private static final int ITERATIONS = 20;

    private AnalysisResults results;
    
    public CentralController(Grid grid, double basePower, double baseVoltage)
    {
    	super(grid, basePower, baseVoltage);
    }
    
    /**
     * Perform power flow analysis and obtain network state.
     */
    public void analyseNetwork()
    {
        this.results = analyse(this.grid);
    }

    public AnalysisResults analyse(Grid grid)
    {
        LoadFlowAnalyser lfa = new LoadFlowAnalyser(grid);
        lfa.setBasePower(basePower);
        lfa.setBaseVoltage(baseVoltage);
        lfa.setIterations(100);
        lfa.setTargetError(1e-6);
        
        AnalysisResults results = lfa.analyse();
        if(!results.getDidConverge())
        	throw new RuntimeException("Did not converge!");
        return results;
    }
    
    public Collection<LocalController> defineSubnetworks(List<Forecast> forecasts)
    {
        Collection<LocalController> subnets = new LinkedHashSet<LocalController>();
        int slackIndex = results.getSlackIndex();
        Map<String, Integer> busNumbers = results.getBusNumbers();
        int size = busses.size()-1; // -1 for slack bus
        
        // 1. Inverse Jacobian:
        RealMatrix jacobianInverse = jacobianInverse();

        // 2. Calculate expected maximum shifts in power at each bus:
        double deltaPQMinArray[] = new double[size*2];
        double deltaPQMaxArray[] = new double[size*2];
        
        calculatePowerShifts(forecasts.get(0), size, deltaPQMinArray, deltaPQMaxArray);
        
//        ArrayRealVector deltaPQMin = new ArrayRealVector(deltaPQMinArray, false);
        ArrayRealVector deltaPQMax = new ArrayRealVector(deltaPQMaxArray, false);
        
        // 3. Calculated corresponding expected maximum shift in voltage:
//        RealVector deltaVMin = jacobianInverse.operate(deltaPQMin);
        RealVector deltaVMax = jacobianInverse.operate(deltaPQMax);
        
        // Get relationship between slack bus and adjacent busses:
        FieldVector<Complex> slackAdmittances = results.getAdmittanceMatrix().getAdmittances().getRowVector(slackIndex);

        for (Bus bus : busses)
        {
            // Skip slack bus and skip non-generator busses:
            if(bus.getSlackVoltage().abs() > 0 || (bus.getSources().isEmpty() && bus.getControllableDemands().isEmpty()))
                continue;

            // Get index, allowing for slack bus:
            int i = busNumbers.get(bus.getName());
            if(slackIndex < i)
                i -= 1;
            
            // 4. Approximate voltage shift by removing low impact busses:
            double[] vArgRow = jacobianInverse.getRow(i);
            double[] vAbsRow = jacobianInverse.getRow(i+size);
            
            double[] vArgTerms = new double[size*2];
            double[] vAbsTerms = new double[size*2];
            
            // Get all the arguments:
            for(int j = 0; j < size; ++j)
            {
                vArgTerms[j] =                                       // arg(v_i) = sum(
                        vArgRow[j]*deltaPQMax.getEntry(j) +          // d(arg(v_i))/dP_j*(Delta P_j) +
                        vArgRow[j+size]*deltaPQMax.getEntry(j+size); // d(arg(v_i))/dQ_j*(Delta Q_j) )
                
                vAbsTerms[j] =                                       // |v_i| = sum(
                        vAbsRow[j]*deltaPQMax.getEntry(j) +          // d|v_i|/dP_j*(Delta P_j) +
                        vAbsRow[j+size]*deltaPQMax.getEntry(j+size); // d|v_i|/dQ_j*(Delta Q_j) )
            }
            
            double vArg = deltaVMax.getEntry(i);
            double vAbs = deltaVMax.getEntry(i+size);
            
            @SuppressWarnings("unused")
			double approxVArg = vArg;
            @SuppressWarnings("unused")
			double approxVAbs = vAbs;
//            double argError = 0;
//            double absError = 0;
            HashSet<Integer> removedTerms = new HashSet<Integer>();
            
            PartitioningLogger.instance().logTerms(busNumbers, slackIndex, bus.getName(), vArgRow, vAbsRow, deltaPQMaxArray, vArgTerms, vAbsTerms, CSVLogger.LOG_LEVEL_DEBUG);
            
            // Remove terms until max error is reached:
            HashSet<Bus> subnet = null;
            int count = size;
            while(countControllables(subnet) > MAX_PARTITION_UNITS || count > MAX_PARTITION_SIZE)//argError < MAX_V_ARG_ERROR && absError < MAX_V_ABS_ERROR) 
            {
                int minAbsIndex = -1;
                double minArgTerm = Integer.MAX_VALUE;
                double minAbsTerm = Integer.MAX_VALUE;
                
                // Find min terms:
                for (int j = 0; j < size; j++)
                {
                    if(removedTerms.contains(j))
                        continue;
                    
                    double term = vAbsTerms[j]+vAbsTerms[j+size]; // v shift due to P + Q shifts
                    if(Math.abs(term) < Math.abs(minAbsTerm))
                    {
                        minArgTerm = vArgTerms[j]+vArgTerms[j+size];
                        minAbsTerm = term;
                        minAbsIndex = j;
                    }
                }
                
                // Remove min term:
                if(minAbsIndex != -1)
                {
                	--count;
                	
                    removedTerms.add(minAbsIndex);
                    subnet = createSubnet(slackIndex, busNumbers, bus, removedTerms);
                    
                    // Update error:
                    approxVArg -= minArgTerm;
                    approxVAbs -= minAbsTerm;
//                    argError = Math.abs(vArg - approxVArg);
//                    absError = Math.abs(vAbs - approxVAbs);
                }
                else
                {
                    break;
                }
            }

            // 5. Construct subnetwork:
//            subnet = createSubnet(slackIndex, busNumbers, bus, removedTerms);
            
            // Add virtual slack bus:
//            Bus vsb = addVirtualSlackBus(slackIndex, busNumbers, size, jacobianInverse, subnet);

            // Get copies of all subnetwork units:
            Grid subgrid = grid.subnet(subnet);
            PartitioningLogger.instance().logSubnet(subnet, CSVLogger.LOG_LEVEL_DEBUG);
            
            // Setup virtual bus (remove children and add slack bus):
//            vsb = subgrid.getBus(vsb.getName()); // replace with copy from subgrid
//            for (Iterator<Unit> iterator = vsb.getChildren().iterator(); iterator.hasNext();)
//            {
//                Unit unit = iterator.next();
//                if(!(unit instanceof Line)) // must leave line connecting vsb to subnetwork
//                    iterator.remove();
//            }
//            SlackSource slack = new SlackSource();
//            slack.setName("Virtual slack source");
//            Complex slackVoltage = results.getBusVoltage(vsb.getName()).multiply(baseVoltage);
//            slack.setVoltage(slackVoltage);
//            vsb.addChild(slack);
            
            // Create subnetwork controller:
			LocalController localController = new LocalController(bus.getName(), basePower, baseVoltage, subgrid, jacobianInverse, results, busNumbers);
//			localController.setGlobalSlackBusVoltage(grid.getSlackBus().getSlackVoltage());
			localController.setSlackAdmittances(slackAdmittances);
			localController.setSlackAdjacentJacobianInverse(slackAdjacentJacobianInverse(localController, jacobianInverse));
			subnets.add(localController);
            
            // Find the max removed term to see what the upper bound of \epsilon is (i.e. \Delta u_{b,~B} = \epsilon R for all R_i < 1):
            double epsilon = 0;
            for (int j = 0; j < size; j++)
            {
                if(!removedTerms.contains(j))
                    continue;
                
                double term = vAbsTerms[j]+vAbsTerms[j+size]; // v shift due to P + Q shifts
                if(Math.abs(term) > epsilon)
                {
                    epsilon = Math.abs(term);
                }
            }
            System.out.println("Epsilon for "+localController+" is "+epsilon);

	        //analyseEigenSystem(jacobianInverse, localController);
        }
        
        return subnets;
    }

	private int countControllables(HashSet<Bus> subnet)
	{
		if(subnet == null)
			return Integer.MAX_VALUE;
		
		int count = 0;
		for (Bus bus : subnet)
		{
			for (Unit unit : bus.getChildren())
			{
				if(unit instanceof ControllableDemand || unit instanceof SolarPanel || unit instanceof WindTurbine)
					++count;
			}
		}
		
		return count;
	}

	private HashSet<Bus> createSubnet(int slackIndex,
			Map<String, Integer> busNumbers, Bus bus,
			HashSet<Integer> removedTerms)
	{
		HashSet<Bus> subnet = new HashSet<Bus>();
		subnet.add(bus);
		for (Bus b : busses)
		{
		    int index = busNumbers.get(b.getName());
		    if(slackIndex == index)
		        continue;
		    if(slackIndex < index)
		        index -= 1;
		    if(!removedTerms.contains(index))
		    {
		        subnet.add(b);
		    }
		}
		return subnet;
	}

	private RealMatrix slackAdjacentJacobianInverse(LocalController localController, RealMatrix jacobianInverse)
	{
		// Get slack adjacent busses:
		Set<Bus> adjacents = new HashSet<Bus>();
		Bus slackBus = grid.getSlackBus();
		for (Line line : slackBus.getLines())
		{
			Bus other = line.getFromBus();
			if(other == slackBus)
				other = line.getToBus();
			adjacents.add(other);
		}
		
		// Zero out irrelevant entries:
		jacobianInverse = jacobianInverse.copy();
		int size = jacobianInverse.getRowDimension()/2;
		Collection<Bus> localBusses = localController.getGrid().getBusses();
	    for (Bus bus : busses)
	    {
	        int i = indexOf(bus);
	        
	        if(!localBusses.contains(bus))
	        {
	        	for (int j = 0; j < jacobianInverse.getRowDimension(); j++)
				{
	        		jacobianInverse.setEntry(j, i, 0);
	        		jacobianInverse.setEntry(j, i+size, 0);
				}
	        }
	        if(!adjacents.contains(bus))
	        {
	        	for (int j = 0; j < jacobianInverse.getColumnDimension(); j++)
				{
	        		jacobianInverse.setEntry(i, j, 0);
	        		jacobianInverse.setEntry(i+size, j, 0);
				}
	        }
	    }
	    
	    return jacobianInverse;
	}

	/**
	 * Sets any rows not belonging to the local network and any columns not
	 * belonging to the external network to zero.
	 * @param localController
	 * @param A
	 */
	private void externalInfluenceOnly(LocalController localController, RealMatrix A)
	{
		int size = A.getRowDimension()/2;
		Collection<Bus> localBusses = localController.getGrid().getBusses();
        for (Bus busi : busses)
        {
            int n = results.getBusNumbers().get(busi.getName());
            int i;
            if(n > results.getSlackIndex())
                i = n - 1;
            else
                i = n;
            
            if(localBusses.contains(busi))
            {
            	for (int j = 0; j < A.getRowDimension(); j++)
				{
					A.setEntry(j, i, 0);
					A.setEntry(j, i+size, 0);
				}
            }
            else
            {
            	for (int j = 0; j < A.getColumnDimension(); j++)
				{
					A.setEntry(i, j, 0);
					A.setEntry(i+size, j, 0);
				}
            }
            
            /*for (Bus busj : busses)
            {
                n = results.getBusNumbers().get(busj.getName());
                
                // Zero out any partial derivatives not connecting external power changes to internal voltages:
                if(!localBusses.contains(busi) || localBusses.contains(busj))
                {
                    int j;
                    if(n > results.getSlackIndex())
                        j = n - 1;
                    else
                        j = n;
                    
                    A.setEntry(i+size, j,        0);
                    A.setEntry(i+size, j+size, 0);
                    A.setEntry(i,      j,        0);
                    A.setEntry(i,      j+size, 0);
                }
            }*/
        }
	}

    public RealMatrix jacobianInverse() {
		DecompositionSolver solver = new LUDecomposition(results.getJacobian()).getSolver();
        RealMatrix jacobianInverse = solver.getInverse();
		return jacobianInverse;
	}

    /**
     * Calculates the minimum as max load + min dg, and the maximum as min load with max dg.
     * All in per unit.
     * @param nextTimePeriod
     * @param size
     * @param deltaPQMinArray
     * @param deltaPQMaxArray
     */
    private void calculatePowerShifts(Forecast nextTimePeriod, int size, double[] deltaPQMinArray, double[] deltaPQMaxArray)
    {
        for (Bus bus : busses)
        {
            // Skip slack bus:
            if(bus.getSlackVoltage().abs() > 0)
                continue;
            
            // Get index, allowing for slack bus:
            int i = results.getBusNumbers().get(bus.getName());
            if(results.getSlackIndex() < i)
                i -= 1;
            
            // Min and max expected DG output:
            Complex minDG = Complex.ZERO;
            Complex maxDG = Complex.ZERO;
            
            for (Source source : bus.getSources())
            {
                if(source instanceof RenewableGenerator)
                {
                    RenewableGenerator rg = (RenewableGenerator)source;
                    Complex min_i = rg.getExpectedMinimumPower(nextTimePeriod);
                    Complex max_i = rg.getExpectedMaximumPower(nextTimePeriod);
                    
                    minDG = minDG.add(min_i);
                    maxDG = maxDG.add(max_i);
                }
            }
            
            // Min and max expected loads:
            Complex minLoad = Complex.ZERO;
            Complex maxLoad = Complex.ZERO;
            
            for (Load load : bus.getLoads())
            {
                if(load instanceof Capacitor)
                    continue;
                
                Complex min_i = load.getExpectedMinimumPower(nextTimePeriod);
                Complex max_i = load.getExpectedMaximumPower(nextTimePeriod);
                
                minLoad = minLoad.add(min_i);
                maxLoad = maxLoad.add(max_i);
            }
            
            // Convert all powers to per unit:
            minDG = minDG.divide(basePower);
            maxDG = maxDG.divide(basePower);
            minLoad = minLoad.divide(basePower);
            maxLoad = maxLoad.divide(basePower);
            Complex presentDGPower = bus.getDistributedSourcePower().divide(basePower);
            Complex presentLoadPower = bus.getLoadPower().divide(basePower);
            
            // Get total present power output at this bus:
            double presentDGActive = presentDGPower.getReal();
            double presentDGReactive = presentDGPower.getImaginary();
            double presentLoadActive = presentLoadPower.getReal();
            double presentLoadReactive = presentLoadPower.getImaginary();
            
            double presentActive = presentDGActive - presentLoadActive;
            double presentReactive = presentDGReactive - presentLoadReactive;
            
            // Calculate power shifts at this bus:
            deltaPQMinArray[i] = (minDG.getReal() - maxLoad.getReal()) - presentActive; // Delta P
            deltaPQMinArray[i+size] = (minDG.getImaginary() - maxLoad.getImaginary()) - presentReactive; // Delta Q
            
            deltaPQMaxArray[i] = (maxDG.getReal() - minLoad.getReal()) - presentActive; // Delta P
            deltaPQMaxArray[i+size] = (maxDG.getImaginary() - minLoad.getImaginary()) - presentReactive; // Delta Q
        }
    }
    
    
    //// Local/Global Optimisation ////

public static LocalController[] sortedSubnets(Collection<LocalController> subnets) {
	LocalController[] aSubnets = subnets.toArray(new LocalController[subnets.size()]);
    Arrays.sort(aSubnets, new Comparator<LocalController>() {
		@Override
		public int compare(LocalController o1, LocalController o2) {
			return o1.getControlledBus().getName().compareTo(o2.getControlledBus().getName());
		}
	});
	return aSubnets;
}

	public void optimiseLocals(Collection<LocalController> subnets, final List<Forecast> forecasts) 
	{
	    startTimer("optimiseLocals");
		CentralControllerLogger costAnalysisLogger = CentralControllerLogger.instance(CentralControllerLogger.LOGGER_COST_AND_ERROR_ANALYSIS);
		final CentralControllerLogger costToGoAnalysisLogger = CentralControllerLogger.instance(CentralControllerLogger.LOGGER_COST_TO_GO_ANALYSIS);
		costToGoAnalysisLogger.startCentralOptimisation(CSVLogger.LOG_LEVEL_DEBUG);
		final ScenarioLogger scenarioLogger = ScenarioLogger.instance();
		
		Map<ControllableDemand, Double> capacities = new HashMap<ControllableDemand, Double>();
		Map<Unit, Complex> powers = new HashMap<Unit, Complex>();
		Map<SolarPanel, Double> irradianceMap = new HashMap<SolarPanel, Double>();
		Map<Load, Complex> loadMap = new HashMap<Load, Complex>();
		storeState(capacities, powers, irradianceMap, loadMap);
		
		// Get the initial powers:
		analyseNetwork();
		RealVector initialPowers = globalPowers();
		RealVector initialVoltages = globalVoltages();
		scenarioLogger.initialState(results, CSVLogger.LOG_LEVEL_DEBUG);
		RealMatrix initialJacobianInverse = jacobianInverse();
		
		// Iterate:
		for(int k = 0; k < ITERATIONS; ++k)
		{
		    timer.start(""+k);
		    System.out.println("Central Iteration "+k);
		    
		    this.centralIteration = k;
		    this.viableSolution = true;

		    costToGoAnalysisLogger.startCentralIteration(k, CSVLogger.LOG_LEVEL_DEBUG);

		    LocalController[] aSubnets = sortedSubnets(subnets);

			// Get all the local controllers' suggested controls:
			for (final LocalController localController : aSubnets)
			{
ThreadPool.getInstance().setDisallowed(true);
				ThreadPool.getInstance().queueTask(new Runnable()
				{
					@Override
					public void run()
					{
						String timerName = localController.getGrid().getName()+":"+centralIteration;
						timer.start(timerName);
						scenarioLogger.startLocalOptimisation(centralIteration, localController, CSVLogger.LOG_LEVEL_DEBUG, forecasts.size());
						
						localController.setCentralIteration(centralIteration);
						localController.optimiseControls(forecasts);
						if(!localController.isViableSolution())
							CentralController.this.viableSolution = false;
						
						scenarioLogger.endLocalOptimisation(CSVLogger.LOG_LEVEL_DEBUG);
						timer.stop(timerName);
						
						costToGoAnalysisLogger.startLocalOptimisation(localController, CSVLogger.LOG_LEVEL_DEBUG);
					}
				});
			}

			ThreadPool.getInstance().waitForAll();
			
			// Update baseline (PHD-17):
			//if(this.viableSolution)
			{
				restoreState(capacities, powers, irradianceMap, loadMap);
				applyScenarioIfViable(subnets, forecasts, 0); // Set state to after variables and controls have been applied for time t=0
				storeState(capacities, powers, irradianceMap, loadMap);
				analyseNetwork();
				initialPowers = globalPowers();
				initialVoltages = globalVoltages();
				initialJacobianInverse = jacobianInverse();
				
				for (LocalController localController : aSubnets)
				{
					localController.update(initialJacobianInverse, results);
				}
			}
			ValidSolutionLogger.instance().logCentralValid(this.viableSolution, CSVLogger.LOG_LEVEL_DEBUG);
			
			// 2014-03-24 - Need to update controls from neighbouring controllers:
			previousCostToGoSchedule = new ArrayList<>();
			boolean first = true;
			for (LocalController localController : aSubnets)
			{
				List<CostToGo> schedule = localController.getCostToGoSchedule(); // The schedule is only from local conroller's controlled bus's units.
				int t = 0;
				for (CostToGo ctg : schedule)
				{
					CostToGo newCtg;
					if(first)
					{
						newCtg = new CostToGo();
						newCtg.controls = new HashMap<>();
						previousCostToGoSchedule.add(newCtg);
					}
					else
					{
						newCtg = previousCostToGoSchedule.get(t);
						++t;
					}
					
					for (ControllableDemand cd : localController.getControlledBus().getControllableDemands())
					{
						if(newCtg.controls.containsKey(cd))
							throw new RuntimeException("Two local controllers are controlling the same unit!");
						
						Complex value = ctg.controls.get(cd);
						newCtg.controls.put(cd, value);
					}
					
					for (DistributedSource dg : localController.getControlledBus().getRenewableSources())
					{
						if(newCtg.controls.containsKey(dg))
							throw new RuntimeException("Two local controllers are controlling the same unit!");
						
						Complex value = ctg.controls.get(dg);
						newCtg.controls.put(dg, value);
					}
				}
				
				first = false;
			}
			for (LocalController localController : aSubnets)
			{
				localController.previousCostToGoSchedule = previousCostToGoSchedule;
			}
			
			// Get global delta S:
			int forecastCount = forecasts.size();
			RealVector fullDeltaS[] = new RealVector[forecastCount];
			int size = busses.size()-1;
			calculateGlobalDeltaS(forecasts, costAnalysisLogger, capacities,
					powers, irradianceMap, loadMap, initialPowers,
					initialVoltages, initialJacobianInverse, k, aSubnets,
					forecastCount, fullDeltaS, size);
			
			// Update local errors:
			for (LocalController localController : aSubnets)
			{
				RealMatrix jacobianInverse = localController.getJacobianInverse();
			
				RealVector errors[] = new RealVector[forecastCount];

			    for (int t = 0; t < forecastCount; ++t)
				{
				    // Find error due to external influence over local controller:
					RealVector extDeltaS = fullDeltaS[t].copy();
					externalOnly(localController, extDeltaS);

					// Iterate over external error terms:
					size = extDeltaS.getDimension();
					RealMatrix intVExtS = jacobianInverse.copy();
					externalInfluenceOnly(localController, intVExtS);
					errors[t] = intVExtS.operate(extDeltaS);
					
					FeederVoltageLogger.instance().logRow(centralIteration, t, localController.getControlledBus(), results, errors[t], CSVLogger.LOG_LEVEL_DEBUG);
					costAnalysisLogger.logError(localController.getControlledBus(), results, errors[t], CSVLogger.LOG_LEVEL_DEBUG);
				} // END forecast ITERATIONS

				localController.setVoltageError(errors);
			} // END LocalController ITERATIONS

			// Finalise iteration:
			timer.stop(""+k);
			System.out.println("\tfinished in "+timer.getTime(""+k)+"ms");

            double costToGo = costToGo(subnets, forecasts, capacities, powers, irradianceMap, loadMap);
			costToGoAnalysisLogger.endCentralIteration(costToGo, CSVLogger.LOG_LEVEL_DEBUG);
		} // END ITERATIONS (k)
		
		// Logging:
		costToGoAnalysisLogger.printDataAndReset(LOG_LEVEL_DEBUG);
		costToGoAnalysisLogger.logOptions(
				LocalGlobalOptimiser.ADP_ITERATIONS,
				LocalGlobalOptimiser.GAMMA,
				LocalGlobalOptimiser.BANDWIDTH,
				LocalGlobalOptimiser.FORGETTING_FACTOR,
				LocalGlobalOptimiser.HISTORY,
				LocalGlobalOptimiser.CENTRAL_ITERATIONS_WITH_MEMORY,
				LocalGlobalOptimiser.INDIVIDUAL_EXPECTATION_APPROXIMATIONS,
				LocalGlobalOptimiser.RANDOM,
				LocalGlobalOptimiser.ALPHA,
				LocalGlobalOptimiser.RESET_EXPECTATION_APPROXIMATIONS,
				LocalGlobalOptimiser.LOCAL_ITERATIONS_WITH_MEMORY,
				LocalGlobalOptimiser.ERRORS_WITH_MEMORY,
				LocalGlobalOptimiser.BETA,
				LOG_LEVEL_DEBUG);
		costAnalysisLogger.printDataAndReset(LOG_LEVEL_DEBUG);
		ScheduleLogger.instance().printDataAndReset(LOG_LEVEL_DEBUG);
		FeederVoltageLogger.instance().printData(LOG_LEVEL_DEBUG);
		ErrorLogger.instance().printData(LOG_LEVEL_DEBUG);
		ValidSolutionLogger.instance().printData(LOG_LEVEL_DEBUG);
		
		logTime("Total time for optimiseLocals()", true);
	}

	private void calculateGlobalDeltaS(final List<Forecast> forecasts,
			CentralControllerLogger costAnalysisLogger,
			Map<ControllableDemand, Double> capacities,
			Map<Unit, Complex> powers, Map<SolarPanel, Double> irradianceMap,
			Map<Load, Complex> loadMap, RealVector initialPowers,
			RealVector initialVoltages, RealMatrix initialJacobianInverse,
			int k, LocalController[] aSubnets, int forecastCount,
			RealVector[] fullDeltaS, int size)
	{
		{	
			// Initialise to t = t_0:
			restoreState(capacities, powers, irradianceMap, loadMap);
			RealVector fullS1 = new ArrayRealVector(size*2);

			// Calculate delta Ses:
			for (int t = 0; t < forecastCount; ++t)
			{
			    costAnalysisLogger.logIterationAndTime(k, t, CSVLogger.LOG_LEVEL_DEBUG);

			    // Calculate delta S at current time:
			    for (Bus bus : busses)
				{
					int i = indexOf(bus);
					Complex busPower = bus.netPower().divide(basePower);
					fullS1.setEntry(i, busPower.getReal());
					fullS1.setEntry(i+size, busPower.getImaginary());
				}
				fullDeltaS[t] = fullS1.subtract(initialPowers);
				
				FeederVoltageLogger.instance().logApproximateCentralVoltages(t, results, approximateVoltages(initialJacobianInverse, fullDeltaS[t], initialVoltages), CSVLogger.LOG_LEVEL_DEBUG);
				
				// Apply controls:
				for (LocalController localController : aSubnets)
				{
					costAnalysisLogger.logControls(localController.getControlledBus(), localController.getCostToGoSchedule(), t, CSVLogger.LOG_LEVEL_DEBUG);
					
					localController.applyControls(this.grid, t);
				}
				
				// Apply forecast loads:
				Forecast forecast = forecasts.get(t);
				double forecastLoad = forecast.getLoadMean();
				for (Load load : grid.getLoads())
				{
					if(load instanceof ControllableDemand.CDLoad || load instanceof Capacitor)
						continue;
					Complex l = load.getLoad();
					double abs = forecastLoad;
					double arg = l.getArgument();
					l = ComplexUtils.polar2Complex(abs, arg);
					load.setLoad(ComplexUtils.polar2Complex(abs, arg));
				}
				
				// Update capacities:
				passTime();
			}
		}
	}

    private RealVector approximateVoltages(RealMatrix initialJacobianInverse, RealVector deltaS, RealVector initialVoltages)
	{
		return initialVoltages.add(initialJacobianInverse.operate(deltaS));
	}

	private double costToGo(Collection<LocalController> subnets,
			final List<Forecast> forecasts, Map<ControllableDemand, Double> capacities, Map<Unit, Complex> powers, Map<SolarPanel, Double> irradianceMap, Map<Load, Complex> loadMap)
	{
		// Initialise to state at t = 0:
    	restoreState(capacities, powers, irradianceMap, loadMap);
		
		// Sum costs up to horizon:
		int t;
		Stack<Double> costs = new Stack<Double>();
		for (t = 0; t < forecasts.size(); t++)
		{
			// Cost at time t:
			costs.push(cost(t));
			
			// Apply cotrols and stochastic variables:
			applyScenario(subnets, forecasts, t);
			
			// Update ready for next time period:
			passTime();
		}
		
		costs.push(cost(t)); // terminating cost
		
		double costToGo = 0;
		while(!costs.isEmpty())
		{
			costToGo = costs.pop() + GAMMA*costToGo;
		}
		
		return costToGo;
	}

	private void applyScenario(Collection<LocalController> subnets, List<Forecast> forecasts, int t)
	{
		applyScenario(subnets, forecasts, t, false);
	}
	private void applyScenarioIfViable(Collection<LocalController> subnets, List<Forecast> forecasts, int t)
	{
		applyScenario(subnets, forecasts, t, true);
	}
	private void applyScenario(Collection<LocalController> subnets, List<Forecast> forecasts, int t, boolean onlyIfViable)
	{
		// Set forecast loads for time t: 
		Forecast forecast = forecasts.get(t);
		double forecastLoad = forecast.getLoadMean();
		for (Load load : grid.getLoads())
		{
			if(load instanceof ControllableDemand.CDLoad || load instanceof Capacitor)
				continue;
			Complex l = load.getLoad();
			double abs = forecastLoad;
			double arg = l.getArgument();
			load.setLoad(ComplexUtils.polar2Complex(abs, arg));
		}
		
		// Set forecast irradiance:
		double irradiance = forecast.getIrradianceMean();
		for (SolarPanel sp : grid.get(SolarPanel.class))
		{
			sp.setIrradiance(irradiance);
		}
		
		// Apply controls for time t:
		for (LocalController localController : subnets)
		{
			if(onlyIfViable && !localController.isViableSolution())
				continue;
			
			localController.applyControls(grid, t);
		}
	}

	private int indexOf(Bus bus)
	{
		return indexOf(bus, results);
	}

	/**
	 * Sets any values in the given vector that are in the local grid to zero.
	 * @param localController
	 * @param v
	 */
	private void externalOnly(LocalController localController, RealVector v)
	{
		int size = v.getDimension()/2;
		Collection<Bus> localBusses = localController.getGrid().getBusses();
        for (Bus bus : localBusses)
        {
	        int i = indexOf(bus);
	        v.setEntry(i, 0);      // arg
	        v.setEntry(i+size, 0); // abs
        }
	}

	private RealVector globalPowers()
	{
		int size = busses.size()-1;
		RealVector initialPowers = new ArrayRealVector(size*2);
		for (Bus bus : busses)
		{
			String busName = bus.getName();
			int i = indexOf(bus);
            
			Complex busPower = results.getBusPower(busName);
            initialPowers.setEntry(i, busPower.getReal());           // P
			initialPowers.setEntry(i+size, busPower.getImaginary()); // Q
		}
		return initialPowers;
	}

	private RealVector globalVoltages()
	{
		int size = busses.size()-1;
		RealVector initialVoltages = new ArrayRealVector(size*2);
		for (Bus bus : busses)
		{
			String busName = bus.getName();
			int i = indexOf(bus);
            
			Complex busPower = results.getBusVoltage(busName);
			initialVoltages.setEntry(i, busPower.getArgument()); // arg(v)
			initialVoltages.setEntry(i+size, busPower.abs());    // |v|
		}
		return initialVoltages;
	}

//	private void reset(Map<Unit, Complex> baseLine, Map<Load, Complex> loads, Map<ControllableDemand, Double> capacities)
//    {
//        for (Unit unit : baseLine.keySet())
//        {
//            Complex control = baseLine.get(unit);
//            
//            if(unit instanceof Source)
//                ((Source)unit).setPowerOutput(control);
//            else if (unit instanceof ControllableDemand)
//                ((ControllableDemand)unit).setChargeRate(control);
//        }
//        
//        for (Load load : loads.keySet())
//		{
//			Complex power = loads.get(load);
//			load.setLoad(power);
//		}
//        
//        for (ControllableDemand cd : capacities.keySet())
//        {
//			Double capacity = capacities.get(cd);
//			cd.setCapacity(capacity);
//		}
//    }

    //// ACCESSORS ////

	public Grid getGrid()
    {
        return grid;
    }
    public void setGrid(Grid grid)
    {
        this.grid = grid;
    }

    public AnalysisResults getResults()
    {
        return results;
    }
    
    
    //// Logging ////

	@SuppressWarnings("all")
	private void logDeltas(LocalController localController, Bus slackBus,
			double approximateCost, double trueCost, double costError,
			Complex trueSlackVoltage) {
		Complex approximateSlackVoltage = localController.getGlobalResults().getBusVoltage(slackBus.getName());
		Complex voltageError = trueSlackVoltage.subtract(approximateSlackVoltage);
		System.out.println(localController.getGrid().getName()+": \n" +
				"\t|S("+slackBus+")| = "+trueCost+"\n" +
				"\t|~S("+slackBus+")| = "+approximateCost+"\n" +
				"\tDelta |S| = "+costError+"\n" +
				"\tv("+slackBus+") = "+trueSlackVoltage+"\n" +
				"\t~v("+slackBus+") = "+approximateSlackVoltage+"\n" +
				"\tDelta v("+slackBus+") = "+voltageError);
	}

    @SuppressWarnings("all")
	private HashMap<Integer, String> busNames(Map<String, Integer> busNumbers) 
	{
		HashMap<Integer, String> busNames = new HashMap<Integer, String>();
		for (String name : busNumbers.keySet()) 
		{
			busNames.put(busNumbers.get(name), name);
		}
		
		System.out.print("Index,");
		for (Bus bus2 : busses)
		{
		    if(bus2.getSlackVoltage().abs() > 0)
		        continue;
		    System.out.print(busNumbers.get(bus2.getName()));
		    System.out.print(',');
		}
		for (Bus bus2 : busses)
		{
		    if(bus2.getSlackVoltage().abs() > 0)
		        continue;
		    System.out.print(busNumbers.get(bus2.getName()));
		    System.out.print(',');
		}
		System.out.println();
		System.out.print("Bus,");
		for (Bus bus2 : busses)
		{
		    if(bus2.getSlackVoltage().abs() > 0)
		        continue;
		    System.out.print(bus2);
		    System.out.print(',');
		}
		for (Bus bus2 : busses)
		{
		    if(bus2.getSlackVoltage().abs() > 0)
		        continue;
		    System.out.print(bus2);
		    System.out.print(',');
		}
		System.out.println();
		return busNames;
	}
	
	protected void logJacobian(int slackIndex, int size,
			HashMap<Integer, String> busNames, RealMatrix jacobianInverse) 
	{
		System.out.println("Jacobian inverse:");
		for (int i = 0; i < size; i++) 
		{
			System.out.print(i >= slackIndex ? busNames.get(i+1) : busNames.get(i) );
			System.out.print(',');
			for (int j = 0; j < size; j++) 
			{
				System.out.print(jacobianInverse.getEntry(i, j));
				System.out.print(',');
			}
			for (int j = 0; j < size; j++) 
			{
				System.out.print(jacobianInverse.getEntry(i, j+size));
				System.out.print(',');
			}
			System.out.println();
		}
		for (int i = 0; i < size; i++) 
		{
			System.out.print(i >= slackIndex ? busNames.get(i+1) : busNames.get(i) );
			System.out.print(',');
			for (int j = 0; j < size; j++) 
			{
				System.out.print(jacobianInverse.getEntry(i+size, j));
				System.out.print(',');
			}
			for (int j = 0; j < size; j++) 
			{
				System.out.print(jacobianInverse.getEntry(i+size, j+size));
				System.out.print(',');
			}
			System.out.println();
		}
	}

	protected void logTerms(int slackIndex, int size,
			HashMap<Integer, String> busNames, double[] vArgTerms,
			double[] vAbsTerms, HashSet<Integer> removedTerms) 
	{
		System.out.print("\t,");
		for(int j = 0; j < size; ++j)
		{
		    int index = j+1 > slackIndex ? j+1 : j;
		    System.out.print(busNames.get(index));
		    System.out.print(',');
		}
		System.out.println();
		System.out.print("\tArg Terms:,");
		for(int j = 0; j < size; ++j)
		{
		    System.out.print(vArgTerms[j]);
		    System.out.print(',');
		}
		System.out.println();
		System.out.print("\tAbs Terms:,");
		for(int j = 0; j < size; ++j)
		{
		    System.out.print(vAbsTerms[j]);
		    System.out.print(',');
		}
		System.out.println();
		System.out.print("\tRemoved?,");
		for(int j = 0; j < size; ++j)
		{
		    System.out.print(!removedTerms.contains(j));
		    System.out.print(',');
		}
		System.out.println();
	}
	
	protected void logSubnet(HashSet<Bus> subnet) {
		System.out.print("\tBusses:,");
		for (Bus bus2 : subnet)
		{
		    System.out.print(bus2);
		    System.out.print(',');
		}
		System.out.println();
	}
}
