package ellipsis.energy.smartgrid.global_local;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexUtils;
import org.apache.commons.math3.linear.ArrayFieldVector;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.FieldVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import com.mls.util.Util;

import ellipsis.energy.calculation.AnalysisResults;
import ellipsis.energy.grid.Bus;
import ellipsis.energy.grid.Grid;
import ellipsis.energy.grid.Source;
import ellipsis.energy.grid.Unit;
import ellipsis.energy.grid.res.RenewableGenerator;
import ellipsis.energy.smartgrid.ControllableDemand;
import ellipsis.energy.smartgrid.global_local.log.ControlLogger;
import ellipsis.energy.smartgrid.global_local.log.ErrorLogger;
import ellipsis.energy.smartgrid.global_local.log.FeederVoltageLogger;
import ellipsis.energy.smartgrid.global_local.log.ScenarioLogger;
import ellipsis.util.CSVLogger;
import ellipsis.util.TreeCache.LazyCallback;

public class LocalController extends LocalGlobalOptimiser
{
    private Bus slackBus;
	private RealMatrix jacobianInverse;
	private Map<String, Integer> globalBusNumbers;
	private Bus controlledBus;
	private AnalysisResults globalResults;
	private RealVector globalVoltages;
	private RealVector globalPowers;
//	private Map<Unit, Complex> controls;
//	private Map<Unit, Complex> means;
//	private Map<Unit, Complex> standardDeviations;
//	private Complex globalSlackVoltage;
	private FieldVector<Complex> slackAdmittances;
	private RealMatrix slackAdjacentJacobianInverse;
	
	// Error term values:
// Some of these terms are probably not needed any more.
//	private RealVector[] errorMean;
//	private RealVector[] errorSD;
	private RealVector[] voltageError;
	private RealVector[] voltageErrorDelta;
//	private RealVector[] voltageMargin;

	public LocalController(String controlledBus, double basePower, double baseVoltage, Grid subnetwork, RealMatrix jacobianInverse, AnalysisResults results, Map<String, Integer> busNumbers)
    {
		super(subnetwork, basePower, baseVoltage);
		
		slackBus = subnetwork.getSlackBus();
		
		subnetwork.setName("Local ("+controlledBus+")");
		
		this.controlledBus = subnetwork.getBus(controlledBus);
        this.globalBusNumbers = busNumbers;
        
        update(jacobianInverse, results);
        
//        int slackIndex = results.getSlackIndex();
//        findSlackBus(slackIndex, results);

		resetControllableUnits();
    }
	
	
	//// FIXME only controlling pilot bus
	
	@Override
	protected void resetControllableUnits()
	{
		if(controlledBus == null) // This will happen when called via super()
			return;
		
		controllableUnits = new ArrayList<Unit>();
		controllableUnits.addAll(controlledBus.getControllableDemands());
		controllableUnits.addAll(controlledBus.getRenewableSources());
	}

	@Override
	protected void applyControls(CostToGo ctg, int time)
	{
		// Apply for all busses from previous controls:
		if(previousCostToGoSchedule != null)
		{
			CostToGo prevCtg = previousCostToGoSchedule.get(time);
			for (Unit unit : prevCtg.controls.keySet())
			{
	            if(unit instanceof RenewableGenerator)
	            	((RenewableGenerator)unit).setPowerOutput(prevCtg.controls.get(unit));
	            else if(unit instanceof ControllableDemand)
	            	((ControllableDemand)unit).setChargeRate(prevCtg.controls.get(unit));
			}
		}
		
		// Apply for controllable units (this is only controlled bus's units):
		if(time == 0)
			Util.nullop();
		super.applyControls(ctg, time);
	}
	
	////
	

	public void update(
			RealMatrix jacobianInverse, 
			AnalysisResults globalResults) 
	{
		this.globalResults = globalResults;
        this.jacobianInverse = jacobianInverse;
        
		int size = jacobianInverse.getColumnDimension();
		
		globalVoltages = new ArrayRealVector(size);
        globalPowers = new ArrayRealVector(size);
        fillInVoltagesAndPowers(globalResults, size, globalVoltages, globalPowers);
	}

	void fillInVoltagesAndPowers(
			AnalysisResults results, 
			int size, 
			RealVector voltages, 
			RealVector powers)
	{
		int globalSlackIndex = results.getSlackIndex();
        for (Bus bus : busses) 
		{
			String busName = bus.getName();
			Integer busIndex = globalBusNumbers.get(busName);
			
        	if(busIndex == globalSlackIndex)
				continue;
			
			if(busIndex > globalSlackIndex)
				--busIndex;
			
			Complex busVoltage = results.getBusVoltage(busName);
			Complex busPower = results.getBusPower(busName);
			voltages.setEntry(busIndex, busVoltage.getArgument());
			voltages.setEntry(busIndex+size/2, busVoltage.abs());
			powers.setEntry(busIndex, busPower.getReal());
			powers.setEntry(busIndex+size/2, busPower.getImaginary());
		}
	}
    
    
    //// Local Optimisation ////
    
    public Map<Unit, Complex> getControls() 
    {
		return getCostToGoSchedule().get(0).controls;
	}
    
    /**
     * Applies all controls to controllable units in the local controller's grid.
     * @see #applyAllControls(Grid)
     */
    public void applyAllControls(int time)
    {
        applyAllControls(this.grid, time);
    }
    
    public void applyAllControls(Grid grid, int time)
    {
        for (Unit unit : controllableUnits)
        {
            applyControls(grid, unit, time);
        }
    }

    /**
     * Applies the controls for the given unit in the given grid.
     * The given unit may not belong to the given grid, but a unit
     * of the same type with the same name must exist in the given grid.
     * @param grid
     * @param unit
     * @return
     */
	private void applyControls(Grid grid, Unit unit, int time)
	{
		List<CostToGo> costToGoSchedule = getCostToGoSchedule();
		if(time >= costToGoSchedule.size())
			return;
		
		Map<Unit, Complex> controls = costToGoSchedule.get(time).controls;
		if(controls == null)
			return;
		Complex control = controls.get(unit);
		if(control == null)
			return;
		
		// Apply control to unit from given grid:
		unit = grid.getUnit(unit.getName());
		if(unit instanceof Source)
		    ((Source)unit).setPowerOutput(control);
		else if (unit instanceof ControllableDemand)
		    ((ControllableDemand)unit).setChargeRate(control);
	}

	/**
	 * Applies the controls to the units in the controlled bus, but applies them to 
	 * the equivalent units in the given grid.
	 * @param grid
	 */
	public void applyControls(Grid grid)
	{
		applyControls(grid, 0);
	}
	
	public void applyControls(Grid grid, int time)
	{
		for (Unit unit : controlledBus.getChildren()) 
		{
			applyControls(grid, grid.getUnit(unit.getName()), time);
		}
	}

    @Override
    public double cost(final int timeIndex)
    {
    	final ScenarioLogger scenarioLogger = ScenarioLogger.instance();
    	
    	if(STATS) timer.start(TIME_COST);

		LazyCallback<Double> callback = new LazyCallback<Double>() {
			@Override
			public Double value() {

				Map<String, Integer> busNumbers = new HashMap<String, Integer>();
				int n = 0;
				for (Bus bus : busses) 
				{
					busNumbers.put(bus.getName(), n);
					++n;
				}

				// Calculate change in power:
				int size = jacobianInverse.getColumnDimension()/2;
				RealVector deltaS = deltaS();
				
				// Calculate change in voltage:
				RealVector deltaV = deltaV(timeIndex, deltaS);

				// Apply voltage changes:
				RealVector newV = globalVoltages.add(deltaV);
				
				// Convert to complex voltages:
				Complex voltages[] = new Complex[busses.size()];
				int globalSlackIndex = globalResults.getSlackIndex();
				double highestVoltage = -1;

				double breachCost = 0;

				for (Bus bus : busses)
				{
					String busName = bus.getName();
					Integer busIndex_arg = globalBusNumbers.get(busName);
					
		        	if(busIndex_arg == globalSlackIndex)
					{
		        		throw new RuntimeException("Subnetwork '"+grid.getName()+"' contained global slack bus");
					}
					
					if(busIndex_arg > globalSlackIndex)
						--busIndex_arg;
					
					int busIndex_abs = busIndex_arg + size;
					
					double abs = newV.getEntry(busIndex_abs);
					double arg = newV.getEntry(busIndex_arg);
					
					int localIndex = busNumbers.get(busName);
					voltages[localIndex] = ComplexUtils.polar2Complex(abs, arg);

					if(finalADPCostToGo && bus.equals(controlledBus))
						FeederVoltageLogger.instance().logLocalVoltage(timeIndex, controlledBus, voltages[localIndex], CSVLogger.LOG_LEVEL_DEBUG);
					
					// Check voltage constraints [PHD-17]:
					double max, min;
					if(voltageErrorDelta == null)
						viableSolution = false;
					if(!viableSolution || timeIndex >= voltageErrorDelta.length)
					{
						max = maxVoltage;
						min = minVoltage;
					}
					else
					{
						double delta = Math.abs(voltageErrorDelta[timeIndex].getEntry(busIndex_abs));
						max = maxVoltage - delta;
						min = minVoltage + delta;
					}
					
					if(abs < min || abs > max) // Invalid state, invalid controls
						breachCost = 1;
//						throw new InvalidStateException();
//					if(abs < min)
//						breachCost += min-abs;
//					else if(abs > max)
//						breachCost += abs-max;
					
					highestVoltage = Math.max(highestVoltage, abs);
					
					// Check target vs actual controls:
//					sd
				}
				
				ControlLogger.instance().setHighestVoltage(highestVoltage, timeIndex, CSVLogger.LOG_LEVEL_DEBUG);

				// Log:
				scenarioLogger.logDeltaS(deltaS, busses, globalResults, CSVLogger.LOG_LEVEL_DEBUG);
				scenarioLogger.logDeltaV(deltaV, busses, globalResults, CSVLogger.LOG_LEVEL_DEBUG);
				scenarioLogger.logNewState(voltages, busses, busNumbers, CSVLogger.LOG_LEVEL_DEBUG);
				
				// Find contribution to global slack bus power:
				// Y(slack)[dV_A/dS_B] Delta S_B
				// Where A is the set of busses adjacent to the global slack bus,
				// and B is the set of busses in the local sub-network.
				RealVector deltaV_A = slackAdjacentJacobianInverse.operate(deltaS);
				Complex aComplexDeltaV_A[] = new Complex[size+1]; // +1: Delta V_A must include the slack bus before being multiplied by Y
				for (int i = 0; i < size; i++)
				{
					double arg = deltaV_A.getEntry(i);
					double abs = deltaV_A.getEntry(i + size);
					if(abs < 0)
					{
						arg += Math.PI; // +180 deg
						abs = -abs;
					}

					int i2 = i >= globalSlackIndex ? i+1 : i;
					aComplexDeltaV_A[i2] = ComplexUtils.polar2Complex(abs, arg);
				}
				aComplexDeltaV_A[globalSlackIndex] = Complex.ZERO;
				FieldVector<Complex> complexDeltaV_A = new ArrayFieldVector<Complex>(aComplexDeltaV_A);
				Complex slackPower = slackAdmittances.dotProduct(complexDeltaV_A);
				
				double cost = slackPower.abs()*Math.signum(slackPower.getReal());
				return cost + breachCost;
			}
		};
		
		if(CACHE)
		{
			Object[] state = new Object[busses.size()];
			int i = 0;
			for (Bus bus : busses) 
			{
				state[i++] = bus.getNetPower(false);
			}
			
			Double slackPower = costCache.lazyGet(callback, state);

			if(STATS) timer.stop(TIME_COST);

			return slackPower;
		}
		else
		{
			if(STATS) timer.stop(TIME_COST);
			
			Double value = callback.value();
			if(timeIndex == 0)
				Util.nullop();
			return value;
		}
	}

	public RealVector deltaS()
	{
		int size = jacobianInverse.getColumnDimension()/2;
		int globalSlackIndex = globalResults.getSlackIndex();
		
		RealVector deltaS = new ArrayRealVector(size*2);
		for (Bus bus : busses) 
		{
			String busName = bus.getName();
			Integer busIndex = globalBusNumbers.get(busName);
			
        	if(busIndex == globalSlackIndex)
				continue;
			
			if(busIndex > globalSlackIndex)
				--busIndex;
			
			if(bus.getSlackVoltage().equals(Complex.ZERO))
			{
				Complex busPower = bus.netPower().divide(basePower);
				Complex previousPower = new Complex(
						globalPowers.getEntry(busIndex),
						globalPowers.getEntry(busIndex+size));
				Complex delta = busPower.subtract(previousPower);
				deltaS.setEntry(busIndex, delta.getReal());           // P
				deltaS.setEntry(busIndex+size, delta.getImaginary()); // Q
			}
		}
		return deltaS;
	}

	public RealVector deltaV(int timeIndex, RealVector deltaS)
	{
if(timeIndex == 1 && centralIteration >= 8 && controlledBus.getName().equals("Bus 6") && this.finalADPCostToGo)
	Util.nullop();
		RealVector deltaV = jacobianInverse.operate(deltaS);
		
		if(voltageError != null && timeIndex < voltageError.length && voltageError[timeIndex] != null)
			deltaV = deltaV.add(voltageError[timeIndex]); // attempt to correct approximation

		return deltaV;
	}
	
    
    //// Accessors ////

	public Bus getSlackBus() {
		return slackBus;
	}

	public Bus getControlledBus() {
		return controlledBus;
	}
	
	public AnalysisResults getGlobalResults()
    {
        return globalResults;
    }

	public RealVector getGlobalVoltages()
	{
		return globalVoltages;
	}

	public RealVector getGlobalPowers()
	{
		return globalPowers;
	}

	public void setVoltageError(RealVector[] error)
	{
		RealVector[] oldError;
		if(voltageError != null)
		{
			oldError = new RealVector[voltageError.length];
			for (int t = 0; t < oldError.length; t++)
			{
				oldError[t] = new ArrayRealVector(voltageError[t]);
			}
		}
		else
		{
			oldError = null;
		}
		
		// Update voltage errors:
		if(ERRORS_WITH_MEMORY)
		{
			if (this.voltageError != null)
			{
				for (int t = 0; t < error.length; t++)
				{
					this.voltageError[t] = this.voltageError[t].mapMultiply(1-BETA).add(error[t].mapMultiply(BETA));
				}
			} 
			else
			{
				this.voltageError = error;
			}
		} 
		else
		{
			this.voltageError = error;
		}

		// Update voltageErrorDelta (PHD-17):
		if(oldError != null)
		{
			if(voltageErrorDelta == null)
				voltageErrorDelta = new RealVector[error.length];
			
			for (int t = 0; t < error.length; t++)
			{
				RealVector delta = voltageError[t].subtract(oldError[t]);
				
				if(ERROR_DELTAS_WITH_MEMORY)
				{
					if(voltageErrorDelta[t] == null)
					{
						voltageErrorDelta[t] = delta;
					}
					else
					{
						voltageErrorDelta[t] = voltageErrorDelta[t].mapMultiply(1-BETA).add(delta.mapMultiply(BETA));
					}
				}
				else
				{
					voltageErrorDelta[t] = delta;
				}
				
				ErrorLogger.instance().logVoltageError(centralIteration, grid.getName(), t, voltageError[t], voltageErrorDelta[t], CSVLogger.LOG_LEVEL_DEBUG);
			}
		}
	}

	public RealMatrix getJacobianInverse()
	{
		return this.jacobianInverse;
	}

//	public void setGlobalSlackBusVoltage(Complex slackVoltage)
//	{
//		this.globalSlackVoltage = slackVoltage;
//	}

	public void setSlackAdmittances(FieldVector<Complex> slackAdmittances)
	{
		this.slackAdmittances = slackAdmittances;
	}

	public void setSlackAdjacentJacobianInverse(RealMatrix slackAdjacentJacobianInverse)
	{
		this.slackAdjacentJacobianInverse = slackAdjacentJacobianInverse;
	}
}
