package ellipsis.energy.smartgrid;

import ellipsis.energy.grid.Grid;
import ellipsis.energy.grid.Unit;

/**
 * A {@link ControllableDemand} that can only charge (i.e. it is always a load and never a source).
 * The load will have a target demand in kWh that it will wish to consume periodically (ref. {@link #setDemandPeriod(double)},
 * {@link #setPeriodicDemand(double)} and {@link #passTime(double)}.
 * @author bmillar
 *
 */
public class ControllableLoad extends ControllableDemand
{
    private double demandPeriod;
    private double timePassed = 0;

    public ControllableLoad()
    {
        setMaxDischargeRate(0);
    }
    
    /**
     * Alias for {@link ControllableDemand#setMaxCapacity(double)}.
     * @param demand Demand in kWh.
     */
    public void setPeriodicDemand(double demand)
    {
        setMaxCapacity(demand);
    }
    
    /**
     * Sets the number of hours after which the capacity should be reset such that
     * over this period the kWh consumed should total the periodic demand.
     * @param hours
     */
    public void setDemandPeriod(double hours)
    {
        this.demandPeriod = hours;
    }
    
    /**
     * After the demand period has passed this will reset the capacity
     * such that it can charge again.
     * @see ControllableDemand#passTime(double)
     */
    @Override
    public void passTime(double hours)
    {
        super.passTime(hours);
        timePassed += hours;
        
        if(timePassed > demandPeriod)
        {
            timePassed -= demandPeriod;
            setCapacity(0);
        }
    }
    
    @Override
    public Unit copy(Grid grid) 
    {
    	ControllableLoad copy = new ControllableLoad();
    	copyInto(copy);
    	copy.setDemandPeriod(demandPeriod);
    	return copy;
    }
}
