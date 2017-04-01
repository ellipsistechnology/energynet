package ellipsis.energy.smartgrid;

import org.apache.commons.math3.complex.Complex;

import ellipsis.energy.grid.Bus;
import ellipsis.energy.grid.Child;
import ellipsis.energy.grid.DistributedSource;
import ellipsis.energy.grid.Grid;
import ellipsis.energy.grid.Load;
import ellipsis.energy.grid.Source;
import ellipsis.energy.grid.Unit;

public class ControllableDemand extends Unit implements Child
{
    private Bus parent;
    
    // Constraints:
    private double minCapacity;
    private double maxCapacity;
    private double minChargeRate;
    private double maxChargeRate;
    private double minDischargeRate;
    private double maxDischargeRate;
    
    // State:
    private double capacity;
    private Complex chargeRate;
    private Complex targetChargeRate;
    
    // Grid units:
    public static class CDLoad extends Load
    {
		@Override
		public Unit copy(Grid grid) 
		{
			return null;
		}
	}
    public static class CDSource extends DistributedSource
    {
		@Override
		public Unit copy(Grid grid) 
		{
			return null;
		}
	}
    private Load load = new CDLoad();
    private Source source = new CDSource();
    
    
    //// CONSTRUCTOR/FACTORY ////
    
    public static class Factory
    {
        private ControllableDemand cd = new ControllableDemand();
        
        private double initialCapacity;
        
        private Factory() {}
        
        public Factory withName(String name)
        {
        	cd.setName(name);
            return this;
        }
        
        public Factory withInitialCapacity(double capacity)
        {
        	initialCapacity = capacity;
            return this;
        }
        
        public Factory withMinCapacity(double minCapacity)
        {
            cd.minCapacity = minCapacity;
            return this;
        }
        
        public Factory withMaxCapacity(double maxCapacity)
        {
            cd.maxCapacity = maxCapacity;
            return this;
        }
        
        public Factory withMinChargeRate(double minChargeRate)
        {
            cd.minChargeRate = minChargeRate;
            return this;
        }
        
        public Factory withMaxChargeRate(double maxChargeRate)
        {
            cd.maxChargeRate = maxChargeRate;
            return this;
        }
        
        public Factory withMinDischargeRate(double minDischargeRate)
        {
            cd.minDischargeRate = minDischargeRate;
            return this;
        }
        
        public Factory withMaxDischargeRate(double maxDischargeRate)
        {
            cd.maxDischargeRate = maxDischargeRate;
            return this;
        }
        
        public ControllableDemand done()
        {
        	cd.setCapacity(initialCapacity);
            return cd;
        }
    }
    
    public static Factory create()
    {
        return new Factory();
    }
    
    
    //// ACCESSORS ////
    
    public double getMinCapacity()
    {
        return minCapacity;
    }
    public void setMinCapacity(double minCapacity)
    {
        this.minCapacity = minCapacity;
    }
    public double getMaxCapacity()
    {
        return maxCapacity;
    }
    public void setMaxCapacity(double maxCapacity)
    {
        this.maxCapacity = maxCapacity;
    }
    public double getMinChargeRate()
    {
        return minChargeRate;
    }
    public void setMinChargeRate(double minChargeRate)
    {
        this.minChargeRate = minChargeRate;
    }
    public double getMaxChargeRate()
    {
        return maxChargeRate;
    }
    public void setMaxChargeRate(double maxChargeRate)
    {
        this.maxChargeRate = maxChargeRate;
    }
    public double getMinDischargeRate()
    {
        return minDischargeRate;
    }
    public void setMinDischargeRate(double minDischargeRate)
    {
        this.minDischargeRate = minDischargeRate;
        this.source.setPmin(minDischargeRate);
    }
    public double getMaxDischargeRate()
    {
        return maxDischargeRate;
    }
    public void setMaxDischargeRate(double maxDischargeRate)
    {
        this.maxDischargeRate = maxDischargeRate;
        this.source.setPmax(maxDischargeRate);
    }
    public double getCapacity()
    {
        return capacity;
    }
	public Complex getTargetChargeRate()
	{
		return targetChargeRate;
	}
    
    /**
     * Capacity in kWh. 
     * Must be between relevant min/max.
     * @param capacity
     */
    public void setCapacity(double capacity)
    {
        // Check capacity is within bounds:
        if(capacity < minCapacity)
            throw new RuntimeException("capacity ("+capacity+") must be greater than minCapacity ("+minCapacity+")");
        if(capacity > maxCapacity)
            throw new RuntimeException("capacity ("+capacity+") must be less than maxCapacity ("+maxCapacity+")");
        
        this.capacity = capacity;
    }
    
    public Complex getChargeRate()
    {
        return chargeRate;
    }
    
    /**
     * Charge rate in kW. For discharge, set to a negative number.
     * Must be between relevant min/max.
     * Source and load values will be updated.
     * @param chargeRate
     */
    public void setChargeRate(Complex chargeRate)
    {
    	this.targetChargeRate = chargeRate;
    	
        // Check rate is within bounds:
        if(chargeRate.getReal() >= 0) // charging
        {
            if(chargeRate.abs() < minChargeRate)
                throw new RuntimeException("chargeRate ("+chargeRate+") must be greater than minChargeRate ("+minChargeRate+")");
            if(chargeRate.abs() > maxChargeRate)
                throw new RuntimeException("chargeRate ("+chargeRate+") must be less than maxChargeRate ("+maxChargeRate+")");
        }
        else // discharging
        {
            if(chargeRate.abs() < minDischargeRate)
                throw new RuntimeException("chargeRate ("+chargeRate+") must be greater than minDischargeRate ("+minDischargeRate+")");
            if(chargeRate.abs() > maxDischargeRate)
                throw new RuntimeException("chargeRate ("+chargeRate+") must be less than maxDischargeRate ("+maxDischargeRate+")");
        }
        
        // Update load and source:
        if(chargeRate.getReal() >= 0) // charging
        {
        	if(getCapacity() >= getMaxCapacity())
        		chargeRate = Complex.ZERO;
        	load.setLoad(chargeRate);
            source.setPowerOutput(Complex.ZERO);
        }
        else // discharging
        {
            load.setLoad(Complex.ZERO);
            if(getCapacity() <= 0)
            	chargeRate = Complex.ZERO;
            source.setPowerOutput(chargeRate.negate());
        }
        
        this.chargeRate = chargeRate;
    }
    
	public void setChargeRate(double active, double reactive)
	{
		setChargeRate(new Complex(active, reactive));
	}
    
	@Override
	public Complex netPower() 
	{
		return Complex.ZERO;
//		if(getChargeRate().getReal() >= 0) // charging
//			return load.netPower();
//		else // discharging
//			return source.netPower();
	}

    
    //// CHILD/UNIT IMLEMENTATION ////

    @Override
    public Bus getParent()
    {
        return parent;
    }

    @Override
    public void setParent(Bus bus)
    {
        this.parent = bus;
        
        // Add load and source to parent:
        if(bus != null)
        {
	        bus.addChild(load);
	        bus.addChild(source);
        }
    }
    
    @Override
    public void wasRemovedFromBus(Bus parent)
    {
        super.wasRemovedFromBus(parent);
        
        // Remove load and source:
        parent.removeChild(load);
        parent.removeChild(source);
    }
    
    @Override
    public void setName(String name)
    {
        super.setName(name);
        load.setName(name+":load");
        source.setName(name+":source");
    }
    
    @Override
    public Unit copy(Grid grid) 
    {
    	ControllableDemand copy = new ControllableDemand();
    	copyInto(copy);
    	
    	return copy;
    }

	private void copyInto(ControllableDemand copy) 
	{
		super.copyInto(copy);
    	copy.setMinCapacity(getMinCapacity());
    	copy.setMaxCapacity(getMaxCapacity());
    	copy.setMinChargeRate(getMinChargeRate());
    	copy.setMaxChargeRate(getMaxChargeRate());
    	copy.setMinDischargeRate(getMinDischargeRate());
    	copy.setMaxDischargeRate(getMaxDischargeRate());
    	copy.setCapacity(getCapacity());
    	copy.setChargeRate(getChargeRate());
	}
    
    
    //// SIMULATION ////
    
    /**
     * Updates the capacity based on the current charge/discharge rate,
     * within the min/max capacity values.
     * If min capacity has been reached (fully discharged) then the charge rate
     * will be set to {@link #minDischargeRate}.
     * If max capacity has been reached (fully charged) then the charge rate
     * will be set to {@link #minChargeRate}.
     * Source and load values will be updated.
     * @see #setChargeRate(double)
     * @param hours Hours that have passed since the last call to {@link #passTime(double)}.
     */
    public void passTime(double hours)
    {
		double chargeDelta = hours*chargeRate.abs();
        if(chargeRate.getReal() < 0)
        	chargeDelta = -chargeDelta;
        double newCapacity = getCapacity()+chargeDelta;
        
        Complex target = this.targetChargeRate;
        if(newCapacity >= maxCapacity)
        {
            newCapacity = maxCapacity;
            setChargeRate(new Complex(minChargeRate));
        }
        else if(newCapacity <= minCapacity)
        {
            newCapacity = minCapacity;
            setChargeRate(new Complex(minDischargeRate));
        }
        this.targetChargeRate = target;
        
        setCapacity(newCapacity);
    }
}
