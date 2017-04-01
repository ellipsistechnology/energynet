package ellipsis.energy.grid;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedHashSet;

import org.apache.commons.math3.complex.Complex;

import ellipsis.energy.grid.res.RenewableGenerator;
import ellipsis.energy.smartgrid.ControllableDemand;

/**
 * Everything is a child of a bus.
 * @author bmillar
 *
 */
public class Bus extends Unit
{
    private Collection<Unit> units = new LinkedHashSet<Unit>();
    private Grid grid;
    
    public Bus(Grid grid)
    {
        this.grid = grid;
    }
    
    @Override
    public void setName(String name)
    {
        if(grid != null)
            grid.nameWillChange(this, name);
        
        super.setName(name);
    }

    /**
     * If the bus does not already contain the given unit then:
     * <ul>
     *     <li>If the unit is a child it will have {@link Child#setParent(Bus)} called,</li>
     *     <li>The given unit will be added to the bus,</li>
     *     <li>And if the bus' grid is not null, {@link Grid#add(Unit)} will be called with the given unit.</li>
     * </ul>
     * @param child
     */
    public void addChild(Unit child)
    {
        // TODO make sure that only one slack source is added
        if(!units.contains(child))
        {
            if(child instanceof Child)
                ((Child)child).setParent(this);
            
            units.add(child);
            
            if(grid != null)
                grid.add(child);
        }
    }

    /**
     * <ul>
     *     <li>If the unit is a child it will have {@link Child#setParent(Bus)} called setting the parent to null,</li>
     *     <li>The given unit will be removed from the bus,</li>
     *     <li>{@link Unit#wasRemovedFromBus(Bus)} will be called,</li>
     *     <li>And if the bus' grid is not null, {@link Grid#remove(Unit)} will be called with the given unit.</li>
     * </ul>
     * @param child
     */
    public void removeChild(String name)
    {
        for (Iterator<Unit> iterator = units.iterator(); iterator.hasNext();)
        {
            Unit unit = iterator.next();
            if(unit.getName().equals(name))
            {
                iterator.remove();
                
                if(unit instanceof Child)
                    ((Child)unit).setParent(null);
                
                unit.wasRemovedFromBus(this);
                
                if(grid != null)
                    grid.remove(unit);
                
                return;
            }
        }
    }
    
    public void removeChild(Unit unit)
    {
        units.remove(unit);
        
        if(unit instanceof Child)
            ((Child)unit).setParent(null);
        unit.wasRemovedFromBus(this);
        
        if(grid != null)
            grid.remove(unit);
    }

	public void removeAllChildren(boolean includeLines) 
	{
		for (Unit unit : new ArrayList<Unit>(units)) 
		{
			if(!includeLines && unit instanceof Line)
				continue;
			
			units.remove(unit);
			
	        if(unit instanceof Child)
	            ((Child)unit).setParent(null);
	        unit.wasRemovedFromBus(this);
            
            if(grid != null)
                grid.remove(unit);
		}
	}
    
    public Collection<Unit> getChildren()
    {
        return units;
    }

    public Complex getLoadPower()
    {
        Complex loadValue = Complex.ZERO;
        for (Load load : getLoads())
        {
            loadValue = loadValue.add(load.getLoad());
        }
        return loadValue;
    }

    public Complex getDistributedSourcePower()
    {
        Complex sourceValue = Complex.ZERO;
        for (Unit unit : units)
        {
            if(unit instanceof DistributedSource)
                sourceValue = sourceValue.add(((Source)unit).getPowerOutput());
        }
        return sourceValue;
    }

    public Collection<Load> getLoads()
    {
        ArrayList<Load> loads = new ArrayList<Load>();
        for (Unit unit : units)
        {
            if(unit instanceof Load)
                loads.add((Load)unit);
        }
        return loads;
    }

    public Collection<Load> getNonCapacitorLoads()
    {
        ArrayList<Load> loads = new ArrayList<Load>();
        for (Unit unit : units)
        {
            if(unit instanceof Load && !(unit instanceof Capacitor))
                loads.add((Load)unit);
        }
        return loads;
    }

    public Complex getGeneratedPower()
    {
        Complex sourceValue = Complex.ZERO;
        for (Source source : getSources())
        {
            sourceValue = sourceValue.add(source.getPowerOutput());
        }
        return sourceValue;
    }

    public double getAvailablePower()
    {
        double available = 0;
        for (Source source : getSources())
        {
            available += source.getAvailablePower();
        }
        return available;
    }
    
    public Complex getNetPower(boolean includeSlack)
    {
        if(includeSlack)
            return getGeneratedPower().subtract(getLoadPower());
        else
            return getDistributedSourcePower().subtract(getLoadPower());
    }
    
    @Override
    public Complex netPower()
    {
    	Complex power = Complex.ZERO;
    	for (Unit unit : units) 
    	{
			power = power.add(unit.netPower());
		}
    	return power;
    }

    public Collection<Source> getSources()
    {
        ArrayList<Source> sources = new ArrayList<Source>();
        for (Unit unit : units)
        {
            if(unit instanceof Source)
                sources.add((Source)unit);
        }
        return sources;
    }

	public Collection<DistributedSource> getRenewableSources()
	{
		ArrayList<DistributedSource> sources = new ArrayList<DistributedSource>();
        for (Unit unit : units)
        {
            if(unit instanceof RenewableGenerator)
                sources.add((DistributedSource)unit);
        }
        return sources;
	}

	public Collection<ControllableDemand> getControllableDemands() 
	{
        ArrayList<ControllableDemand> demands = new ArrayList<ControllableDemand>();
        for (Unit unit : units)
        {
            if(unit instanceof ControllableDemand)
                demands.add((ControllableDemand)unit);
        }
        return demands;
	}

    public Collection<Line> getLines()
    {
        ArrayList<Line> lines = new ArrayList<Line>();
        for (Unit unit : units)
        {
            if(unit instanceof Line)
                lines.add((Line)unit);
        }
        
        return lines;
    }

    public Complex getSlackVoltage()
    {
        for (Unit unit : units)
        {
            if(unit instanceof SlackSource)
                return ((SlackSource)unit).getVoltage();
        }
        return Complex.ZERO;
    }

    public void setSlackVoltage(Complex v)
    {
        for (Unit unit : units)
        {
            if(unit instanceof SlackSource)
                ((SlackSource)unit).setVoltage(v);
        }
    }

    public Collection<ShuntAdmittance> getShuntAdmittances()
    {
        ArrayList<ShuntAdmittance> shuntAdmittances = new ArrayList<ShuntAdmittance>();
        for (Unit unit : units)
        {
            if(unit instanceof ShuntAdmittance)
                shuntAdmittances.add((ShuntAdmittance)unit);
        }
        return shuntAdmittances;
    }
    
    @Override
    public String toString()
    {
        return getName();
    }

    /**
     * Finds the {@link Line} that is a child of this bus and the bus
     * with the given name. It does not matter which bus is the To and
     * which bus is the From.
     * @param busName The name of the bus this bus is connected to via the line.
     * @return The line between this bus and the bus with the given name.
     */
    public Line lineTo(String busName)
    {
        for (Unit unit : units)
        {
            if(unit instanceof Line)
            {
                Line line = (Line)unit;
                if(this.equals(line.getFromBus()) && line.getToBus().getName().equals(busName) ||
                   this.equals(line.getToBus()) && line.getFromBus().getName().equals(busName))
                {
                    return line;
                }
            }
        }
        
        return null;
    }
    
    protected static class BusPlaceHolder extends Bus {
		public BusPlaceHolder(Grid grid) {
			super(grid);
		}
	}

    @Override
	public Unit copy(Grid grid)
	{
		Bus newBus = new Bus(grid);
		newBus.setName(getName());
		for (Unit unit: units)
		{
			Unit copy = unit.copy(grid);
			if(copy != null)
				newBus.addChild(copy);
		}
		return newBus;
	}
}
