package ellipsis.energy.grid;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Stack;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexUtils;

import ellipsis.energy.smartgrid.ControllableDemand;
import ellipsis.energy.smartgrid.ControllableDemand.CDLoad;
import ellipsis.energy.smartgrid.ControllableDemand.CDSource;

public class Grid
{
    Map<String, Unit> units = new LinkedHashMap<String, Unit>();
    private String name;
    private double basePower;
    private double baseVoltage;
    
    public Grid()
    {
    }
    
    public int getUnitCount()
    {
    	return units.size();
    }
    
    /**
     * Adds the given unit to the grid.
     * If the given unit is a {@link Bus} then each child will be added to the grid
     * through a recursive call to {@link #add(Unit)} if that unit is not already a
     * unit if the grid.
     * @param unit
     */
    public void add(Unit unit)
    {
        if(unit.getName() == null)
            throw new NullPointerException("Unit name cannot be null: "+unit.getClass().getName());
        units.put(unit.getName(), unit);
        
        if(unit instanceof Bus)
        {
            Bus bus = (Bus)unit;
            for (Unit child : bus.getChildren())
            {
                if(!units.containsKey(child))
                    add(child);
            }
        }
    }
    
    /**
     * Removes the given unit from the grid.
     * If the given unit is a {@link Bus} then each child will be removed from the grid
     * through a recursive call to {@link #remove(Unit)}.
     * @param unit
     */
    public void remove(Unit unit)
    {
        units.remove(unit.getName());
        
        if(unit instanceof Bus)
        {
            Bus bus = (Bus)unit;
            for (Unit child : bus.getChildren())
            {
                remove(child);
            }
        }
    }

    /**
     * Notifies the grid that the bus name is about to change.
     * This must be called BEFORE the name is changed.
     * This will not change the name of the bus.
     * @param bus
     * @param newName
     */
    public void nameWillChange(Bus bus, String newName)
    {
        units.remove(bus.getName());
        units.put(newName, bus);
    }
    
    public Iterable<Unit> unitIterator()
    {
        return units.values();
    }
    
    public Unit getUnit(String name)
    {
        return units.get(name);
    }

    @SuppressWarnings("unchecked")
    public <T extends Unit> T getUnit(Class<T> class1, String name)
    {
        return (T)getUnit(name);
    }

    public Transformer getTransformer(String name)
    {
        return (Transformer)units.get(name);
    }

    public Source getSource(String name)
    {
        return (Source)units.get(name);
    }

    public SlackSource getSlackSource(String name)
    {
        return (SlackSource)units.get(name);
    }
    
    public DistributedSource getDistributedSource(String name)
    {
        return (DistributedSource)units.get(name);
    }

    public Bus getBus(String name)
    {
        return (Bus)units.get(name);
    }

    public Bus getOrCreateBus(String name)
    {
        Bus bus = (Bus)units.get(name);
        if(bus == null)
        {
        	bus = new Bus(this);
        	bus.setName(name);
        	add(bus);
        }
		return bus;
    }

    public Line getLine(String name)
    {
        return (Line)units.get(name);
    }
    
    public Load getLoad(String name)
    {
        return (Load)units.get(name);
    }
    
    public ControllableDemand getControllableDemand(String name)
    {
        return (ControllableDemand)units.get(name);
    }

    public Capacitor getCapacitor(String name)
    {
        return (Capacitor)units.get(name);
    }

    public LowVoltageGroup getLowVoltageGroup(String name)
    {
        return (LowVoltageGroup)units.get(name);
    }
    
    public Collection<Source> getSources()
    {
        HashSet<Source> sources = new LinkedHashSet<Source>();
        for (Unit unit : units.values())
        {
            if(unit instanceof Source)
                sources.add((Source)unit);
        }
        return sources;
    }
    
    public Collection<DistributedSource> getDistributedSources(boolean includeCDSources)
    {
        HashSet<DistributedSource> sources = new LinkedHashSet<>();
        for (Unit unit : units.values())
        {
            if(
            		unit instanceof DistributedSource && 
            		(includeCDSources || !(unit instanceof CDSource)))
            {
                sources.add((DistributedSource)unit);
            }
        }
        return sources;
    }
    
    public Collection<Load> getLoads()
    {
    	return getLoads(true);
    }

    public Collection<Load> getLoads(boolean includeCDLoads)
    {
        HashSet<Load> loads = new LinkedHashSet<Load>();
        for (Unit unit : units.values())
        {
            if(
            		unit instanceof Load &&
            		(includeCDLoads || !(unit instanceof CDLoad)))
            {
                loads.add((Load)unit);
            }
        }
        return loads;
    }

    public Collection<Line> getLines()
    {
        HashSet<Line> lines = new LinkedHashSet<Line>();
        for (Unit unit : units.values())
        {
            if(unit instanceof Line)
                lines.add((Line)unit);
        }
        return lines;
    }

    public Collection<Bus> getBusses()
    {
        HashSet<Bus> busses = new LinkedHashSet<Bus>();
        for (Unit unit : units.values())
        {
            if(unit instanceof Bus)
                busses.add((Bus)unit);
        }
        return busses;
    }

    @SuppressWarnings("unchecked")
    public <T extends Unit> Collection<T> get(Class<T> class1)
    {
        HashSet<T> ts = new LinkedHashSet<T>();
        for (Unit unit : units.values())
        {
            if(class1.isAssignableFrom(unit.getClass()))
                ts.add((T)unit);
        }
        return ts;
    }
    
    /**
     * Constructs a new Grid based on the given busses.
     * The constructed subnetwork will include all busses given,
     * and also any busses required to ensure that the given
     * busses are connected. The intermediate busses are chosen
     * according to the shortest path between the given busses.
     * The returned grid's units are copies of the given grid's.
     * @param busses
     * @return
     */
    public Grid subnet(Collection<Bus> busses)
    {
        if(busses.isEmpty())
            return new Grid();
            
        HashSet<Bus> unorphaned = new HashSet<Bus>();
        HashSet<Bus> orphaned = new HashSet<Bus>(busses);
        
        // Move from orphaned to unorphaned:
        search(orphaned.iterator().next(), orphaned, unorphaned);
        
        // For each orphaned bus find the shortest path:
        for (Iterator<Bus> orphanIterator = orphaned.iterator(); orphanIterator.hasNext();)
        {
            Bus orphan = (Bus) orphanIterator.next();
            
            // Assign to every node a tentative distance value: set it to zero for our 
            // initial node and to infinity for all other nodes.
            HashMap<Bus, Double> distances = new HashMap<Bus, Double>();
            for (Bus bus : getBusses())
            {
                distances.put(bus, Double.MAX_VALUE);
            }
            distances.put(orphan, Double.valueOf(0));
            
            // Mark all nodes unvisited. Set the initial node as current. Create a set 
            // of the unvisited nodes called the unvisited set consisting of all the 
            // nodes except the initial node.
            HashSet<Bus> unvisited = new HashSet<Bus>(getBusses());
            unvisited.remove(orphan);
            HashMap<Bus, Bus> previous = new HashMap<Bus,Bus>();
            Bus current = orphan;
            
            while(!unvisited.isEmpty()) // NOTE: Loop should break below since we expect to find a valid path.
            {
                // If the destination node has been marked visited (when planning a route between 
                // two specific nodes) or if the smallest tentative distance among the nodes in the 
                // unvisited set is infinity (when planning a complete traversal), then stop. 
                // The algorithm has finished.
                if(unorphaned.contains(current))
                {
                    // Find all busses along the path (working backwards) and add them to the unorphaned group:
                    Bus previousBus = current;
                    while(previousBus != null)
                    {
                        unorphaned.add(previousBus);
                        previousBus = previous.get(previousBus);
                    }
                    
                    // We've found our way back to the original group, so the orphaned bus is no longer orphaned:
                    orphanIterator.remove();
                    
                    // Don't need to keep searching, so go to the next orphan:
                    break;
                }
                
                // For the current node, consider all of its unvisited neighbors and calculate 
                // their tentative distances. For example, if the current node A is marked 
                // with a tentative distance of 6, and the edge connecting it with a neighbor B 
                // has length 2, then the distance to B (through A) will be 6+2=8. If this 
                // distance is less than the previously recorded tentative distance of B, then 
                // overwrite that distance. Even though a neighbor has been examined, it is not 
                // marked as "visited" at this time, and it remains in the unvisited set.
                Double currentDistance = distances.get(current);
                for (Line line : current.getLines())
                {
                    Bus from = line.getFromBus();
                    Bus to = line.getToBus();
                    Bus neighbour = (from == current) ? to : from;
                    if(unvisited.contains(neighbour))
                    {
                        double newDistance = currentDistance + line.impedance().abs();
                        double neighbourDistance = distances.get(neighbour);
                        if(newDistance < neighbourDistance)
                        {
                            distances.put(neighbour, newDistance);
                            previous.put(neighbour, current);
                        }
                    }
                }
                
                // When we are done considering all of the neighbors of the current node, mark 
                // the current node as visited and remove it from the unvisited set. A visited 
                // node will never be checked again; its distance recorded now is final and minimal.
                unvisited.remove(current);
                
                // Set the unvisited node marked with the smallest tentative distance 
                // as the next "current node".
                double minDistance = Double.MAX_VALUE;
                for (Bus bus : unvisited)
                {
                    double distance = distances.get(bus);
                    if(distance < minDistance)
                    {
                        minDistance = distance;
                        current = bus;
                    }
                }
            }
        }
        
        // Create a grid from the unorphaned busses:
        Grid grid = new Grid();
        for (Bus bus : unorphaned)
        {
            Unit copy = bus.copy(grid);
			grid.add(copy);
        }
        
        // Replace placeholder busses:
        for (Iterator<Line> iterator = grid.getLines().iterator(); iterator.hasNext();) 
        {
			Line line = iterator.next();
			
        	Bus from = line.getFromBus();
        	Bus to = line.getToBus();
        	
			if(from instanceof Bus.BusPlaceHolder)
        	{
        		Bus bus = grid.getBus(from.getName());
        		if(bus == null)
        		{
        			grid.getBus(to.getName()).removeChild(line);
        			grid.remove(line);
        			continue;
        		}
        		line.setFromBus(bus);
        	}

        	if(to instanceof Bus.BusPlaceHolder)
        	{
        		Bus bus = grid.getBus(to.getName());
        		if(bus == null)
        		{
        			line.getFromBus().removeChild(line);
        			grid.remove(line);
        			continue;
        		}
        		line.setToBus(bus);
        	}
		}
        
//        Grid grid = new Grid();
//        for (Bus bus : busses)
//        {
//            Unit copy = bus.copy(grid);
//			grid.add(copy);
//        }
        
        return grid;
    }

    private void search(Bus current, HashSet<Bus> orphaned, HashSet<Bus> unorphaned)
    {
        if(orphaned.remove(current))
        {
            unorphaned.add(current);
            
            for (Line line : current.getLines())
            {
                Bus from = line.getFromBus();
                Bus to = line.getToBus();
                if(from == current)
                    search(to, orphaned, unorphaned);
                else // to == current
                    search(from, orphaned, unorphaned);
            }
        }
    }

    public static GridFactory grid()
    {
        return new GridFactory();
    }
    
    public static class GridFactory
    {
        private Grid grid = new Grid();
        private Stack<Unit> parents = new Stack<Unit>();
        private Collection<Runnable> finalTasks = new LinkedHashSet<Runnable>();
        
        public GridFactory()
        {
        }
        
        private void add(String name, Unit unit, boolean isParent)
        {
            unit.setName(name);
            grid.add(unit);
            
            if(!parents.isEmpty() && parents.peek() instanceof Bus)
            {
                ((Bus)parents.peek()).addChild(unit);
            }
            
            if(isParent)
                parents.push(unit);
        }
        
        public GridFactory terminate()
        {
            parents.pop();
            return this;
        }

        public Grid grid()
        {
            for (Runnable task : finalTasks)
            {
                task.run();
            }
            return grid;
        }

        private void assertParentIsBus(String type)
        {
            if(parents.isEmpty() || !(parents.peek() instanceof Bus))
                throw new RuntimeException("Can only add a "+type+" to a parent of type Bus.");
        }
        
        public GridFactory Bus(String name)
        {
            Bus bus = new Bus(grid);
            if(!parents.isEmpty())
            {
                if(!(parents.peek() instanceof Line))
                    throw new RuntimeException("Can only add a bus to a parent of type Line.");
            
                Line line = (Line)parents.peek();
                line.setToBus(bus);
            }
            
            add(name, bus, true);
            
            return this;
        }
        
        public GridFactory Line(String name, double length, double resistancePerMetre, double inductancePerMetre)
        {
            assertParentIsBus("Line");
            Bus parent = (Bus) parents.peek();
            
            Line line = new Line();
            add(name, line, true);
            
            line.setLength(length);
            line.setResistancePerMetre(resistancePerMetre);
            line.setInductancePerMetre(inductancePerMetre);
            line.setFromBus(parent);
            
            return this;
        }
        
        public GridFactory Line(String name, final String toBusName, double length, double resistancePerMetre, double inductancePerMetre)
        {
            // Create the line which will set the parent to the created line:
            Line(name, length, resistancePerMetre, inductancePerMetre);
            final Line line = (Line)parents.pop();
            finalTasks.add(new Runnable()
            {
                @Override
                public void run()
                {
                    Bus toBus = grid.getBus(toBusName);
                    line.setToBus(toBus);
                }
            });
            
            return this;
        }
        
        public GridFactory Switch(String name, double resistance, double inductance)
        {
            assertParentIsBus("Switch");
            Bus parent = (Bus) parents.peek();
            
            Switch s = new Switch();
            add(name, s, true);
            
            s.setResistance(resistance);
            s.setInductance(inductance);
            s.setFromBus(parent);
            
            return this;
        }
        
        public GridFactory Transformer(String name, double baseVoltageRatio, double resistance, double inductance)
        {
            assertParentIsBus("Transformer");
            Bus parent = (Bus) parents.peek();
            
            Transformer t = new Transformer();
            add(name, t, true);
            
            t.setBaseVoltageRatio(baseVoltageRatio);
            t.setImpedance(new Complex(resistance, inductance));
            t.setFromBus(parent);
            
            return this;
        }

        public GridFactory swapEndpoints()
        {
            if(!(parents.peek() instanceof Transformer))
                throw new RuntimeException("Can only swapEndpoints on a Transformer.");
            
            final Transformer txr = (Transformer)parents.peek();
            finalTasks.add(new Runnable()
            {
                @Override
                public void run()
                {
                    txr.swapEndpoints();
                }
            });
            
            return this;
        }

        public GridFactory SlackSource(String name, double voltage, double phase, double resistance, double inductance)
        {
            SlackSource hvSource = new SlackSource();
            hvSource.setVoltage(ComplexUtils.polar2Complex(voltage, phase));
            hvSource.setResistance(resistance);
            hvSource.setInductance(inductance);
            
            add(name, hvSource, false);
            
            return this;
        }

        public GridFactory Capacitor(String name, double var)
        {
            Capacitor c = new Capacitor();
            c.setVAR(var);
            add(name, c, false);
            return this;
        }

        public GridFactory Load(String name)
        {
            add(name, new Load(), false);
            return this;
        }

        public GridFactory DistributedSource(String name)
        {
            add(name, new DistributedSource(), false);
            return this;
        }
        
        @Deprecated
		public GridFactory Storage(String name)
        {
            add(name, new Storage(), false);
            return this;
        }
        
        public GridFactory ControllableDemand(String name)
        {
        	add(name, new ControllableDemand(), false);
        	return this;
        }
    }

	public Bus getSlackBus() 
	{
		for (Bus bus : getBusses()) 
		{
			if(!bus.getSlackVoltage().equals(Complex.ZERO))
				return bus;
		}
		
		return null;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}
	
	@Override
	public String toString() 
	{
		StringBuffer sb = new StringBuffer();
		sb.append(getName());
		sb.append(":\n");
		for (Unit unit : units.values()) 
		{
			sb.append(unit.getName());
			sb.append("\n");
		}
		return sb.toString();
	}

	public double getBasePower()
	{
		return basePower;
	}

	public void setBasePower(double basePower)
	{
		this.basePower = basePower;
	}

	public double getBaseVoltage()
	{
		return baseVoltage;
	}

	public void setBaseVoltage(double baseVoltage)
	{
		this.baseVoltage = baseVoltage;
	}
}
