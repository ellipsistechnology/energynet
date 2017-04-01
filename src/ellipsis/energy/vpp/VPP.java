package ellipsis.energy.vpp;

import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import org.apache.commons.math3.complex.Complex;

import ellipsis.energy.calculation.AnalysisResults;
import ellipsis.energy.calculation.BasicPowerGenerationOptimiser;
import ellipsis.energy.grid.Bus;
import ellipsis.energy.grid.Load;
import ellipsis.energy.grid.Source;
import ellipsis.energy.scenario.Regulator;

/**
 * A Virtual Power Plant.
 * @author bmillar
 *
 */
public class VPP implements Regulator
{
//    private static final double TARGET_PF = 0.9;
    
	private Set<Source> sources = new HashSet<Source>();
	private Set<Load> loads = new HashSet<Load>();
	private Set<Bus> busses = new HashSet<Bus>();
    private double basePower;
    private String name;
    
    
    /* Control/Regulation */

    @Override
    public void regulate(AnalysisResults results)
    {
        if(basePower <= 0)
            throw new RuntimeException("basePower must be positive for regulate");
        
        // Regulate all sources:
        Complex loadPower = loadPower();
        int sourceCount = sources.size();
        if(sourceCount > 0)
        {
            double targetQ = loadPower.getImaginary()/sourceCount;
            double targetP = loadPower.getReal()/sourceCount;
            for (Source source : sources)
            {
                double availablePower = source.getAvailablePower();
                if(targetP < availablePower)
                {
                    double availableQ = Math.sqrt(availablePower*availablePower - targetP*targetP);
                    double Q = Math.min(targetQ, availableQ);
                    source.setReactivePowerOutput(Q);
                }
            }
        }
    }
    
    public void optimise(double targetTotalPowerOutput)
    {
        BasicPowerGenerationOptimiser optimiser = new BasicPowerGenerationOptimiser();
        optimiser.addAllSources(getSources());
        optimiser.optimise(targetTotalPowerOutput);
    }
    
    public void balancePower()
    {   
        // Assign power equally:
        Complex load = loadPower();
                //new Complex(loadPower().getReal());
        Set<Source> sources = new HashSet<Source>(getSources());
        Complex loadPerSource = load.divide(sources.size());
        assignSourcePower(loadPerSource, sources);
        
        // Assign remaining power:
        Complex remainder = load.subtract(sourcePower());
        final double PRECISION = 0.001;
        while(remainder.abs() > PRECISION && sources.size() > 0)
        {
            loadPerSource = remainder.divide(sources.size());
            assignSourcePower(loadPerSource, sources);
            remainder = sourcePower().subtract(load);
        }
    }

    public void assignSourcePower(Complex loadPerSource, Set<Source> sources)
    {
        for (Iterator<Source> iterator = sources.iterator(); iterator.hasNext();)
        {
            Source source = iterator.next();

            source.setPowerOutput(loadPerSource);
            if(source.getPowerOutput().abs() < loadPerSource.abs())
                iterator.remove();
        }
    }
    
    
    /* Organisation Operations */
    

    /**
     * Adds all the busses from the given VPP bus but does not remove them from that VPP.
     * @param virtualBus
     */
    public void merge(VPP vpp2)
    {
        addBusses(vpp2.busses);
    }
    
    
    /* Information */

    public Complex sourcePower()
    {
        Complex s = Complex.ZERO;
        for (Source source : sources)
        {
            s = s.add(source.getPowerOutput());
        }
        return s;
    }

    public Complex loadPower()
    {
        Complex l = Complex.ZERO;
        for (Load load : loads)
        {
            l = l.add(load.getLoad());
        }
        return l;
    }
    
    public Complex netPower()
    {
        return sourcePower().subtract(loadPower());
    }
    
    public Complex availablePower()
    {
        Complex s = Complex.ZERO;
        for (Source source : sources)
        {
            s = s.add(source.getAvailablePower());
        }
        return s;
    }
    
    
    /* Accessors */


    public boolean add(Source source)
    {
        busses.add(source.getParent());
        
        return sources.add(source);
    }

    public void addSources(Collection<Source> sources2)
    {
        for (Source source : sources2)
        {
            busses.add(source.getParent());
        }
        
        sources.addAll(sources2);
    }

    public boolean contains(Source source)
    {
        return sources.contains(source);
    }

    public Set<Source> getSources()
    {
        return sources;
    }

    public void setSources(Set<Source> sources)
    {
        // Remove current source busses:
        for (Source source : this.sources)
        {
            busses.remove(source.getParent());
        }
        
        // Add new source busses:
        for (Source source : sources)
        {
            busses.add(source.getParent());
        }
        
        this.sources = sources;
    }

    public double getBasePower()
    {
        return basePower;
    }

    public void setBasePower(double basePower)
    {
        this.basePower = basePower;
    }

    public Set<Load> getLoads()
    {
        return loads;
    }

    public void setLoads(Set<Load> loads)
    {
        // Remove current load busses:
        for (Load load : this.loads)
        {
            busses.remove(load.getParent());
        }
        
        // Add new load busses:
        for (Load load : loads)
        {
            busses.add(load.getParent());
        }
        
        this.loads = loads;
    }

    public boolean add(Load load)
    {
        busses.add(load.getParent());
        
        return loads.add(load);
    }

    public void addLoads(Collection<Load> loads2)
    {
        for (Load load : loads2)
        {
            busses.add(load.getParent());
        }
        
        loads.addAll(loads2);
    }

    public boolean contains(Load load)
    {
        return loads.contains(load);
    }
    
    public Set<Bus> getBusses()
    {
        return busses;
    }

    public String getName()
    {
        return name;
    }

    public void setName(String name)
    {
        this.name = name;
    }
    
    @Override
    public String toString()
    {
        return name;
    }

    public void remove(Bus bus)
    {
        if(busses.remove(bus))
        {
            for (Source source : bus.getSources())
            {
                sources.remove(source);
            }
            
            for (Load load : bus.getLoads())
            {
                loads.remove(load);
            }
        }
    }

    public void addBus(Bus bus)
    {
        busses.add(bus);
        sources.addAll(bus.getSources());
        loads.addAll(bus.getLoads());
    }
    
    public void addBusses(Collection<Bus> busses)
    {
        for (Bus bus : busses)
        {
            addBus(bus);
        }
    }

    public int busCount()
    {
        return busses.size();
    }
}
