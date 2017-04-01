package ellipsis.energy.vpp;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.apache.commons.math3.complex.Complex;

import ellipsis.energy.calculation.AnalysisResults;
import ellipsis.energy.calculation.LoadFlowAnalyser;
import ellipsis.energy.grid.Bus;
import ellipsis.energy.grid.Grid;
import ellipsis.util.SAProblem;
import ellipsis.util.SASolver;

public class VPPPartitionFactory
{
    private static final double COST_DIVISOR = 1e16;
    private static final String EMPTY_VPP_NAME = "VPP0";
    public boolean singleMovePerIteration = false;
    private NetworkPartitionUtil util;
    
    
    //// Members ////
    
    private AnalysisResults results;
    private Grid grid;
    private double max_radius;
    private int maxSources;
    private boolean approximate;
    
    
    //// Constructors ////
    
    public VPPPartitionFactory(Grid grid)
    {
        this.grid = grid;
    }
    
    
    //// Partitioning ////

    public Collection<VPP> partition(double basePower, double baseVoltage)
    {
        recalculateResults(basePower, baseVoltage);
        Collection<VPP> vpps = partition();
        for (VPP vpp : vpps)
        {
            vpp.setBasePower(basePower);
        }
        return vpps;
    }


    public void recalculateResults(double basePower, double baseVoltage)
    {
        LoadFlowAnalyser lfa = new LoadFlowAnalyser(grid);
        lfa.setIterations(10000);
        lfa.setTargetError(0.00000000001);
        lfa.setBasePower(basePower);
        lfa.setBaseVoltage(baseVoltage);
        AnalysisResults results = lfa.analyse();
        setResults(results);
        if(!results.getDidConverge())
            throw new RuntimeException("Did not converge!");
    }

    /**
     * Ensure {@link #recalculateResults(double, double)} or {@link #setResults(AnalysisResults)} is called prior to partitioning.
     * @return
     */
    @SuppressWarnings("null")
    public Collection<VPP> partition()
    {
        List<VPP> vpps = createVPPs(grid, results);
        
        VPP vpp_A = null;
        VPP vpp_B = null;
        double minDistance = Double.MAX_VALUE;
        while(true)
        {
            for(int i = 0; i < vpps.size(); ++i)
            {
                VPP vpp_i = vpps.get(i);
                for(int j = i+1; j < vpps.size(); ++j) // j = i+1: Ignore busses visited.
                {
                    VPP vpp_j = vpps.get(j);
                    if(vpp_i != vpp_j)
                    {
                        double d = util.distance(vpp_i, vpp_j);
                        int sourceCount = vpp_i.busCount()+vpp_j.busCount();
                        if(d < minDistance && sourceCount <= maxSources)
                        {
                            minDistance = d;
                            vpp_A = vpp_i;
                            vpp_B = vpp_j;
                        }
                    }
                }
            }
            
            if(minDistance <= max_radius)
            {
                vpp_A.merge(vpp_B);
                vpps.remove(vpp_B);
                
                vpp_A = null;
                vpp_B = null;
                minDistance = Double.MAX_VALUE;
            }
            else
            {
                break;
            }
        }
        
    	return vpps;
    } // END partition

    private List<VPP> createVPPs(Grid grid, AnalysisResults results) 
    {
        List<VPP> vpps = new ArrayList<VPP>();
    	for (Bus bus : grid.getBusses())
    	{
    	    // Don't create VPPs for slack busses or busses without loads or sources:
    	    int sourceCount = bus.getSources().size();
            int loadCount = bus.getNonCapacitorLoads().size();
            boolean noSlackVoltage = bus.getSlackVoltage().equals(Complex.ZERO);
            if(noSlackVoltage && (sourceCount > 0 || loadCount > 0))
    	    {
        	    VPP vpp = new VPP();
        	    vpp.addBus(bus);
        	    vpps.add(vpp);
    	    }
    	}
    	return vpps;
    }

    /**
     * 
     * @param vpps
     * @return true if changes were made and {@link #singleMovePerIteration} is true, false otherwise.
     */
    public boolean revisePartitions(Collection<VPP> vpps)
    {
        // Attempt to move each bus to each other VPP:
        double minCost;
        Bus busToMove = null;
        VPP fromVPP = null;
        VPP toVPP = null;
        do
        {
            minCost = Double.MAX_VALUE;
            
            for (VPP vpp_A : vpps)
            {
                for (Bus bus_i : vpp_A.getBusses())
                {
                    for (VPP vpp_B : vpps)
                    {
                        if(vpp_A == vpp_B /*|| vpp_B.getBusses().size() >= maxSources*/) // NOTE max_sources or max_busses?
                            continue;
                        
                        // Current system cost:
                        Set<Bus> busses_A = new HashSet<Bus>(vpp_A.getBusses());
                        Set<Bus> busses_B = new HashSet<Bus>(vpp_B.getBusses());
                        double cost_before = cost(busses_A) + cost(busses_B);
                        
                        // Move bus_i from A to B:
                        busses_A.remove(bus_i);
                        busses_B.add(bus_i);
                        
                        // New system cost:
                        double cost_after = cost(busses_A) + cost(busses_B);

                        // Change in cost:
                        double cost = cost_after - cost_before;
                        
                        if(cost < minCost)
                        {
                            minCost = cost;
                            busToMove = bus_i;
                            fromVPP = vpp_A;
                            toVPP = vpp_B;
                        }
                    }
                }
            }
            
            // If there is a benefit to moving the bus then move it:
            if(minCost < 0)
            {
                fromVPP.remove(busToMove);
                toVPP.addBus(busToMove);
                
                if(fromVPP.getBusses().isEmpty())
                    vpps.remove(fromVPP);

                if(singleMovePerIteration)
                    return true;
            }
        } while(minCost < 0);
        
        return false;
    }
    
    public Collection<VPP> partitionV2(final double basePower, final double baseVoltage, final int maxVPPs)
    {
        // Initialise:
        recalculateResults(basePower, baseVoltage);
        
        // Define Simulated Annealing problem:
        SAProblem<List<VPP>> problem = new SAProblem<List<VPP>>()
        {
            Random rand = new Random(0);
            int vppIndex = 1; // Used for unique naming of VPPs
            HashMap<String, Double> vppCosts = new HashMap<String, Double>();
            
            @Override
            public void prepareNextSolution(double T)
            {
                // Nothing to do here since each iteration is the same.
            }
            
            @Override
            public List<VPP> nextSolution(List<VPP> vpps, int attempt)
            {
                // Randomly choose a bus to move:
                int vppCount = vpps.size();
                VPP fromVPP = null;
                int fromVPPBusCount = 0;
                while(fromVPPBusCount == 0)
                {
                    int fromVPPNumber = rand.nextInt(vppCount);
                    fromVPP = vpps.get(fromVPPNumber);
                    fromVPPBusCount= fromVPP.busCount();
                }
                int busNumber = rand.nextInt(fromVPPBusCount);
                
                // Choose VPP to move to:
                // Note: Can't move to an empty bus if we already have out limit of VPPs.
                //       Empty VPP does not count towards limit.
                VPP toVPP = null;
                do
                {
                    int toVPPNumber = rand.nextInt(vppCount);
                    toVPP = vpps.get(toVPPNumber);
                } while(toVPP == fromVPP);// || vppCount > maxVPPs && toVPP.getBusses().isEmpty());
                
                // Create new list of VPPs:
                List<VPP> newList = new ArrayList<VPP>();
                newList.addAll(vpps);
                newList.remove(fromVPP);
                newList.remove(toVPP);
                
                // Copy from VPP and remove chosen bus:
                VPP newFromVPP = new VPP();
                newFromVPP.setName(nextVPPName());
                newFromVPP.setBasePower(basePower);
                List<Bus> fromBusses = new ArrayList<Bus>(fromVPP.getBusses());
                Bus bus = fromBusses.remove(busNumber);
                newFromVPP.addBusses(fromBusses);
                if(!newFromVPP.getBusses().isEmpty())
                    newList.add(newFromVPP);
                
                // Copy to VPP and add chosen bus:
                VPP newToVPP = new VPP();
                newToVPP.setName(nextVPPName());
                newToVPP.setBasePower(basePower);
                newToVPP.addBusses(toVPP.getBusses());
                newToVPP.addBus(bus);
                newList.add(newToVPP);
                
                // If "To VPP" was empty add new empty VPP (unless we already have max number of VPPs)
                // or if new "From VPP" is now empty (it won't have been added):
                if((toVPP.getBusses().isEmpty() || newFromVPP.getBusses().isEmpty()) && newList.size() < maxVPPs)
                {
                    addEmptyVPP(basePower, newList);
                }
                
// Verify new list:
for (VPP vpp : vpps) // has all busses from original list?
{
    for (Bus b : vpp.getBusses())
    {
        boolean contains = false;
        for (VPP newVPP : newList)
        {
            if(newVPP.getBusses().contains(b))
            {
                contains = true;
                break;
            }
        }
        if(!contains)
            throw new RuntimeException();
    }
}

for (Bus b : grid.getBusses()) // has all grid busses?
{
    if(b.getName().equals("Bus 01") || b.getName().equals("Bus 02") || b.getName().equals("OLTC Bus"))
        continue;
    
    boolean contains = false;
    for (VPP vpp : newList)
    {
        if(vpp.getBusses().contains(b))
        {
            contains = true;
            break;
        }
    }
    if(!contains)
        throw new RuntimeException();
}

for (VPP vpp1 : newList) // no duplicates?
{
    for (Bus b1 : vpp1.getBusses())
    {
        for (VPP vpp2 : newList)
        {
            if(vpp1 != vpp2 && vpp2.getBusses().contains(b1))
                throw new RuntimeException();
        }
    }
}

boolean contains = false;
for (VPP vpp : newList)
{
    if(vpp.getName().equals(EMPTY_VPP_NAME))
    {
        contains = true;
        break;
    }
}
if(!contains && newList.size() != maxVPPs)
    throw new RuntimeException();

                return newList;
            }

            public String nextVPPName()
            {
                return "VPP"+(vppIndex++);
            }

            public void addEmptyVPP(final double basePower, List<VPP> newList)
            {
                VPP emptyVPP = new VPP();
                emptyVPP.setName(EMPTY_VPP_NAME);
                emptyVPP.setBasePower(basePower);
                newList.add(emptyVPP);
            }
            
            @Override
            public List<VPP> getInitialSolution()
            {
                List<VPP> vpps = new ArrayList<VPP>();
                
                VPP initialVPP = new VPP();
                initialVPP.setName(nextVPPName());
                initialVPP.setBasePower(basePower);
                for (Bus bus : grid.getBusses())
                {
                    // Don't add slack busses or busses without loads or sources:
                    int sourceCount = bus.getSources().size();
                    int loadCount = bus.getNonCapacitorLoads().size();
                    boolean noSlackVoltage = bus.getSlackVoltage().equals(Complex.ZERO);
                    if(noSlackVoltage && (sourceCount > 0 || loadCount > 0))
                    {
                        initialVPP.addBus(bus);
                    }
                }
                vpps.add(initialVPP);
                
                addEmptyVPP(basePower, vpps);

                return vpps;
            }
            
            @Override
            public double cost(List<VPP> vpps)
            {
                double cost = 0;
                for (VPP vpp : vpps)
                {
                    String name = vpp.getName();
                    
                    if(EMPTY_VPP_NAME.equals(name))
                        continue;
                    
                    // Lazily add cost to cache:
                    if(!vppCosts.containsKey(name))
                    {
                        Set<Bus> busses = vpp.getBusses();
                        double vppCost = VPPPartitionFactory.this.cost(busses);
                        vppCosts.put(name, vppCost/COST_DIVISOR);
                    }
                    
                    // Add cached cost:
                    cost += vppCosts.get(name);
                }
                return cost;
            }

            @Override
            public List<VPP> finish(List<VPP> vpps)
            {
                // Remove empty vpp:
                for (Iterator<VPP> iterator = vpps.iterator(); iterator.hasNext();)
                {
                    VPP vpp = (VPP) iterator.next();
                    if(vpp.getBusses().isEmpty())
                    {
                        iterator.remove();
                        break;
                    }
                }
                
                return vpps;
            }
        };
        
        // Solve:
        SASolver solver = new SASolver();
        solver.setMaxAttempts(grid.getBusses().size());
        solver.setMinT(0.0000001);
        return solver.solve(problem);
    }

    private double square(double d)
    {
        return d*d;
    }

    /**
     * Sums the power imbalance of each bus scaled by the square of the 
     * deviation from the target radius.
     * @param busses_A
     * @return
     */
    public double cost(Collection<Bus> busses_A)
    {
        double S_A = 0;
        double radius = radius(busses_A);
        
        for (Bus bus_i : busses_A)
        {
            double l = bus_i.getLoadPower().abs();
            double s = bus_i.getAvailablePower();
            S_A += s - l;
        }
        
        double deviation = radius - max_radius;
        if(deviation < 0)
            deviation = 0;

        return square(S_A) * (1 + square(deviation));
    }
    
    public double systemCost(Collection<VPP> vpps)
    {
        double cost = 0;
        for (VPP vpp : vpps)
        {
            cost += cost(vpp.getBusses());
        }
        return cost;
    }

    public double radius(Collection<Bus> busses)
    {
        switch (busses.size())
        {
        case 0:
            throw new RuntimeException("VPP of 0 busses cannot have it's radius calculated.");
        case 1:
            return 0;
        case 2:
            Iterator<Bus> iterator = busses.iterator();
            return 0.5*util.distance(iterator.next(), iterator.next());
        }
        
        // For 3 or more busses find the maximum radius:
        double radius = 0;
        for (Bus bus : busses)
        {
            radius = Math.max(radius, util.distance(bus, busses));
        }
        return radius;
    }

    
    //// Accessors ////

    public AnalysisResults getResults()
    {
        return results;
    }

    public void setResults(AnalysisResults results)
    {
        this.results = results;
        util = new NetworkPartitionUtil(grid, results, approximate, 0);
    }

    public double getRadius()
    {
        return max_radius;
    }

    public void setRadius(double radius)
    {
        this.max_radius = radius;
    }

    public int getMaxSources()
    {
        return maxSources;
    }

    public void setMaxSources(int maxSources)
    {
        this.maxSources = maxSources;
    }

    public boolean isApproximate()
    {
        return approximate;
    }

    public void setApproximate(boolean approximate)
    {
        this.approximate = approximate;
    }
}
