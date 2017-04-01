package ellipsis.energy.vpp;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.apache.commons.math3.complex.Complex;

import ellipsis.energy.calculation.AnalysisResults;
import ellipsis.energy.calculation.LoadFlowAnalyser;
import ellipsis.energy.grid.Bus;
import ellipsis.energy.grid.Grid;
import ellipsis.util.SAProblem;
import ellipsis.util.SASolver;

@Deprecated
public class SAPartitionFactory
{
    private static final boolean TEST = false;
    
    private static final double COST_DIVISOR = 1e13; // TODO this needs to be determined programatically/statistically
    private static final String EMPTY_VPP_NAME = "VPP0";

// TODO this should be private (testing at the moment)
    public NetworkPartitionUtil util;

    
    //// Members ////
    
    private AnalysisResults results;
    private Grid grid;
    private double max_radius;
    private int maxSources;
    private boolean approximate;
    private List<Map<Bus, Double>> forecast;
    
    
    //// Constructors ////
    
    public SAPartitionFactory(Grid grid)
    {
        this.grid = grid;
    }
    
    
    //// Partitioning ////

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
    
    public Collection<VPP> partition(final double basePower, final double baseVoltage, final int maxVPPs)
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
                VPP toVPP = null;
                do
                {
                    int toVPPNumber = rand.nextInt(vppCount);
                    toVPP = vpps.get(toVPPNumber);
                } while(toVPP == fromVPP); // make sure we don't select the same as the from VPP
                
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
                
                if(TEST)
                    verify(maxVPPs, vpps, newList);

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
                        double vppCost = SAPartitionFactory.this.cost(busses);
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
        double radius = util.radius(busses_A);
        
        // If forecast is null just use bus powers:
        if(forecast == null)
        {
            S_A = powerImbalance(busses_A);
        }
        // Otherwise sum over forecast values and get the average cost:
        else
        {
            S_A = averageForecastPowerImbalance(busses_A);
        }
        
        double deviation = radius - max_radius;
        if(deviation < 0)
            deviation = 0;

        return square(S_A) * (1 + square(deviation));
//        return (Math.pow(2, S_A/1e6)-1) * (1 + square(deviation)); // possible alternative
    }

 // TODO this should be private (testing at the moment)
    public double averageForecastPowerImbalance(Collection<Bus> busses_A)
    {
        double imbalance = 0;
        for (Map<Bus, Double> netPowers : forecast)
        {
            double imbalance_t = 0;
            for (Bus bus : busses_A)
            {
                imbalance_t += netPowers.get(bus);
            }
            imbalance += Math.abs(imbalance_t);
        }

        int sampleSize = forecast.size();
        return imbalance/sampleSize;
    }

    public double powerImbalance(Collection<Bus> busses_A)
    {
        double imbalance = 0;
        for (Bus bus_i : busses_A)
        {
            double l = bus_i.getLoadPower().abs();
            double s = bus_i.getAvailablePower();
            imbalance += s - l;
        }
        return imbalance;
    }

    
    //// Accessors ////

    public AnalysisResults getResults()
    {
        return results;
    }

    public void setResults(AnalysisResults results)
    {
        if(results == null)
            throw new NullPointerException("results must not be null");
        
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

    public List<Map<Bus, Double>> getForecast()
    {
        return forecast;
    }

    /**
     * Sets the forecast. Note that when using a forecast the cost function will
     * integrate cost over all forecast values, and distance/radius calculations
     * must be approximate (i.e. NetworkPartitionUtil is set to approximate mode
     * and uses the admittance matrix imaginary values only for distance calculations
     * rather than using the full Jacobian matrix).
     * @param forecast
     * @param baseImpedance
     */
    public void setForecast(List<Map<Bus, Double>> forecast, double baseImpedance)
    {
        if(forecast == null)
            throw new NullPointerException("forecast must not be null");
        
        approximate = true;
        util = new NetworkPartitionUtil(grid, null, true, baseImpedance);
        
        this.forecast = forecast;
    }
    
    
    //// TEST ////
    
    public void verify(int maxVPPs, List<VPP> vpps, List<VPP> newList)
    {
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
    }
}
