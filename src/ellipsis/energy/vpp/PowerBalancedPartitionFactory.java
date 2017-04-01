package ellipsis.energy.vpp;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ellipsis.energy.grid.Bus;
import ellipsis.energy.grid.Grid;

public class PowerBalancedPartitionFactory
{
    private Grid grid;
    private Map<VPP, Double> currentRadii = new HashMap<VPP, Double>(); // r_A
    private Map<VPP, Double> currentNetPowers = new HashMap<VPP, Double>(); // \Delta S_A
//    private Map<String, Double> radii = new HashMap<String, Double>();
//    private Map<String, Double> netPower = new HashMap<String, Double>();
    private Map<Bus, Double> powerBalances = new HashMap<Bus, Double>();
    private NetworkPartitionUtil util;
    private double maxRadius;
    
    public PowerBalancedPartitionFactory(Grid grid, double baseImpedance)
    {
        this.grid = grid;
        util =  new NetworkPartitionUtil(grid, null, true, baseImpedance); // must be approximate since we don't have any results at this time
    }
    
    /**
     * NOTE Because vpps can't be merged, this algorithm will potentially create many new VPPs if it encounters
     * mostly loads before DGs.
     * NOTE Change this to start with a single partition and to divide bus-by-bus in the same manner as the original partitioning algorithm.
     * @param forecast Mapping of busses to expected net power output/demand at each time period.
     * @return A set of partitioned busses (i.e. a set of VPPs).
     */
    public Set<VPP> partition(List<Map<Bus, Double>> forecast)
    {
        Set<VPP> vpps = new HashSet<VPP>();
        VPP emptyVPP = addEmptyVPP(vpps); // One empty VPP is always needed (will be removed at the end)
        
        Collection<Bus> busses = grid.getBusses();
        Bus[] aBusses = busses.toArray(new Bus[busses.size()]);
        Arrays.sort(aBusses);
        
        // Add busses one by one to VPPs:
        for (int i = 0; i < aBusses.length; i++)
        {
            Bus bus = aBusses[i];
            if(bus.getSlackVoltage().abs() > 0)
                continue;
            
            VPP minCostVPP = null;
            double minCost = Double.MAX_VALUE;
            double minCostRadius = 0;
            double minCostPowerBalance = 0;
            
            // Find the best VPP to merge with:
            for (VPP vpp : vpps)
            {
                // Find value of adding this bus to this VPP:
                double radius = radius(vpp, bus);
                double powerBalance = currentNetPowers.get(vpp) + integratedPowerBalance(forecast, bus);
                double cost = radiusCost(radius)*powerBalance*powerBalance;
                if(cost < minCost)
                {
                    minCost = cost;
                    minCostVPP = vpp;
                    minCostRadius = radius;
                    minCostPowerBalance = powerBalance;
                }
            }
            
            // If the selected VPP is currently empty, create a new empty VPP before adding the bus to it:
            if(minCostVPP == emptyVPP)
                emptyVPP = addEmptyVPP(vpps);
            
            // Add the bus to the best case:
            currentNetPowers.put(minCostVPP, minCostPowerBalance);
            currentRadii.put(minCostVPP, minCostRadius);
            minCostVPP.addBus(bus);
            minCostVPP.setName(bus.getName()+","+minCostVPP.getName());
        }
        
        // Remove the empty VPP and return the set:
        vpps.remove(emptyVPP);
        return vpps;
    }

    public double integratedPowerBalance(List<Map<Bus, Double>> forecast, Bus bus)
    {
        if(powerBalances.containsKey(bus))
            return powerBalances.get(bus);
        
        double sum = 0;
        
        for (Map<Bus, Double> balances : forecast)
        {
            sum += balances.get(bus);
        }
        
        powerBalances.put(bus, sum);
        
        return sum;
    }

    private double radiusCost(double radius)
    {
        if(radius < maxRadius)
        {
            return 1;
        }
        else
        {
            double deviation = maxRadius - radius;
            return 1 + deviation*deviation;
        }
    }

    public VPP addEmptyVPP(Set<VPP> vpps)
    {
        VPP emptyVPP = new VPP();
        emptyVPP.setName("");
        currentNetPowers.put(emptyVPP, 0.0);
        currentRadii.put(emptyVPP, 0.0);
        vpps.add(emptyVPP);
        
        return emptyVPP;
    }

    /**
     * 
     * @param vpp
     * @param bus
     * @return The radius of the given VPP after adding the given bus.
     */
    private double radius(VPP vpp, Bus bus)
    {
        Set<Bus> busses = vpp.getBusses();
        switch (busses.size())
        {
        case 0:
            return 0;
        case 1:
            Iterator<Bus> iterator = busses.iterator();
            return 0.5*util.distance(iterator.next(), bus);
        }
        
        double radius = currentRadii.get(vpp);
        double distance = util.distance(bus, busses);
        return Math.max(radius, distance);
    }
    
    
    //// ACCESSORS ////

    public double getMaxRadius()
    {
        return maxRadius;
    }

    public void setMaxRadius(double maxRadius)
    {
        this.maxRadius = maxRadius;
    }

//    private Set<VPP> initialVPPs(Collection<Bus> busses)
//    {
//        HashSet<VPP> vpps = new HashSet<VPP>();
//        
//        // Add a VPP for each bus:
//        for (Bus bus : busses)
//        {
//            VPP vpp = new VPP();
//            vpp.addBus(bus);
//            vpps.add(vpp);
//        }
//        
//        // Add an empty VPP:
//        vpps.add(new VPP());
//        
//        return vpps;
//    }
}