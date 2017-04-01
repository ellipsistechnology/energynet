package ellipsis.energy.calculation;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import ellipsis.energy.grid.Bus;
import ellipsis.energy.grid.Grid;
import ellipsis.energy.grid.Line;

public class AnalysisResults
{
    private static final int busResultTypeCount = 2;
    private static final int V = 0;
    private static final int S = 1;
    private Map<String, Complex[]> busResults;
    
    private static final int lineResultTypeCount = 4;
    private static final int I = 0;
    private static final int Sij = 1;
    private static final int Sji = 2;
    private static final int SL = 3;
    private Map<String, Complex[]> lineResults;

    Map<String, Integer> busNumbers;
    private int slackIndex;
    private Bus slackBus;
    AdmittanceMatrix admittances; // bus to bus p.u. admittances
    
    boolean didConverge;
    private RealMatrix jacobian;
    
    AnalysisResults() {}
    
    public static AnalysisResults initWith(Map<String, Integer> busNumbers, AdmittanceMatrix admittances)
    {
    	AnalysisResults r = new AnalysisResults();
    	r.admittances = admittances;
    	r.busNumbers = busNumbers;
    	return r;
    }
    
    public AdmittanceMatrix getAdmittanceMatrix()
    {
        return admittances;
    }

    public boolean getDidConverge()
    {
        return didConverge;
    }

    public Complex getBusVoltage(String busName)
    {
        return busResults.get(busName)[V];
    }
    
    public Complex getBusPower(String busName)
    {
        return busResults.get(busName)[S];
    }
    
    public Complex getLineCurrent(String lineName)
    {
        return lineResults.get(lineName)[I];
    }
    
    public Complex getLinePower(String lineName)
    {
        return lineResults.get(lineName)[Sij];
    }
    
    public Complex getLinePowerReverse(String lineName)
    {
        return lineResults.get(lineName)[Sji];
    }
    
    public Complex getLineLoss(String lineName)
    {
        return lineResults.get(lineName)[SL];
    }
    
    public Complex[] pathVoltages(String[] pathBusNames)
    {
        Complex[] voltages = new Complex[pathBusNames.length];
        for (int i = 0; i < voltages.length; i++)
        {
            Complex[] results = this.busResults.get(pathBusNames[i]);
            voltages[i] = results[V];
        }
        return voltages;
    }
    
    public Complex[] pathPowers(String[] pathBusNames)
    {
        Complex[] powers = new Complex[pathBusNames.length];
        for (int i = 0; i < powers.length; i++)
        {
            Complex[] results = this.busResults.get(pathBusNames[i]);
            powers[i] = results[S];
        }
        return powers;
    }
    
    public Complex busCurrent(Bus bus)
    {
        Complex current = Complex.ZERO;
        for (Line line : bus.getLines())
        {
            if(line.getFromBus().equals(bus))
                current = current.add(lineResults.get(line.getName())[I]);
            else
                current = current.add(lineResults.get(line.getName())[I].negate());
        }
        return current;
    }

    public Map<String, Integer> getBusNumbers()
    {
        return busNumbers;
    }

    public void applyVoltages(Complex[] Vk, Grid grid)
    {
        // Calculate line currents and flows,
        // and set up line->results map:
        lineResults = new HashMap<String, Complex[]>();
        for (Line line : grid.getLines())
        {
            Bus from = line.getFromBus();
            Bus to = line.getToBus();
            int fromIndex = busNumbers.get(from.getName());
            int toIndex = busNumbers.get(to.getName());
            Complex Vi = Vk[fromIndex];
            Complex Vj = Vk[toIndex];
            Complex Yij = admittances.get(fromIndex, toIndex).negate();
            Complex Iij = Vi.subtract(Vj).multiply(Yij);

            Complex Iji = Iij.negate();
            Complex Sij_value = Vi.multiply(Iij.conjugate());
            Complex Sji_value = Vj.multiply(Iji.conjugate());
            
            Complex[] results_array = new Complex[lineResultTypeCount];
            results_array[I] = Iij;
            results_array[Sij] = Sij_value;
            results_array[Sji] = Sji_value;
            results_array[SL] = Sij_value.add(Sji_value);
            lineResults.put(line.getName(), results_array);
        }
        
        // Set up bus->results map:
        busResults = new HashMap<String, Complex[]>();
        for (String name : busNumbers.keySet())
        {
            Bus bus = grid.getBus(name);
            if(bus == null) // account for grid being a subnetwork
            	continue;
            
            Integer index = busNumbers.get(name);
            Complex voltage = Vk[index];
            Complex current = busCurrent(bus);
            
            Complex[] results_array = new Complex[busResultTypeCount];
            results_array[V] = voltage;
            results_array[S] = voltage.multiply(current.conjugate());
            busResults.put(name, results_array);
            
            if(bus.getSlackVoltage().abs() > 0)
            {
            	slackBus = bus;
                slackIndex = index;
            }
        }
    }
    
    public RealMatrix sensitivities()
    {
		DecompositionSolver solver = new LUDecomposition(jacobian).getSolver();
        RealMatrix sensitivities = solver.getInverse();
        return sensitivities;
    }

	public RealVector voltages()
	{
		RealVector v = new ArrayRealVector(2*(busNumbers.size()-1));
		int offset = v.getDimension()/2;
		for (String bus : busNumbers.keySet())
		{
			int i = busNumbers.get(bus);
			if(i == slackIndex)
				continue;
			if(i > slackIndex)
				--i;
			
			Complex v_i = getBusVoltage(bus);
			
			v.setEntry(i, v_i.getArgument());
			v.setEntry(i+offset, v_i.abs());
		}
		return v;
	}
    
    public Bus getSlackBus()
    {
    	return slackBus;
    }

    public void setJacobian(RealMatrix jacobian)
    {
        this.jacobian = jacobian;
    }
    
    public RealMatrix getJacobian()
    {
        return jacobian;
    }

    public int getSlackIndex()
    {
        return slackIndex;
    }

    public void setSlackIndex(int slackIndex)
    {
        this.slackIndex = slackIndex;
    }
	
	public int indexWithoutSlack(String bus)
	{
		int i = busNumbers.get(bus);
		if(i > slackIndex)
			--i;
		return i;
	}
}