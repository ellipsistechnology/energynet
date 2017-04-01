package ellipsis.energy.calculation;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.FieldLUDecomposition;
import org.apache.commons.math3.linear.FieldMatrix;

import ellipsis.energy.grid.Bus;
import ellipsis.energy.grid.Grid;
import ellipsis.energy.grid.Load;

public class FaultAnalyser
{
    private Grid grid;
    private double basePower;
    private double baseVoltage;
    private int iterations;
    private double targetError;
    private Complex faultImpedance = Complex.ZERO;
    
    
    /* Constructors */
    
    public FaultAnalyser(Grid grid)
    {
        this.grid = grid;
    }
    
    
    /* Analysis */
    
    /**
     * Performs the fault analysis based on the results of a load flow analysis.
     * The load flow analysis will use the base power, base voltage, iterations and
     * target error (if set) specified in this FaultAnalyser.
     */
    public Map<String,AnalysisResults> analyse()
    {
        LoadFlowAnalyser lfa = new LoadFlowAnalyser(grid);
        lfa.setBasePower(basePower);
        lfa.setBaseVoltage(baseVoltage);
        if(iterations != 0)
            lfa.setIterations(iterations);
        if(targetError != 0)
            lfa.setTargetError(targetError);
        
        return analyse(lfa.analyse());
    }
    
    /**
     * Performs the fault analysis based on the assumption that all voltages
     * from the given result are valid the instant before the fault occurs. 
     */
    public Map<String,AnalysisResults> analyse(AnalysisResults results)
    {
        FieldMatrix<Complex> admittances = results.getAdmittanceMatrix().getAdmittances();
for (int i = 0; i < admittances.getRowDimension(); i++)
{
    for (int j = 0; j < admittances.getColumnDimension(); j++)
    {
        System.out.print(admittances.getEntry(i, j));
        System.out.print(", ");
    }
    System.out.println();
}

        FieldLUDecomposition<Complex> yDecomp = new FieldLUDecomposition<Complex>(admittances);
        FieldMatrix<Complex> Z = yDecomp.getSolver().getInverse();
        Map<String, Integer> busNumbers = results.getBusNumbers();
        Collection<Bus> busses = grid.getBusses();
        
        // Replace bus loads with impedances:
        for (Bus bus : busses)
        {
            for (Load load : bus.getLoads())
            {
                int i = busNumbers.get(bus.getName()); // Bus index in admittance/impedance matrix
                
                Complex loadPower = load.getLoad().divide(basePower); // Bus load p.u.
                Complex v = results.getBusVoltage(bus.getName());     // Bus voltage p.u.
                
                // V = IZ
                // S = VI*
                // :. S = VV*/Z*
                // :. Z = VV*/S*
                Complex loadImpedance = v.conjugate().multiply(v).divide(loadPower.conjugate());
                Z.addToEntry(i, i, loadImpedance);
            }
        }
        
        // Apply a fault at each bus:
        Map<String,AnalysisResults> faultResults = new HashMap<String, AnalysisResults>();
        for (Bus bus_k : busses)
        {
//System.out.println("\nFault on bus '"+bus_k+"'");
            int k = busNumbers.get(bus_k.getName());
            Complex Vs[] = new Complex[busses.size()];
            
            // Pre-fault voltage at bus k:
            Complex Vk0 = results.getBusVoltage(bus_k.getName());
            
            // Fault current and voltage:
            Complex Ik = Vk0.divide(Z.getEntry(k, k).add(faultImpedance));
//System.out.println("\tIk = "+Ik);
//System.out.println("\tVk = "+Vk);
//System.out.println();
            Vs[k] = faultImpedance.multiply(Ik);
            for (Bus bus_i : busses)
            {
                int i = busNumbers.get(bus_i.getName());
                if(i == k)
                    continue;
                
                // Pre-fault voltage on bus i:
                Complex Vi0 = results.getBusVoltage(bus_i.getName());
                
                // Fault voltage and current on bus i:
                Vs[i] = Vi0.subtract(Z.getEntry(i, k).multiply(Ik));
//System.out.println();
            }
            
            // Setup results:
            AnalysisResults faultResult = new AnalysisResults();
            FieldLUDecomposition<Complex> zDecomp = new FieldLUDecomposition<Complex>(Z);
            faultResult.admittances = new AdmittanceMatrix(zDecomp.getSolver().getInverse());
            faultResult.busNumbers = busNumbers;
            faultResult.didConverge = true;
            faultResult.applyVoltages(Vs, grid);
            faultResults.put(bus_k.getName(), faultResult);
        }
        
        return faultResults;
    }

    
    /* Accessors */
    
    public Grid getGrid()
    {
        return grid;
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

    public int getIterations()
    {
        return iterations;
    }

    public void setIterations(int iterations)
    {
        this.iterations = iterations;
    }

    public double getTargetError()
    {
        return targetError;
    }

    public void setTargetError(double targetError)
    {
        this.targetError = targetError;
    }


    public Complex getFaultImpedance()
    {
        return faultImpedance;
    }


    public void setFaultImpedance(Complex faultImpedance)
    {
        this.faultImpedance = faultImpedance;
    }
}
