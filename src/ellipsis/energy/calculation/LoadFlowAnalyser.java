package ellipsis.energy.calculation;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.FieldMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import ellipsis.energy.grid.Bus;
import ellipsis.energy.grid.Grid;
import ellipsis.energy.grid.Line;

public class LoadFlowAnalyser
{
    /* Analysis parameters */

    private Grid grid;
    private double basePower;
    private double baseVoltage;
    private int iterations = 20;
    private double targetError = 0.00001;
    
    
    /* Analysis variables */
    
    private boolean[] isSlackSource;
    private double baseImpedance;
    private Complex[] power; // p.u.
    private Complex[] voltages; // p.u.
    
    
    /* Analysis Results */

    private AnalysisResults results;
    
    
    /* Constructors */
    
    public LoadFlowAnalyser(Grid grid)
    {
        this.grid = grid;
    }
    
    
    /* Analysis */
    
    private void init()
    {
        // Reinitialise results:
        results = new AnalysisResults();
        
        // Get a sorted list of busses:
        List<Bus> busses = new ArrayList<Bus>(grid.getBusses());
        
        // Set bus loads and slack sources (per unit):
        int nBusses = busses.size();
        power = new Complex[nBusses];
        isSlackSource = new boolean[nBusses];
        voltages = new Complex[nBusses];
        for (int i = 0; i < nBusses; i++)
        {
            Bus bus = busses.get(i);

            // Set bus load and power output:
            power[i] = bus.getLoadPower().divide(basePower).negate();
            power[i] = power[i].add(bus.getDistributedSourcePower().divide(basePower));
            
            // Set default bus voltage:
            Complex slackVoltage = bus.getSlackVoltage().divide(baseVoltage);
            isSlackSource[i] = !slackVoltage.equals(Complex.ZERO);
            if(isSlackSource[i])
                voltages[i] = slackVoltage;
            else
                voltages[i] = Complex.ONE; // Default for iteration one
        }
        
        // Assign numbers to busses:
        results.busNumbers = new HashMap<String, Integer>();
        for (int i = 0; i < nBusses; ++i)
        {
            Bus bus = busses.get(i);
            results.busNumbers.put(bus.getName(), i);
        }
        
        // Set up admittance matrix:
        List<Line> lines = new ArrayList<Line>(grid.getLines());
        baseImpedance = baseVoltage*baseVoltage/basePower;

        // Create admittance matrix:
        results.admittances = new AdmittanceMatrix(busses, results.busNumbers, lines, baseImpedance);
    }
    
    /**
     * Analyses the grid but makes no changes to either the grid or the analyser.
     * @return A mapping of busses with their respective voltages.
     */
    public AnalysisResults analyse()
    {
    	return analyse(null);
    }
    
    public AnalysisResults analyse(Complex[] initialVoltages)
    {
        init();
        
        if(initialVoltages != null)
        	voltages = initialVoltages;
        
        results.didConverge = false;

        // Iterate to find voltages:
        boolean hasError = true;
        PowerFlowSolver solver = new NewtonRaphsonPowerFlowSolver(power, isSlackSource, results.admittances);
        for(int k = 0; k < iterations && hasError; ++k)
        {
            voltages = solver.solve(voltages);
            hasError = solver.getError() > targetError;
        }
        
        results.applyVoltages(voltages, grid);
        results.setJacobian(((NewtonRaphsonPowerFlowSolver)solver).getJacobian());
        
        results.didConverge = !hasError;
        
        return results;
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
    
    
    /* Debugging */
    
    public void printBusNumbers()
    {
        System.out.println("\nBus Mapping:");
        for (String bus : results.busNumbers.keySet())
        {
            System.out.println(bus+":\t"+results.busNumbers.get(bus));
        }
        System.out.println();
    }
    
    public void printAdmittanceMatrix()
    {
        System.out.println("\nAdmittance matrix:");
		printMatrix(results.admittances.Y);
    }

	protected void printMatrix(@SuppressWarnings("rawtypes") FieldMatrix matrix)
	{
		for (int i = 0; i < matrix.getRowDimension(); i++)
        {
            for (int j = 0; j < matrix.getColumnDimension(); j++)
            {
                System.out.print(matrix.getEntry(i, j));
                System.out.print("\t");
            }
            System.out.println();
        }
	}

	protected void printMatrix(RealMatrix matrix)
	{
		for (int i = 0; i < matrix.getRowDimension(); i++)
        {
			System.out.print("{");
            for (int j = 0; j < matrix.getColumnDimension(); j++)
            {
                System.out.print(matrix.getEntry(i, j));
                System.out.print(",\t");
            }
            System.out.println("},");
        }
	}


	public void printJacobianMatrix()
	{
		System.out.println("\nJacobian Matrix:");
		printMatrix(results.getJacobian());
	}


	public void printSensitivityMatrix()
	{
		System.out.println("\nSensitivity Matrix:");
		RealMatrix jacobian = results.getJacobian();
		DecompositionSolver solver = new LUDecomposition(jacobian).getSolver();
        RealMatrix sensitivities = solver.getInverse();
		printMatrix(sensitivities);
	}
}
