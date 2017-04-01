package ellipsis.energy.sandbox;

import java.util.Map;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.FieldMatrix;

import ellipsis.energy.calculation.AnalysisResults;
import ellipsis.energy.calculation.LoadFlowAnalyser;
import ellipsis.energy.test.IEEE13BusGrid;

public class Sandbox010
{
	public static void main(String[] args)
	{
		new Sandbox010().run();
	}
	
	public void run()
	{
		IEEE13BusGrid grid = new IEEE13BusGrid();
		
		LoadFlowAnalyser lfa = new LoadFlowAnalyser(grid);
        lfa.setBasePower(grid.getBasePower());
        lfa.setBaseVoltage(grid.getBaseVoltage());
        AnalysisResults results = lfa.analyse();
        
        double minVarg = Double.MAX_VALUE;
        double minYarg = Double.MAX_VALUE;
        double maxVarg = Double.MIN_VALUE;
        double maxYarg = Double.MIN_VALUE;
        double avVarg = 0;
        double avYarg = 0;
        
        double minVabs = Double.MAX_VALUE;
        double minYabs = Double.MAX_VALUE;
        double maxVabs = Double.MIN_VALUE;
        double maxYabs = Double.MIN_VALUE;
        double avVabs = 0;
        double avYabs = 0;
        
        double minCombo = Double.MAX_VALUE;
        double maxCombo = Double.MIN_VALUE;
        double avCombo = 0;
        
        Map<String, Integer> numbers = results.getBusNumbers();
        for (String busName : numbers.keySet())
		{
			Complex v_i = results.getBusVoltage(busName);
			double varg = v_i.getArgument();
			double vabs = v_i.abs();
			
			minVarg = Math.min(minVarg, varg);
			maxVarg = Math.max(maxVarg, varg);
			avVarg += varg;
			
			minVabs = Math.min(minVabs, vabs);
			maxVabs = Math.max(maxVabs, vabs);
			avVabs += vabs;
		}
        avVarg /= numbers.size();
        avVabs /= numbers.size();
        
        FieldMatrix<Complex> Y = results.getAdmittanceMatrix().Y;
		int dimension = Y.getColumnDimension();
		int count = 0;
		for(int i = 0; i < dimension; ++i)
		{
			for(int j = 0; j < dimension; ++j)
			{
				if(i == j)
					continue;
				
				Complex Y_ij = Y.getEntry(i, j);
				double Yarg = Y_ij.getArgument();
				double Yabs = Y_ij.abs();
				
				minYarg = Math.min(minYarg, Yarg);
				maxYarg = Math.max(maxYarg, Yarg);
				avYarg += Yarg;
				
				minYabs = Math.min(minYabs, Yabs);
				maxYabs = Math.max(maxYabs, Yabs);
				avYabs += Yabs;
				
				 ++count;
			}
		}
		avYabs /= count;
		avYarg /= count;
		
		count = 0;
		for (String busName_i : numbers.keySet())
		{
			int i = numbers.get(busName_i);
			for (String busName_j : numbers.keySet())
			{
				int j = numbers.get(busName_j);
				
				if(i == j)
					continue;
				
				Complex Y_ij = Y.getEntry(i, j);
				double Yarg = Y_ij.getArgument();
				double varg_i = results.getBusVoltage(busName_i).getArgument();
				double varg_j = results.getBusVoltage(busName_j).getArgument();
				double combo = Yarg - varg_i + varg_j;
				
				minCombo = Math.min(minCombo, combo);
				maxCombo = Math.max(maxCombo, combo);
				avCombo += combo;
				
				 ++count;
			}
		}
		avCombo /= count;
        
        System.out.println("Min Vabs = "+minVabs);
        System.out.println("Max Vabs = "+maxVabs);
        System.out.println("Av Vabs = "+avVabs);

        System.out.println("Min Yabs = "+minYabs);
        System.out.println("Max Yabs = "+maxYabs);
        System.out.println("Av Yabs = "+avYabs);

        System.out.println("pi/2 = "+Math.PI/2);
        System.out.println("Min Varg = "+minVarg);
        System.out.println("Max Varg = "+maxVarg);
        System.out.println("Av Varg = "+avVarg);

        System.out.println("Min Yarg = "+minYarg);
        System.out.println("Max Yarg = "+maxYarg);
        System.out.println("Av Yarg = "+avYarg);

        System.out.println("Min Combo = "+minCombo);
        System.out.println("Max Combo = "+maxCombo);
        System.out.println("Av Combo = "+avCombo);
        
        showY(Y, false);
	}

	private void showY(FieldMatrix<Complex> Y, boolean abs)
	{
		int dimension = Y.getColumnDimension();
		for(int i = 0; i < dimension; ++i)
		{
			for(int j = 0; j < dimension; ++j)
			{
				if(abs)
					System.out.print(Y.getEntry(i, j).abs());
				else
					System.out.print(Y.getEntry(i, j).getArgument());
				System.out.print(',');
			}
			System.out.println();
		}
	}
}