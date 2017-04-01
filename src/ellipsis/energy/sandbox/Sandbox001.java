package ellipsis.energy.sandbox;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import ellipsis.energy.calculation.AnalysisResults;
import ellipsis.energy.calculation.LoadFlowAnalyser;
import ellipsis.energy.grid.DistributedSource;
import ellipsis.energy.grid.Grid;
import ellipsis.energy.test.IEEE13BusGrid;
import ellipsis.util.MatrixHelper;

public class Sandbox001
{

	/**
	 * @param args
	 */
	public static void main(String[] args)
	{
		IEEE13BusGrid grid = new IEEE13BusGrid();
		grid.bus680.addChild(makeDG("DG680"));
		grid.bus675.addChild(makeDG("DG675"));		
		grid.bus684.addChild(makeDG("DG684"));
		grid.bus645.addChild(makeDG("DG645"));
		grid.bus611.addChild(makeDG("DG611"));
		grid.bus646.addChild(makeDG("DG646"));
		Set<String> names = new HashSet<>();
		names.add("DG680");
		names.add("DG675");
		names.add("DG684");
		names.add("DG645");
		names.add("DG611");
		names.add("DG646");
		
		System.out.println("Initial voltages:");
		AnalysisResults results = displayVoltages(grid);
		
		int n = results.getBusNumbers().size();
		RealVector v0 = voltages(results);
		RealVector s0 = powers(results);
		RealVector vplus = new ArrayRealVector(n, 1.05);
		RealMatrix L = dgSensitivities(results, names);
		RealMatrix L_inv = MatrixHelper.invert(L);
		RealVector sstar = s0.add(L_inv.operate(vplus).subtract(v0));
		
		System.out.println("Final voltages:");
		results = displayVoltages(grid);
	}

	private static RealMatrix dgSensitivities(AnalysisResults results, Set<String> names)
	{
		RealMatrix L = results.sensitivities();
		Map<String, Integer> busNumbers = results.getBusNumbers();
		for (String name_i : busNumbers.keySet())
		{
			int i = busNumbers.get(name_i);
			for (String name_j : busNumbers.keySet())
			{
				int j = busNumbers.get(name_j);
				if(!names.contains(name_j))
					L.setEntry(i, j, 0.0);
			}
		}
		return L;
	}

	private static RealVector powers(AnalysisResults results)
	{
		Map<String, Integer> busNumbers = results.getBusNumbers();
		int n = busNumbers.size();
		RealVector s = new ArrayRealVector(n*2);
		for (String name : busNumbers.keySet())
		{
			int i = busNumbers.get(name);
			Complex s_i = results.getBusPower(name);
			s.setEntry(i, s_i.getArgument());
			s.setEntry(i+n, s_i.abs());
		}
		return s;
	}

	private static RealVector voltages(AnalysisResults results)
	{
		Map<String, Integer> busNumbers = results.getBusNumbers();
		int n = busNumbers.size();
		RealVector v = new ArrayRealVector(n*2);
		for (String name : busNumbers.keySet())
		{
			int i = busNumbers.get(name);
			Complex v_i = results.getBusVoltage(name);
			v.setEntry(i, v_i.getArgument());
			v.setEntry(i+n, v_i.abs());
		}
		return v;
	}

	public static AnalysisResults displayVoltages(IEEE13BusGrid grid)
	{
		AnalysisResults results = analyse(grid);
		for (String name : results.getBusNumbers().keySet())
		{
			System.out.println(name+",\t"+results.getBusVoltage(name));
		}
		return results;
	}

	public static DistributedSource makeDG(String name)
	{
		DistributedSource dg = new DistributedSource();
		dg.setName(name);
		dg.setPmax(10e6);
		dg.setPowerOutput(1400e3, 0.0);
		return dg;
	}

	public static AnalysisResults analyse(Grid grid)
	{
		LoadFlowAnalyser lfa = new LoadFlowAnalyser(grid);
		lfa.setBasePower(IEEE13BusGrid.BASE_POWER);
		lfa.setBaseVoltage(IEEE13BusGrid.BASE_VOLTAGE);
		lfa.setIterations(100);
		lfa.setTargetError(1e-6);
		
		AnalysisResults results = lfa.analyse();
		if(!results.getDidConverge())
			System.err.println("Did not converge!");
		return results;
	}
}